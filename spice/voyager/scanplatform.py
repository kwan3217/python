"""
We have a scan platform CK file with discrete pointing instances at each
"""

import spiceypy as cspice
import spiceypy.utils.support_types
import numpy as np

data_path="/home/jeppesen/workspace/Data/spice/"

cspice.furnsh(data_path+"generic/lsk/naif0012.tls")
cspice.furnsh(data_path+"Voyager/fk/vg2_v02.tf")
cspice.furnsh(data_path+"Voyager/sclk/vg200022.tsc")
ScanPlatformCK=data_path+"Voyager/ck/vg2_jup_version1_type1_iss_sedr.bc"
cspice.furnsh(data_path+"Voyager/ck/vgr2_super.bc")
cspice.furnsh(ScanPlatformCK)
print((1/3).hex())
print(float.fromhex('0x1.5555555555555p-2'))
print(float.fromhex('0x1.5555555555555p-2')-(1/3))

#eventually we will figure out how to get ckcov to work
#ckids=spiceypy.ckobj(ScanPlatformCK)
#print(ckids[0])
#cover=cspice.ckcov(ScanPlatformCK,ckids[0],False,"INTERVAL",0.0,"SCLK")

csv=ScanPlatformCK[:-3]+".csv"
print(csv)
header=None
sclks=[]
ets=[]
with open(csv) as inf:
    for line in inf:
        parts=line.split(",")
        if len(parts)<10:
            continue
        if header is None:
            header=parts
            print(header)
            continue
        sclks.append(float.fromhex(parts[3]))
        ets.append(float.fromhex(parts[6]))


# Checks if a matrix is a valid rotation matrix.
def isRotationMatrix(R):
    Rt = np.transpose(R)
    shouldBeIdentity = np.dot(Rt, R)
    I = np.identity(3, dtype=R.dtype)
    n = np.linalg.norm(I - shouldBeIdentity)
    return n < 1e-6


def Rx(theta=None, c=None, s=None):
    if c is None:
        c = np.cos(theta)
        s = np.sin(theta)
    return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])

def Ry(theta=None, c=None, s=None):
    if c is None:
        c = np.cos(theta)
        s = np.sin(theta)
    return np.array([[c,0,s],[0,1,0],[-s,0,c]])

def Rz(theta=None, c=None, s=None):
    if c is None:
        c = np.cos(theta)
        s = np.sin(theta)
    return np.array([[c,-s,0],[s,c,0],[0,0,1]])

def mtx_to_euler(R,verbose=False):
    """
    Calculate the Euler angles in the Proper ZYZ Extrinsic form, compatible
    with co-elevation, co-azimuth and twist for the Voyager spacecraft.
    See https://omoikane.kwansystems.org/wiki/index.php/Matrix_to_Euler_Angle
    :param R: rotation matrix. If this is not a proper rotation matrix, you
              will silently get a wrong answer.
    :return: a tuple (twist,co-elevation,co-azimuth) of Euler angles

    The proper ZYZ extrinsic form means that from the home position (~90deg pitch up
    from squared-up Von Karman position) you:

    * rotate around reference Z axis by twist
    * rotate around reference Y axis by co-elevation
    * rotate around reference Z axis again by co-azimuth
    """
    ce= R[2,2]
    se=np.sqrt(1-ce**2)
    ca= R[0,2]/se
    sa= R[1,2]/se
    ct=-R[2,0]/se
    st= R[2,1]/se

    #Check Pythagorean identity
    if verbose:
        print(ca**2+sa**2," Pythag(a) should be 1")
        print(ct**2+st**2," Pythag(t) should be 1")

    #Check other components
    if verbose:
        print(R[0,0]," r00 should be ",ca*ce*ct-sa*st)
        print(R[0,1]," r01 should be ",-ca*ce*st-ct*sa)
        print(R[1,0]," r10 should be ",ca*st+ce*ct*sa)
        print(R[1,1]," r11 should be ",ca*ct-ce*sa*st)
    twist=np.arctan2(st,ct)
    coelevation=np.arctan2(se,ce)
    coazimuth=np.arctan2(sa,ca)

    #Check that the result is a rotation matrix with no reflection
    R=Rz(c=ca,s=sa) @ Ry(c=ce,s=se) @ Rz(c=ct, s=st)
    if verbose:
        print(np.linalg.det(R)," det should be 1")
    if np.linalg.det(R)<0:
        raise ValueError("We hit one")
    return (twist,coelevation,coazimuth)

et0=cspice.str2et("1979-07-07 11:39:45 TDB")
et1=cspice.str2et("1979-07-11 05:11:11 TDB")

import datetime
import matplotlib.dates

tt=[] #Twists
ee=[] #Coelevations
aa=[] #Coazimuths
eet=[] #Teph's
j2000=matplotlib.dates.date2num(datetime.datetime(2000,1,1,12,0,0))
tts=[] #Timestamps
for (i,(et,sclk)) in enumerate(zip(ets,sclks)):
    if et<et0 or et>et1:
        continue
    print(i,cspice.etcal(et))
    try:
        M = cspice.pxform("VG2_SCAN_PLATFORM", "VG2_AZ_EL", et)
    except:
        print("Spice error (probably no super_ck)")
        continue
    t,e,a=mtx_to_euler(M)
    tt.append(t)
    ee.append(e)
    aa.append(a)
    eet.append(et)
    tts.append(et / 86400 + j2000)
tt=np.array(tt)
ee=np.array(ee)
aa=np.array(aa)

import matplotlib.pyplot as plt
plt.plot_date(tts,np.degrees(tt),'r+',label='twist')
plt.plot_date(tts,np.degrees(ee),'g+',label='coelevation')
plt.plot_date(tts,np.degrees(aa),'b+',label='coazimuth')
plt.legend()
plt.show()

msopck_txt="""
Voyager 2 Jupiter encounter scan platform kernel

This kernel represents the orientation of the scan platform
relative to the az/el frame. It is based upon two kernels:

* One is derived from the ISS SEDR, whose raw data gives the
  position of the scan platform at each NAC image, but no
  data in between. This has been processed using optical 
  navigation, and therefore is relative to the star catalog
  and not the orientation of the spacecraft. This gives pointing
  only at discrete points in time.
* The other is the orientation of the spacecraft relative
  to the stars. This has the measured position at a period
  of as low as 48s over almost the entire mission, including
  the Jupiter near encounter.
  
The present kernel is constructed by figuring the orientation
of the scan platform relative to the az/el frame, and then
resolving that orientation into azimuth, elevation, and twist
(around the instrument boresight) angles. Since Voyager doesn't
have a twist actuator, the twist angle should be zero. Experience
shows that it is near zero almost all of the time.

Once these angles are obtained, we do the following:

* Force the twist angle to zero. This will not change the 
direction of the boresight, but will change the orientation
of the image around the boresight.
* Figure out the slew. The spacecraft always used one axis 
of motion at a time, in any of three slew rates. It always 
hits its mark in time to take the image. We will simulate this by:
  - Figure the distance to move in azimuth and elevation
  - Assume that the scan platform hits its mark X seconds before
    the exposure (time of data point in NAC kernel). We will try
    x=5sec for now.
  - Before the scan platform reaches its mark, it slews in one axis
    at a time at the slowest rate allowed by the distance and time
    between two NAC kernel data points. We add up the distance to
    slew in each direction, and divide by each slew rate to get
    the time required to do the slew at that rate. The longest
    slew time which fits in between the NAC kernel points is used --
    this would put minimum wear on the mechanisms in real life,
    and generates slower slews which are more likely to be seen in
    an animation. We call the slew time Y, with Y_el the elevation
    portion and Y_az the azimuth portion.
  - Generate a zero-motion position from the previous NAC point to 
    X+Y_el+Y_az from the current NAC time
  - Generate an elevation slew starting at the previous NAC point
    with a constant rate around the Y axis, starting at X+Y_el+Y_az
    before the current NAC time
  - Generate an azimuth slew starting at the azimuth of the previous
    NAC point and the elevation of the current point, with a constant
    rate around the Z axis, starting at X+Y_az before the current 
    NAC time. This will result in hitting the mark X seconds before
    the NAC time.
  - Generate a zero-motion position from X seconds before the current
    NAC time to the NAC time.
    
Since all of this is done with Euler azimuth/elevation angles, the result
will be in the az/el frame of the spacecraft.

We call this instrument -32200, so that we can use both it and the nominal
scan platform -32100 at the same time.
    
Times are in ET. I would have preferred to specify in SCLK, but
I am not sure how to interpolate. At all of the ETs of the photos,
the ET matches exactly the SCLK in the source kernel given the 
SCLK kernel listed below, raw from 

\\begindata
INPUT_DATA_TYPE = 'EULER ANGLES'
INPUT_TIME_TYPE = 'ET'
ANGULAR_RATE_PRESENT= 'YES'
CK_TYPE = 2
INSTRUMENT_ID=-32500
REFERENCE_FRAME_NAME='VG2_AZ_EL'
PRODUCER_ID='C. Jeppesen, Kwan Systems'
FRAMES_FILE_NAME='%s/fk/vg2_v02.tf'
SCLK_FILE_NAME='%s/sclk/vg200022.tsc'
LSK_FILE_NAME='%s/lsk/naif0012.tls'
EULER_ANGLE_UNITS='DEGREES'
EULER_ROTATIONS_ORDER=('Z','Y','Z')
EULER_ROTATIONS_TYPE='SPACE'
ANGULAR_RATE_FRAME='REFERENCE'
""" % ('../../../Data/spice/Voyager','../../../Data/spice/Voyager','../../../Data/spice/generic')

with open('v2j_msopck_header.txt','w') as ouf:
    print(msopck_txt,file=ouf)

allowed_slew_rates=[0.08,0.33,1.0]
slew_rates=[]
with open('v2j_msopck_data.txt','w') as ouf:
    for i_nac in range(1,len(eet)):
        M0 = cspice.pxform("VG2_SCAN_PLATFORM", "VG2_AZ_EL", eet[i_nac-1])
        (twist0,coel0,coaz0)=mtx_to_euler(M0)
        twist0=np.degrees(twist0)
        coel0=np.degrees(coel0)
        coaz0=np.degrees(coaz0)
        M1 = cspice.pxform("VG2_SCAN_PLATFORM", "VG2_AZ_EL", eet[i_nac])
        (twist1,coel1,coaz1)=mtx_to_euler(M1)
        twist1=np.degrees(twist1)
        coel1=np.degrees(coel1)
        coaz1=np.degrees(coaz1)
        d_el=np.abs(coel1-coel0)
        d_az=np.abs(coaz1-coaz0)
        d_t=np.abs(eet[i_nac]-eet[i_nac-1])
        d=d_el+d_az
        x=5
        time_available=d_t-x
        req_spd=d/time_available
        i_slew=None
        found=False
        for i_slew,slew_rate in enumerate(allowed_slew_rates):
            if slew_rate>=req_spd:
                found=True
                break
        if not found:
            print(i_nac,"Warning: Slew is too fast")
            slew_rate=req_spd
        slew_rates.append(slew_rate)
        print(i_nac,req_spd)
        y_el=d_el/slew_rate
        y_az=d_az/slew_rate
        #Set up first segment - non-motion from previous time to current time-(x+y_el+y_az)
        t00=eet[i_nac-1]
        t01=eet[i_nac]-(x+y_el+y_az)
        if t00>t01:
            raise ValueError("Time travel?")
        elif t00<t01:
            print(("%25.15e %25.15e"+
                 " 0 %25.15e %25.15e"+
                 " 0 0 0")%(t00,t01,
                           coel0,coaz0)
                  ,file=ouf)
        #Set up second segment - motion around elevation axis
        t10=eet[i_nac]-(x+y_el+y_az)
        t11=eet[i_nac]-(x+y_az)
        if t10<t01:
            raise ValueError("Time travel?")
        elif t10>t11:
            raise ValueError("Time travel?")
        elif t10<t11:
            print(("%25.15e %25.15e"+
                 " 0 %25.15e %25.15e"+
                 " 0 %25.15e 0")%(t10,t11,
                           coel0,coaz0,slew_rate),file=ouf)
        #Set up third segment - motion around azimuth axis
        t20=eet[i_nac]-(x+y_az)
        t21=eet[i_nac]-(x)
        if t20<t11:
            raise ValueError("Time travel?")
        elif t20>t21:
            raise ValueError("Time travel?")
        elif t20<t21:
            print(("%25.15e %25.15e"+
                 " 0 %25.15e %25.15e"+
                 " 0 0 %25.15e")%(t20,t21,
                           coel1,coaz0,slew_rate),file=ouf)
        #Set up fourth segment - non-motion from current time-x to current time
        t30=eet[i_nac]-x
        t31=eet[i_nac]
        if t30<t21:
            raise ValueError("Time travel?")
        elif t30>t31:
            raise ValueError("Time travel?")
        elif t30<t31:
            print(("%25.15e %25.15e"+
                 " 0 %25.15e %25.15e"+
                 " 0 0 0")%(t30,t31,
                           coel1,coaz1),file=ouf)
plt.plot(slew_rates)
plt.show()
import os
import subprocess

try:
    os.remove(data_path+"Voyager/ck/v2j_slew.bc")
except FileNotFoundError:
    pass #no error, file is already not present
subprocess.call("~/bin/msopck v2j_msopck_header.txt v2j_msopck_data.txt "+data_path+"Voyager/ck/v2j_slew.bc",shell=True)
