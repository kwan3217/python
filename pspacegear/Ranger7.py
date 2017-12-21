'''
Created on Dec 6, 2017

@author: chrisj

Ranger was before the era of leap seconds, and Spice uses an incorrect version of UTC prior to the 
beginning of the leap second table. "Absurd Accuracy is our Obsession" 

These are instants which have the same *name* as the UTC (actually GMT) times given in the Ranger report,
but are on a uniform time scale which must be converted. TDT is ahead of TAI by 32.184s, so to get the
ET of the TAI time with the same name, add 32.184s to the ET of the TDT time. TAI is ahead of GMT by deltaAT. 
Ranger 7 was launched on mjd 38505, when the following row in tai-utc.dat was valid

 1964 APR  1 =JD 2438486.5  TAI-UTC=   3.3401300 S + (MJD - 38761.) X 0.001296 S

so deltaAT on that day was 3.1379540

'''

import csv
from collections import namedtuple
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize as opt
import spiceypy as cspice
import os
import bmw

old=os.getcwd()
os.chdir('../../Data/spice/Ranger/')
cspice.furnsh('Ranger7Background.tm')
os.chdir(old)

r_moon=1735.46        #Reported radius of Moon at impact point. This happens to be about 1700m below the LRO/LROC reference sphere
#mu_moon=4904.8695     #Value from Vallado of gravitational parameter of Moon in km and s
#mu_earth=398600.4415  #Value from Vallado of gravitational parameter of Earth in km and s
mu_moon =cspice.gdpool("BODY301_GM",0,1)[0]  #DE431 value of gravitational parameter of Moon in km and s
mu_earth=cspice.gdpool("BODY399_GM",0,1)[0] #DE431 value of gravitational parameter of Earth in km and s
#Documented Ranger 7 impact points
#
def gmt_to_et(gmt):
    """
Calculate ET from given GMT, bypassing spice leap second kernels. Spice does not properly handle the "rubber second"
era. Ranger 7 was launched during this era, so we have to deal with it.

Since the times have a rated accuracy of 5ms (even though they appear to have a sub-millisecond moment-to-moment
consistency) we don't need to worry about such things as the change in MJD over the approach time. However, "Absurd
Accuracy is Our Obsession", so since we *can* do it, we *will* do it.
    
Ranger 7 was flown when the following row in tai-utc.dat was valid

1964 APR  1 =JD 2438486.5  TAI-UTC=   3.3401300 S + (MJD - 38761.) X 0.001296 S

:param str gmt: GMT to convert. Any time acceptable to cspice.str2et is acceptable, except that
                there must not be a time scale tag like UTC on it. Suggested is 1964-Jul-31 13:25:12.345
:return: Spice ET of given time
:rtype: float64
"""
    #Tell Spice that the incoming times are TDT. This is a lie, but we will correct it piece by piece.
    #Tell Spice that it is TDT rather than TDB. Spice itself will return an ET value, which 
    #includes the TDT to TDB (Spice ET) correction
    gmt_num=cspice.str2et(gmt+" TDT")
    #Calculate the MJD. This doesn't specify what timescale (TAI, TDB, UTC, etc) is used, but
    #at our required level of accuracy, it doesn't matter.
    jd=cspice.timout(gmt_num,"JULIAND.#########")
    mjd=float(jd)-2400000.5
    #Calculate TAI-UTC from the formula row above 
    tai_utc=3.3401300+(mjd-38761.0)*0.001296
    #Add the correction to get TAI
    tai=gmt_num+tai_utc
    #Add the ET-TAI correction to get Spice ET
    et=tai+32.184
    return et

def floatN(x):
    """
    Convert the input to a floating point number if possible, but return NaN if it's not
    :param x: Value to convert
    :return: floating-point value of input, or NaN if conversion fails
    """
    try:
        return float(x)
    except ValueError:
        return 0.0 #float('NaN')

image_a_tuple=namedtuple('image_a_tuple',['PhotoNum','GMT',
                                          'sc_alt','sc_lat','sc_lon',
                                          'p2_lat','p2_lon','p2_srange',
                                          'v','pth','az',
                                          'p1_lat','p1_lon','p1_srange',
                                          'azn'])

def readImageA(lonofs=0.089):
    """
    Read the Image A table.
    
    :param float lonofs: offset in longitude to subtract from all longitudes in the file. This is used to make the trajectory
                         match a given impact longitude from another source.
    :rtype: list of namedtuple
    """
    result=[]
    with open('Ranger 7 Trajectory - Image A Parameters.csv','r') as inf:
        reader=csv.reader(inf)
        #Read past the first header line
        for row in reader:
            break
        #Read past the second header line
        for row in reader:
            break
        #Read the rows
        for row in reader:
            if not row[0].isdigit():
                if row[0]!="Impact":
                    break
            result.append(image_a_tuple(floatN(row[ 0]), #PhotoNum
                                               row[ 1] , #GMT
                                        floatN(row[ 2]), #sc_alt
                                        floatN(row[ 3]), #sc_lat
                                        floatN(row[ 4])-lonofs, #sc_lon
                                        floatN(row[ 5]), #p2_lat
                                        floatN(row[ 6])-lonofs, #p2_lon
                                        floatN(row[ 7]), #p2_srange
                                        floatN(row[ 8]), #v
                                        floatN(row[ 9]), #pth
                                        floatN(row[10]), #az
                                        floatN(row[11]), #p1_lat
                                        floatN(row[12])-lonofs, #p1_lon
                                        floatN(row[13]), #p1_srange
                                        floatN(row[14]))) #azn
    return result

def ray_sphere_intersect(r0,v,re):
    """
    Calculate the point of intersection between a given ray r(t)=r0+v*t and a sphere dot(r,r)=re**2
    
    :param numpy vector r0: initial point of ray
    :param numpy vector v:  direction of ray. If this is a unit vector, then the units of t returned will be
                            the same as the units of r0
    :param float re: radius of sphere
    :return: tuple, first element is ray parameter, second is intersect point.
    
    Parametric equation for a ray
      r(t)=r0+v*t
    Equation for a sphere of radius r_m
      dot(r,r)=re**2
    solve simultaneous equations
      dot(r0+v*t,r0+v*t)=re**2
      (rx+vx*t)**2+                           #expand dot product to components
      (ry+vy*t)**2+
      (rz+vz*t)**2=r_m**2
      (r_+v_*t)**2=r_**2+2*r_*v_*t+v_**2*t**2 #square each term
      rx**2+2*rx*vx*t+vx**2*t**2+             #substitute terms
      ry**2+2*ry*vy*t+vy**2*t**2+
      rz**2+2*rz*vz*t+vz**2*t**2=r_m**2
      t**2*(vx**2  +vy**2  +vz**2  )+         #gather into coefficients of P(t)
      t**1*(2*rx*vx+2*ry*vy+2*rz*vz)+
      t**0*(rx**2  +ry**2  +rz**2  )=re**2
      t**2*dot(v,v)+                          #recognize dot products
      t**1*2*dot(r0,v)+
      t**0*dot(r,r)=re**2
    So, the quadratic coefficients are:
      * A=dot(v,v)
      * B=2*dot(r0,v)
      * C=dot(r0,r0)-re**2
    """
    A=np.dot(v,v)
    B=2*np.dot(r0,v)
    C=np.dot(r0,r0)-re**2
    D=B**2-4*A*C
    #Since A is positive, it is always the case that using the negative sign will give the
    #lower root, which is what we want. If this root is negative, then the spacecraft is
    #inside the sphere or the sphere is behind the spacecraft.
    t=(-B-np.sqrt(D))/2*A 
    #Finish using the ray equation to find the coordinates of p1 from the spacecraft pos/vel
    return (t,r0+v*t)

def processImageA(image_a,plot_check_values=False):
    """
    Convert Image A table to usable state vectors, and calculate the check values

    The tables include the spacecraft position in lat/lon/alt coordinates. Altitude is defined to be zero at impact, so the
    reference surface is a sphere centered on the Moon's center of mass and has radius equal to the radius at impact. The
    latitude and longitude are referenced to the Mean-Earth Polar system, as realized by this report. This gives a final
    impact point about 2700m away from where LRO found the impact crater.

    All reticle marks are numbered, with point 2 being the center mark. Point 1 is the velocity vector as described above.
    For all such marks, the table includes the latitude and longitude on the reference surface, the distance from the
    spacecraft to that point, several azimuths including the azimuth between the image vertical and true North.

    Table "Point 1" is the point on the Lunar surface that the spacecraft is moving directly towards, based on
    its instantaneous velocity vector in a lunar body-fixed frame. This is the point that doesn't move as the
    spacecraft approaches the moon. The image appears to zoom in centered around this point.

    Table "Point 2" is the point on the Lunar surface covered by the center reticle mark. It also has latitude, longitude,
    and slant range

    Much of the data which was entered is redundant - point 1 is completely determined by the spacecraft position and velocity
    vector. Slant ranges are always calculable from the point latitudes and longitudes and the spacecraft position.
    These values were transcribed anyway, to validate the position and velocity transcription. Ideally, the calculated
    values will be exactly equal to the table values, but due to limited precision in the table, particularly the latitudes
    and longitudes, the values will be inconsistent on the few-meter level. Any discrepancy of more than 10m drew my
    attention, and any larger than 20m indicated a transcription error, which was corrected.
    
    :param array of namedtuple image_a: Rows from original table
    :param bool plot_check_values: If true, plot the check value residuals
    :rtype: tuple
    :return: First element is numpy array of position vectors, one row for each row in the table, 3 columns
             Second element is numpy array of velocity vectors, one row for each row in the table, 3 columns
             Third element
    """
    rs=[]
    vs=[]
    #dsrange2 is the difference between the table slant range to point 2 and that calculated from the spacecraft lat/lon/alt
    #and point 2 lat/lon
    dsrange2s=[]
    #dsrange1a is the difference between the table slant range to point 1 and that calculated from the spacecraft lat/lon/alt
    dsrange1as=[]
    #dsrange1b is the difference between the table slant range to point 1 and that calculated from the spacecraft pos/vel
    #and the quadratic method (quadratic parameter t is the calculated distance between the ray origin and the ray/sphere
    #intersect point)
    dsrange1bs=[]
    #dsrange1c is the distance between the table point 1 calculated from lat/lon and that calculated by the quadratic
    #method. This isn't a difference in srange like the others are, but it is measured in the same units. However, dsrange1c
    #will always be positive, while the other measures can be positive or negative.
    dsrange1cs=[]
    #dv is the difference between the table velocity and that calculated by dividing the distance from the previous row's
    #position to this row's position by the difference in time. Keep track of last row position in r_last, use constant
    #5.12s as dt. Note that this will be biased from zero because it doesn't take into account the acceleration of gravity
    #over the time step.
    dvs=[float('NaN')]
    r_last=None
    t_last=None
    #Size of 1 millidegree of latitude at the current spacecraft altitude. This is an idea of the precision we can expect
    #in using vectors with latitudes and longitudes specified in millidegree precision.
    mds=[]
    ts=[]
    for row in image_a:
        #Zenith vector
        rbar=np.array([np.cos(np.radians(row.sc_lat))*np.cos(np.radians(row.sc_lon)),
                       np.cos(np.radians(row.sc_lat))*np.sin(np.radians(row.sc_lon)),
                       np.sin(np.radians(row.sc_lat))])
        #Position in selenocentric moon-fixed mean-earth/pole coordinates
        r=rbar*(row.sc_alt+r_moon)
        mds.append((row.sc_alt+r_moon)*np.pi*2.0/360000.0)
        rs.append(r)
        #East vector
        e=np.cross(np.array([0,0,1]),rbar)
        ebar=e/np.linalg.norm(e)
        #North vector
        nbar=np.cross(rbar,ebar)
        #Velocity in selenocentric moon-fixed mean-earth/pole coordinates
        vbarr=np.sin(np.radians(row.pth))
        vbare=np.cos(np.radians(row.pth))*np.sin(np.radians(row.az))
        vbarn=np.cos(np.radians(row.pth))*np.cos(np.radians(row.az))
        vbar=vbarr*rbar+vbare*ebar+vbarn*nbar
        v=vbar*row.v
        vs.append(v)
        #p2 position
        p2=r_moon*np.array([np.cos(np.radians(row.p2_lat))*np.cos(np.radians(row.p2_lon)),
                            np.cos(np.radians(row.p2_lat))*np.sin(np.radians(row.p2_lon)),
                            np.sin(np.radians(row.p2_lat))])
        p2_srange_calc=np.sqrt((r[0]-p2[0])**2+(r[1]-p2[1])**2+(r[2]-p2[2])**2)
        dsrange2=row.p2_srange-p2_srange_calc
        dsrange2s.append(dsrange2)
        #p1 position from table
        p1a=r_moon*np.array([np.cos(np.radians(row.p1_lat))*np.cos(np.radians(row.p1_lon)),
                             np.cos(np.radians(row.p1_lat))*np.sin(np.radians(row.p1_lon)),
                             np.sin(np.radians(row.p1_lat))])
        p1a_srange_calc=np.sqrt((r[0]-p1a[0])**2+(r[1]-p1a[1])**2+(r[2]-p1a[2])**2)
        dsrange1a=row.p1_srange-p1a_srange_calc
        dsrange1as.append(dsrange1a)
        #p1 position from velocity vector
        (t,p1c)=ray_sphere_intersect(r, vbar, r_moon)
        #Since vbar is a unit vector, and r is measured in units of km, t has units of km itself,
        #and is therefore directly comparable to p1_srange.
        dsrange1b=row.p1_srange-t
        dsrange1bs.append(dsrange1b)
        dsrange1c=np.sqrt(np.sum((p1c-p1a)**2))
        dsrange1cs.append(dsrange1c)
        #Calculate velocity from last row
        ts.append(gmt_to_et(row.GMT))
        if r_last is not None:
            dt=ts[-1]-ts[-2]
            dr=np.sqrt(np.sum((r-r_last)**2))
            dv=dr/dt-row.v
            dvs.append(dv)
        r_last=r
        #print(row.GMT, gmt, cspice.etcal(gmt), mjd,tai_utc,tai,et, cspice.etcal(et))
    
    if plot_check_values:
        #This plot is meant to duplicate the residual plot on the spreadsheet
        fig,ax1=plt.subplots()
        ax2=ax1.twinx()
        ax1.plot(np.array(ts)-ts[-1],dsrange2s,'bo')
        ax1.plot(np.array(ts)-ts[-1],dsrange1as,'ro')
        ax1.plot(np.array(ts)-ts[-1],dsrange1bs,'yo')
        ax1.plot(np.array(ts)-ts[-1],dsrange1cs,'go')
        ax2.plot(np.array(ts)-ts[-1],dvs,'m+')
        ax1.plot(np.array(ts)-ts[-1],mds,'k--')
        ax1.plot(np.array(ts)-ts[-1],-np.array(mds),'k--')
        plt.show()
    return (rs,vs,ts)

def convertImageACanonical(rs,vs,ts):
    rcus=np.zeros((len(ts),3))
    vcus=np.zeros((len(ts),3))
    tcus=np.zeros( len(ts)   )
    for (i,(r,v,t)) in enumerate(zip(rs,vs,ts)):
        M=cspice.sxform("IAU_MOON","ECI_TOD",t)
        s=np.concatenate((r,v))
        Ms=np.matmul(M,s)
        rcus[i,:]=bmw.su_to_cu(Ms[0:3],r_moon,mu_moon,1, 0)
        vcus[i,:]=bmw.su_to_cu(Ms[3:6],r_moon,mu_moon,1,-1)
        tcus[i  ]=bmw.su_to_cu(t-ts[0],r_moon,mu_moon,0, 1)
    return (rcus,vcus,tcus)

def wrap_kepler(r0,v0,ts):
    """
    Use bmw.kepler to evaluate the trajectory at many times.

    :param numpy vector r0: Start position in canonical units
    :param numpy vector v0: Start velocity in canonical units
    :param numpy array  ts: Times to propagate to. One RK4 time step (4 function evaluations) between
                            each consecutive pair of times. First time should be zero.
    :rtype tuple:
    :return: First element is numpy array of position vectors, second element is numpy array of velocities 
    """
    result=np.zeros((ts.size,6))
    yi=np.concatenate((r0,v0))
    result[0,:]=yi
    for i in range(ts.size-1):
        (r,v)=bmw.kepler(r0,v0,ts[i+1])
        result[i+1,0:3]=r
        result[i+1,3:6]=v
    return result[:,0:3],result[:,3:6]

def threeBodyRK4(r0,v0,ts):
    """
    Three-body propagation at a list of discrete times. Intended to match the interface of bmw.kepler, but
    uses direct numerical integration, and takes into account the gravity of the Earth.
    
    :param numpy vector r0: Start position in canonical units
    :param numpy vector v0: Start velocity in canonical units
    :param numpy array  ts: Times to propagate to. One RK4 time step (4 function evaluations) between
                            each consecutive pair of times. First time should be zero.
    :rtype tuple:
    :return: First element is numpy array of position vectors, second element is numpy array of velocities 

    """
    def f(t,y):
        """
        Derivative function, following the form in the Wikipedia article on Runge-Kutta
        
        :param float t: Time to evaluate at, in canonical time units
        :param numpy vector y: Six element state vector, containing position and velocity at current time in canonical units
        :rtype numpy array:
        :return: Derivative of state vector with respect to time in canonical units

        Global variables used:
        * tas - first element used as Spice time of start, for calculating ephemeris of Earth
        * r_moon, mu_moon - used for converting ephemeris of Earth to canonical units
        * mu_earth - used for calculating accelerations towards Earth
        """
        dvdtMoon=-y[0:3]/np.linalg.norm(y[0:3])**3
        et=bmw.su_to_cu(t,r_moon,mu_moon,0,1,inverse=True)+tas[0]
        (xEarth,ltime)=cspice.spkezr('399',+et,'ECI_TOD','NONE','301')
        rEarth=bmw.su_to_cu(xEarth[0:3],r_moon,mu_moon,1,0)
        drEarth=y[0:3]-rEarth
        #Acceleration of the *probe* from the gravity of the Earth
        dvdtEarth=-(mu_earth/mu_moon)*drEarth/np.linalg.norm(drEarth)**3
        #acceleration of the *moon*  from the gravity of the Earth. This isn't quite right,
        #as it doesn't take into account the non-negligible mass of the moon.
        dvdtEM   = (mu_earth/mu_moon)* rEarth/np.linalg.norm( rEarth)**3
        return np.concatenate((y[3:6],dvdtMoon+dvdtEarth-dvdtEM))
    result=np.zeros((ts.size,6))
    yi=np.concatenate((r0,v0))
    result[0,:]=yi
    for i in range(ts.size-1):
        ti0=ts[i]
        ti3=ts[i+1]
        h=ti3-ti0
        ti1=ti0+h/2
        ti2=ti1
        ki1=f(ti0,yi)
        ki2=f(ti1,yi+h/2*ki1)
        ki3=f(ti2,yi+h/2*ki2)
        ki4=f(ti3,yi+h*ki3)
        yi+=(ki1+2*ki2+2*ki3+ki4)*h/6
        result[i+1]=yi
    return (result[:,0:3],result[:,3:6])

#Calculate a trajectory from the current start position to the final impact point, and return the difference in position at table rows

def cost(r0,rs,ts,dr=None,propagate=wrap_kepler):
    """
    Calculate the cost function using biased Gauss targeting and three-body propagation. Compatibile with scipy.optimize.minimize
    :param numpy vector r0: Initial position for Gauss targeting, in canonical units.
    :param numpy array of vectors rs: Positions to fit, in canonical units. Last vector is used as position for Gauss targeting to aim at.
    :param numpy array ts: Times for each position in the fit, in canonical units.
    :param float      et0: Spice time of first position
    :param numpy array dr: Target bias vector, in canonical units. This vector is subtracted from the last position and used as the target
                           for Gauss targeting. A correct bias will result in nearly hitting the actual last position.
    :rtype float:
    :return: Cost function, square of distance from each given position to the calculated position on the targeted trajectory at the 
             corresponding time. Value is in square canonical distance units.
    """
    #unpack *args rs,ts,dr=None,propagate=wrap_kepler
    result=0
    target=rs[-1]
    if dr is not None:
        target-=dr
    (v0,v1)=bmw.gauss(r0,target,ts[-1]-ts[0])
    (rcalcs,vcalcs)=propagate(r0,v0,ts)
    for i in range(ts.size):
        result+=np.sum((rs[i]-rcalcs[i,:])**2)
    return result

def trajectory_to_su(rs,vs):
    rsus=bmw.su_to_cu(rs,r_moon,mu_moon,1, 0,inverse=True)
    vsus=bmw.su_to_cu(vs,r_moon,mu_moon,1,-1,inverse=True)
    return (rsus,vsus)

def plot_residuals(rcalcs,vcalcs,rs,vs,ts,subplot=411, title=''):
    #Use the fit trajectory, use
    #Kepler propagation to evaluate at each image time, and graph the difference
    drxfs=[]
    dryfs=[]
    drzfs=[]
    dvxfs=[]
    dvyfs=[]
    dvzfs=[]
    for (i,(r,v,t)) in enumerate(zip(rs,vs,ts)):
        drxfs.append(r[0]-rcalcs[i,0])
        dryfs.append(r[1]-rcalcs[i,1])
        drzfs.append(r[2]-rcalcs[i,2])
        dvxfs.append(v[0]-vcalcs[i,0])
        dvyfs.append(v[1]-vcalcs[i,1])
        dvzfs.append(v[2]-vcalcs[i,2])        
    plt.subplot(subplot)
    plt.title(title+', pos residuals')
    plt.ylabel('pos residual/(m)')
    plt.xlabel('Time from impact/s')
    tsus=bmw.su_to_cu(ts-ts[-1],r_moon,mu_moon,0,1,inverse=True)
    xh,=plt.plot(tsus,bmw.su_to_cu(np.array(drxfs),r_moon,mu_moon,1,0,inverse=True)*1000,'rx',label='dx')
    yh,=plt.plot(tsus,bmw.su_to_cu(np.array(dryfs),r_moon,mu_moon,1,0,inverse=True)*1000,'gx',label='dy')
    zh,=plt.plot(tsus,bmw.su_to_cu(np.array(drzfs),r_moon,mu_moon,1,0,inverse=True)*1000,'bx',label='dz')
    plt.legend(handles=(xh,yh,zh))
    plt.subplot(subplot+1)
    plt.title(title+', vel residuals')
    plt.ylabel('vel residual/(m/s)')
    plt.xlabel('Time from impact/s')
    xh,=plt.plot(tsus,bmw.su_to_cu(np.array(dvxfs),r_moon,mu_moon,1,-1,inverse=True),'r+',label='dvx')
    yh,=plt.plot(tsus,bmw.su_to_cu(np.array(dvyfs),r_moon,mu_moon,1,-1,inverse=True),'g+',label='dvy')
    zh,=plt.plot(tsus,bmw.su_to_cu(np.array(dvzfs),r_moon,mu_moon,1,-1,inverse=True),'b+',label='dvz')
    plt.legend(handles=(xh,yh,zh))

image_a=readImageA() #Read table A
(recias,vecias,tas)=processImageA(image_a)
#Convert the coordinates from moon body-fixed to moon-centered inertial canonical
(racus,vacus,tacus)=convertImageACanonical(recias,vecias,tas)

#Use Gauss targeting to get a trajectory from the initial to final positions, without any target bias
(v0_gauss,v1_gauss)=bmw.gauss(racus[0],racus[-1],tacus[-1],Type=1)

#Use three-body propagation to calculate the target bias
(rs_for_dr,vs_for_dr)=threeBodyRK4(racus[0],v0_gauss,tacus)
#How close do we get to table A?
plt.figure(1)
plot_residuals(rs_for_dr,vs_for_dr,racus,vacus,tacus,subplot=211,title='Unbiased Gauss targeting')
plt.show()
dr=rs_for_dr[-1,:]-racus[-1]

#Fit the observations using biased Gauss targeting and three-body propagation
plt.figure(2)
(rcu0_fit,vcu0_fit)=opt.minimize(cost,racus[0,:],args=(racus,tacus,dr,threeBodyRK4))
(rcus_fit,vcus_fit)=threeBodyRK4(rcu0_fit,vcu0_fit,tacus)
plot_residuals(rcus_fit,vcus_fit,racus,vacus,tacus,subplot=211)

#Now bias the aimpoint for Gauss and re-find the position and velocity
dr1=r1_rk-rcus[-1]
print("State before tide correction (km): ",bmw.su_to_cu(r0_fit,r_moon,mu_moon,1,0,inverse=True),bmw.su_to_cu(v0_fit,r_moon,mu_moon,1,-1,inverse=True))
print("Error before tide correction (km): ",bmw.su_to_cu(dr1,r_moon,mu_moon,1,0,inverse=True))
(v0_fit2,v1_fit2)=bmw.gauss(r0_fit,rcus[-1]-dr1,tcus[-1])

done_fit=False
r0_fit2=r0_fit
costm=0.0
costs2=[]
while not done_fit:
    cost0=fit_cost2(r0_fit2,rcus,dr1,tcus)
    costs2.append(cost0)
    rstep=1e-7
    dcostdrx=(fit_cost2(r0_fit2+np.array([rstep,0.000,0.000]),rcus,dr1,tcus)-cost0)/rstep        
    dcostdry=(fit_cost2(r0_fit2+np.array([  0.0,rstep,0.000]),rcus,dr1,tcus)-cost0)/rstep        
    dcostdrz=(fit_cost2(r0_fit2+np.array([  0.0,0.000,rstep]),rcus,dr1,tcus)-cost0)/rstep        
    print("Old cost: %e, new cost: %e" %(costm,cost0))
    if costm==0.0:
        done_fit=False
    else:
        done_fit=abs(cost0-costm)<1e-4*cost0
    if done_fit:
        break
    costm=cost0
    dr=np.array([dcostdrx,dcostdry,dcostdrz])
    gamma=rstep/np.linalg.norm(dr) #Move 1m closer to correct position
    dr*=gamma
    r0_fit2-=dr

(v0_fit2,v1_fit2)=bmw.gauss(r0_fit2,rcus[-1]-dr1,tcus[-1])
(r1_rk2,v1_rk2)=threeBodyRK4(r0_fit2,v0_fit2,tcus)

#Now bias the aimpoint for Gauss and re-find the position and velocity
dr2=r1_rk2[-1,:]-rcus[-1]
print("State after tide correction (km): ",bmw.su_to_cu(r0_fit2,r_moon,mu_moon,1,0,inverse=True),bmw.su_to_cu(v0_fit2,r_moon,mu_moon,1,-1,inverse=True))
print("Error after tide correction (km): ",bmw.su_to_cu(dr2,r_moon,mu_moon,1,0,inverse=True))

if False:    
    plt.plot(costs)
    plt.plot(costs2)
    plt.show()

#Use the fit trajectory, use
#three-body propagation to evaluate at each image time, and graph the difference
drxfs2=[]
dryfs2=[]
drzfs2=[]
dvxfs2=[]
dvyfs2=[]
dvzfs2=[]
recis=np.zeros([len(ts),3])
vecis=np.zeros([len(ts),3])
print("elorb0: ",bmw.elorb(r0_fit,v0_fit))
for (i,(rcu,vcu,tcu)) in enumerate(zip(rcus,vcus,tcus)):
    (rcalc,vcalc)=(r1_rk2[i,:],v1_rk2[i,:])
    recis[i,:]=bmw.su_to_cu(rcalc,r_moon,mu_moon, 1, 0,inverse=True)
    vecis[i,:]=bmw.su_to_cu(vcalc,r_moon,mu_moon, 1,-1,inverse=True)
    drxfs2.append(rcu[0]-rcalc[0])
    dryfs2.append(rcu[1]-rcalc[1])
    drzfs2.append(rcu[2]-rcalc[2])
    dvxfs2.append(vcu[0]-vcalc[0])
    dvyfs2.append(vcu[1]-vcalc[1])
    dvzfs2.append(vcu[2]-vcalc[2])
    
if True:
    plt.subplot(413)
    plt.plot(np.array(ts)-ts[-1],bmw.su_to_cu(np.array(drxfs2),r_moon,mu_moon,1,0,inverse=True),'rx')
    plt.plot(np.array(ts)-ts[-1],bmw.su_to_cu(np.array(dryfs2),r_moon,mu_moon,1,0,inverse=True),'gx')
    plt.plot(np.array(ts)-ts[-1],bmw.su_to_cu(np.array(drzfs2),r_moon,mu_moon,1,0,inverse=True),'bx')
    plt.subplot(414)
    plt.plot(np.array(ts)-ts[-1],bmw.su_to_cu(np.array(dvxfs2),r_moon,mu_moon,1,-1,inverse=True),'r+')
    plt.plot(np.array(ts)-ts[-1],bmw.su_to_cu(np.array(dvyfs2),r_moon,mu_moon,1,-1,inverse=True),'g+')
    plt.plot(np.array(ts)-ts[-1],bmw.su_to_cu(np.array(dvzfs2),r_moon,mu_moon,1,-1,inverse=True),'b+')
    plt.show()

trajtuple=namedtuple('trajtuple',['GMT',
                                  'GeoRX','GeoRY','GeoRZ',
                                  'GeoVX','GeoVY','GeoVZ',
                                  'SelenoRX','SelenoRY','SelenoRZ',
                                  'SelenoVX','SelenoVY','SelenoVZ'])
traj=[]
with open('Ranger 7 Trajectory - Trajectory.csv','r') as inf:
    reader=csv.reader(inf)
    #Read past the first header line
    for row in reader:
        break
    #Read the rows
    for row in reader:
        traj.append(trajtuple(      row[ 0] , #GMT
                             floatN(row[ 1]), #GeoRX
                             floatN(row[ 2]), #GeoRY
                             floatN(row[ 3]), #GeoRZ
                             floatN(row[ 4]), #GeoVX
                             floatN(row[ 5]), #GeoVY
                             floatN(row[ 6]), #GeoVZ
                             floatN(row[ 7]), #SelenoRX
                             floatN(row[ 8]), #SelenoRY
                             floatN(row[ 9]), #SelenoRZ
                             floatN(row[10]), #SelenoVX
                             floatN(row[11]), #SelenoVY
                             floatN(row[12])))#SelenoVZ
        print(traj[-1])
        
t=np.zeros(len(traj))
GeoState=np.zeros((len(traj),6))
SelenoState=np.zeros((len(traj),6))
for i,row in enumerate(traj):
    t[i]=gmt_to_et(row.GMT)
    GeoState[i,:]=np.array(row[1:7])
    SelenoState[i,:]=np.array(row[7:13])

Ranger7Geo_txt="""
Ranger 7 - first completely successful Ranger lunar impact mission. Data from
'Ranger 7 flight path and its determination from tracking data', 15 Dec 1964
available at https://archive.org/details/nasa_techdoc_19650003678 (but I got
it from NTRS)

The report had a table of geocentric state vectors in the True of Date system
starting at injection and continuing every hour on the hour until impact. The
report also had selenocentric state vectors starting at 6:00 GMT 31 Jul 1964,
including one at time of impact, 1964 Jul 31 13:25:48.724 GMT

These vectors were manually entered into a spreadsheet and verified by checking
the magnitude of the radius and velocity vectors, which always matched to within
1 unit in the last place (several errors were caught this way).

\\begindata
INPUT_DATA_TYPE = 'STATES'
OUTPUT_SPK_TYPE = 5
OBJECT_ID=-1007
OBJECT_NAME='RANGER 7'
CENTER_ID=399
CENTER_NAME='EARTH'
REF_FRAME_NAME='ECI_TOD'
PRODUCER_ID='C. Jeppesen, Kwan Systems'
DATA_ORDER='EPOCH X Y Z VX VY VZ'
TIME_WRAPPER='# ETSECONDS'
INPUT_DATA_UNITS = ('ANGLES=DEGREES' 'DISTANCES=km')
DATA_DELIMITER=';'
LINES_PER_RECORD=1
CENTER_GM=%11.4f
FRAME_DEF_FILE='%s/fk/eci_tod.tf'
LEAPSECONDS_FILE='%s/lsk/naif0011.tls'
""" % (mu_earth,'../../Data/spice/generic','../../Data/spice/generic')

Ranger7Seleno_txt="""
Ranger 7 - first completely successful Ranger lunar impact mission. Data from
'Ranger 7 flight path and its determination from tracking data', 15 Dec 1964
available at https://archive.org/details/nasa_techdoc_19650003678 (but I got
it from NTRS)

The report had a table of geocentric state vectors in the True of Date system
starting at injection and continuing every hour on the hour until impact. The
report also had selenocentric state vectors starting at 6:00 GMT 31 Jul 1964,
including one at time of impact, 1964 Jul 31 13:25:48.724 GMT

These vectors were manually entered into a spreadsheet and verified by checking
the magnitude of the radius and velocity vectors, which always matched to within
1 unit in the last place (several errors were caught this way).

This file contains the selenocentric states converted to elements. I have not checked
the accuracy of the lunar ephemeris (which is in the Ranger report, but would require
typing the same number of numbers again) so any error in that ephemeris is propagated
here as a difference between this kernel and the geocentric kernel covering the same
time. Again, the elements should be completely accurate at the original report points.  

\\begindata
INPUT_DATA_TYPE = 'STATES'
OUTPUT_SPK_TYPE = 5
OBJECT_ID=-1007
OBJECT_NAME='RANGER 7'
CENTER_ID=301
CENTER_NAME='MOON'
REF_FRAME_NAME='ECI_TOD'
PRODUCER_ID='C. Jeppesen, Kwan Systems'
DATA_ORDER='EPOCH X Y Z VX VY VZ'
TIME_WRAPPER='# ETSECONDS'
INPUT_DATA_UNITS = ('ANGLES=DEGREES' 'DISTANCES=km')
DATA_DELIMITER=';'
LINES_PER_RECORD=1
CENTER_GM=%9.4f
FRAME_DEF_FILE='%s/fk/eci_tod.tf'
LEAPSECONDS_FILE='%s/lsk/naif0011.tls'
""" % (mu_moon,'../../Data/spice/generic','../../Data/spice/generic')

Ranger7Seleno2_txt="""
Ranger 7 - first completely successful Ranger lunar impact mission. Data from
'Ranger VII Photographic Parameters', JPL Technical Report No. 32-964, 1 Nov 1966
available at NTRS as document number 19670002488 

The report had a table of camera parameters including the position of the spacecraft
in a selenocentric body-fixed (Mean-earth-polar) frame for each picture published in
the photo atlases. This spice segment uses the positions from the A camera, starting 
about 15min before impact at %s, and ending with the last A
camera image 2.5s before impact at %s. There is also a 
calculated state at impact, at %s. 

The table seems to have had a bias in longitude, as it matches neither the previous
segments nor the actual location of the crater as found by LRO/LROC. This bias is 
corrected in this spice segment, such that the final longitude matches that of LRO, and
matches the selenocentric trajectory from the other segments much better.

These table values were manually entered into a spreadsheet and verified by checking
the slant ranges, which were consistent with the table values to the precision allowed
by the table latitude and longitude (several errors were caught this way).

\\begindata
INPUT_DATA_TYPE = 'STATES'
OUTPUT_SPK_TYPE = 5
OBJECT_ID=-1007
OBJECT_NAME='RANGER 7'
CENTER_ID=301
CENTER_NAME='MOON'
REF_FRAME_NAME='ECI_TOD'
PRODUCER_ID='C. Jeppesen, Kwan Systems'
DATA_ORDER='EPOCH X Y Z VX VY VZ'
TIME_WRAPPER='# ETSECONDS'
INPUT_DATA_UNITS = ('ANGLES=DEGREES' 'DISTANCES=km')
DATA_DELIMITER=';'
LINES_PER_RECORD=1
CENTER_GM=%9.4f
FRAME_DEF_FILE='%s/fk/eci_tod.tf'
LEAPSECONDS_FILE='%s/lsk/naif0011.tls'
""" % (rows[0].GMT,rows[-2].GMT,rows[-1].GMT,mu_moon,'../../Data/spice/generic','../../Data/spice/generic')

print(t)
tofs=None
with open('geo.txt','w') as ouf_geo:
    with open('seleno.txt','w') as ouf_seleno:
        for i in range(len(t)):
            this_state=np.array(GeoState[i,:])
            this_pos=this_state[0:3]
            this_vel=this_state[3:6]
            print("%23.6f;%23.15e;%23.15e;%23.15e;%23.15e;%23.15e;%23.15e"%(t[i],
                this_pos[0],
                this_pos[1],
                this_pos[2],
                this_vel[0],
                this_vel[1],
                this_vel[2]),file=ouf_geo)
            this_state=np.array(SelenoState[i,:])
            if np.isfinite(this_state[0]):
                if tofs is None:
                    tofs=i
                print("%23.6f;%23.15e;%23.15e;%23.15e;%23.15e;%23.15e;%23.15e"%(t[i],
                this_state[0],
                this_state[1],
                this_state[2],
                this_state[3],
                this_state[4],
                this_state[5]),file=ouf_seleno)

with open('seleno2.txt','w') as ouf:
    for i in range(len(ts)):
        print("%23.6f;%23.15e;%23.15e;%23.15e;%23.15e;%23.15e;%23.15e"%(ts[i],
            recis[i,0],
            recis[i,1],
            recis[i,2],
            vecis[i,0],
            vecis[i,1],
            vecis[i,2]),file=ouf)

with open('Ranger7Geo_mkspk.txt','w') as ouf:
    print(Ranger7Geo_txt,file=ouf)

with open('Ranger7Seleno_mkspk.txt','w') as ouf:
    print(Ranger7Seleno_txt,file=ouf)

with open('Ranger7Seleno2_mkspk.txt','w') as ouf:
    print(Ranger7Seleno2_txt,file=ouf)

import subprocess
import os
try:
    os.remove("Ranger7.bsp")
except FileNotFoundError:
    pass #no error, file is already not present
subprocess.call("~/bin/mkspk -setup Ranger7Geo_mkspk.txt     -input geo.txt     -output Ranger7.bsp        ",shell=True)
subprocess.call("~/bin/mkspk -setup Ranger7Seleno_mkspk.txt  -input seleno.txt  -output Ranger7.bsp -append",shell=True)
subprocess.call("~/bin/mkspk -setup Ranger7Seleno2_mkspk.txt -input seleno2.txt -output Ranger7.bsp -append",shell=True)

cspice.furnsh("Ranger7.bsp")

selenostatepos_x=[0.0]*SelenoState.shape[0]
selenostatepos_y=[0.0]*SelenoState.shape[0]

for i in range(SelenoState.shape[0]):
    this_state = SelenoState[i,:]
    if np.isfinite(this_state[0]):
        this_pos=this_state[0:3]
        selenostatepos_x[i]=this_pos[0]
        selenostatepos_y[i]=this_pos[1]

n_step=1000
step=sorted(list(t[tofs:])+list(ts)+list(np.linspace(t[tofs],t[-1],n_step)))
n_step=len(step)
spicepos_x=[0.0]*n_step
spicepos_y=[0.0]*n_step

for i,tt in enumerate(step):
    (spice_state,ltime)=cspice.spkezr('-1007',tt,'ECI_TOD','NONE','301')
    spice_pos=spice_state[0:3]
    spicepos_x[i]=spice_pos[0]
    spicepos_y[i]=spice_pos[1]
    #print(spice_pos)
    spice_vel=spice_state[3:6]
    #print(spice_vel)

if True:
    plt.plot(selenostatepos_x,selenostatepos_y,'b*')
    plt.plot(spicepos_x,spicepos_y,'g-*')
    plt.plot(recis[:,0],recis[:,1],'r+')
    plt.axis('equal')
    plt.show()


