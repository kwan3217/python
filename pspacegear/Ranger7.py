"""
Created on Dec 6, 2017

@author: chrisj
"""

import csv
from collections import namedtuple
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
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
        return float('NaN')

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

def processImageA(image_a,plot=False):
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
    :param bool plot: If true, plot the check value residuals
    :rtype: tuple
    :return: First element is numpy array of position vectors, one row for each row in the table, 3 columns
             Second element is numpy array of velocity vectors, one row for each row in the table, 3 columns
             Third element
    """
    rs=np.zeros((len(image_a),3),np.float64)
    vs=np.zeros((len(image_a),3),np.float64)
    #dsrange2 is the difference between the table slant range to point 2 and that calculated from the spacecraft lat/lon/alt
    #and point 2 lat/lon
    dsrange2s=np.zeros(len(image_a),np.float64)
    #dsrange1a is the difference between the table slant range to point 1 and that calculated from the spacecraft lat/lon/alt
    dsrange1as=np.zeros(len(image_a),np.float64)
    #dsrange1b is the difference between the table slant range to point 1 and that calculated from the spacecraft pos/vel
    #and the quadratic method (quadratic parameter t is the calculated distance between the ray origin and the ray/sphere
    #intersect point)
    dsrange1bs=np.zeros(len(image_a),np.float64)
    #dsrange1c is the distance between the table point 1 calculated from lat/lon and that calculated by the quadratic
    #method. This isn't a difference in srange like the others are, but it is measured in the same units. However, dsrange1c
    #will always be positive, while the other measures can be positive or negative.
    dsrange1cs=np.zeros(len(image_a),np.float64)
    #dv is the difference between the table velocity and that calculated by dividing the distance from the previous row's
    #position to this row's position by the difference in time. Keep track of last row position in r_last, use constant
    #5.12s as dt. Note that this will be biased from zero because it doesn't take into account the acceleration of gravity
    #over the time step.
    dvs=np.zeros(len(image_a),np.float64)
    dvs[0]=float('NaN')
    r_last=None
    #Size of 1 millidegree of latitude at the current spacecraft altitude. This is an idea of the precision we can expect
    #in using vectors with latitudes and longitudes specified in millidegree precision.
    mds=np.zeros(len(image_a),dtype=np.float64)
    ts=np.zeros(len(image_a),dtype=np.float64)
    for i,row in enumerate(image_a):
        #Zenith vector
        rbar=np.array([np.cos(np.radians(row.sc_lat))*np.cos(np.radians(row.sc_lon)),
                       np.cos(np.radians(row.sc_lat))*np.sin(np.radians(row.sc_lon)),
                       np.sin(np.radians(row.sc_lat))])
        #Position in selenocentric moon-fixed mean-earth/pole coordinates
        r=rbar*(row.sc_alt+r_moon)
        mds[i]=(row.sc_alt+r_moon)*np.pi*2.0/360000.0
        rs[i,:]=r
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
        vs[i,:]=v
        #p2 position
        p2=r_moon*np.array([np.cos(np.radians(row.p2_lat))*np.cos(np.radians(row.p2_lon)),
                            np.cos(np.radians(row.p2_lat))*np.sin(np.radians(row.p2_lon)),
                            np.sin(np.radians(row.p2_lat))])
        p2_srange_calc=np.sqrt((r[0]-p2[0])**2+(r[1]-p2[1])**2+(r[2]-p2[2])**2)
        dsrange2=row.p2_srange-p2_srange_calc
        dsrange2s[i]=dsrange2
        #p1 position from table
        p1a=r_moon*np.array([np.cos(np.radians(row.p1_lat))*np.cos(np.radians(row.p1_lon)),
                             np.cos(np.radians(row.p1_lat))*np.sin(np.radians(row.p1_lon)),
                             np.sin(np.radians(row.p1_lat))])
        p1a_srange_calc=np.sqrt((r[0]-p1a[0])**2+(r[1]-p1a[1])**2+(r[2]-p1a[2])**2)
        dsrange1a=row.p1_srange-p1a_srange_calc
        dsrange1as[i]=dsrange1a
        #p1 position from velocity vector
        (t,p1c)=ray_sphere_intersect(r, vbar, r_moon)
        #Since vbar is a unit vector, and r is measured in units of km, t has units of km itself,
        #and is therefore directly comparable to p1_srange.
        dsrange1b=row.p1_srange-t
        dsrange1bs[i]=dsrange1b
        dsrange1c=np.sqrt(np.sum((p1c-p1a)**2))
        dsrange1cs[i]=dsrange1c
        #Calculate velocity from last row
        ts[i]=gmt_to_et(row.GMT)
        if r_last is not None:
            dt=ts[-1]-ts[-2]
            if dt==0:
                pass
            dr=np.sqrt(np.sum((r-r_last)**2))
            dv=dr/dt-row.v
            dvs[i]=dv
        r_last=r
        #print(row.GMT, gmt, cspice.etcal(gmt), mjd,tai_utc,tai,et, cspice.etcal(et))
    
    if plot:
        #This plot is meant to duplicate the residual plot on the spreadsheet
        fig,ax1=plt.subplots()
        ax2=ax1.twinx()
        ax1.plot(ts-ts[-1],dsrange2s,'bo')
        ax1.plot(ts-ts[-1],dsrange1as,'ro')
        ax1.plot(ts-ts[-1],dsrange1bs,'yo')
        ax1.plot(ts-ts[-1],dsrange1cs,'go')
        ax2.plot(ts-ts[-1],dvs,'m+')
        ax1.plot(ts-ts[-1],mds,'k--')
        ax1.plot(ts-ts[-1],-np.array(mds),'k--')
        plt.show()
    return rs, vs, ts

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
    def f(t,y,i):
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
        drEarth=y[0:3]-rEarth[i,:] #Use cached Earth position
        #Acceleration of the *probe* from the gravity of the Earth
        dvdtEarth=-(mu_earth/mu_moon)*drEarth/np.linalg.norm(drEarth)**3
        return np.concatenate((y[3:6],dvdtMoon+dvdtEarth-dvdtEM[i,:])) #Use cached EM acceleration
    result=np.zeros((ts.size,6))
    yi=np.concatenate((r0,v0))
    result[0,:]=yi
    for i in range(ts.size-1):
        ti0=ts[i]
        ti3=ts[i+1]
        h=ti3-ti0
        ti1=ti0+h/2
        ti2=ti1
        ki1=f(ti0,yi        ,i*2  )
        ki2=f(ti1,yi+h/2*ki1,i*2+1)
        ki3=f(ti2,yi+h/2*ki2,i*2+1)
        ki4=f(ti3,yi+h  *ki3,i*2+2)
        yi+=(ki1+2*ki2+2*ki3+ki4)*h/6
        result[i+1,:]=yi
    return (result[:,0:3],result[:,3:6])

def cost(r0,rs,ts,bias=None,propagate=wrap_kepler):
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
    target=rs[-1,:]
    if bias is not None:
        target-=bias
    (v0,v1)=bmw.gauss(r0,target,ts[-1]-ts[0])
    (rcalcs,vcalcs)=propagate(r0,v0,ts)
    result=np.sum((rs-rcalcs)**2)
    print(r0,result)
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
    mdegs=[]
    for (i,(r,v,t)) in enumerate(zip(rs,vs,ts)):
        drxfs.append(r[0]-rcalcs[i,0])
        dryfs.append(r[1]-rcalcs[i,1])
        drzfs.append(r[2]-rcalcs[i,2])
        dvxfs.append(v[0]-vcalcs[i,0])
        dvyfs.append(v[1]-vcalcs[i,1])
        dvzfs.append(v[2]-vcalcs[i,2])
        mdegs.append(bmw.su_to_cu(np.linalg.norm(r)*2*np.pi/360000.0,r_moon,mu_moon,1,0,inverse=True))        
    plt.subplot(subplot)
    plt.title(title+', pos residuals')
    plt.ylabel('pos residual/(m)')
    plt.xlabel('Time from impact/s')
    tsus=bmw.su_to_cu(ts-ts[-1],r_moon,mu_moon,0,1,inverse=True)
    xh,=plt.plot(tsus,bmw.su_to_cu(np.array(drxfs),r_moon,mu_moon,1,0,inverse=True)*1000,'rx',label='dx')
    yh,=plt.plot(tsus,bmw.su_to_cu(np.array(dryfs),r_moon,mu_moon,1,0,inverse=True)*1000,'gx',label='dy')
    zh,=plt.plot(tsus,bmw.su_to_cu(np.array(drzfs),r_moon,mu_moon,1,0,inverse=True)*1000,'bx',label='dz')
    plt.plot(tsus,np.array(mdegs)*500,'k--')
    plt.plot(tsus,np.array(mdegs)*-500,'k--')
    plt.legend(handles=(xh,yh,zh))
    plt.subplot(subplot+1)
    plt.title(title+', vel residuals')
    plt.ylabel('vel residual/(m/s)')
    plt.xlabel('Time from impact/s')
    xh,=plt.plot(tsus,bmw.su_to_cu(np.array(dvxfs),r_moon,mu_moon,1,-1,inverse=True),'r+',label='dvx')
    yh,=plt.plot(tsus,bmw.su_to_cu(np.array(dvyfs),r_moon,mu_moon,1,-1,inverse=True),'g+',label='dvy')
    zh,=plt.plot(tsus,bmw.su_to_cu(np.array(dvzfs),r_moon,mu_moon,1,-1,inverse=True),'b+',label='dvz')
    plt.legend(handles=(xh,yh,zh))

def cache_earth(ts):
    """"
    cache the position vector of the Earth, and the EM acceleration, since
    we always use threeBodyRK4 with the same ts.

    :param numpy array ts: array of Spice times to calculate the positions at
    :rtype tuple:
    :return: First elemet is rEarth, a numpy array [ts.size*2-1,3]. Each row is the
               position of the Earth relative to the Moon in selenocentric ECI_TOD frame,
               in lunar canonical units. Every even row is the position at one of
               the requested times in ts, and every odd row is the position exactly
               in between two consecutive times in ts.
             Second element is dvdtEarth, a numpy array [ts.size*2-1,3]. Each row is the
               acceleration of the moon towards the earth in the same frame and units
               as above. Same even/odd breakdown too.
    """
    rEarth=np.zeros((ts.size*2-1,3),np.float64)
    dvdtEM=np.zeros((ts.size*2-1,3),np.float64)
    for i in range(ts.size):
        (xEarth, _) = cspice.spkezr('399', ts[i], 'ECI_TOD', 'NONE', '301')
        rEarth[i*2,:] = bmw.su_to_cu(xEarth[0:3], r_moon, mu_moon, 1, 0)
        # acceleration of the *moon*  from the gravity of the Earth. This isn't quite right,
        # as it doesn't take into account the non-negligible mass of the moon.
        dvdtEM[i*2,:] = (mu_earth/mu_moon)*rEarth[i*2,:]/np.linalg.norm(rEarth[i*2,:]) ** 3
        if i<tas.size-1:
            (xEarth, _) = cspice.spkezr('399', (ts[i]+ts[i+1])/2, 'ECI_TOD', 'NONE', '301')
            rEarth[i*2+1,:] = bmw.su_to_cu(xEarth[0:3], r_moon, mu_moon, 1, 0)
            dvdtEM[i*2+1,:] = (mu_earth/mu_moon)*rEarth[i*2,:]/np.linalg.norm(rEarth[i*2,:]) ** 3
    return rEarth, dvdtEM

def gradient_descent(F,x0,args=(),delta=1e-14,gamma0=1e-12,adapt=False,plot=False):
    """
    Find a local minimum of a scalar-valued function of a vector (scalar field)
    by the Gradient Descent method.

    https://en.wikipedia.org/wiki/Gradient_descent

    :param function F: Function to minimize. Initial use case is a cost
                       function computed as a sum of squares of difference
                       between points on a calculated trajectory and observed
                       points on the same trajectory
    :param numpy vector x0: Initial guess as to minimum of function. This method
                            will "probably" find the local minimum nearest
                            to this point
    :param float delta: Amount to perturb the function parameters to calculate
                        the gradient
    :param float gamma0: Initial step size. The step size is adaptive after
                         the first step, so this controls the size of that
                         first step.

    Gradient descent minimizes the function by steps. At each step, it figures
    out the gradient of the function. This is a vector which points in the
    direction of maximum increase of the function. Since we want to minimize,
    we take a step in the opposite direction. Gamma is a scale factor which
    determines how far to step in that direction.

    \delta_x=\gamma_n*grad(F(x_n)
    x_{n+1}=x_{n+1}-\delta_x

    Now the only hard part is to figure out the right gamma. A simple
    non-working optimizer would use the gradient function to project to zero,
    and take a step of that size. Unfortunately, that is only sound if the
    minimum is exactly zero. The wiki (boo!) gives a formula for gamma as:

    \gamma_n=\frac{(x_n-x_{n-1})^T[grad(F(x_n))-grad(F(x_{n-1}))]}
                  {norm(grad(F(x_n))-grad(F(x_{n-1})))**2}

    Since this requires the gradient both at the current step and the previous one,
    we need an initial step size in order to take that first step.
    """
    def grad(F,x0,delta,*args):
        """
        Estimate the gradient of a scalar field by finite differences
        :param function F: Function to calculate the value of the field at a given point
        :param numpy vector x0: Point to take the gradient at
        :param float delta: Step size to take
        :return: Estimate of gradient vector
        """
        F0=F(x0,*args)
        kd=np.zeros(x0.size)
        result=np.zeros(x0.size)
        min=True
        for i in range(x0.size):
            kd[i]=delta
            Fp=F(x0+kd,*args)
            Fm=F(x0-kd,*args)
            if Fp<F0 or Fm<F0:
                min=False
            result[i]=(Fp-Fm)/(2*delta)
            kd[i]=0.0
        return (result,F0,min)
    #Take the first step
    xnm1=x0
    gFxnm1,Fxnm1,min=grad(F,x0,delta,*args)
    if min:
        #won't be able to find an improvement, delta is too big
        return x0
    xn=xnm1-gamma0*gFxnm1
    #Now take steps with dynamic step size
    done=False
    Fs=[Fxnm1]
    xs=[0.0]
    ys=[0.0]
    zs=[0.0]
    while not done:
        gFxn,Fxn,min=grad(F,xn,delta,*args)
        Fs.append(Fxn)
        xs.append(xn[0]-x0[0])
        ys.append(xn[1]-x0[1])
        zs.append(xn[2]-x0[2])
        if adapt:
            #Adaptive size (doesn't work for my problem)
            dx=xn-xnm1
            dF=gFxn-gFxnm1
            gamma=np.dot(dx,dF)/np.dot(dF,dF)
        else:
            gamma=gamma0
        xnm1=xn
        gFxnm1=gFxn
        xn=xnm1-gamma*gFxnm1/np.linalg.norm(gFxnm1)
        done=min
    if plot:
        mpl.rcParams['legend.fontsize'] = 10
        fig = plt.figure(3)
        ax = fig.gca(projection='3d')
        ax.plot(ys, zs, Fs)
        plt.show()
    return xn

image_a=readImageA() #Read table A
(recias,vecias,tas)=processImageA(image_a)

#Convert the coordinates from moon body-fixed to moon-centered inertial canonical
(racus,vacus,tacus)=convertImageACanonical(recias,vecias,tas)

#Cache the earth positions and accelerations
(rEarth,dvdtEM)=cache_earth(tas)

#Use Gauss targeting to get a trajectory from the initial to final positions,
#without any target bias
(v0_gauss,v1_gauss)=bmw.gauss(racus[0],racus[-1],tacus[-1],Type=1)

#Use three-body propagation to calculate the target bias
(rs_for_bias,vs_for_bias)=threeBodyRK4(racus[0],v0_gauss,tacus)
plt.figure(1)
plot_residuals(rs_for_bias,vs_for_bias,racus,vacus,tacus,subplot=211,title='Unbiased Gauss targeting')
bias=rs_for_bias[-1,:]-racus[-1]

#Fit the observations using biased Gauss targeting and three-body propagation
plt.figure(2)
(v0_gaussb,_)=bmw.gauss(racus[0],racus[-1,:]-bias,tacus[-1])
(rcus_gaussb,vcus_gaussb)=threeBodyRK4(racus[0],v0_gaussb,tacus)
plot_residuals(rcus_gaussb,vcus_gaussb,racus,vacus,tacus,subplot=211,title='Biased Gauss targeting')

rcus_for_spice=rcus_gaussb
vcus_for_spice=vcus_gaussb

#Previous results work from the first table position, which has limited precision. The code below is
#an attempt to find a start point which fits all of the points better, but unfortunately it is failing.
#It looks like the initial guess is so close to the minimum that numerical precision is getting in the
#way of finding a better initial guess.
if False:
    rcu0_fit=gradient_descent(cost,racus[0,:],args=(racus,tacus,bias,threeBodyRK4),plot=False)
    (v0_fit,_)=bmw.gauss(rcu0_fit,racus[-1,:]-bias,tacus[-1])
    (rcus_fit,vcus_fit)=threeBodyRK4(racus[0],v0_gaussb,tacus)
    plt.figure(3)
    plot_residuals(rcus_fit,vcus_fit,racus,vacus,tacus,subplot=211,title='Fit Gauss targeting')
    rcus_for_spice=rcus_fit
    vcus_for_spice=vcus_fit
    
(rs_for_spice,vs_for_spice)=trajectory_to_su(rcus_for_spice, vcus_for_spice)

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
""" % (image_a[0].GMT,image_a[-2].GMT,image_a[-1].GMT,mu_moon,'../../Data/spice/generic','../../Data/spice/generic')

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

with open('seleno2.txt','w') as ouf_seleno2:
    for i in range(tas.size):
        print("%23.6f;%23.15e;%23.15e;%23.15e;%23.15e;%23.15e;%23.15e"%(tas[i],
            rs_for_spice[i,0],
            rs_for_spice[i,1],
            rs_for_spice[i,2],
            vs_for_spice[i,0],
            vs_for_spice[i,1],
            vs_for_spice[i,2]),file=ouf_seleno2)

with open('Ranger7Geo_mkspk.txt','w') as ouf:
    print(Ranger7Geo_txt,file=ouf)

with open('Ranger7Seleno_mkspk.txt','w') as ouf:
    print(Ranger7Seleno_txt,file=ouf)

with open('Ranger7Seleno2_mkspk.txt','w') as ouf:
    print(Ranger7Seleno2_txt,file=ouf)

import subprocess
try:
    os.remove("Ranger7.bsp")
except FileNotFoundError:
    pass #no error, file is already not present
subprocess.call("~/bin/mkspk -setup Ranger7Geo_mkspk.txt     -input geo.txt     -output Ranger7.bsp        ",shell=True)
subprocess.call("~/bin/mkspk -setup Ranger7Seleno_mkspk.txt  -input seleno.txt  -output Ranger7.bsp -append",shell=True)
subprocess.call("~/bin/mkspk -setup Ranger7Seleno2_mkspk.txt -input seleno2.txt -output Ranger7.bsp -append",shell=True)

cspice.furnsh("Ranger7.bsp")

selenostatepos_x=[]
selenostatepos_y=[]
selenostatepos_z=[]

for i in range(SelenoState.shape[0]):
    this_state = SelenoState[i,:]
    if np.isfinite(this_state[0]):
        this_pos=this_state[0:3]
        selenostatepos_x.append(this_pos[0])
        selenostatepos_y.append(this_pos[1])
        selenostatepos_z.append(this_pos[2])

n_step=1000
step=sorted(list(t[tofs:])+list(tas)+list(np.linspace(t[tofs],t[-1],n_step)))
n_step=len(step)
spicepos_x=np.zeros(n_step)
spicepos_y=np.zeros(n_step)
spicepos_z=np.zeros(n_step)

for i,tt in enumerate(step):
    (spice_state,ltime)=cspice.spkezr('-1007',tt,'ECI_TOD','NONE','301')
    spice_pos=spice_state[0:3]
    spicepos_x[i]=spice_pos[0]
    spicepos_y[i]=spice_pos[1]
    spicepos_z[i]=spice_pos[2]
    #print(spice_pos)
    spice_vel=spice_state[3:6]
    print(cspice.etcal(tt),tt,spice_state)

if True:
    plt.figure(4)
    plt.subplot(211)
    plt.xlabel('x selenocentric/km')
    plt.ylabel('y selenocentric/km')
    plt.plot(selenostatepos_x,selenostatepos_y,'b*',label='selenostate')
    plt.plot(spicepos_x,spicepos_y,'g-*',label='spice output')
    plt.plot(rs_for_spice[:,0],rs_for_spice[:,1],'r+',label='spice input')
    plt.legend()
    plt.axis('equal')
    plt.subplot(212)
    plt.xlabel('x selenocentric/km')
    plt.ylabel('z selenocentric/km')
    plt.plot(selenostatepos_x,selenostatepos_z,'b*',label='selenostate')
    plt.plot(spicepos_x,spicepos_z,'g-*',label='spice output')
    plt.plot(rs_for_spice[:,0],rs_for_spice[:,2],'r+',label='spice input')
    plt.legend()
    plt.axis('equal')
    plt.show()


