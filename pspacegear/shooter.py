'''
Created on Jan 17, 2017

@author: chrisj
'''

import numpy as np
from numpy.linalg import norm as vlength
from numpy import cross, dot, cos, sin, arccos as acos, sqrt, pi
import matplotlib.pyplot as plt
import collections

def vangle(a,b):
    return acos(dot(a,b)/(vlength(a)*vlength(b)))

#Earth
#mu=398600.4415e9  #Gravitational parameter of Earth
#w=2*pi/86164.09 #Rotation speed of Earth
#Re=6378137         #Equatorial radius of Earth

#Kerbin
mu=3.5316000e12       #Gravitational parameter, m^3/s^2
#w=2*pi/(21549.425) #Rotation speed, rad/s (calculated from Sidereal period)
w=0
Re=600000             #Equatorial radius, m
pole=w*np.array([0,0,1]) #Direction is direction of North Pole, magnitude is spin
                         #rate in radians/sec. It happens that the wind speed
                         #at any point is the cross product of this vector with
                         #the position vector of the point 
g0=9.82               #m/s^2, standard acceleration of gravity used for Isp-ve
                      #conversion
def gf(t,x,extra=None):
    """
    Calculate gravitational acceleration on vehicle
    
    Parameters
    ----------
    t : real
        range time in time units. (SI - second) 
    r : numpy array
        Position vector of vessel relative to center of planet
        in distance units. (SI - meters)
    v : numpy array
        Inertial velocity vector in meters relative to center of planet
        in the same units as t and r (SI - m/s)
    extra : dictionary
        The passed dictionary must include the following elements (none yet 
        but this will include any parameters for the steering model)
        
    Returns
    -------
    numpy array
        Acceleration in units consistent with t, r, and v (SI - m/s^2)
    
    Notes
    -----
        This is directly an acceleration, not a force. It doesn't need
        to be multiplied by mass.
    """
    r=np.array(x[:3])
    return -r*mu/(np.linalg.norm(r)**3)

def istage(t,extra=None):
    """
    Figure out which stage is currently active. 

    Parameters
    ----------
    t : real
        range time in time units. (SI - second) 
    extra : dictionary
        The passed dictionary must include the following elements (none yet 
        but this will include any parameters for the steering model)

    Returns
    -------
    A tuple
      [0] Which stage is active. Stages below this should not include either
          their propellant or inert mass
      [1] Amount of fuel left in the current stage in kg. Zero if empty, never
          negative 
    
    For now, any current stage lower than the top stage has fuel. If it 
    doesn't, it's not the current stage.
    """
    ttot=0 #total amount of time that all previous stages have burned
    for i in range(len(extra["mp"])):
        stageTime=extra["mp"][i]/extra["mdot"][i]
        if stageTime+ttot>=t:
            tused=t-ttot #Amount of time that this stage has been running
            return (i,extra["mp"][i]-extra["mdot"][i]*tused)
        ttot=ttot+stageTime
    return (len(extra["mp"])-1,0) #All stages are empty

def mf(t,x,extra=None):
    """
    Calculate mass of vehicle
    """
    i,mp=istage(t,extra)
    m=mp+extra["m0"][i] #remaining prop and inert mass of current stage
    for j in range(i+1,len(extra["mp"])):
        m=m+extra["m0"][j]+extra["mp"][j] #All prop and inert mass of each higher stage
    return m    

Ffextra=collections.namedtuple("Ffextra",
        ["pitch","vvert","vhorz","vverthat","vhorzhat","Fv","mode","alpha","pitchgrav","pitchpoly"])
def Ff(t,x,extra=None):
    """
    Calculate thrust vector for vessel
    Parameters
    ----------
    t : real
        range time in time units. (SI - second) 
    r : numpy array
        Position vector of vessel in meters relative to center of planet
        in implied distance units. (SI - meters)
    v : numpy array
        Inertial velocity vector in meters relative to center of planet
        in the same units as t and r (SI - m/s)
    extra : dictionary
        The passed dictionary must include the following elements (none yet 
        but this will include any parameters for the steering model)

    Returns
    -------
    numpy array
        Force in units consistent with t, r, v and m() (SI - N)
    
    Notes
    -----
        This is a force, which must be multiplied by mass.
    """
    #Calculate the direction of the force
    r=np.array(x[:3])
    vorb=np.array(x[3:6])
    wind=cross(pole,r)
    vsur=vorb-wind
    #Resolve the relative velocity into vertical and horizontal projections
    vvert=dot(vsur,r)/dot(r,r)*r #projection of vsur in vertical
    vhorz=vsur-vvert             #rejection of vsur from vertical
    vverthat=r/vlength(r)        #Since the rocket may travel down, don't
                                 #use the vertical component as the basis
                                 #vector, instead use the position vector. 
    vhorzhat=vhorz/vlength(vhorz)
    pitchpoly=np.deg2rad(90-vlength(vsur)/extra["pitchover"]) #pitch down 1 degree for each 12m/s of velocity
    pitchgrav=vangle(vhorzhat,vsur)
    if vlength(vsur)>300 and pitchpoly<pitchgrav:
        pitch=pitchgrav
        mode=1
    else:
        #Velocity-dependent pitch
        pitch=pitchpoly
        mode=0
    alpha=pitch-pitchgrav
    Fv=vverthat*sin(pitch)+vhorzhat*cos(pitch)
    #Calculate the force magnitude, which depends on specific impulse, 
    #throttle, and presence of sufficient propellant
    Fm=0
    i,prop=istage(t,extra)    
    if prop>0:
        #Only add thrust if we have propellant left in this stage
        Fm=extra["ve"][i]*extra["mdot"][i]
    return Fv*Fm,Ffextra(np.rad2deg(pitch),vvert,vhorz,vverthat,vhorzhat,
                         Fv,mode,np.rad2deg(alpha),np.rad2deg(pitchgrav),np.rad2deg(pitchpoly))
    
def Df(t,x,extra=None):
    """
    Calculate aero force on vehicle

    Parameters
    ----------
    t : real
        range time in time units. (SI - second) 
    r : numpy array
        Position vector of vessel in meters relative to center of planet
        in implied distance units. (SI - meters)
    v : numpy array
        Inertial velocity vector in meters relative to center of planet
        in the same units as t and r (SI - m/s)
    extra : dictionary
        The passed dictionary must include the following elements (none yet 
        but this will include any parameters for the steering model)
        
    Returns
    -------
    numpy array
        Force in units consistent with t, r, v and m() (SI - N)
    
    Notes
    -----
        This is a force, which must be multiplied by mass. It is letter
        D for drag, but includes the total aerodynamic force on the vehicle
        in all three dimensions.
        This function may call F() to get the thrust vector, the direction
        of which implies the vehicle attitude. This may be needed for
        angle-of-attack calculation.
    """
    return np.array([0,0,0])

xdotextra=collections.namedtuple("xdotextra",["m","g","F","D","Fextra"])
def xdot(t,x,extra=None):
    m=mf(t,x,extra)
    g=gf(t,x,extra)
    F,Fextra=Ff(t,x,extra)
    Fa=F/m
    D=Df(t,x,extra)
    Da=D/m
    a=Fa+Da+g
    return np.concatenate((np.array(x[3:]),a)), xdotextra(m,g,F,D,Fextra)

def RK4(t,x,dt,extra=None):
    k1,extout=xdot(t     ,x        ,extra)
    k2       =xdot(t+dt/2,x+dt*k1/2,extra)[0]
    k3       =xdot(t+dt/2,x+dt*k2/2,extra)[0]
    k4       =xdot(t+dt  ,x+dt*k3  ,extra)[0]
    return x+dt*(k1+2*k2+2*k3+k4)/6,extout
    
def shoot(x0,t0,t1,dt=0.125,extra=None):
    """
    Shoot a trajectory.
    Parameters
    ----------
    t : real
        range time in time units. (SI - second) 
    r0 : numpy array
        Initial position vector of vessel in meters relative to center of planet
        in implied distance units. (SI - meters)
    v0 : numpy array
        Initial inertial velocity vector in meters relative to center of planet
        in the same units as t and r (SI - m/s)
    t0 : real
        initial range time in time units
    t1 : real
        final range time in time units
    extra : dictionary
        The passed dictionary will be passed to the individual force functions
    
    Returns
    -------
    A tuple:
      # final position, numpy vector
      # final velocity, numpy vector
      # list of times
      # list of positions for each time
      # list of velocities for each time
    
    Notes
    -----
    For an integration step (or sub-step for higher-order numerical integrator
    like Runge-Kutta) the acceleration is calculated for the conditions
    at the beginning of the time step. The initial velocity is presumed constant
    across the timestep, and is used to calculate the change in position. Then,
    the acceleration (which is presumed constant across the time step) is used
    to calculate the change in velocity.
    """
    tlist=[]
    xlist=[]
    extoutlist=[]
    n=0
    t=t0
    x=x0
    rl=vlength(x[:3])
    vl=vlength(x[3:6])
    vcirc=sqrt(mu/rl)
    h=rl-Re
    extout=RK4(t=t,x=x,dt=dt,extra=extra)[1]
    print(t,vl,vcirc)
    while t<t1 and h>=0 and vlength(extout.F)>0 and vl<vcirc:
        #Calculate the forces and accelerations
        tlist.append(t)
        xlist.append(x)
        x,extout=RK4(t=t,x=x,dt=dt,extra=extra)
        extoutlist.append(extout)
        n=n+1
        h=vlength(x[:3])-Re
        t=t0+dt*n
        rl=vlength(x[:3])
        vl=vlength(x[3:6])
        vcirc=sqrt(mu/rl)
    return (x,tlist,xlist,extoutlist)

if __name__ == '__main__':
    x0=np.array([Re,0,0,1,0.001,0])
    extra={'ve':[320*g0,345*g0], #Specific impulse for each stage
           'mdot':[215000/(320*g0),60000/(345*g0)], #Mass flow rate for each stage
           'mp':[4*1000,4*1000],   #Propellant mass of each FL-T800
           'm0':[(0.4    #TR16A decoupler (stage 1/2)
                 +0.5    #Empty mass of FL-T800 tank
                 +1.5    #LV-T45 engine
                )*1000
                ,(0.8    #Mk1 command pod with no monoprop
                 +0.3    #Heat shield with ablator
                 +0.1    #Mk16 Parachute
                 +0.4    #TR16A decoupler (stage 2/capsule)
                 +0.5    #Empty mass of FL-T800 tank
                 +0.5    #LV-909 engine
                )*1000], #convert tons to kg, mass for each stage
           'pitchover':10
           }    
    print(istage(60,extra))
    x1,tlist,xlist,extoutlist=shoot(x0,0,1200,extra=extra)
    print(x1)
    hlist=[]
    Xlist=[]
    Ylist=[]
    vlist=[]
    pitchlist=[]
    mlist=[]
    modelist=[]
    alphalist=[]
    ppolylist=[]
    pgravlist=[]
    for i in range(len(tlist)):
        Xlist.append(xlist[i][0])
        Ylist.append(xlist[i][1])
        vlist.append(vlength(xlist[i][3:6]))
        hlist.append(vlength(xlist[i][0:3])-Re)
        pitchlist.append(extoutlist[i].Fextra.pitch)
        mlist.append(extoutlist[i].m)
        modelist.append(extoutlist[i].Fextra.mode)
        alphalist.append(extoutlist[i].Fextra.alpha)
        ppolylist.append(extoutlist[i].Fextra.pitchpoly)
        pgravlist.append(extoutlist[i].Fextra.pitchgrav)
        print(tlist[i],xlist[i],extoutlist[i])
    rl=vlength(x1[:3])
    vl=vlength(x1[3:6])
    vcirc=sqrt(mu/rl)
    print(vl,vcirc,pitchlist[-1])
    Xlist=np.array(Xlist)
    Ylist=np.array(Ylist)
    plt.figure(7)
    surf=sqrt(Re**2-Ylist**2)-Re
    plt.plot(Ylist,Xlist-Re,Ylist,surf)
    plt.axis('equal')
    plt.xlabel("Y/m")
    plt.ylabel("X/m")
    plt.figure(8)
    plt.plot(tlist,ppolylist,tlist,pgravlist,tlist,pitchlist)
    plt.xlabel("t/s")
    plt.ylabel("pitch/deg")
    plt.show()
