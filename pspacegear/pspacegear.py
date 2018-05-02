'''
Created on Jan 17, 2017

@author: chrisj
'''

import numpy as np
from numpy.linalg import norm as vlength
from numpy import cross, dot, cos, sin, arccos as acos, sqrt, pi
import matplotlib.pyplot as plt
import collections
from atmosphere.earth import lower_atmosphere as atm
from scipy.interpolate import interp1d

class shooter:
    pass

def vangle(a, b):
    return acos(dot(a, b) / (vlength(a) * vlength(b)))


# Earth
mu=398600.4415e9  #Gravitational parameter of Earth
#w=2*pi/86164.09 #Rotation speed of Earth
Re=6378137         #Equatorial radius of Earth
g0=9.80665         #Standard 1-G magnitude in m/s^2

# Kerbin
#mu = 3.5316000e12  # Gravitational parameter, m^3/s^2
#w=2*pi/(21549.425) #Rotation speed, rad/s (calculated from Sidereal period)
#Re = 600000  # Equatorial radius, m

#Rotation
w = 0
pole = w * np.array([0, 0, 1])  # Direction is direction of North Pole, magnitude is spin
                                # rate in radians/sec. It happens that the wind speed
                                # at any point is the cross product of this vector with
                                # the position vector of the point
#g0 = 9.81  # m/s^2, standard acceleration of gravity used for Isp-ve

#Sea-level atmosphere
a1=atm(0)

# conversion
def gf(t, x, extra=None):
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
    r = np.array(x[:3])
    return -r * mu / (np.linalg.norm(r) ** 3)

def istage(t, extra=None):
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
    ttot = 0  # total amount of time that all previous stages have burned
    for i in range(len(extra["mp"])):
        stageTime = extra["mp"][i] / extra["mdot"][i]
        if stageTime + ttot >= t:
            tused = t - ttot  # Amount of time that this stage has been running
            return (i, extra["mp"][i] - extra["mdot"][i] * tused)
        ttot = ttot + stageTime
    return (len(extra["mp"]) - 1, 0)  # All stages are empty


def mf(t, x, extra=None):
    """
    Calculate mass of vehicle
    """
    i, mp = istage(t, extra)
    m = mp + extra["m0"][i]  # remaining prop and inert mass of current stage
    for j in range(i + 1, len(extra["mp"])):
        m = m + extra["m0"][j] + extra["mp"][j]  # All prop and inert mass of each higher stage
    for tdrop in extra["mdrop"]:
        if t<tdrop:
            m+=extra["mdrop"][tdrop]
    return m


Ffextra = collections.namedtuple("Ffextra",
                                 ["pitch", "vvert", "vhorz", "vverthat", "vhorzhat", "Fv", "mode", "alpha", "pitchgrav",
                                  "pitchpoly"])


def Ff(t, x, a, m, extra=None):
    """

    :param t: Range time in time units (SI - second)
    :param x: State vector in length and time units (SI - m and m/s)
    :param a: Atmosphere properties at given altitude
    :param m: Mass of vehicle (SI - kg)
    :param extra: Dictionary with whatever other parameters are needed
    :return: A tuple with three elements:
        Force vector in units consistent with t, x, and m (SI - N)
        mdot for each stage in units consistent with t and m (SI - kg/s)
        Extra output structure (Ffextra named tuple)
    Notes
    -----
        This is a force, which must be multiplied by mass.
    """
    # Calculate the direction of the force
    r = np.array(x[:3])
    vorb = np.array(x[3:6])
    props=np.array(x[6:]) #Propellant masses remaining in kg
    wind = cross(pole, r)
    vsur = vorb - wind
    # Resolve the relative velocity into vertical and horizontal projections
    vvert = dot(vsur, r) / dot(r, r) * r  # projection of vsur in vertical
    vhorz = vsur - vvert  # rejection of vsur from vertical
    vverthat = r / vlength(r)  # Since the rocket may travel down, don't
                               # use the vertical component as the basis
                               # vector, instead use the position vector.
    vhorzhat = vhorz / vlength(vhorz)
    pitchpoly = np.deg2rad(90 - vlength(vsur) / extra["pitchover"])  # pitch down 1 degree for each 12m/s of velocity
    pitchgrav = vangle(vhorzhat, vsur)
    if vlength(vsur) > 300 and pitchpoly < pitchgrav:
        pitch = pitchgrav
        mode = 1
    else:
        # Velocity-dependent pitch
        pitch = pitchpoly
        mode = 0
    alpha = pitch - pitchgrav
    Fv = vverthat * sin(pitch) + vhorzhat * cos(pitch)
    # Calculate the force magnitude, which depends on specific impulse,
    # throttle, and presence of sufficient propellant
    Fm = 0
    Fm_max=50*m #Maximum thrust, which will limit to maximum acceleration of 5g
    i=0
    for prop in props:
        if prop>0:
            break
        i+=1
    prop=props[i]
    if prop > 0:
        # Only add thrust if we have propellant left in this stage
        ve=extra["ve0"][i]*(1-a.P/a1.P)+extra["ve1"][i]*(a.P/a1.P) #Weighted average of sea-level and vacuum ve
        Fm = ve * extra["mdot"][i]
        mdot=extra["mdot"][i]
    return Fv * Fm, Ffextra(np.rad2deg(pitch), vvert, vhorz, vverthat, vhorzhat,
                            Fv, mode, np.rad2deg(alpha), np.rad2deg(pitchgrav), np.rad2deg(pitchpoly))

#Axial force for Atlas SLV3. This can be taken as drag coefficient as long as the rocket has zero angle of attack.
drag_M =np.array((0.0  ,0.25 ,0.5  ,0.6  ,0.7  ,0.8  ,0.85 ,0.9  ,0.93 ,0.95 ,1.0  ,1.05 ,1.1  ,1.15 ,1.25 ,1.4  ,1.5  ,1.75 ,2.0  ,2.5  ,3.5  ,4.5  ,6.0  ,8.0  ,10.0  ))
drag_Ca=np.array((0.373,0.347,0.345,0.350,0.365,0.391,0.425,0.481,0.565,0.610,0.725,0.760,0.773,0.770,0.740,0.665,0.622,0.530,0.459,0.374,0.303,0.273,0.259,0.267, 0.289))
f_Ca=interp1d(drag_M,drag_Ca)

#drag_r=(5*0.3048) #Effective radius of rocket body, m. This is used to describe a circle which is the drag reference area
drag_r=0.625
drag_S=drag_r**2*np.pi

Dfextra = collections.namedtuple("Dfextra",
                                 ["spd","M", "q", "Ca", "Fa"])
def Df(t, x, a, m, F, extra=None):
    """
    Calculate aero force on vehicle

    :param t: Range time in time units (SI - second)
    :param x: State vector in length and time units (SI - m and m/s)
    :param a: Atmosphere properties at given altitude
    :param m: Mass of vehicle (SI - kg)
    :param F: Thrust vector in force units consistent with t, x, and m (SI - N)
              Passed so that this function can calculate angle of attack
    :param extra: Dictionary with whatever other parameters are needed
    :return: Force in units consistent with t, x, and m (SI - N)

    Notes
    -----
        This is a force, which must be multiplied by mass. It is letter
        D for drag, but includes the total aerodynamic force on the vehicle
        in all three dimensions.
    """
    r = np.array(x[:3])
    vorb = np.array(x[3:6])
    wind = cross(pole, r)
    vsur = vorb - wind
    spd=vlength(vsur)
    M=spd/a.Cs #Mach number
    q=a.rho*spd**2/2
    if M>drag_M[-1]:
        Ca=drag_Ca[-1]
    elif M<drag_M[0]:
        Ca=drag_Ca[0]
    else:
        Ca=f_Ca(M)
    Fa=q*drag_S*Ca
    if spd==0:
        return np.zeros(3),Dfextra(spd,M,q,Ca,Fa)
    if a.rho==0:
        return np.zeros(3),Dfextra(spd,M,q,Ca,Fa)
    return -Fa*vsur/spd, Dfextra(spd,M,q,Ca,Fa)

xdotextra = collections.namedtuple("xdotextra", ["m", "g", "Z", "atm", "F", "D", "Fextra", "Dextra"])

def xdot(t, x, extra=None):
    Z=vlength(x[:3])-Re
    at=atm(Z)
    m = mf(t=t, x=x, extra=extra)
    g = gf(t=t, x=x, extra=extra)
    F, Fextra = Ff(t=t, x=x, a=at, m=m, extra=extra)
    D, Dextra = Df(t=t, x=x, a=at, m=m, F=F, extra=extra)
    a = (F + D)/m + g
    return np.concatenate((np.array(x[3:]), a)), xdotextra(m, g, Z, at, F, D, Fextra, Dextra)

def RK4(t, x, dt, extra=None):
    k1, extout = xdot(t, x, extra)
    k2 = xdot(t + dt / 2, x + dt * k1 / 2, extra)[0]
    k3 = xdot(t + dt / 2, x + dt * k2 / 2, extra)[0]
    k4 = xdot(t + dt, x + dt * k3, extra)[0]
    return x + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6, extout

def shoot(x0, t0, t1, dt=0.125, extra=None):
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
    tlist = []
    xlist = []
    extoutlist = []
    n = 0
    t = t0
    x = x0
    rl = vlength(x[:3])
    vl = vlength(x[3:6])
    vcirc = sqrt(mu / rl)
    h = rl - Re
    extout = RK4(t=t, x=x, dt=dt, extra=extra)[1]
    print(t, vl, vcirc)
    while t < t1 and h >= 0 and vlength(extout.F) > 0 and vl < vcirc:
        # Calculate the forces and accelerations
        tlist.append(t)
        xlist.append(x)
        x2, extout = RK4(t=t, x=x, dt=dt, extra=extra)
        vl = vlength(x2[3:6])
        if not np.isfinite(vl):
            print("Something happened!")
        extoutlist.append(extout)
        h = vlength(x[:3]) - Re
        rl = vlength(x2[:3])
        vcirc = sqrt(mu / rl)
        n = n + 1
        x=x2
        t = t0 + dt * n
    if not (t<t1):
        term="Time expired"
    elif not (h>=0):
        term="Crashed into ground"
    elif not (vlength(extout.F)>0):
        term="Out of fuel"
    elif not (vl<vcirc):
        term="Vcirc achieved"
    else:
        term="Huh? None of the termination conditions tripped"
    return (x, tlist, xlist, extoutlist,term)


if __name__ == '__main__':
    x0 = np.array([Re, 0, 0, 1, 0.001, 0])
    Sandstone={'ve0': [320 * g0, 345 * g0],  # Vacuum specific impulse for each stage
               've1': [320 * g0, 345 * g0],  #Sea-level specific impulse
             'mdot': [215000 / (320 * g0), 60000 / (345 * g0)],  # Mass flow rate for each stage
             'mp': [4 * 1000, 4 * 1000],  # Propellant mass of each FL-T800
             'm0': [(0.4  # TR16A decoupler (stage 1/2)
                     + 0.5  # Empty mass of FL-T800 tank
                     + 1.5  # LV-T45 engine
                     ) * 1000
                 , (0.8  # Mk1 command pod with no monoprop
                    + 0.3  # Heat shield with ablator
                    + 0.1  # Mk16 Parachute
                    + 0.4  # TR16A decoupler (stage 2/capsule)
                    + 0.5  # Empty mass of FL-T800 tank
                    + 0.5  # LV-909 engine
                    ) * 1000],  # convert tons to kg, mass for each stage
             'mdrop':{}, #No fairings etc to drop
             'pitchover': 15
             }
    Atlas401={'ve0': [337.8*g0,450.5*g0],  #Vacuum specific impulse for each stage
              've1': [311.3*g0,450.5*g0],  #Sea-level specific impulse
             'mdot': [4.152e6/(337.8*g0), 99.2e3 / (450.5*g0)],  # Mass flow rate derived from vacuum thrust
             'mp': [284089, 20830],        # Propellant mass of each FL-T800
             'm0': [(21054    # Inert mass of booster
                     + 947    # Interstage
                     + 181.7  # Stub adapter
                     )
                 , (2243  # Inert mass of Centaur
                    +39.5 #D1666 Payload Sep Ring
                    +32.2 #C22 Launch Vehicle Adapter, 0.120" wall thickness
                    )],
             'mdrop':{250:2127}, #LPF
             'pitchover': 15
             }
    extra=Atlas401
    print(istage(60, extra))
    x1, tlist, xlist, extoutlist, term = shoot(x0, 0, 1200, extra=extra)
    print(term)
    print(x1)
    hlist = []
    Xlist = []
    Ylist = []
    vlist = []
    pitchlist = []
    mlist = []
    modelist = []
    alphalist = []
    ppolylist = []
    pgravlist = []
    Zlist=[]
    spdlist=[]
    qlist=[]
    Falist=[]
    Calist=[]
    Mlist=[]
    Flist=[]
    for i in range(len(tlist)):
        Xlist.append(xlist[i][0])
        Ylist.append(xlist[i][1])
        vlist.append(vlength(xlist[i][3:6]))
        Flist.append(vlength(extoutlist[i].F)/extoutlist[i].m)
        hlist.append(vlength(xlist[i][0:3]) - Re)
        pitchlist.append(extoutlist[i].Fextra.pitch)
        mlist.append(extoutlist[i].m)
        modelist.append(extoutlist[i].Fextra.mode)
        alphalist.append(extoutlist[i].Fextra.alpha)
        ppolylist.append(extoutlist[i].Fextra.pitchpoly)
        pgravlist.append(extoutlist[i].Fextra.pitchgrav)
        Zlist.append(extoutlist[i].Z)
        spdlist.append(extoutlist[i].Dextra.spd)
        qlist.append(extoutlist[i].Dextra.q/extoutlist[i].m)
        Falist.append(extoutlist[i].Dextra.Fa/extoutlist[i].m)
        Calist.append(extoutlist[i].Dextra.Ca)
        Mlist.append(extoutlist[i].Dextra.M)
        print(tlist[i], xlist[i], extoutlist[i])
    print(term)
    rl = vlength(x1[:3])
    vl = vlength(x1[3:6])
    vcirc = sqrt(mu / rl)
    print(vl, vcirc, pitchlist[-1])
    Xlist = np.array(Xlist)
    Ylist = np.array(Ylist)
    plt.figure(7)
    surf = sqrt(Re ** 2 - Ylist ** 2) - Re
    plt.plot(Ylist, Xlist - Re,'b-', Ylist, surf,'g-',Ylist[0::80],Xlist[0::80]-Re,'bx')
    plt.axis('equal')
    plt.xlabel("Y/m")
    plt.ylabel("X/m")
    plt.figure(8)
    plt.plot(spdlist, ppolylist, 'b-',spdlist, pgravlist, 'g-',spdlist, pitchlist,'r--')
    plt.xlabel("spd/(m/s)")
    plt.ylabel("pitch/deg")
    plt.figure(9)
    plt.plot(tlist,Zlist)
    plt.xlabel("t/s")
    plt.ylabel("Altitude/m")
    plt.figure(10)
    plt.plot(tlist,spdlist)
    plt.xlabel("t/s")
    plt.ylabel("spd/(m/s)")
    plt.figure(11)
    plt.plot(tlist,qlist,'b-',tlist,Falist,'g-',tlist,Flist,'r-')
    plt.xlabel("t/s")
    plt.ylabel("acc/(m/s^2)")
    plt.figure(12)
    plt.plot(tlist,Calist,'b-',tlist,Mlist,'g-')
    plt.xlabel("t/s")
    plt.ylabel("Ca,Mach")
    plt.figure(13)
    plt.plot(tlist,mlist,'b-')
    plt.xlabel("t/s")
    plt.ylabel("Mass/kg")
    plt.show()
