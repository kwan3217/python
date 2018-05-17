'''
Created on Jan 17, 2017

@author: chrisj
'''

import numpy as np
from numpy.linalg import norm as vlength
import matplotlib.pyplot as plt
import collections
from atmosphere.earth import lower_atmosphere as atmf
from scipy.interpolate import interp1d
import guidance.peg
import copy

def vangle(a, b):
    return np.arccos(np.dot(a, b) / (vlength(a) * vlength(b)))

# Earth
mu=398600.4415e9  #Gravitational parameter of Earth
w=2*np.pi/86164.09 #Rotation speed of Earth
Re=6378137         #Equatorial radius of Earth
g0=9.80665         #Standard 1-G magnitude in m/s^2

# Kerbin
#mu = 3.5316000e12  # Gravitational parameter, m^3/s^2
#w=2*pi/(21549.425) #Rotation speed, rad/s (calculated from Sidereal period)
#Re = 600000  # Equatorial radius, m

#Rotation
#w = 0
pole = w * np.array([0, 0, 1])  # Direction is direction of North Pole, magnitude is spin
                                # rate in radians/sec. It happens that the wind speed
                                # at any point is the cross product of this vector with
                                # the position vector of the point
#g0 = 9.81  # m/s^2, standard acceleration of gravity used for Isp-ve

#Sea-level atmosphere
atm_sl=atmf(0)

# conversion
def gf(t, x):
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

xdotextra = collections.namedtuple("xdotextra", ["m", "g", "Z", "atm", "F", "D", "Fextra", "Dextra"])
def xdot(t=None, x=None, discrete=None, dt=None, extra=None):
    """
    Differential equation to integrate. This is specific to vessel motion.
    :param t: Range time in time units (SI - second)
    :param x: State vector. This is now specialized, so we can say for sure that:
              *Elements [0:3] are position
              *Elements [3:6] are velocity
              *Remainder of elements are fuel levels in tanks in mass units
    :param discrete: State of vessel that *isn't* integrated. Passed along to mf, through Ff, and to Df.
                     Only Ff should change it.
    :param dt: Time step in same units as t
    :param extra: Passed to mf, Ff, and Df. Intended to be read-only.
    :return: A tuple
      * Derivative of state with respect to time, so:
        *Elements [0:3] are velocity
        *Elements [3:6] are acceleration
        *Remainder are mass flow rate of each fuel tank
      * Discrete as modified by the mf, Ff, and Df functions
    """
    Z=vlength(x[:3])-Re
    atm=atmf(Z)
    g                = gf(t=t, x=x)
    m                = extra.mf(t=t, x=x, discrete=discrete)
    F, mdot  ,Fextra = extra.Ff(t=t, x=x, discrete=discrete, atm=atm, m=m, dt=dt)
    D,        Dextra = extra.Df(t=t, x=x, discrete=discrete, atm=atm, m=m, F=F)
    a=(F+D)/m + g
    return np.concatenate((np.array(x[3:6]), a, np.array(mdot))),xdotextra(m, g, Z, atm, F, D, Fextra, Dextra)

def RK4(t=None, x=None, discrete=None, dt=None, extra=None):
    """
    Numerically integrate a differential equation using the Runge-Kutta 4th order method

    :param t: Range time in time units (SI - second)
    :param x: State vector
    :param discrete: All of the state for a vessel that *isn't* integrated.
                     For instance, things like if a fairing is attached
                     or if an engine is running or shut down.
    :param dt: Time step in same time units as t
    :param extra: Passed to xdot. Intended for things like vessel description.
    :return: A tuple
       * New state vector at time t+dt
       * New discrete value(s)
       * extout from first step
    """
    k1,extout = xdot(t=t     , x=x        , discrete=discrete, dt=dt/2, extra=extra)
    k2,_      = xdot(t=t+dt/2, x=x+dt*k1/2, discrete=discrete, dt=dt/2, extra=extra)
    k3,_      = xdot(t=t+dt/2, x=x+dt*k2/2, discrete=discrete, dt=dt/2, extra=extra)
    k4,_      = xdot(t=t+dt  , x=x+dt*k3  , discrete=discrete, dt=dt/2, extra=extra)
    return x+dt*(k1+2*k2+2*k3+k4)/6, extout

def shoot(x0, t0, t1, dt=0.125, discrete=None, extra=None):
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
    discretelist=[]
    n = 0
    t = t0
    x = x0
    rl = vlength(x[:3])
    vl = vlength(x[3:6])
    vcirc = np.sqrt(mu / rl)
    h = rl - Re
    print(t, vl, vcirc)
    while t < t1 and h >= 0 and not discrete.Done:
        # Calculate the forces and accelerations
        tlist.append(t)
        xlist.append(copy.deepcopy(x))
        discretelist.append(copy.deepcopy(discrete))
        x, extout = RK4(t=t, x=x, dt=dt, discrete=discrete, extra=extra)
        extoutlist.append(extout)
        vl = vlength(x[3:6])
        if not np.isfinite(vl):
            print("Something happened!")
        h = vlength(x[:3]) - Re
        rl = vlength(x[:3])
        vcirc = np.sqrt(mu / rl)
        n = n + 1
        t = t0 + dt * n
    if not (t<t1):
        term="Time expired"
    elif not (h>=0):
        term="Crashed into ground"
    elif not (vl<vcirc):
        term="Vcirc achieved"
    else:
        term="Huh? None of the termination conditions tripped"
    return (x, tlist, xlist, extoutlist,term, discretelist)

class Discrete:
    def __init__(self,PropAttached,InertAttached,EngineOn,Done):
        self.PropAttached=PropAttached
        self.InertAttached=InertAttached
        self.EngineOn=EngineOn
        self.Done=Done

class Vessel:
    # Axial force for Atlas SLV3. This can be taken as drag coefficient as long as the rocket has zero angle of attack.
    drag_M = np.array((0.0, 0.25, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.93, 0.95, 1.0, 1.05, 1.1, 1.15, 1.25, 1.4, 1.5, 1.75,
                       2.0, 2.5, 3.5, 4.5, 6.0, 8.0, 10.0))
    drag_Ca = np.array((
                       0.373, 0.347, 0.345, 0.350, 0.365, 0.391, 0.425, 0.481, 0.565, 0.610, 0.725, 0.760, 0.773, 0.770,
                       0.740, 0.665, 0.622, 0.530, 0.459, 0.374, 0.303, 0.273, 0.259, 0.267, 0.289))
    f_Ca = interp1d(drag_M, drag_Ca)
    def __init__(self):
        pass
    def mf(self, t=None, x=None, discrete=None):
        """
        Calculate mass of vehicle
        """
        m=0
        props = x[6:]
        for i,prop in enumerate(props):
            if discrete.PropAttached[i]:
                m+=prop
        for i,m0 in enumerate(self.m0):
            if discrete.InertAttached[i]:
                m+=m0
        return m
    Dfextra = collections.namedtuple("Dfextra", ["spd","M", "q", "Ca", "Fa"])
    def Df(self, t=None, x=None, discrete=None, atm=None, m=None, F=None):
        """
        Calculate aero force on vehicle

        :param t: Range time in time units (SI - second)
        :param x: State vector in length and time units (SI - m and m/s)
        :param discrete: All of the state for a vessel that *isn't* integrated.
                         For instance, things like if a fairing is attached
                         or if an engine is running or shut down.
        :param a: Atmosphere properties at given altitude
        :param m: Mass of vehicle (SI - kg)
        :param F: Thrust vector in force units consistent with t, x, and m (SI - N)
                  Passed so that this function can calculate angle of attack
        :param extra: Dictionary with whatever other parameters are needed
        :return: A tuple
          * Force in units consistent with t, x, and m (SI - N)
          * Extra output, a namedtuple of type Dfextra
        Notes
        -----
            This is a force, which must be multiplied by mass. It is letter
            D for drag, but includes the total aerodynamic force on the vehicle
            in all three dimensions.
        """
        r = np.array(x[:3])
        vorb = np.array(x[3:6])
        wind = np.cross(pole, r)
        vsur = vorb - wind
        spd=vlength(vsur)
        M=spd/atm.Cs #Mach number
        q=atm.rho*spd**2/2
        if M>self.drag_M[-1]:
            Ca=self.drag_Ca[-1]
        elif M<self.drag_M[0]:
            Ca=self.drag_Ca[0]
        else:
            Ca=self.f_Ca(M)
        Fa=q*self.drag_S*Ca
        if spd==0: #Avoid 0/0 error (denominator is 0 in -Fa*vsur/spd below)
            return np.zeros(3),self.Dfextra(spd,M,q,Ca,Fa)
        if atm.rho==0: #Avoid 1/Nan error (M is spd/atm.Cs, and atm.Cs is NaN in vacuum)
            return np.zeros(3),self.Dfextra(spd,M,q,Ca,Fa)
        return -Fa*vsur/spd, self.Dfextra(spd,M,q,Ca,Fa)
    Ffextra = collections.namedtuple("Ffextra",
                                     ["pitch", "vvert", "vhorz", "vverthat", "vhorzhat", "Fv", "mode", "alpha", "pitchgrav",
                                      "pitchpoly"])
    def Ff(self,t=None, x=None, discrete=None, atm=None, m=None, dt=None):
        """
        Calculate thrust on vessel. This is the right place to apply the guidance system.

        :param t: Range time in time units (SI - second)
        :param x: State vector in length and time units (SI - m and m/s)
        :param discrete:
        :param atm: Atmosphere properties at current altitude
        :param m: Mass of vehicle (SI - kg)
        :param dt: Time step length in same units as t. This is the size of a Runge-Kutta substep.
                   Ideally, this wouldn't be necessary, but it is because we have to know the time step
                   in order to figure out if we will run out of fuel in the next time step.
        :return: A tuple
            *Force vector in units consistent with t, x, and m (SI - N)
            *discrete, as modified by the guidance system
            *mdot for each stage in units consistent with t and m (SI - kg/s)
            *Extra output, a namedtuple of type Ffextra
        """
        props=np.array(x[6:]) #Propellant masses remaining in kg

        Fmax=[0.0]*len(self.ve0)
        mdot=[0.0]*len(self.ve0)
        ve=[0.0]*len(self.ve0)
        for i,this_mdot in enumerate(self.mdot):
            if props[i]>0:
                ve[i]=self.ve0[i]*(1-atm.P/atm_sl.P)+self.ve1[i]*(atm.P/atm_sl.P) #Weighted average of sea-level and vacuum ve
                if self.mdot[i] * dt > props[i]:
                    #Limit thrust if engine is about to run out of propellant
                    this_mdot=props[i]/dt
                Fmax[i]=ve[i]*this_mdot
            else:
                #Engine is out of propellant
                Fmax[i]=0
        Fhat,throttle,guidextra=self.guide(t=t,x=x,discrete=discrete,Fmax=Fmax, m=m,props=props)
        # Calculate the force magnitude, which depends on specific impulse,
        # throttle, and presence of sufficient propellant
        Fmag=0
        for i in range(len(Fmax)):
            if(ve[i]>0):
                Fmag+=Fmax[i]*throttle[i]
                mdot[i]=-(Fmax[i]*throttle[i])/ve[i]
            else:
                mdot[i]=0
        return Fhat * Fmag, mdot, guidextra

class Sandstone(Vessel):
    ve0=[320 * g0, 345 * g0]  # Vacuum specific impulse for each stage
    ve1=[320 * g0, 345 * g0]  #Sea-level specific impulse
    mdot=[215000 / ve0[0], 60000 /ve0[1]]  # Mass flow rate for each stage
    mp=[4 * 1000, 4 * 1000]  # Propellant mass of each FL-T800
    m0=[(  0.4  # TR16A decoupler (stage 1/2)
         + 0.5  # Empty mass of FL-T800 tank
         + 1.5  # LV-T45 engine
        ) * 1000, (
           0.8  # Mk1 command pod with no monoprop
         + 0.3  # Heat shield with ablator
         + 0.1  # Mk16 Parachute
         + 0.4  # TR16A decoupler (stage 2/capsule)
         + 0.5  # Empty mass of FL-T800 tank
         + 0.5  # LV-909 engine
        ) * 1000]  # convert tons to kg, mass for each stage
    mdrop={}
    pitchover=15
    discrete0=None
    drag_r=0.625 #Radius of FT800 tank and Mk1 command pod
    drag_S=drag_r**2*np.pi #Reference drag area

resolved=collections.namedtuple("resolved",["zhat","ehat","nhat","hihat","hrelhat","vrel","zv","nv","eiv","erelv","hiv","hrelv"])
def resolve(x,pole):
    r=x[0:3]
    vi=x[3:6]
    #Component of velocity in vertical direction
    zhat=r/vlength(r)
    wind=np.cross(pole, r)
    ehat=wind/vlength(wind)
    nhat=np.cross(zhat,ehat)
    vrel=vi-wind
    ziv=np.dot(zhat,vi)
    zrelv=np.dot(zhat,vrel) #Should be the same as ziv
    niv=np.dot(nhat,vi)
    nrelv=np.dot(nhat,vrel) #Should be the same as niv
    eiv=np.dot(ehat,vi)
    erelv=np.dot(ehat,vrel) #Not the same as eiv
    hi=vi-zhat*ziv #Subtract off vertical component of velocity to get horizontal inertial velocity
    hiv=vlength(hi)
    hihat=hi/hiv
    hrel=vrel-zhat*ziv #Subtract off vertical component of velocity to get horizontal relative velocity
    hrelv=vlength(hrel)
    hrelhat=hrel/hrelv
    return resolved(zhat=zhat,ehat=ehat,nhat=nhat,hihat=hihat,hrelhat=hrelhat,vrel=vrel,zv=ziv,nv=niv,eiv=eiv,erelv=erelv,hiv=hiv,hrelv=hrelv)

class AtlasDiscrete(Discrete):
    def __init__(self):
        super().__init__([True,True],[True,True,True],[True,False],False)
        self.tAtlasDrop=250
        self.tCentaurStart=260
        self.tFairingDrop=268
        self.pegA=-2
        self.pegT=400
        self.pegB=-self.pegA/self.pegT
        self.pitchover=16
        self.pegLastt=float('nan')
        self.warn=True

class Atlas401(Vessel):
    ve0=[337.8*g0,450.5*g0]  #Vacuum specific impulse for each stage
    ve1=[311.3*g0,450.5*g0]  #Sea-level specific impulse
    mdot=[4.152e6/ve0[0], 99.2e3 / ve0[1]]  # Mass flow rate derived from vacuum thrust
    mp=[284089, 20830]        # Propellant mass of booster and centaur
    m0=[(21054    # Inert mass of booster
         + 947    # Interstage
         + 181.7  # Stub adapter
        ),(
          2243  # Inert mass of Centaur
           +39.5 #D1666 Payload Sep Ring
           +32.2 #C22 Launch Vehicle Adapter, 0.120" wall thickness
           +721  #Insight total payload mass
        ),(
          2127  #LPF
        )]
    discrete0=AtlasDiscrete()
    drag_r = 2  # Atlas with 4m fairing
    drag_S=drag_r**2*np.pi #Reference drag area
    GuidExtra=collections.namedtuple("GuidExtra",["fdotr","fdoth","pegextra"])
    def guide(self,t=None,x=None,discrete=None,Fmax=None,m=None,props=None):
        """
        Decide what direction the vehicle is going to point

        We pass a lot of things in that are calculated outside. These things are all calculable from x (that's how Ff
        does it) but we want to focus on doing actual guidance calculations here, and separate out as much as possible
        the resolving of vectors into components in common code.

        :param t: Range time in time units (SI - second)
        :param x: State vector in length and time units (SI - m and m/s)
        :param discrete:
        :param zv: vertical surface-relative speed in same units as state vector
        :param hvi: horizontal inertial speed in same units as state vector
        :param hvrel: horizontal surface-relative speed in same units as state vector
        :param ev: Eastbound inertial speed in same units as state vector
        :param evrel: Eastbound surface-relative speed in same units as state vector
        :param nv: Northbound surface-relative speed in same units as state vector
        :param zvhat: normalized local vertical
        :param evhat: normalized east direction
        :param nvhat: normalized north direction
        :param hvihat: normalized horizontal vector
        :param hvrelhat: normalized horizontal vector
        :param Fmax: Maximum available thrust (SI - N)
        :param m: Current mass (SI - kg)
        :return: A tuple
          * Normalized force vector
          * Fractional throttle for each engine - 0 means no thrust, 1.0 means full thrust
          * Extra guidance output
        """

        #Decide which engines are on
        Fn=[0.0]*2
        if discrete.EngineOn[0]:
            #Atlas RD-180 is on. Do thrust limiting.
            if t<100:
                Fn_max=1.0
            else:
                Fn_max=0.94 #Throttle back to 94% after 100s
            if t<230:
                a_max=5*g0 #Maximum acceleration
            else:
                a_max=4.6*g0 #Reduced maximum acceleration
            #Throttle which gives maximum acceleration
            Fn_amax=a_max*m/Fmax[0]
            if Fn_max>Fn_amax:
                Fn[0]=Fn_amax
            else:
                Fn[0]=Fn_max
            if props[0]<500:
                discrete.EngineOn[0]=False
                discrete.tAtlasDrop=t+8
                discrete.tCentaurStart=discrete.tAtlasDrop+10
                discrete.tFairingDrop=discrete.tCentaurStart+8
        elif discrete.PropAttached[0] and t>discrete.tAtlasDrop:
            discrete.PropAttached[0] = False
            discrete.InertAttached[0] = False
        if discrete.EngineOn[1]:
            #Centaur engine is on
            Fn[1]=1
        elif t>discrete.tCentaurStart:
            discrete.EngineOn[1]=True
        if discrete.InertAttached[2] and t>discrete.tFairingDrop:
            discrete.InertAttached[2]=False

        res=resolve(x,pole)
        if discrete.InertAttached[2]: #Run pitch program/gravity turn through fairing drop
            pitchpoly = np.deg2rad(90 - vlength(res.vrel) / discrete.pitchover)  # pitch down 1 degree for each 12m/s of velocity
            pitchgrav = vangle(res.hrelhat, res.vrel)
            if vlength(res.vrel) > 300 and pitchpoly < pitchgrav:
                pitch = pitchgrav
                mode = 1
            else:
                # Velocity-dependent pitch
                pitch = pitchpoly
                mode = 0
            alpha = pitch - pitchgrav
            fdotr=np.sin(pitch)
            fdoth=np.cos(pitch)
            Fv = res.zhat * fdotr + res.hrelhat * fdoth
            pegextra=None
        else:
            hT=200000
            rT=Re+hT
            vT=np.sqrt(mu/rT)
            rdotT=-100
            vqT=np.sqrt(vT**2-rdotT**2)
            #Run PEG after fairing drop
            if not np.isfinite(discrete.pegLastt):
                dt=0
                n=10
            else:
                dt=t-discrete.pegLastt
                n=1
            discrete.pegLastt=t
            (discrete.pegA,discrete.pegB,discrete.pegT,fdotr,fdoth,pegextra,discrete.warn)=guidance.peg.peg(
                A=discrete.pegA, B=discrete.pegB, T=discrete.pegT, dt=dt,
                a=Fmax[1]/m, ve=self.ve0[1], r=vlength(x[0:3]), rdot=res.zv, vq=res.hiv,
                rT=rT, rdotT=rdotT, vqT=vqT, mu=mu, n=n,warn=discrete.warn
            )
            discrete.Done=discrete.pegT<=0
            Fv=res.zhat*fdotr+res.hihat*fdoth
        return Fv,Fn,self.GuidExtra(fdotr=fdotr,fdoth=fdoth,pegextra=pegextra)

def plot_shot(x1,tlist,xlist,extoutlist,term,discretelist,fignum):
    fdotrlist = []
    for i in range(len(tlist)):
        fdotrlist.append(extoutlist[i].Fextra.fdotr)
        #print(tlist[i], xlist[i], extoutlist[i])
    print(term,discretelist[0].pitchover,tlist[-1],fdotrlist[-1])
#    plt.figure(fignum)
#    plt.plot(tlist, fdotrlist, 'b-')
#    plt.xlabel("t/s")
#    plt.ylabel("fdotr")
#    plt.show()

def plot_shot_all(x1,tlist,xlist,extoutlist,term,discretelist):
    print(term)
    print(x1)
    hlist = []
    Xlist = []
    Ylist = []
    Zlist = []
    vlist = []
    fdotrlist = []
    mlist = []
    modelist = []
    alphalist = []
    ppolylist = []
    pgravlist = []
    altlist = []
    spdlist = []
    qlist = []
    Falist = []
    Calist = []
    Mlist = []
    Flist = []
    aFlist = []
    aqlist = []
    adlist = []
    Aproplist = []
    Alist = []
    Blist = []
    Tlist = []
    for i in range(len(tlist)):
        Xlist.append(xlist[i][0])
        Ylist.append(xlist[i][1])
        Zlist.append(xlist[i][2])
        vlist.append(vlength(xlist[i][3:6]))
        Flist.append(vlength(extoutlist[i].F))
        aFlist.append(vlength(extoutlist[i].F) / extoutlist[i].m)
        hlist.append(vlength(xlist[i][0:3]) - Re)
        #       pitchlist.append(extoutlist[i].Fextra.pitch)
        fdotrlist.append(extoutlist[i].Fextra.fdotr)
        mlist.append(extoutlist[i].m)
        #       modelist.append(extoutlist[i].Fextra.mode)
        #       alphalist.append(extoutlist[i].Fextra.alpha)
        #        ppolylist.append(extoutlist[i].Fextra.pitchpoly)
        #        pgravlist.append(extoutlist[i].Fextra.pitchgrav)
        altlist.append(extoutlist[i].Z)
        spdlist.append(extoutlist[i].Dextra.spd)
        qlist.append(extoutlist[i].Dextra.q)
        aqlist.append(extoutlist[i].Dextra.q / extoutlist[i].m)
        Falist.append(extoutlist[i].Dextra.Fa)
        adlist.append(extoutlist[i].Dextra.Fa / extoutlist[i].m)
        Calist.append(extoutlist[i].Dextra.Ca)
        Mlist.append(extoutlist[i].Dextra.M)
        Aproplist.append(xlist[i][6])
        if extoutlist[i].Fextra.pegextra is not None:
            Alist.append(extoutlist[i].Fextra.pegextra.A)
            Blist.append(extoutlist[i].Fextra.pegextra.B * 100)
            Tlist.append(extoutlist[i].Fextra.pegextra.T / 100)
        else:
            Alist.append(float('NaN'))
            Blist.append(float('NaN'))
            Tlist.append(float('NaN'))
        print(tlist[i], xlist[i], extoutlist[i])
    print(term)
    rl = vlength(x1[:3])
    vl = vlength(x1[3:6])
    vcirc = np.sqrt(mu / rl)
    print(vl, vcirc)  # , pitchlist[-1])
    Xlist = np.array(Xlist)
    Ylist = np.array(Ylist)
    Zlist = np.array(Zlist)
    plt.figure(7)
    surf = np.sqrt(Re ** 2 - Zlist ** 2) - Re
    plt.plot(Zlist, Xlist - Re, 'b-', Zlist, surf, 'g-', Zlist[0::80], Xlist[0::80] - Re, 'bx')
    plt.axis('equal')
    plt.xlabel("Z/m")
    plt.ylabel("X/m")
    #    plt.figure(8)
    #    plt.plot(spdlist, ppolylist, 'b-',spdlist, pgravlist, 'g-',spdlist, pitchlist,'r--')
    #    plt.xlabel("spd/(m/s)")
    #    plt.ylabel("pitch/deg")
    plt.figure(9)
    plt.plot(tlist, altlist)
    plt.xlabel("t/s")
    plt.ylabel("Altitude/m")
    plt.figure(10)
    plt.plot(tlist, spdlist)
    plt.xlabel("t/s")
    plt.ylabel("spd/(m/s)")
    plt.figure(11)
    plt.plot(tlist, aqlist, 'b-', tlist, adlist, 'g-', tlist, aFlist, 'r-')
    plt.xlabel("t/s")
    plt.ylabel("acc/(m/s^2)")
    plt.figure(12)
    plt.plot(tlist, Calist, 'b-', tlist, Mlist, 'g-')
    plt.xlabel("t/s")
    plt.ylabel("Ca,Mach")
    plt.figure(13)
    plt.plot(tlist, mlist, 'b-')
    plt.xlabel("t/s")
    plt.ylabel("Mass/kg")
    plt.figure(14)
    plt.plot(tlist, qlist, 'b-', tlist, Falist, 'g-', tlist, Flist, 'r-')
    plt.xlabel("t/s")
    plt.ylabel("F/N")
    plt.figure(15)
    plt.plot(tlist, fdotrlist, 'b-')
    plt.xlabel("t/s")
    plt.ylabel("fdotr")
    plt.figure(16)
    plt.plot(tlist, Alist, 'b-', tlist, Blist, 'g-', tlist, Tlist, 'r-')
    plt.xlabel("t/s")
    plt.ylabel("A,B,T/100")
    plt.show()

if __name__ == '__main__':
    for i in range(40):
        extra=Atlas401()
        lat0=np.radians(34.5)
        r0=np.array([Re*np.cos(lat0),0,Re*np.sin(lat0)])
        wind = np.cross(pole, r0)
        azimuth=np.radians(148.9)
        v0=np.array([1,0.001*np.sin(azimuth),0.001*np.cos(azimuth)])+wind
        x0 = np.hstack((r0,v0,np.array(extra.mp)))
        discrete=copy.deepcopy(extra.discrete0)
        discrete.pitchover=14+i/10
        print(i,discrete.pitchover)
        x1, tlist, xlist, extoutlist, term, discretelist = shoot(x0=x0, t0=0, t1=900, discrete=discrete, extra=extra)
        plot_shot(x1,tlist,xlist,extoutlist,term,discretelist,i)
        if i==24:
            plot_shot_all(x1,tlist,xlist,extoutlist,term,discretelist)