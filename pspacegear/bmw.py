"""
Algorithms for solving the Kepler and Gauss problems. the two fundamental problems in
two-body mechanics

The <i>Kepler Problem</i>, or prediction problem, is as follows: Given the GM of the
central body, and the position and velocity of a test particle of negligible mass at 
one time, find the position and velocity at any other time.

The <i>Gauss Problem</i>, or targeting problem, is as follows: Given the GM of the 
central body, the position at which a test particle is now, the target position where
it will be, and the time between the postions, find the velocity of the particle at 
both points such that it will travel from the initial position to the target position
in the given time.

This basically follows the Universal Variable formation found in Bate, Muller, and White
chapters 4 and 5.

These algorithms make use of canonical units. Canonical units are distance and time units
relating to a particular central body, such that the GM of that body is 1. In canonical 
units, an object in a circular orbit of radius one Distance Unit (DU) has a speed of
one DU per Time Unit (TU), and therefore an angular velocity of one radian per TU. The 
length of a DU is arbitrary, but is customarily the radius of the central body for planets
and moons, and 1 AU for the Sun.   
"""

import numpy as np
from collections import namedtuple
import math

def su_to_cu(x,a,mu,LL,TT,inverse=False):
    """
    Convert standard units to canonical units
    input
  x - measurement to convert, may be scalar or any size or shape of vector
  a - length of canonical distance unit in standard units. Standard length unit is implied by this. 
  mu - gravitational constant in standard distance and time units. Standard length unit same as 
       above, standard time unit implied by this.
  LL - power of length dimension used in x
  TT - power of time dimension used in x
  /inverse - convert x from canonical units back to standard units
return
  a scalar or array the same size as x, in canonical units (or standard if /inv was set)
Example
 An object orbiting Earth has a position of <1131340,-2282343,6672423> m 
 and a speed of <-5643.05,4303.33,2428.79> m/s. Convert this to canonical units
 Earth radius used as distance unit length: 6378137m
 Earth gravitational constant: 398600.4415d9 m,s
 print,su_to_cu([1131340d,-2282343d,6672423d],6378137d,398600.4415d9,1,0)
      0.17737781     -0.35783850       1.0461398
 print,su_to_cu([-5643.05d,4303.33d,2428.79d],6378137d,398600.4415d9,1,-1)
     -0.71382529      0.54435559      0.30723310
 We are going to solve the Kepler problem over a time of 40min=2400s. How many canonical time units?
 print,su_to_cu(2400d,6378137d,398600.4415d9,0,1)
       2.9746739
 The answer is r_t=<-0.6616125, 0.6840739,-0.6206809> and
               v_t=< 0.4667380,-0.2424455,-0.7732126>. What is this in SI units?
 print,su_to_cu([-0.6616125d, 0.6840739d,-0.6206809d],6378137d,398600.4415d9,1,0,/inv)
      -4219855.2       4363117.1      -3958787.8
 print,su_to_cu([ 0.4667380d,-0.2424455d,-0.7732126d],6378137d,398600.4415d9,1,-1,/inv)
       3689.7346      -1916.6203      -6112.5284
To convert:  Standard Unit to       Canonical Unit          Multiply by
Distance     m                      DU                      1/a
Time         s                      TU                      sqrt(mu/a**3)
For derived units, raise the base unit for each dimension to the power of the dimension needed, then multiply. For example
Speed        m/s                    DU/TU (LL=1,TT=-1)      1/a*sqrt(a**3/mu)
For inverse conversion, divide instead of multiply
"""
    DU=(1.0/a)**LL
    TU=np.sqrt(mu/a**3)**TT
    DUTU=DU*TU
    if inverse:
        DUTU=1.0/DUTU
    return x*DUTU

def vangle(a,b):
   return np.arccos(np.dot(a,b)/np.linalg.norm(a)/np.linalg.norm(b))

Elorb=namedtuple('elorb',['p','a','e','i','an','ap','ta','tp','rp','MM','n','t'])
def elorb(r,v,l_DU=1,mu=1,t0=0.0):
    """
    Given state vector, calculate orbital elements.

    input
      r - position vector, in distance units implied by mu, any inertial frame is fine, but center of attraction is assumed to be
          at the origin of the frame
      v - inertial velocity vector, distance and time units implied by mu
          must be in same frame as r
      l_du - length of a distance unit, used for conversion to canonical units internally
      mu - gravity parameter, implies distance and time units
    return
      a named tuple
        p:  semi-parameter, distance from focus to orbit at TA=+-90deg, in original distance units, always positive for any eccentricity
        a:  semimajor axis, in original distance units
        e:  eccentricity
        i:  inclination, radians
        an: longitude of ascending node, angle between x axis and line of intersection between orbit plane and xy plane, radians
        ap: argument of periapse, angle between xy plane and periapse along orbit plane, radians
        ta: true anomaly, angle between periapse and object, radians
        tp: time to next periapse, in original time units. Negative if only one periapse and in the past
        rp: radius of periapse in original distance units
    """
    rv=su_to_cu(r,l_DU,mu,1,0)
    vv=su_to_cu(v,l_DU,mu,1,-1)
    r=np.linalg.norm(rv)
    v=np.linalg.norm(vv)
        
    hv=np.cross(rv,vv)
    h=np.linalg.norm(hv)
    nv=np.cross([0,0,1],hv)
    n=np.linalg.norm(nv)
    ev=(v**2-1.0/r)*rv-np.dot(rv,vv)*vv
    e=np.linalg.norm(ev)
    
    xi=v**2/2.0-1/r
    if e!=1.0:
        a=-1/(2*xi)
        p=a*(1-e**2)
    else:
        #parabolic case
        p=h**2
        a=math.inf
    (hx,hy,hz)=hv
    (nx,ny,nz)=nv
    i=np.arccos(hz/h)
    an=np.arccos(nx/n)
    if ny<0:
        an=2*np.pi-an
    ap=vangle(ev,nv)
    (ex,ey,ez)=ev
    if ez<0:
        ap=2*np.pi-ap
    ta=vangle(ev,rv)
    if np.dot(rv,vv)<0:
        ta=2*np.pi-ta
    if not np.isfinite(a):
        EE=np.tan(ta)/2
        MM=EE**3/3+EE
        n=np.sqrt(2/(p**3))
        rp=p/2
    elif a>0:
        EE=2*np.arctan(np.sqrt((1-e)/(1+e))*np.tan(ta/2))
        MM=EE-e*np.sin(EE)
        n=np.sqrt(1/(a**3))
        rp=a*(1-e)
    elif a<0:
        EE=np.arcsinh(np.sin(ta)*np.sqrt(e**2-1)/(1+e*np.cos(ta)))
        MM=e*np.sinh(EE)-EE
        n=np.sqrt(-1/(a**3))
        rp=a*(1-e)
    tp=-MM/n
    return Elorb(p =su_to_cu(p,l_DU,mu,1,0,inverse=True),   
                 a =su_to_cu(a,l_DU,mu,1,0,inverse=True),
                 e =e,                    
                 i =i,
                 an=an,
                 ap=ap,
                 ta=ta,
                 tp=su_to_cu(tp,l_DU,mu,0,1,inverse=True)+t0,
                 rp=su_to_cu(rp,l_DU,mu,1,0,inverse=True),
                 MM=MM, 
                 n =su_to_cu(n,l_DU,mu,0,-1,inverse=True),
                 t =su_to_cu(2*math.pi*np.sqrt(a**3),l_DU,mu,0,1,inverse=True))

def CC(z):
    """
Calculate the Universal Variable C(z) function
    """
    if z>0:
        return (1-np.cos(np.sqrt(z)))/z
    elif z==0:
        return 0.5
    else:
        return (1-np.cosh(np.sqrt(-z)))/z

def SS(z):
    """
Calculate the Universal Variable S(z) function
    """
    if z>0:
        sz=np.sqrt(z)
        return (sz-np.sin(sz))/sz**3
    elif z==0:
        return 1.0/6.0
    else:
        sz=np.sqrt(-z)
        return (np.sinh(sz)-sz)/np.sqrt((-z)**3)

def kepler(rv0_,vv0_,t_,l_DU=1,mu=1,eps=1e-9):
    """
Given a state vector and time interval, calculate the state after the time interval elapses
input
  rv0  - position vector relative to central body
  vv0  - velocity vector relative to central body
  t    - time interval from given state to requested state - may be negative
  l_DU=- optional - used to convert rv0 and vv0 to canonical units. Length of distance canonical
         unit in standard units implied by rv0. If not set, input state and time are already
         presumed to be in canonical units
  mu=  - optional, but required if l_DU= is set - used to convert to canonical units. Gravitational
         constant in standard units
  eps= - optional - loop termination criterion, default to 1e-9
output
  rv_t= - Position vector after passage of time t, in same units as rv0 
  vv_t= - Velocity vector after passage of time t, in same units as vv0
    """
    if t_==0.0:
        #Shortcut if we ask for zero time interval
        return (rv0_,vv0_)
    tau=np.pi*2.0
    rv0=su_to_cu(rv0_,l_DU,mu,1, 0)
    vv0=su_to_cu(vv0_,l_DU,mu,1,-1)
    t  =su_to_cu(t_  ,l_DU,mu,0, 1)

    r0=np.linalg.norm(rv0)
    v0=np.linalg.norm(vv0)
  
    r0dv0=np.dot(rv0,vv0)
  
    #Determine specific energy, and thereby the orbit shape
    #These are dependent on the initial state and therefore scalar
    E=v0**2/2-1/r0
    alpha=-2*E #alpha is 1/a, reciprocal of semimajor axis (in case of parabola. a can be infinite, but never zero
  
    #Starting guess for x
    if alpha > 0:
        #elliptical
        x0=t*alpha
    elif alpha==0:
        #parabolic (this will never really happen)
        hv=np.cross(rv0,vv0)
        p=np.dot(hv,hv)
        #acot(x)=tau/4-atan(x)
        s=(tau/4.0-np.atan(3.0*t*np.sqrt(1.0/p**2)))/2.0
        w=np.atan(np.tan(s)**(1.0/3.0))
        x0=np.sqrt(p)/np.tan(2*w) #cot(x)=1/tan(x)
    else:
        #hyperbolic
        sa=np.sqrt(-1.0/alpha)
        st=1 if t>0 else -1
        x0_a=st*sa
        x0_n=-2*alpha*t
        x0_d=r0dv0+st*sa*(1-r0*alpha)
        x0=x0_a*np.log(x0_n/x0_d)

    done=False
    xn=x0
    while not done:
        z=xn**2*alpha
        C=CC(z)
        S=SS(z)
        r=xn**2*C+r0dv0*xn*(1-z*S)+r0*(1.0-z*C)
        tn=xn**3*S+r0dv0*xn**2*C+r0*xn*(1.0-z*S)
        xnp1=xn+(t-tn)/r
        done=np.abs(xn-xnp1)<eps
        xn=xnp1
    x=xn
    f=1-x**2*C/r0
    g=t-x**3*S
    fdot=x*(z*S-1)/(r*r0)
    gdot=1-x**2*C/r
    rv_t=f   *rv0+g   *vv0
    vv_t=fdot*rv0+gdot*vv0
    if l_DU!=1:
        rv_t=su_to_cu(rv_t,l_DU,mu,1, 0,inverse=True)
        vv_t=su_to_cu(vv_t,l_DU,mu,1,-1,inverse=True)
    return(rv_t,vv_t)

def gauss(rv1_,rv2_,t_,Type=-1,l_DU=1,mu=1,eps=1e-9):
    def FindTTrialCore(A,S,X,Y):
        return (X**3)*S+A*np.sqrt(Y)
    def FindTTrial(A,r1,r2,Z):
        S=SS(Z)
        C=CC(Z)
        if C==0:
            return float('inf')
        Y=r1+r2-A*(1-Z*S)/np.sqrt(C)
        X=np.sqrt(Y/C)
        return FindTTrialCore(A,S,X,Y)
    def FindZLo(A,r1,r2):
        eps=1e-9
        #Find the Z which results in a Y of exactly zero, by bisection
        Zhi=0.0;
        Y=1.0;
        Zlo=-1.0;

        while Y>0:
            Zlo*=2.0;
            Y=r1+r2-A*(1-Zlo*SS(Zlo))/np.sqrt(CC(Zlo))
  
        while True: #Emulate a repeat/until loop
            Z=(Zlo+Zhi)/2
            Y=r1+r2-A*(1-Z*SS(Z))/np.sqrt(CC(Z))
            if Y*Zlo>0:
                Zlo=Z
            else:
                Zhi=Z
            if np.abs(Zlo-Zhi)<eps: #repeat until this condition is true
                break
        return Z+1e-5

    def FindZLo2(A,r1,r2,T):
        #Find the Z which results in a TTrial of less than T
        Z=-1.0
        TTrial=FindTTrial(A,r1,r2,Z)-T
        while TTrial>0:
            Z*=2
            TTrial=FindTTrial(A,r1,r2,Z)-T
        return Z

    tau=2*np.pi
    rv1=su_to_cu(rv1_,l_DU,mu,1,0)
    rv2=su_to_cu(rv2_,l_DU,mu,1,0)
    t  =su_to_cu(t_  ,l_DU,mu,0,1)

    if(Type<0):
        pole=np.cross(rv1,rv2)
        if(pole[2])>0:
            #prograde is short way
            Type=(-Type-1)*2+1
        else:
            #prograde is long way way
            Type=(-Type-1)*2+2

    if t<0:
        raise ValueError("Time to intercept is negative. Time travel is not allowed in this universe!")
    r1=np.linalg.norm(rv1)
    r2=np.linalg.norm(rv2)
    r1dr2=np.dot(rv1,rv2);
    DeltaNu=vangle(rv1,rv2);
    Revs=(Type-1)/2;
    #short-way and long-way are reversed for odd-numbers of complete revs
    if ((Revs%2)==1) ^ ((Type%2)==1):
        #Short way
        DM=1.0
    else:
        #Long way
        DM=-1.0
        DeltaNu=tau-DeltaNu
    if Revs>0:
        minA=r1/2.0
        minT=Revs*np.sqrt(tau*minA**3)
        if minT>t:
            raise ValueError("Can't do it! Minimum trip time for %d revs is %fTU, more than requested %fTU" %(Revs,minT,t))
    A=DM*np.sqrt(r1*r2*(1+np.cos(DeltaNu)))
    if(Revs<1):
        #less than one rev
        if Type==1:
            Zlo=FindZLo(A,r1,r2)
        else:
            Zlo=FindZLo2(A,r1,r2,t)
        Zhi=tau**2;
    else:
        #more than one rev
        #Use Zeno's method
        Zlo=((2*Revs+1)*tau/2.0)^2# Z that gives the lowest TIME, not necessarily lowest Z
        #Zbound is the value of Z which gives an infinite T
        if (Type%2)==1:
            Zbound=(Revs*tau)**2
        else:
            Zbound=((Revs+1)*tau)**2
        Zhi=(Zbound+Zlo)/2.0 #Z that gives the highest TIME, not necessarily highest Z
        while True: #emulate repeat/until 
            Thi=FindTTrial(A,r1,r2,Zhi)
            Zhi=(Zbound+Zhi)/2; #Split the difference between current Zhi and bound
            if Thi>=t:
                break

    #Solve it by bisection
    tnlo=FindTTrial(A,r1,r2,Zlo)
    tnhi=FindTTrial(A,r1,r2,Zhi)
    while True: #emulate repeat/until
        Z=(Zlo+Zhi)/2.0
        tn=FindTTrial(A,r1,r2,Z)
        if (t-tn)*tnlo>0:
            Zlo=Z
        else:
            Zhi=Z
        if abs(Zlo-Zhi)<=eps:
            break
    
    S=(SS(Z))
    C=(CC(Z))
    Y=r1+r2-A*(1.0-Z*S)/np.sqrt(C)
    f=1.0-Y/r1
    g=A*np.sqrt(Y)
    gdot=1.0-Y/r2
    vv1=su_to_cu((rv2     -rv1*f)/g,l_DU,mu,1,-1,inverse=True)
    vv2=su_to_cu((rv2*gdot-rv1  )/g,l_DU,mu,1,-1,inverse=True)
    return (vv1,vv2)

def herrick_gibbs(rr1,rr2,rr3,t1,t2,t3,mu=1):
    """
    Given three closely-spaced position observations, calculate the
    velocity at the middle observation. From Vallado p444, algorithm 52

    This is much less computationally intensive than the Gauss method,
    and has a better chance of numerical stability.

    :param rr1: Position vector at t1
    :param rr2: Position vector at t2
    :param rr3: Position vector at t3
    :param t1:  Time of first position vector
    :param t2:  Time of second position vector
    :param t3:  Time of third position vector
    :param mu:  Gravitational parameter
    :return:    Velocity at t2
    """
    dt31=t3-t1
    dt32=t3-t2
    dt21=t2-t1
    r1=np.linalg.norm(rr1)
    r2=np.linalg.norm(rr2)
    r3=np.linalg.norm(rr3)
    #Coplanarity (not strictly needed)
    Z23=np.cross(rr2,rr3)
    alpha_cop=np.pi/2-np.arccos(np.dot(Z23,rr1)/(np.linalg.norm(Z23)*np.linalg.norm(rr1)))
    #position spread (not strictly needed)
    cosalpha12=np.dot(rr1,rr2)/(r1*r2)
    cosalpha23=np.dot(rr2,rr3)/(r2*r3)
    vv2=      -dt32 *(1/(dt21*dt31)+mu/(12*r1**3))*rr1+\
         (dt32-dt21)*(1/(dt21*dt32)+mu/(12*r2**3))*rr2+\
               dt21 *(1/(dt32*dt31)+mu/(12*r3**3))*rr3
    return vv2

def test_su_to_cu():
    """
Example
 An object orbiting Earth has a position of <1131340,-2282343,6672423> m
 and a speed of <-5643.05,4303.33,2428.79> m/s. Convert this to canonical units
 Earth radius used as distance unit length: 6378137m
 Earth gravitational constant: 398600.4415d9 m,s
 print,su_to_cu([1131340d,-2282343d,6672423d],6378137d,398600.4415d9,1,0)
      0.17737781     -0.35783850       1.0461398
 print,su_to_cu([-5643.05d,4303.33d,2428.79d],6378137d,398600.4415d9,1,-1)
     -0.71382529      0.54435559      0.30723310
 We are going to solve the Kepler problem over a time of 40min=2400s. How many canonical time units?
 print,su_to_cu(2400d,6378137d,398600.4415d9,0,1)
       2.9746739
 The answer is r_t=<-0.6616125, 0.6840739,-0.6206809> and
               v_t=< 0.4667380,-0.2424455,-0.7732126>. What is this in SI units?
 print,su_to_cu([-0.6616125d, 0.6840739d,-0.6206809d],6378137d,398600.4415d9,1,0,/inv)
      -4219855.2       4363117.1      -3958787.8
 print,su_to_cu([ 0.4667380d,-0.2424455d,-0.7732126d],6378137d,398600.4415d9,1,-1,/inv)
       3689.7346      -1916.6203      -6112.5284
    """
    re=6378137.0
    mu=398600.4415e9
    print(su_to_cu(np.array([1131340.0 ,-2282343.0 ,6672423.0 ]),re,mu,1, 0))
    print(np.array([0.17737781,     -0.35783850,      1.0461398]))
    print(su_to_cu(np.array([  -5643.05,    4303.33,   2428.79]),re,mu,1,-1))
    print(np.array([-0.71382529,     0.54435559,     0.30723310]))
    print(su_to_cu(2400.0,6378137.0,398600.4415e9,0,1))
    print(2.9746739)
    print(su_to_cu(np.array([-0.6616125, 0.6840739,-0.6206809]),re,mu,1,0,inverse=True))
    print(np.array([-4219855.2,      4363117.1,     -3958787.8]))
    print(su_to_cu(np.array([ 0.4667380,-0.2424455,-0.7732126]),re,mu,1,-1,inverse=True))
    print(np.array([  3689.7346,      -1916.6203,      -6112.5284]))

def test_kepler():
    #From canonical values in Vallado, example 2-4, p102-103
    r0_cu=np.array([ 0.17738,-0.35784,1.04614])
    v0_cu=np.array([-0.71383, 0.54436,0.30723])
    r1_cu_standard=np.array([-0.6616125, +0.6840739, -0.6206809])
    v1_cu_standard=np.array([ 0.4667380, -0.2424455, -0.7732126])
    dt=2.974674
    (r1_cu,v1_cu)=kepler(r0_cu,v0_cu,dt)
    print("test_kepler() Calculated: ",r1_cu,v1_cu)
    print("test_kepler() Documented: ",r1_cu_standard,v1_cu_standard)

def plot_kepler():
    import matplotlib.pyplot as plt
    r0_cu=np.array([ 0.17738,-0.35784,1.04614])
    v0_cu=np.array([-0.71383, 0.54436,0.30723])
    print("plot_kepler elorb: ",elorb(r0_cu, v0_cu))
    r_cu=np.zeros([50,3])
    for i,dt in enumerate(np.linspace(0,2.974674)):
        (r1_cu,v1_cu)=kepler(r0_cu,v0_cu,dt)
        r_cu[i,:]=r1_cu
    plt.subplot(221)
    plt.plot(r_cu[:,0],r_cu[:,1])
    plt.subplot(222)
    plt.plot(r_cu[:,0],r_cu[:,2])
    plt.subplot(223)
    plt.plot(r_cu[:,1],r_cu[:,2])
    plt.show()

def test_gauss1():
    #Canonical unit numbers from Vallado example 7-5, p467
    #and time between, and find velocities
    r0_cu=np.array([2.5,0.0,0.0])
    v0_cu_standard=np.array([0.2604450,0.3688589,0.0])
    r1_cu=np.array([1.915111,1.606969,0.0])
    v1_cu_standard=np.array([-0.4366104,0.1151515,0.0])
    dt=5.6519
    (v0_cu,v1_cu)=gauss(r0_cu,r1_cu,dt)
    print("test_gauss1() Calculated: ",v0_cu,v1_cu)
    print("test_gauss1() Documented: ",v0_cu_standard,v1_cu_standard)

def test_gauss2():
    #Same numbers as from test_kepler, but we use begin and end positions
    #and time between, and find velocities
    r0_cu         =np.array([ 0.17738,-0.35784, 1.04614])
    v0_cu_standard=np.array([-0.71383, 0.54436, 0.30723])
    r1_cu         =np.array([-0.6616125, +0.6840739, -0.6206809])
    v1_cu_standard=np.array([ 0.4667380, -0.2424455, -0.7732126])
    dt=2.974674
    (v0_cu,v1_cu)=gauss(r0_cu,r1_cu,dt,Type=1)
    print("test_gauss2() Calculated: ",v0_cu,v1_cu)
    print("test_gauss2() Documented: ",v0_cu_standard,v1_cu_standard)

def test_herrick_gibbs():
    rr1=np.array([3419.85564,6019.82602,2784.60022])
    rr2=np.array([2935.91195,6326.18324,2660.59584])
    rr3=np.array([2434.95205,6597.38674,2521.52311])
    t1=0
    t2=1*60+16.48
    t3=2*60+33.04
    re=6378.1363
    mu=398600.4415
    rr1_cu=su_to_cu(rr1,re,mu,1, 0)
    rr2_cu=su_to_cu(rr2,re,mu,1, 0)
    rr3_cu=su_to_cu(rr3,re,mu,1, 0)
    t1_cu =su_to_cu(t1 ,re,mu,0, 1)
    t2_cu =su_to_cu(t2 ,re,mu,0, 1)
    t3_cu =su_to_cu(t3 ,re,mu,0, 1)
    vv2_cu=herrick_gibbs(rr1_cu,rr2_cu,rr3_cu,t1_cu,t2_cu,t3_cu)
    vv2=su_to_cu(vv2_cu,re,mu,1,-1,inverse=True)
    print(vv2)

if __name__=="__main__":
    test_herrick_gibbs()
    #test_kepler()
    #test_gauss1()
    #test_gauss2()
    #print(elorb(np.array([1,0,0]),np.array([0,0,1])))

