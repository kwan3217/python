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

def test_su_to_cu():
   print(su_to_cu(np.array([1131340.0,-2282343.0,6672423.0]),6378137.0,398600.4415e9,1,0))
   print(np.array([0.17737781,     -0.35783850,      1.0461398]))
   print(su_to_cu(np.array([-5643.05,4303.33,2428.79]),6378137,398600.4415e9,1,-1))
   print(np.array([-0.71382529,     0.54435559,     0.30723310]))
   print(su_to_cu(2400.0,6378137.0,398600.4415e9,0,1))
   print(2.9746739)
   print(su_to_cu(np.array([-0.6616125, 0.6840739,-0.6206809]),6378137.0,398600.4415e9,1,0,inverse=True))
   print(np.array([-4219855.2,      4363117.1,     -3958787.8]))
   print(su_to_cu(np.array([ 0.4667380,-0.2424455,-0.7732126]),6378137.0,398600.4415e9,1,-1,inverse=True))
   print(np.array([  3689.7346,      -1916.6203,      -6112.5284]))

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
            
if __name__=="__main__":
   test_su_to_cu()
   print(elorb(np.array([1,0,0]),np.array([0,0,1])))

