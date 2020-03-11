"""
Stuff for Solar System integration
"""
import spiceypy as cspice
from collections import namedtuple
import numpy as np

bodytype = namedtuple("body", ["spice_id", "frame", "gm", "j2", "re"])

def getbody(bodyid,frame):
    try:
        gm=cspice.gdpool("BODY%d_GM"%bodyid,0,1)[0] #Not finding a GM is fatal
    except:
        print("Could not find value for BODY%d_GM"%bodyid)
        raise
    try:
        j2=cspice.gdpool("BODY%d_J2"%bodyid,0,1)[0]
    except:
        j2=0 #Assume these objects are perfectly spherical
    try:
        re=cspice.gdpool("BODY%d_RADII" % bodyid, 0, 1)[0]
    except:
        re=0  #Assume these bodies are perfect point masses
    return bodytype(spice_id=bodyid,frame=frame,gm=gm,j2=j2,re=re)

def getbodies(bodies,planetframe):
    result=[]
    for i,body in enumerate(bodies):
        result.append(getbody(body,planetframe[i]))
    return result

def j2(planetrel_state,et,body,ssframe=None):
#    sc_state=cspice.spkezr("-189",et,"J2000","NONE","399")[0][0:3] #State of spacecraft relative to Earth
    rotate=cspice.pxform(ssframe,body.frame,et)
    framestate=cspice.mxv(rotate,planetrel_state)
    x=framestate[0]
    y=framestate[1]
    z=framestate[2]
    r=np.sqrt(x**2+y**2+z**2)
    c=-3*body.j2*body.gm*body.re**2/2
    ax=c*x/r**5*(1-5*z**2/r**2)
    ay=c*y/r**5*(1-5*z**2/r**2)
    az=c*z/r**5*(3-5*z**2/r**2)
    frameacc=c*framestate[0:3]/r**5*(np.array([1,1,3])-5*z**2/r**2)
    rotate=cspice.pxform(body.frame,ssframe,et)
    acc=cspice.mxv(rotate,frameacc)
    return acc

def srp(sunrel_state,magnitude_1au=2e-10):
    """
    Generate an acceleration directly away from the Sun, with the given magnitude.
    :param sunrel_state: Position relative to the sun in km
    :param magnitude: Acceleration magnitude in km/s**2
    :return:
    """
    x=sunrel_state[0]
    y=sunrel_state[1]
    z=sunrel_state[2]
    r=np.sqrt(x**2+y**2+z**2)
    return sunrel_state/r*magnitude_1au

def xdot(t=None,x=None,bodies=None,ssframe="J2000",srp_mag=2.95e-11,center=0,do_j2=True):
    acc = np.zeros(3)
    for body in bodies:
        body_state = cspice.spkezr(str(body.spice_id), t, ssframe, "NONE", str(center))[0][0:3]
        if body.spice_id==10:
            sun_state=body_state
        rel = x[0:3]-body_state
        body_acc = -body.gm * rel / np.sqrt(sum(rel ** 2)) ** 3
        acc += body_acc
        if do_j2 and body.j2>0:
            acc+=j2(rel,t,body,ssframe=ssframe)
    #acc += srp(x[0:3]-sun_state,magnitude_1au=2.95e-11)
    return np.array(tuple(x[3:6])+tuple(acc))

def RK4(t=None, x=None, discrete=None, dt=None, xdot=xdot, **kwargs):
    """
    Numerically integrate a differential equation using the Runge-Kutta 4th order method

    :param t: Range time in time units (SI - second)
    :param x: State vector
    :param dt: Time step in same time units as t
    :param kwargs: Passed to xdot. Intended for things like vessel description.
    :return: A tuple
       * New state vector at time t+dt
       * New discrete value(s)
       * extout from first step
    """
    k1 = xdot(t=t     , x=x        ,**kwargs)
    k2 = xdot(t=t+dt/2, x=x+dt*k1/2,**kwargs)
    k3 = xdot(t=t+dt/2, x=x+dt*k2/2,**kwargs)
    k4 = xdot(t=t+dt  , x=x+dt*k3  ,**kwargs)
    return x+dt*(k1+2*k2+2*k3+k4)/6
