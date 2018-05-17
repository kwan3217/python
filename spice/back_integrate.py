import spiceypy as cspice
import daf
import numpy as np
import copy

cspice.furnsh("../../Data/spice/generic/spk/planets/de430s.bsp")
cspice.furnsh("../../Data/spice/generic/spk/satellites/mar097s.bsp")
cspice.furnsh("../../Data/spice/generic/lsk/naif0012.tls")
cspice.furnsh("../../Data/spice/insight/nsyt_spk_cruise_od013_v1.bsp")
cspice.furnsh("../../Data/spice/generic/pck/pck00010.tpc")
cspice.furnsh("../../Data/spice/generic/pck/gm_de431.tpc")
et0=cspice.str2et("2018 MAY 05 21:09:46 TDB")
bodies=(10,1,2,399,301,4,5,6,7,8,9)
gm=[]
for body in bodies:
    gm.append(cspice.gdpool("BODY%d_GM"%body,0,1)[0])
gm=tuple(gm)

def j2(earthrel_state):
#    sc_state=cspice.spkezr("-189",et,"J2000","NONE","399")[0][0:3] #State of spacecraft relative to Earth
    x=earthrel_state[0]
    y=earthrel_state[1]
    z=earthrel_state[2]
    r=np.sqrt(x**2+y**2+z**2)
    re=6.3781363000000000E+03
    j2=1.0826254500000000E-03
    mu=398600.435436
    c=j2*mu*re**2/2
    return 3*c*earthrel_state[0:3]/r**5*(np.array([1,1,3])-5*z**2/r**2)

def srp(sunrel_state):
#    sc_state=cspice.spkezr("-189",et,"J2000","NONE","10")[0][0:3] #State of spacecraft relative to Earth
    x=sunrel_state[0]
    y=sunrel_state[1]
    z=sunrel_state[2]
    r=np.sqrt(x**2+y**2+z**2)
    return sunrel_state/r*(2e-10)

et1=et0+1
n=3600
def xdot(t=None,x=None):
    acc = np.zeros(3)
    for j, body in enumerate(bodies):
        body_state = cspice.spkezr(str(body), t, "J2000", "NONE", "0")[0][0:3]
        if body==10:
            sun_state=body_state
        if body==399:
            earth_state=body_state
        rel = body_state - x[0:3]
        body_acc = gm[j] * rel / np.sqrt(sum(rel ** 2)) ** 3
        acc += body_acc
    acc += j2(x[0:3]-earth_state)
    acc += srp(x[0:3]-sun_state)
    return np.array(tuple(x[3:6])+tuple(acc))

#with open("nsyt2.csv","w") as ouf:
#    for i in range(n):
#        et=et0+i*(et1-et0)/n
#        sc_state = cspice.spkezr("-189", et, "J2000", "NONE", "0")[0]
#        acc=xdot(et,sc_state)[3:6]
#        print("%04d,%50.20f,%20.16e,%20.16e,%20.16e,%20.16e,%20.16e,%20.16e,%20.16e,%20.16e,%20.16e"%((i,et)+tuple(sc_state)+tuple(acc)),file=ouf)

#Actually do the back integration

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
    k1 = xdot(t=t     , x=x        )
    k2 = xdot(t=t+dt/2, x=x+dt*k1/2)
    k3 = xdot(t=t+dt/2, x=x+dt*k2/2)
    k4 = xdot(t=t+dt  , x=x+dt*k3  )
    return x+dt*(k1+2*k2+2*k3+k4)/6

x0=cspice.spkezr("-189", et0, "J2000", "NONE", "0")[0]

#Backpropagation
do_backprop=False
if do_backprop:
    t0=cspice.str2et("2018 MAY 05 11:05:00 UTC") #Best available value of T0 for launch
    print(cspice.etcal(t0))
    tsep=t0+5600
    print(cspice.etcal(tsep))
    x=copy.copy(x0)
    et=et0
    dt=-1

    with open("nsyt_backprop.csv","w") as ouf:
        i=0
        while et>tsep:
            earth_state = cspice.spkezr("399", et, "J2000", "NONE", "0")[0]
            acc=xdot(t=et,x=x)[3:6]
            print("%d,%50.20f,%20.16e,%20.16e,%20.16e,%20.16e,%20.16e,%20.16e,%20.16e,%20.16e,%20.16e"%((i,et)+tuple(x-earth_state)+tuple(acc)),file=ouf)
            x=RK4(et,x,dt=dt)
            et+=dt
            i+=1
            if i%1000==0:
                print(i)

# Fwd propagation, to check accuracy and effect of rcs
do_fwdprop = True
if do_fwdprop:
    x = copy.copy(x0)
    et = et0

    dt=1
    et2=et0+30000
    with open("nsyt_fwdprop.csv", "w") as ouf:
        i = 0
        print("i,et,x,y,z,xd,yd,zd,x_spice,y_spice,z_spice,xd_spice,yd_spice,zd_spice,ax,ay,az",file=ouf)
        while et < et2:
            earth_state = cspice.spkezr("399", et, "J2000", "NONE", "0")[0]
            spice_state=cspice.spkezr("-189",et,"J2000","NONE","0")[0]
            acc = xdot(t=et, x=x)[3:6]
            print("%d,%50.20f,%20.16e,%20.16e,%20.16e,%20.16e,%20.16e,%20.16e,%20.16e,%20.16e,%20.16e,%20.16e,%20.16e,%20.16e,%20.16e,%20.16e,%20.16e" % (
            (i, et) + tuple(x - earth_state)+ tuple(spice_state - earth_state) + tuple(acc)), file=ouf)
            x = RK4(et, x, dt=dt)
            et += dt
            i += 1
            if i % 1000 == 0:
                print(i)

