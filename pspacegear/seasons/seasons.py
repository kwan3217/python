import spiceypy as cspice
import numpy as np
from scipy.optimize import minimize_scalar, bisect

cspice.furnsh('naif0012.tls')
cspice.furnsh('pck00010.tpc')
cspice.furnsh('de432s.bsp')
cspice.furnsh('earth_720101_070426.bpc')
cspice.furnsh('earth_000101_190310_181217.bpc')

et0=0
frame=""
abcor=""
def dec(dt):
    global et0,frame,abcor
    (pos,lt)=cspice.spkezr("10",et0+dt,frame,abcor,"399")
    return np.degrees(np.arctan2(pos[2],np.sqrt(pos[0]**2+pos[1]**2)))

def ls(dt,ofs):
    """

    :param dt: time offset from et0 in seconds
    :param ofs: Number of degrees to subtract from result, for root finder. So
                to find the time for an Ls of 270deg, pass 270 and this will
                return 0 when dt is at the solstice.
    :return: difference between solar longitude and ofs in degrees
    """
    global et0,frame,abcor
    return np.degrees(cspice.lspcn("399",dt+et0,abcor))-ofs

def own_ls(dt,ofs):
    global et0,frame,abcor
    et=dt+et0
    (state,lt)=cspice.spkezr("10",et,frame,abcor,"399")
    pos=np.array(state[:3])
    vel=np.array(state[3:])


et0=cspice.str2et("2018-12-15 01:00:00 UTC")
dt1=cspice.str2et("2018-12-20 00:00:00 UTC")-et0
dt2=cspice.str2et("2018-12-25 00:00:00 UTC")-et0

frames=["ITRF93","IAU_EARTH","J2000"]
abcors=["NONE","LT","LT+S"]
xtols=[1e-10]

print("Minimize declination, use DE432")
for frame in frames:
    for abcor in abcors:
        for xtol in xtols:
            res=minimize_scalar(dec,(0,dt1,dt2),options={"xtol":xtol})
            dt=res.x
            print("%-12s %-12s %s TDT %s UTC %f %f"%(frame,abcor,cspice.etcal(dt+et0),cspice.timout(dt+et0,"YYYY-MM-DD HR:MN:SC.### ::UTC"),dec(dt),ls(dt,0)))

cspice.kclear()
cspice.furnsh('naif0012.tls')
cspice.furnsh('pck00010.tpc')
cspice.furnsh('de200.bsp')
cspice.furnsh('earth_720101_070426.bpc')
cspice.furnsh('earth_000101_190310_181217.bpc')

print("Minimize declination, use DE200")
for frame in frames:
    for abcor in abcors:
        for xtol in xtols:
            res=minimize_scalar(dec,(0,dt1,dt2),options={"xtol":xtol})
            dt=res.x
            print("%-12s %-12s %s TDT %s UTC %f %f"%(frame,abcor,cspice.etcal(dt+et0),cspice.timout(dt+et0,"YYYY-MM-DD HR:MN:SC.### ::UTC"),dec(dt),ls(dt,0)))

cspice.kclear()
cspice.furnsh('naif0012.tls')
cspice.furnsh('pck00010.tpc')
cspice.furnsh('de118.bsp')
cspice.furnsh('earth_720101_070426.bpc')
cspice.furnsh('earth_000101_190310_181217.bpc')

print("Minimize declination, use DE118")
for frame in frames:
    for abcor in abcors:
        for xtol in xtols:
            res=minimize_scalar(dec,(0,dt1,dt2),options={"xtol":xtol})
            dt=res.x
            print("%-12s %-12s %s TDT %s UTC %f %f"%(frame,abcor,cspice.etcal(dt+et0),cspice.timout(dt+et0,"YYYY-MM-DD HR:MN:SC.### ::UTC"),dec(dt),ls(dt,0)))

import collections
lowprecResult=collections.namedtuple("lowprecResult",["lam","alp","delta"])

def lowprec(dt):
    global et0
    et=dt+et0
    n=et/86400.0 #Low precision formula is so low precision that it doesn't care about et/utc
    L0=280.460  #deg
    L1=0.9856474 #deg/day
    L=L0+L1*n   #Mean solar longitude
    while(L<0):
        L+=360
    while(L>360):
        L-=360
    g0=357.528 #deg
    g1=0.9856003 #deg/day
    g=g0+g1*n #Mean anomaly of Sun
    while(g<0):
        g+=360
    while(g>360):
        g-=360
    lam1=1.915 #deg
    lam2=0.020 #deg
    lam=L+lam1*np.sin(np.radians(g))+lam2*np.sin(np.radians(g)) #Ecliptic longitude of Sun, equal to Ls
    bet=0
    eps0=23.429 #deg
    eps1=-0.0000004 #deg/day
    eps=eps0+eps1*n
    alp=np.degrees(np.arctan(np.cos(np.radians(eps))*np.tan(np.radians(lam))))
    delta=np.degrees(np.arcsin(np.sin(np.radians(eps))*np.sin(np.radians(lam))))
    R0=1.00014
    R1=-0.01671
    R2=-0.00014
    R=R0+R1*np.cos(np.radians(g))+R2*np.cos(2*np.radians(g))
    return lowprecResult(lam=lam,alp=alp,delta=delta)

def lowprec_dec(dt):
    return lowprec(dt).delta

def lowprec_ls(dt):
    return lowprec(dt).lam


res = minimize_scalar(lowprec_dec, (0, dt1, dt2), options={"xtol": xtol})
dt = res.x
print("%-12s %-12s %s TDT %s UTC %f %f" % ("lowprec", "lowprec", cspice.etcal(dt + et0), cspice.timout(dt + et0, "YYYY-MM-DD HR:MN:SC.### ::UTC"), lowprec_dec(dt), lowprec_ls(dt)))

cspice.kclear()
cspice.furnsh('naif0012.tls')
cspice.furnsh('pck00010.tpc')
cspice.furnsh('de432s.bsp')
cspice.furnsh('earth_720101_070426.bpc')
cspice.furnsh('earth_000101_190310_181217.bpc')

print("Aim for Ls=270deg, use DE432")
for abcor in abcors:
    for xtol in xtols:
        dt=bisect(ls,dt1,dt2,args=270,)
        print("%-12s %-12s %s TDT %s UTC %f %f" % ("ls", abcor, cspice.etcal(dt + et0), cspice.timout(dt + et0, "YYYY-MM-DD HR:MN:SC.### ::UTC"), dec(dt),ls(dt,0)))

#Try to duplicate example 24.a using DE432
#Time of example
et_meeus=cspice.str2et("JD 2448908.5 TDB") #Exactly 1992-Oct-13 00:00:00 TDB
#Position of Earth-Moon Barycenter in ECLIPJ2000 frame, used to establish dynamical ecliptic
(state,_)=cspice.spkezr("3",et_meeus,"ECLIPJ2000","NONE","10")
#Pole of dynamical ecliptic of date in ECLIPJ2000. Should be within a small fraction of a degree of z-axis of ECLIPJ2000
kecl_pole=np.cross(state[:3],state[3:])
ecl_pole=kecl_pole/np.linalg.norm(kecl_pole)
print(ecl_pole)
print(1-np.linalg.norm(ecl_pole))
print(np.degrees(np.arccos(np.dot(ecl_pole,np.array([0,0,1])))))
#Earth pole in ECLIPJ2000, using best available model, so accounting for precession, nutation, polar motion, DeltaT, etc.
#Should be about 23.5deg away from z-axis of ECLIPJ2000
M=cspice.pxform("ITRF93","ECLIPJ2000",et_meeus)
equ_pole=M @ np.array([0,0,1])
print(equ_pole)
print(1-np.linalg.norm(equ_pole))
print(np.degrees(np.arccos(np.dot(equ_pole,np.array([0,0,1])))))
#Direction of dynamical vernal equinox in ECLIPJ2000, found by cross product of equatorial pole with ecliptic pole.
#Should be near x-axis of ECLIPJ2000
kequinox=np.cross(equ_pole,ecl_pole)
equinox=kequinox/np.linalg.norm(kequinox)
print(equinox)
print(1-np.linalg.norm(equinox))
print(np.degrees(np.arccos(np.dot(equinox,np.array([1,0,0])))))
solstice=np.cross(ecl_pole,equinox)
print(solstice)
print(1-np.linalg.norm(solstice))
print(np.degrees(np.arccos(np.dot(solstice,np.array([0,1,0])))))

M=np.hstack([equinox,solstice,ecl_pole])

print("Done")