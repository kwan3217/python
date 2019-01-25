import spiceypy as cspice
import bmw
import numpy as np
import os
from numpy.linalg import norm
from numpy import dot

def vangle(a,b):
    return dot(a,b)/(norm(a)*norm(b))

utc0_v1='1977-09-05 12:56:00 UTC' #Actual launch time for Voyager 1, from Wikipedia
utc1_v1='1979-03-05 12:05:26 TDB' #Jupiter-centered hyperbolic element epoch, at which mean anomaly is zero (IE periapse)
utc2_v1='1980-11-12 23:46:30 TDB' #Saturn-centered hyperbolic element epoch, at which mean anomaly is zero (IE periapse)
utc0_v2='1977-08-20 14:29:00 UTC' #Actual launch time for Voyager 2, from Wikipedia
utc1_v2='1979-07-09 22:29:51 TDB' #Jupiter-centered hyperbolic element epoch, at which mean anomaly is zero (IE periapse)
utc2_v2='1981-08-26 03:24:57 TDB' #Saturn-centered hyperbolic element epoch, at which mean anomaly is zero (IE periapse)
utc3_v2='1986-01-24 17:59:47 TDB' #Uranus-centered hyperbolic element epoch, at which mean anomaly is zero (IE periapse)
utc4_v2='1989-08-25 03:56:36 TDB' #Neptune-centered hyperbolic element epoch, at which mean anomaly is zero (IE periapse)
os.chdir("../../Data/spice/generic/")
cspice.furnsh("generic.tm")
mu = cspice.gdpool("BODY10_GM", 0, 1)[0]
au = cspice.gdpool("AU", 0, 1)[0]
et0_v1 = cspice.str2et(utc0_v1)
et1_v1 = cspice.str2et(utc1_v1)
et2_v1 = cspice.str2et(utc2_v1)
et0_v2 = cspice.str2et(utc0_v2)
et1_v2 = cspice.str2et(utc1_v2)
et2_v2 = cspice.str2et(utc2_v2)
et3_v2 = cspice.str2et(utc3_v2)
et4_v2 = cspice.str2et(utc4_v2)
frame="ECLIPJ2000"
(rv0_v1, _) = cspice.spkezr("399", et0_v1, frame, "NONE", "SUN")
(rv1_v1, _) = cspice.spkezr("5", et1_v1, frame, "NONE", "SUN")
(rv2_v1, _) = cspice.spkezr("6", et2_v1, frame, "NONE", "SUN")
(rv0_v2, _) = cspice.spkezr("399", et0_v2, frame, "NONE", "SUN")
(rv1_v2, _) = cspice.spkezr("5", et1_v2, frame, "NONE", "SUN")
(rv2_v2, _) = cspice.spkezr("6", et2_v2, frame, "NONE", "SUN")
(rv3_v2, _) = cspice.spkezr("7", et3_v2, frame, "NONE", "SUN")
(rv4_v2, _) = cspice.spkezr("8", et4_v2, frame, "NONE", "SUN")
(dv0_v1, av1_v1) = bmw.gauss(rv0_v1[0:3], rv1_v1[0:3], et1_v1 - et0_v1, l_DU=au, mu=mu, Type=-1)
(dv1_v1, av2_v1) = bmw.gauss(rv1_v1[0:3], rv2_v1[0:3], et2_v1 - et1_v1, l_DU=au, mu=mu, Type=-1)
#Heliocentric orbital elements
el01_v1=bmw.elorb(rv0_v1[0:3],dv0_v1,l_DU=au,mu=mu,t0=et0_v1)
el12_v1=bmw.elorb(rv1_v1[0:3],dv1_v1,l_DU=au,mu=mu,t0=et1_v1)
print("V1 EJ orbit: ",el01_v1)
print("V1 JS orbit: ",el12_v1)
dvinf_v1=dv0_v1-rv0_v1[3:6]
rla_v1=np.degrees(np.arctan2(dvinf_v1[1],dvinf_v1[0]))
rve_v1=np.degrees(np.arctan2(rv0_v1[4],rv0_v1[3]))
dla_v1=np.degrees(np.arcsin(dvinf_v1[2]/norm(dvinf_v1)))
c3_v1=norm(dvinf_v1)**2
print("v1 rla (deg): %f, dla (deg): %f, c3 (km**2/s**2): %f"%(rla_v1,dla_v1,c3_v1))
print("v1 rle (deg): %f, launch heading (deg, positive is inward): %f"%(rve_v1,rla_v1-rve_v1))
avinf_v1=av1_v1-rv1_v1[3:6]
print(avinf_v1,norm(avinf_v1))
(dv0_v2, av1_v2) = bmw.gauss(rv0_v2[0:3], rv1_v2[0:3], et1_v2 - et0_v2, l_DU=au, mu=mu, Type=-1)
(dv1_v2, av2_v2) = bmw.gauss(rv1_v2[0:3], rv2_v2[0:3], et2_v2 - et1_v2, l_DU=au, mu=mu, Type=-1)
(dv2_v2, av3_v2) = bmw.gauss(rv2_v2[0:3], rv3_v2[0:3], et3_v2 - et2_v2, l_DU=au, mu=mu, Type=-1)
(dv3_v2, av4_v2) = bmw.gauss(rv3_v2[0:3], rv4_v2[0:3], et4_v2 - et3_v2, l_DU=au, mu=mu, Type=-1)
el_v2=bmw.elorb(rv0_v2[0:3],dv0_v2,l_DU=au,mu=mu,t0=et0_v2)
print(el_v2)
print(rv0_v2[0:3]/au,norm(rv0_v2[0:3])/au)
print(rv1_v2[0:3]/au,norm(rv1_v2[0:3])/au)
dvinf_v2=dv0_v2-rv0_v2[3:6]
rla_v2=np.degrees(np.arctan2(dvinf_v2[1],dvinf_v2[0]))
rve_v2=np.degrees(np.arctan2(rv0_v2[4],rv0_v2[3]))
dla_v2=np.degrees(np.arcsin(dvinf_v2[2]/norm(dvinf_v2)))
c3_v2=norm(dvinf_v2)**2
print("v2 rla (deg): %f, dla (deg): %f, c3 (km**2/s**2): %f"%(rla_v2,dla_v2,c3_v2))
print("v1 rle (deg): %f, launch heading (deg, positive is inward): %f"%(rve_v2,rla_v2-rve_v2))
print(dvinf_v1,norm(dvinf_v1),np.degrees(vangle(dvinf_v2[0:2],rv0_v2[3:5])))
avinf_v2=av1_v2-rv1_v2[3:6]
print(avinf_v2,norm(avinf_v2))

earth=[]
jupiter=[]
saturn=[]
uranus=[]
neptune=[]
v1=[]
v1_01=[]
v1_12=[]
v2=[]
v2_01=[]
v2_12=[]
v2_23=[]
v2_34=[]
t=et0_v2
while t<et4_v2:
    earth.append(cspice.spkezr("399",t,frame,"NONE","SUN")[0][0:3])
    jupiter.append(cspice.spkezr("5",t,frame,"NONE","SUN")[0][0:3])
    saturn.append(cspice.spkezr("6",t,frame,"NONE","SUN")[0][0:3])
    uranus.append(cspice.spkezr("7",t,frame,"NONE","SUN")[0][0:3])
    neptune.append(cspice.spkezr("8",t,frame,"NONE","SUN")[0][0:3])
    v1_01.append(bmw.kepler(rv0_v1[0:3], dv0_v1, t - et0_v1, l_DU=au, mu=mu)[0])
    if t>=et0_v1 and t<et1_v1:
        v1.append(v1_01[-1])
    v1_12.append(bmw.kepler(rv1_v1[0:3],dv1_v1,t-et1_v1,l_DU=au,mu=mu)[0])
    if t>=et1_v1 and t<et2_v1:
        v1.append(v1_12[-1])

    v2_01.append(bmw.kepler(rv0_v2[0:3], dv0_v2, t - et0_v2, l_DU=au, mu=mu)[0])
    if t>=et0_v2 and t<et1_v2:
        v2.append(v2_01[-1])
    v2_12.append(bmw.kepler(rv1_v2[0:3], dv1_v2, t - et1_v2, l_DU=au, mu=mu)[0])
    if t>=et1_v2 and t<et2_v2:
        v2.append(v2_12[-1])
    v2_23.append(bmw.kepler(rv2_v2[0:3], dv2_v2, t - et2_v2, l_DU=au, mu=mu)[0])
    if t>=et2_v2 and t<et3_v2:
        v2.append(v2_23[-1])
    v2_34.append(bmw.kepler(rv3_v2[0:3], dv3_v2, t - et3_v2, l_DU=au, mu=mu)[0])
    if t>=et3_v2 and t<et4_v2:
        v2.append(v2_34[-1])
    t+=86400
earth=np.array(earth)/au
jupiter=np.array(jupiter)/au
saturn=np.array(saturn)/au
uranus=np.array(uranus)/au
neptune=np.array(neptune)/au
v1=np.array(v1)/au
v2=np.array(v2)/au
v1_01=np.array(v1_01)/au
v1_12=np.array(v1_12)/au
v2_01=np.array(v2_01)/au
v2_12=np.array(v2_12)/au
v2_23=np.array(v2_23)/au
v2_34=np.array(v2_34)/au
import matplotlib.pyplot as plt
plt.plot(earth[:365,0],earth[:365,1],'w',
         jupiter[:,0],jupiter[:,1],'w',
         saturn[:,0],saturn[:,1],'w',
         uranus[:,0],uranus[:,1],'w',
         neptune[:,0],neptune[:,1],'w',
         )
interval=30
plt.plot(jupiter[::interval,0],jupiter[::interval,1],'w*',
         saturn[::interval,0],saturn[::interval,1],'w*',
         uranus[::interval,0], uranus[::interval,1], 'w*',
         neptune[::interval,0], neptune[::interval,1], 'w*',
         )
plt.plot([0],[0],'yo',markersize=5)
plt.plot(v1_01[:,0],v1_01[:,1],'#800000')
plt.plot(v1_01[interval:600:interval,0],v1_01[interval:600:interval,1],color='#800000',marker='*',linestyle='None')
plt.plot(v1_12[:,0],v1_12[:,1],'#800000')
plt.plot(v1_12[::interval,0],v1_12[::interval,1],color='#800000',marker='*',linestyle='None')
plt.plot(v2_01[:,0],v2_01[:,1],'#000080')
plt.plot(v2_01[0:800:interval,0],v2_01[0:800:interval,1],color='#000080',marker='*',linestyle='None')
plt.plot(v2_12[:,0],v2_12[:,1],'#000080')
plt.plot(v2_12[::interval,0],v2_12[::interval,1],color='#000080',marker='*',linestyle='None')
plt.plot(v2_23[:,0],v2_23[:,1],'#000080')
plt.plot(v2_23[::interval,0],v2_23[::interval,1],color='#000080',marker='*',linestyle='None')
plt.plot(v2_34[:,0],v2_34[:,1],'#000080')
plt.plot(v2_34[::interval,0],v2_34[::interval,1],color='#000080',marker='*',linestyle='None')
plt.plot(v1[:,0],v1[:,1],'r',v2[:,0],v2[:,1],'b')
plt.xlabel('ECLIPJ2000 X/au')
plt.ylabel('ECLIPJ2000 Y/au')
plt.axis('equal')
ax = plt.gca()
ax.set_facecolor((0, 0, 0))
plt.figure(2)
plt.plot(norm(v2_01[0:1164,0:2]-v1_01[0:1164,0:2],axis=1))
plt.show()