import spiceypy as cspice
import numpy as np
import matplotlib.pyplot as plt
import integrate

#Solar system positions
cspice.furnsh("../../Data/spice/generic/spk/planets/de430.bsp")
#Satellite positions
cspice.furnsh("../../Data/spice/generic/spk/satellites/jup230l.bsp")
#Planet constants
cspice.furnsh("../../Data/spice/generic/pck/pck00010.tpc")
cspice.furnsh("../../Data/spice/generic/pck/gm_de431.tpc")
cspice.furnsh("../../Data/spice/generic/pck/juno.tpc") #Jupiter gravity field
#Leap seconds
cspice.furnsh("../../Data/spice/generic/lsk/naif0012.tls")
#Spacecraft kernel
#Reference kernel
infn0="../../Data/spice/Voyager/spk/voyager_1.ST+1991_a54418u.merged.bsp" #Supertrajectory for Voyager 1
#Comparison kernel
infn1="../../Data/spice/Voyager/spk/vgr1_jup230.bsp" #Jupiter encounter trajectory, 1979-02-01/03-19

#Get coverage intervals for both kernels and compute intersection
spacecraft=-31
cov0=cspice.stypes.SPICEDOUBLE_CELL(100)
cspice.spkcov(infn0,spacecraft,cov0)
cov1=cspice.stypes.SPICEDOUBLE_CELL(100)
cspice.spkcov(infn1,spacecraft,cov1)

cov=cspice.wnintd(cov0,cov1)
et0=cspice.wnfetd(cov,0)[0]
et1=cspice.wnfetd(cov,cspice.wncard(cov)-1)[1]
#et0=cspice.str2et("1979 MAR 05 12:00:00 TDB")
#et1=cspice.str2et("1979 MAR 05 12:10:00 TDB")
print(cspice.etcal(et0),cspice.etcal(et1))

center=0
#ndt=10000
#dt=(et1-et0)/ndt
dt=60
et0+=dt
et1-=dt
ndt=int((et1-et0)/dt)
print(ndt)
#dt=(et1-et0)/ndt
dtacc=1
pos0=np.zeros((ndt,9))
jstate=np.zeros((ndt,9))
pos1=np.zeros((ndt,9))
bodies=(10,1,2,399,301,4,599,6,7,8,9,501,502,503,504)
planetframe=("IAU_SUN","IAU_MERCURY","IAU_VENUS","J2000","IAU_MOON","IAU_MARS","IAU_JUPITER","IAU_SATURN","IAU_URANUS","IAU_NEPTUNE","IAU_PLUTO","IAU_IO","IAU_EUROPA","IAU_GANYMEDE","IAU_CALLISTO")
spacecraft=-31
bodies=integrate.getbodies(bodies,planetframe)

import datetime
import matplotlib.dates
j2000=matplotlib.dates.date2num(datetime.datetime(2000,1,1,12,0,0))
ssframe="ECLIPB1950"
et=np.zeros(ndt)
cspice.furnsh(infn0)
accj2=np.zeros((ndt,3))
accj0=np.zeros((ndt,3))
for i in range(ndt):
    et[i]=et0+(et1-et0)*i/ndt
    spice_statem=cspice.spkezr(str(spacecraft),et[i]-dtacc,ssframe,"NONE",str(center))[0]
    spice_state0=cspice.spkezr(str(spacecraft),et[i]      ,ssframe,"NONE",str(center))[0]
    spice_statep=cspice.spkezr(str(spacecraft),et[i]+dtacc,ssframe,"NONE",str(center))[0]
    jstatem=cspice.spkezr("599",et[i]-dtacc,ssframe,"NONE",str(center))[0]
    jstate0=cspice.spkezr("599",et[i]      ,ssframe,"NONE",str(center))[0]
    jstatep=cspice.spkezr("599",et[i]-dtacc,ssframe,"NONE",str(center))[0]
    jstate[i,0:6]=jstate0
    jstate[i,6:9]=(jstatep[3:6]-jstatem[3:6])/(2*dtacc)
    accj2[i,:]=integrate.xdot(t=et[i],x=spice_state0,bodies=bodies,ssframe=ssframe,center=center,do_j2=True )[3:6]
    accj0[i,:]=integrate.xdot(t=et[i],x=spice_state0,bodies=bodies,ssframe=ssframe,center=center,do_j2=False)[3:6]
    pos0[i,0:6]=spice_state0
    pos0[i,6:9]=(spice_statep[3:6]-spice_statem[3:6])/(2*dtacc)
    if i%1000==0:
        print(i,cspice.etcal(et[i]))
    i+=1
cspice.unload(infn0)

#cspice.furnsh(infn1)
#for i in range(ndt):
#    et[i]=et0+(et1-et0)*i/ndt
#    spice_statem=cspice.spkezr(str(spacecraft),et[i]-dt,ssframe,"NONE",str(center))[0]
#    spice_state0=cspice.spkezr(str(spacecraft),et[i]   ,ssframe,"NONE",str(center))[0]
#    spice_statep=cspice.spkezr(str(spacecraft),et[i]+dt,ssframe,"NONE",str(center))[0]
#    pos1[i,0:6]=spice_state0
#    pos1[i,6:9]=(spice_statep[3:6]-spice_statem[3:6])/(2*dt)
#    i+=1
#cspice.unload(infn1)

scale=1e12
#red=pos1
#dpos=pos1-pos0
#red=pos0+dpos*scale
tick=ndt//100
plt.figure(1)
#plt.plot(pos0[:,3],pos0[:,4],'b',red[:,3],red[:,4],'r')
#plt.plot(pos0[::tick,3],pos0[::tick,4],'b*',red[::tick,3],red[::tick,4],'r*')
plt.plot(pos0[:,0]-jstate[:,0],pos0[:,1]-jstate[:,1],'b')
plt.plot(pos0[::tick,0]-jstate[::tick,0],pos0[::tick,1]-jstate[::tick,1],'b*')
vtick=ndt//100
for i in range(0,len(et)):
    daccj2=pos0[i,6:9]-accj2[i,:]
    if np.linalg.norm(daccj2)>1e-8:
         plt.plot((pos0[i,0]-jstate[i,0],pos0[i,0]-jstate[i,0]+daccj2[0]*scale),(pos0[i,1]-jstate[i,1],pos0[i,1]-jstate[i,1]+daccj2[1]*scale),'r-')
plt.axis('equal')
#plt.figure(2)
#plt.plot_date(et/86400+j2000,pos[:,3],'b-')
plt.figure(3)
plt.plot_date(et/86400+j2000,np.linalg.norm(pos0[:,6:9]-accj2,axis=1)*1000,'r-')
#plt.plot_date(et/86400+j2000,(pos0[:,6]-accj0[:,0])*1000 ,'g*')
plt.ylabel('dacc/(m/s**2)')
plt.show()
pass