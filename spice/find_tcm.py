import spiceypy as cspice
import numpy as np
import copy
import integrate

#Solar system positions
cspice.furnsh("../../Data/spice/generic/spk/planets/de430.bsp")
#Satellite positions
cspice.furnsh("../../Data/spice/generic/spk/satellites/jup230l.bsp")
#Planet constants
cspice.furnsh("../../Data/spice/generic/pck/pck00010.tpc")
cspice.furnsh("../../Data/spice/generic/pck/de-421-masses.tpc")
cspice.furnsh("../../Data/spice/generic/pck/juno.tpc") #Jupiter gravity field
#Leap seconds
cspice.furnsh("../../Data/spice/generic/lsk/naif0011.tls")
#Spacecraft kernel
cspice.furnsh("../../Data/spice/Voyager/spk/voyager_1.ST+1991_a54418u.merged.bsp") #Supertrajectory for Voyager 1
cspice.furnsh("../../Data/spice/Voyager/spk/vgr1_jup230.bsp") #Jupiter encounter trajectory, 1979-02-01/03-19
cspice.furnsh("../../Data/spice/Voyager/spk/voyager_2.ST+1992_m05208u.merged.bsp") #Supertrajectory for Voyager 2
cspice.furnsh("../../Data/spice/Voyager/spk/vgr2_jup230.bsp") #Jupiter encounter trajectory, 1979-02-01/03-19
#First day -- check J2
#et0=cspice.str2et("1977 SEP 5 14:00:00 TDB")
#et1=cspice.str2et("1977 SEP 6 00:00:00 TDB")
#Look for TCM1
#et0=cspice.str2et("1977 SEP 13 00:00:00 TDB")
#et1=cspice.str2et("1977 SEP 14 00:00:00 TDB")
#Look for cruise TCMs
#Look for Voyager 2 TCM during Jupiter Encounter
et0=cspice.str2et("1979 JUL 09 22:29:51 TDB")-12*3600
et1=et0+24*3600

#et1=cspice.str2et("1979 MAR 31 00:00:00 TDB")
bodies=(10,1,2,399,301,4,599,6,7,8,9,501,502,503,504)
planetframe=("IAU_SUN","IAU_MERCURY","IAU_VENUS","J2000","IAU_MOON","IAU_MARS","IAU_JUPITER","IAU_SATURN","IAU_URANUS","IAU_NEPTUNE","IAU_PLUTO","IAU_IO","IAU_EUROPA","IAU_GANYMEDE","IAU_CALLISTO")
spacecraft=-32
ssframe="J2000"
#dt = 20 * 60
dt=1
et0=et0+dt
bodies=integrate.getbodies(bodies,planetframe)

et_list=[]
dacc_list=[]
do_fwdprop = True
import datetime
import matplotlib.dates
j2000=matplotlib.dates.date2num(datetime.datetime(2000,1,1,12,0,0))
earth_state=cspice.spkezr("599",(et1+et0)/2,ssframe,"NONE","399")
print(earth_state)
if do_fwdprop:
    et = et0
    i1=(et1-et0)/dt
    with open("find_tcm.csv", "w") as ouf:
        i = 0
        print("i,et,etcal,ax_o,ay_o,az_o,ax_c,ay_c,az_c,dax,day,daz,dam",file=ouf)
        while et < et1:
            spice_statem=cspice.spkezr(str(spacecraft),et-dt,ssframe,"NONE","0")[0]
            spice_state0=cspice.spkezr(str(spacecraft),et   ,ssframe,"NONE","0")[0]
            spice_statep=cspice.spkezr(str(spacecraft),et+dt,ssframe,"NONE","0")[0]
            acco=(spice_statep[3:6]-spice_state0[3:6])/(dt)
            accc = integrate.xdot(t=et, x=spice_state0, bodies=bodies)[3:6]
            dacc=acco-accc
            et_list.append(et/86400+j2000)
            dacc_list.append(np.linalg.norm(dacc)*1000)
            print("%d,%50.20f,%s,%20.16e,%20.16e,%20.16e,%20.16e,%20.16e,%20.16e,%20.16e,%20.16e,%20.16e,%20.16e" % (
                    i,     et,cspice.etcal(et),
                                 acco[0],acco[1],acco[2],
                                                         accc[0],accc[1],accc[2],
                                                                                 dacc[0],dacc[1],dacc[2],np.linalg.norm(dacc)*1000),file=ouf)

            et += dt
            i += 1
            if i % 1000 == 0:
                print(i,i1)
    import matplotlib.pyplot as plt
    (fig,ax)=plt.subplots(1,1)
    plt.plot_date(et_list,dacc_list,'b-')
    plt.ylabel('acc/(m/s**2)')
    plt.show()
