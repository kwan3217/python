import spiceypy as cspice
import copy
import integrate

cspice.furnsh("../../Data/spice/generic/spk/planets/de430s.bsp")
cspice.furnsh("../../Data/spice/generic/spk/satellites/mar097s.bsp")
cspice.furnsh("../../Data/spice/generic/lsk/naif0012.tls")
cspice.furnsh("../../Data/spice/insight/nsyt_spk_cruise_od013_v1.bsp")
cspice.furnsh("../../Data/spice/generic/pck/pck00010.tpc")
cspice.furnsh("../../Data/spice/generic/pck/gm_de431.tpc")
#cspice.furnsh("../../Data/spice/generic/pck/juno.tpc")
et0=cspice.str2et("2018 MAY 05 21:09:46 TDB")
bodies=(10,1,2,399,301,4,5,6,7,8,9)
frames=("","","","J2000","","","","","","","")
bodies=integrate.getbodies(bodies,frames)

n=3600

#with open("nsyt2.csv","w") as ouf:
#    for i in range(n):
#        et=et0+i*(et1-et0)/n
#        sc_state = cspice.spkezr("-189", et, "J2000", "NONE", "0")[0]
#        acc=xdot(et,sc_state)[3:6]
#        print("%04d,%50.20f,%20.16e,%20.16e,%20.16e,%20.16e,%20.16e,%20.16e,%20.16e,%20.16e,%20.16e"%((i,et)+tuple(sc_state)+tuple(acc)),file=ouf)

#Actually do the back integration

x0=cspice.spkezr("-189", et0, "J2000", "NONE", "0")[0]

#Backpropagation
do_backprop=True
if do_backprop:
    t0=cspice.str2et("2018 MAY 05 11:05:00 UTC") #Best available value of T0 for launch
    print("tlaunch: ",cspice.etcal(t0))
    tsep=t0+5600
    print("tsep:    ",cspice.etcal(tsep))
    x=copy.copy(x0)
    et=et0
    dt=-1

    with open("nsyt_backprop.csv","w") as ouf:
        i=0
        while et>tsep:
            earth_state = cspice.spkezr("399", et, "J2000", "NONE", "0")[0]
            acc=integrate.xdot(t=et,x=x,bodies=bodies)[3:6]
            print("%d,%50.20f,%20.16e,%20.16e,%20.16e,%20.16e,%20.16e,%20.16e,%20.16e,%20.16e,%20.16e"%((i,et)+tuple(x-earth_state)+tuple(acc)),file=ouf)
            x=integrate.RK4(et,x,dt=dt,bodies=bodies)
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
            spice_state0=cspice.spkezr("-189",et,"J2000","NONE","0")[0]
            spice_statep=cspice.spkezr("-189",et+dt,"J2000","NONE","0")[0]
            acc = integrate.xdot(t=et, x=x,bodies=bodies)[3:6]
            acc_spice=(spice_statep-spice_state0)[3:6]/dt
            print("%d,%50.20f,%20.16e,%20.16e,%20.16e,%20.16e,%20.16e,%20.16e,%20.16e,%20.16e,%20.16e,%20.16e,%20.16e,%20.16e,%20.16e,%20.16e,%20.16e" % (
            (i, et) + tuple(x - earth_state)+ tuple(spice_state0 - earth_state) + tuple(acc-acc_spice)), file=ouf)
            x = integrate.RK4(et, x, dt=dt,bodies=bodies)
            et += dt
            i += 1
            if i % 1000 == 0:
                print(i)

