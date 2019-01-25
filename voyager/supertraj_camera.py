import spiceypy as cspice
import numpy as np
import matplotlib.pyplot as plt

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
cspice.furnsh("../../Data/spice/Voyager/spk/voyager_2.ST+1992_m05208u.merged.bsp") #Supertrajectory for Voyager 1

#From original movie:
#Frame 0-93 is title sequence
#Frame 94 is high above Solar system (about 115AU), immediately begins move down
#Frame 504, latitude stops changing
#Frame 508, distance stops changing
#Frame 708, Voyager 2 Launch 1977-08-20T14:29:45
#Frame 712, Voyager 1 Launch 1977-09-05T12:56
#Frame 862, Voyager 1 at Jupiter, 1979-03-05
#Frame 896, Voyager 2 at Jupiter, 1979-07-09T22:29:00
#Frame 1032, Voyager 1 at Saturn, 1980-11-12
#Frame 1111, Voyager 2 at Saturn, 1981-08-25T03:24:05
#Frame 1553, Voyager 2 at Uranus, 1986-01-24T17:59:47
#Frame 1910, Voyager 2 at Neptune, 1989-08-25T03:56:36
#Frame 1983, start moving back up
#Frame 2171, end of movie

def linterp(x0,y0,x1,y1,x):
    t=(x-x0)/(x1-x0)
    return y0*(1-t)+y1*t

frame0=94
frame1=708
frame2=1910
frame3=2171
et1=cspice.str2et("1977-08-20 14:29:45 UTC")
et_v2=cspice.str2et("1977-08-21 00:00:00 TDB")
et2=cspice.str2et("1989-08-25 03:56:36 TDB")
et0=linterp(frame1,et1,frame2,et2,frame0)
et3=linterp(frame1,et1,frame2,et2,frame3)
frames=np.arange(0,frame3)
ets=linterp(frame1,et1,frame2,et2,frames)
v2_headings=ets*float('nan')
cam_headings=ets*0-35.56
for i,et in enumerate(ets):
    if et>et_v2:
        sc_state=cspice.spkezr("-32",et,"ECLIPB1950","NONE","0")[0]
        v2_headings[i]=np.degrees(np.arctan2(-sc_state[4],-sc_state[3]))

piecewise_x=[708]
piecewise_y=[-35.56]

cam_headings[708:866]=linterp(708,-35.56,866,v2_headings[866],frames[708:866])
f0=866
piecewise_x.append(f0)
piecewise_y.append(v2_headings[f0])
f1=937
cam_headings[f0:f1]=linterp(f0,v2_headings[f0],f1,v2_headings[f1],frames[f0:f1])
f0=f1
piecewise_x.append(f0)
piecewise_y.append(v2_headings[f0])
f1=1040
cam_headings[f0:f1]=v2_headings[f0:f1]
f0=f1
piecewise_x.append(f0)
piecewise_y.append(v2_headings[f0])
f1=1137
cam_headings[f0:f1]=linterp(f0,v2_headings[f0],f1,v2_headings[f1],frames[f0:f1])
f0=f1
piecewise_x.append(f0)
piecewise_y.append(v2_headings[f0])
piecewise_x.append(1200)
piecewise_y.append(v2_headings[1200])
piecewise_x.append(1300)
piecewise_y.append(v2_headings[1300])
piecewise_x.append(1400)
piecewise_y.append(v2_headings[1400])
f1=1525
cam_headings[f0:f1]=v2_headings[f0:f1]
f0=f1
piecewise_x.append(f0)
piecewise_y.append(v2_headings[f0])
f1=1595
cam_headings[f0:f1]=linterp(f0,v2_headings[f0],f1,v2_headings[f1],frames[f0:f1])
f0=f1
piecewise_x.append(f0)
piecewise_y.append(v2_headings[f0])
f1=1627
cam_headings[f0:f1]=v2_headings[f0:f1]
f0=f1
piecewise_x.append(f0)
piecewise_y.append(v2_headings[f0])
f1=1755
cam_headings[f0:f1]=v2_headings[f0:f1]+linterp(f0,0,f1,-15,frames[f0:f1])
f0=f1
piecewise_x.append(f0)
piecewise_y.append(v2_headings[f0]-15)
f1=1823
cam_headings[f0:f1]=v2_headings[f0:f1]+linterp(f0,-15,f1,-5,frames[f0:f1])
f0=f1
piecewise_x.append(f0)
piecewise_y.append(v2_headings[f0]-5)
f1=1900
cam_headings[f0:f1]=v2_headings[f0:f1]+linterp(f0,-5,f1,90,frames[f0:f1])
f0=f1
piecewise_x.append(f0)
piecewise_y.append(v2_headings[f0]+90)
f1=frame3
cam_headings[f0:f1]=v2_headings[f0]+90
import scipy.interpolate
interpol=scipy.interpolate.CubicSpline(piecewise_x,piecewise_y,bc_type='clamped')
xs=np.arange(piecewise_x[0],piecewise_x[-1])
plt.plot(frames,v2_headings,frames,cam_headings,piecewise_x,piecewise_y,'*',xs,interpol(xs))
plt.show()
cam_headings[piecewise_x[0]:piecewise_x[-1]]=interpol(xs)
with open("cam_headings.inc","w") as ouf:
    print("#declare cam_headings=array[%d] {"%len(cam_headings),file=ouf)
    for i,cam_heading in enumerate(cam_headings):
        #print("/*%4d*/ %f,"%(i,cam_heading),file=ouf)
        print("/*%4d*/ %f,"%(i,cam_heading),file=ouf)
    print("}",file=ouf)
