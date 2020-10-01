"""
Scenes for video "Why doesn't a rocket need a heat shield?"

"""

from picturebox import PictureBox
import numpy as np
import matplotlib.pyplot as plt
from actor import *

pb=PictureBox(1280,720,autodraw=False,facecolor='#e0e0ff')

stime=[]
salt=[]
svrel=[]
sdownrange=[]
time=[]
alt=[]
vrel=[]
downrange=[]
with open("STS122Launch.csv","rt") as inf:
    columns=inf.readline()[:-1].split(",")
    for line in inf:
        parts=line[:-1].split(",")
        time.append(float(parts[0]))
        alt.append(float(parts[1]))
        vrel.append(float(parts[2]))
        downrange.append(float(parts[3]))
        stime.append("%.1f"%float(parts[0]))
        salt.append("%.0f"%(float(parts[1])-48)) #Display rocket as starting at alt=0
        svrel.append("%.0f"%float(parts[2]))
        sdownrange.append("%.0f"%float(parts[3]))
time=np.array(time)
alt=np.array(alt)
vrel=np.array(vrel)
downrange=np.array(downrange)

#Kdenlive audio is split this same way. Each time code is the beginning of the indicated phrase.
f=[tc(2,17,20), #  0
   tc(2,31,19), #  1 <END>
]

actors=[]

def expatm(vrel,alt):
    return np.broadcast_arrays(np.exp(-alt/11000),vrel)[0]

def diagfade(vrel,alt):
    return vrel+alt

midAxis=Axis(ts=[f[0],f[0],f[-1],f[-1]],x0=100,y0=100,x1=pb.w-100,y1=pb.h-100,color='k')
actors.append(midAxis)
actors.append(Image(ts=[f[ 0],f[ 0]+48,f[-1]-48,f[-1]],
                   px0=midAxis.x0,dx0=0,px1=midAxis.x1,dx1=120000,nx=1000,
                   py0=midAxis.y1,dy0=120000,py1=midAxis.y0,dy1=0,ny= 500,
                   f    =expatm,ffade=diagfade))

perform(pb,actors,f[0],f[-1],"../../video/atmosphere_factor/atmosphere_factor%04d.png",shadow=False)

print("Done")