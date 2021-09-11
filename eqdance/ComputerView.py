"""
Scenes for video "Why doesn't a rocket need a heat shield?"

"""

from picturebox import PictureBox
import numpy as np
import matplotlib.pyplot as plt
from actor import *

pb=PictureBox(1,1280,720,autodraw=False,facecolor='#e0e0ff')

class AltMark(Actor):
    def __init__(self,ts,x0,y0,y1,d0,d1,data_t,data_y,**kwargs):
        """
        Draw an altitude arrow rising near the vertical axis of the chart

        :param x0: horizontal position
        :param y0: one vertical extreme
        :param y1: other vertical extreme
        :param d0: data value to draw at y0
        :param d1: data value to draw at y1
        :param data_t: data times. Must be monotonically increasing, but no other constraints
        :param data_y: data values at each time data_t, must match length of data_t

        AltMark appears instantly at beginning of entrance, and follows data_t/data_y during
        entrance. Independent variable data_t[0] is mapped to beginning of entrance, and
        data_t[-1] to end. AltMark does default Actor act() and leave().
        """
        super().__init__(ts,**kwargs)
        self.x0=x0
        self.y0=y0
        self.y1=y1
        self.d0=d0
        self.d1=d1
        self.data_t=np.array(data_t)
        self.data_y=np.array(data_y)
    def drawMark(self,pb,x,y,**kwargs):
        pb.line(x,y,x,y-50,**kwargs)
        pb.line(x,y-50,x-10,y-40,**kwargs)
        pb.line(x,y-50,x+10,y-40,**kwargs)
    def enter(self,pb,tt,alpha=1.0,shadow=False):
        this_kwargs=self.kwargs.copy()
        if shadow:
            xofs=5
            yofs=5
            this_kwargs["color"]=shadowcolor
        else:
            xofs=0
            yofs=0
        data_t=linterp(0,self.data_t[0],1,self.data_t[-1],tt)
        data_y=np.interp(data_t,self.data_t,self.data_y)
        self.drawMark(pb,xofs+self.x0,yofs+linterp(self.d0,self.y0,self.d1,self.y1,data_y),alpha=alpha,**this_kwargs)

class DownrangeMark(Actor):
    """
    Should probably be combined with AltMark somehow
    """
    def __init__(self,ts,x0,y0,x1,d0,d1,data_t,data_y,**kwargs):
        super().__init__(ts,**kwargs)
        self.x0=x0
        self.x1=x1
        self.y0=y0
        self.d0=d0
        self.d1=d1
        self.data_t=np.array(data_t)
        self.data_y=np.array(data_y)
    def drawMark(self,pb,x,y,**kwargs):
        pb.line(x,y,x+50,y,**kwargs)
        pb.line(x+50,y,x+40,y-10,**kwargs)
        pb.line(x+50,y,x+40,y+10,**kwargs)
    def enter(self,pb,tt,alpha=1.0,shadow=False):
        this_kwargs=self.kwargs.copy()
        if shadow:
            xofs=5
            yofs=5
            this_kwargs["color"]=shadowcolor
        else:
            xofs=0
            yofs=0
        data_t=linterp(0,self.data_t[0],1,self.data_t[-1],tt)
        data_y=np.interp(data_t,self.data_t,self.data_y)
        self.drawMark(pb,xofs+linterp(self.d0,self.x0,self.d1,self.x1,data_y),yofs+self.y0,alpha=alpha,**this_kwargs)

class TableHighlight(Actor):
    def __init__(self,ts,x0,x1,y0,dy,i,**kwargs):
        super().__init__(ts,**kwargs)
        self.x0=x0
        self.x1=x1
        self.y0=y0+dy*i
        self.y1=y0+dy*(i+1)
    def enter(self,pb,tt,alpha=1.0,shadow=False):
        if shadow:
            return
        this_kwargs=self.kwargs.copy()
        if "alpha" in this_kwargs:
            this_kwargs["alpha"]*=alpha
        else:
            this_kwargs["alpha"]=alpha
        pb.rectangle(self.x0,self.y0,self.x1,self.y1,**this_kwargs)

class SpeedMark(Actor):
    def __init__(self,ts,xc,yc,r,d0,d1,data_t,data_y,**kwargs):
        super().__init__(ts,**kwargs)
        self.xc=xc
        self.yc=yc
        self.r=r
        self.d0=d0
        self.d1=d1
        self.data_t=np.array(data_t)
        self.data_y=np.array(data_y)
    def drawMark(self,pb,speed,shadow=False,**kwargs):
        this_kwargs=kwargs.copy()
        if shadow:
            xofs=5
            yofs=5
            this_kwargs["color"]=shadowcolor
        else:
            xofs=0
            yofs=0
        theta=linterp(self.d0,np.radians(45),self.d1,np.radians(315),speed)
        pb.line(xofs+self.xc-0.5*self.r*np.sin(theta),yofs+self.yc+0.5*self.r*np.cos(theta),
                xofs+self.xc-0.9*self.r*np.sin(theta),yofs+self.yc+0.9*self.r*np.cos(theta),**this_kwargs)
        if not shadow:
            this_kwargs["color"] = 'k'
        for theta in np.radians(np.arange(45, 315)):
            pb.line(xofs+self.xc - self.r * np.sin(theta), yofs+self.yc + self.r * np.cos(theta),
                    xofs+self.xc - self.r * np.sin(theta + np.radians(1)), yofs+self.yc + self.r * np.cos(theta + np.radians(1)),
                    **this_kwargs)
    def enter(self,pb,tt,alpha=1.0,shadow=False):
        data_t=linterp(0,self.data_t[0],1,self.data_t[-1],tt)
        data_y=np.interp(data_t,self.data_t,self.data_y)
        self.drawMark(pb,data_y,alpha=alpha,**self.kwargs,shadow=shadow)

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
f=[tc(0,21, 7), #  0 - Computers love numbers, so they usually use a function or table of numbers.
   tc(0,25, 3), #  1 - This allows the computer to calculate or look up
   tc(0,27,13), #  2 - how high the rocket is,
   tc(0,28,18), #  3 - how far downrange it is,
   tc(0,30, 3), #  4 - and how fast it's going, based on time from launch.
   tc(0,33, 9), #  5 - Looking at all the numbers at once, the computer can figure out where the rocket is at any time.
   tc(0,37,16), #  6 - Humans on the other hand, like pictures.
   tc(0,40, 5), #  7 - Our computer screens are two-dimensional,
   tc(0,42,10), #  8 - so we usually look at two variables at once.
   tc(0,44,22), #  9 - We can look at height,
   tc(0,46, 8), # 10 - downrange,
   tc(0,46,21), # 11 - and so on as a curve versus time.
   tc(0,49, 8), # 12 - Or, we can look at different variables versus each other.
   tc(0,52,22), # 13 - For this problem, maybe the most useful curve has
   tc(0,55,12), # 14 - speed on the horizontal axis,
   tc(0,57, 9), # 15 - and height on the vertical.
   tc(0,59,19), # 16 - <END>
]

actors=[]
rightAxis=Axis(ts=[f[0],f[1],f[5],f[5]+12],x0=pb.w*1/3+100,y0=100,x1=pb.w-100,y1=pb.h-100,color='k')
actors.append(rightAxis)
actors.append(TableGrid     (ts=[f[2]   ,f[2]+20,f[ 6]   ,f[ 6]+12],x0=50,x1=420,yt=80,y0=102,yb=pb.h-70,xs=[105,205,315],color='k'))
actors.append(TableColumn   (ts=[f[2]   ,f[2]+20,f[ 6]   ,f[ 6]+12],header="Time/s"     ,data=stime     [6::15],x=100,y0=100,dy=15,horizontalalignment='right',color='k'))
actors.append(TableColumn   (ts=[f[2]+ 4,f[2]+24,f[ 6]   ,f[ 6]+12],header="Altitude/m" ,data=salt      [6::15],x=200,y0=100,dy=15,horizontalalignment='right',color='b'))
actors.append(AltMark       (ts=[f[2]+ 4,f[2]+24,f[ 5]   ,f[ 5]+12],x0=rightAxis.x0+20,y0=rightAxis.y1,y1=rightAxis.y0,d0=0,d1=110000,data_t=time,data_y=alt,color='b'))
actors.append(TableColumn   (ts=[f[3]   ,f[3]+20,f[ 6]   ,f[ 6]+12],header="Downrange/m",data=sdownrange[6::15],x=310,y0=100,dy=15,horizontalalignment='right',color='g'))
actors.append(DownrangeMark (ts=[f[3]   ,f[3]+20,f[ 5]   ,f[ 5]+12],x0=rightAxis.x0,y0=rightAxis.y1-20,x1=rightAxis.x1,d0=0,d1=1600000,data_t=time,data_y=downrange,color='g'))
actors.append(TableColumn   (ts=[f[4]   ,f[4]+20,f[ 6]   ,f[ 6]+12],header="Speed/(m/s)",data=svrel     [6::15],x=400,y0=100,dy=15,horizontalalignment='right',color='r'))
actors.append(SpeedMark     (ts=[f[4]   ,f[4]+20,f[ 5]   ,f[ 5]+12],xc=(rightAxis.x0+rightAxis.x1)/2,yc=(rightAxis.y0+rightAxis.y1)/2,r=200,d0=0,d1=8000,data_t=time,data_y=vrel,color='r'))
actors.append(TableHighlight(ts=[f[5]   ,f[5]+ 6,f[ 6]- 6,f[ 6]   ],x0=50,x1=105,y0=102,dy=15,i=10,color='k',alpha=0.3,fill=True))
actors.append(TableHighlight(ts=[f[5]   ,f[5]+ 6,f[ 6]- 6,f[ 6]   ],x0=105,x1=205,y0=102,dy=15,i=10,color='b',alpha=0.3,fill=True))
actors.append(TableHighlight(ts=[f[5]   ,f[5]+ 6,f[ 6]- 6,f[ 6]   ],x0=205,x1=315,y0=102,dy=15,i=10,color='g',alpha=0.3,fill=True))
actors.append(TableHighlight(ts=[f[5]   ,f[5]+ 6,f[ 6]- 6,f[ 6]   ],x0=315,x1=420,y0=102,dy=15,i=10,color='r',alpha=0.3,fill=True))
actors.append(Text          (ts=[f[5]   ,f[5]+ 6,f[ 6]- 6,f[ 6]   ],s="At t="+stime[6::15][10]+" s,",x=(rightAxis.x0+rightAxis.x1)/2,y=(rightAxis.y0+rightAxis.y1)/2- 60,fontfamily="serif",fontsize=40,horizontalalignment='center',color='k'))
actors.append(Text          (ts=[f[5]+12,f[5]+18,f[ 6]- 6,f[ 6]   ],s="altitude="+salt[6::15][10]+" m,",x=(rightAxis.x0+rightAxis.x1)/2,y=(rightAxis.y0+rightAxis.y1)/2    ,fontfamily="serif",fontsize=40,horizontalalignment='center',color='b'))
actors.append(Text          (ts=[f[5]+24,f[5]+30,f[ 6]- 6,f[ 6]   ],s="downrange="+sdownrange[6::15][10]+" m,",x=(rightAxis.x0+rightAxis.x1)/2,y=(rightAxis.y0+rightAxis.y1)/2+ 60,fontfamily="serif",fontsize=40,horizontalalignment='center',color='g'))
actors.append(Text          (ts=[f[5]+36,f[5]+42,f[ 6]- 6,f[ 6]   ],s="and speed="+svrel[6::15][10]+" m/s",x=(rightAxis.x0+rightAxis.x1)/2,y=(rightAxis.y0+rightAxis.y1)/2+120,fontfamily="serif",fontsize=40,horizontalalignment='center',color='r'))

midAxis=Axis(ts=[f[8],f[8]+12,f[16],f[16]],x0=100,y0=100,x1=pb.w-100,y1=pb.h-100,color='k')
actors.append(midAxis)
actors.append(Plot(ts=[f[ 9],f[ 9]+12,f[12]-12,f[12]],
                   px0=midAxis.x0,dx0=0,px1=midAxis.x1,dx1=   600,data_x=time,
                   py0=midAxis.y1,dy0=0,py1=midAxis.y0,dy1=120000,data_y=alt,
                   t0=0,t1=time[-1],data_t=time,color='b'))
actors.append(Plot(ts=[f[10],f[10]+12,f[12]-12,f[12]],
                   px0=midAxis.x0,dx0=0,px1=midAxis.x1,dx1=    600,data_x=time,
                   py0=midAxis.y1,dy0=0,py1=midAxis.y0,dy1=1400000,data_y=downrange,
                   t0=0,t1=time[-1],data_t=time,color='g'))
actors.append(Plot(ts=[f[11],f[11]+12,f[12]-12,f[12]],
                   px0=midAxis.x0,dx0=0,px1=midAxis.x1,dx1=    600,data_x=time,
                   py0=midAxis.y1,dy0=0,py1=midAxis.y0,dy1=   8000,data_y=vrel,
                   t0=0,t1=time[-1],data_t=time,color='r'))
actors.append(Plot(ts=[f[13],f[13]+24,f[16],f[16]],
                   px0=midAxis.x0,dx0=0,px1=midAxis.x1,dx1=   8000,data_x=vrel,
                   py0=midAxis.y1,dy0=0,py1=midAxis.y0,dy1= 120000,data_y=alt,
                   t0=0,t1=time[-1],data_t=time,color='y'))

perform(pb,actors,f[0],f[-1],"../../video/ComputerView/ComputerView%04d.png")

