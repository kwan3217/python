import spiceypy as cspice
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.widgets as widgets
from bsc import *
import os

os.chdir('/home/jeppesen/workspace/Data/spice/Voyager/')
cspice.furnsh('vgr1.tm')
cspice.furnsh('vgr2.tm')
cspice.furnsh('../generic/lsk/naif0012.tls')
cspice.furnsh('../generic/spk/planets/de430.bsp')
if True:
    framepath='/home/jeppesen/workspace/pov/Blinn/SuperTrajFrames/crop/'
    framenum=94
    infn=framepath+"frame%04dcrop.png"%framenum
else:
    framepath='/home/jeppesen/workspace/pov/Blinn/'
    infn=framepath+"SuperTrajectory.png"

backimg=mpimg.imread(infn)
#backimg=backimg[:,(1280-960)//2:(1280-960)//2+960]

def cmatrix(loc=None,look=None,sky=None):
    """

    :param loc: Location of camera, equivalent to camera{location...}
    :param look: Look-at point of camera, equivalent to camera{look_at...}
    :param sky: Sky vector of camera, equivalent to camera{sky...}
    :return: Camera matrix which transforms a global vector to a vector in camera space
    We will return a 4x4 matrix which will transform a vector in homogeneous coordinates into
    the camera frame. This matrix is:
     [a b c xt]
     [d e f yt]
     [g h i zt]
     [0 0 0  1]
    Multiply this matrix by a 4-element column vector, with the last element being 1
    if the point is a finite distance from the camera (and thus affected by translation)
    zero if the point is infinitely far away, like a star, and thus unaffected by translation.
    The camera frame is set up such that the look direction is +z, right is +x, and up is +y.

    The coefficients a through i could in theory describe a matrix with arbitrary scaling,
    shearing, mirroring, etc. In practice this code will only ever return a matrix which describes
    a pure rotation.
    """

    #Calculate the relative look direction
    look_rel=look-loc
    assert np.linalg.norm(look_rel)>0,"Camera look_at same as location"
    look_rel=look_rel/np.linalg.norm(look_rel)

    #Calculate the right vector as the cross product of the relative look and sky vector
    right=np.cross(look_rel,sky)
    assert np.linalg.norm(right)>0,"Camera looking at sky"
    right=right/np.linalg.norm(right)

    down=np.cross(look_rel,right) #guaranteed to be unit-length since product of two perpendicular unit-length vectors

    r=np.zeros((4,4))
    r[0:3,0]=right
    r[0:3,1]=down
    r[0:3,2]=look_rel
    #result[0:3,3]=-loc
    r[3,3]=1
    t=np.zeros((4,4))
    t[0,0]=1
    t[1,1]=1
    t[2,2]=1
    t[3,3]=1
    t[0:3,3]=loc
    result=t@r
    result=np.linalg.inv(result)
    return result


def project(up=None,right=None,angle=None,width=None,height=None,target_c=None):
    """
    Project the target into the camera field of view
    :param up: Length of Up vector, equivalent to camera{up y*...}
    :param right: Length of Right vector, equivalent to camera{right -x*...}. Note that a right-handed camera should use a positive value for right
    :param angle: Field-of-view angle in degrees, equivalent to camera{angle...}
    :param width: Width of image in pixels, equivalent to image_width
    :param height: Height of image in pixels, equivalent to image_height
    :param target: 3D position of point to project in camera coordinates. If your point is in world coordinates, use
                   cmatrix(...)@target
    :return: 2D position on camera
    """
    #Convert to normalized screen coordinates. In this frame, the screen is on a plane perpendicular and
    #out along the z axis The edges of the screen are at +-0.5*up and +-0.5*right. Angle determines the distance
    #between the camera and the plane of the screen. Using the image at http://www.povray.org/documentation/view/3.7.0/246/
    #as a reference, tan(angle/2)=0.5*right/direction. We can solve this for direction:
    # tan(angle/2)*direction=0.5*right
    # direction=0.5*right/tan(angle/2)
    direction=0.5*right/np.tan(np.radians(angle)/2)
    #if the z component is negative, we don't want to plot. Do this by setting the z component to NaN if it was negative.
    target_c[2,target_c[2,...]<0]=float('NaN')
    #If the target is at the screen, then the x and y coordinates are already what we want. If it is twice as far,
    #then we need to divide x and y by 2. If half as far, then they need to multiply by two. In general, multiply
    # by direction/z. If we do this right, the z coordinate will become equal to direction, which indicates the other
    # components are normalized screen coordinates
    target_scl=target_c[0:2,...]*direction/target_c[2,...]
    result=np.zeros(target_scl.shape)
    result[0,...]=linterp(-0.5*right,0,0.5*right,width, target_scl[0,...])
    result[1,...]=linterp(-0.5*up,   0,0.5*up   ,height,target_scl[1,...])
    result[:,result[0,...]<0]=float('NaN')
    result[:,result[1,...]<0]=float('NaN')
    result[:,result[0,...]>width]=float('NaN')
    result[:,result[1,...]>height]=float('NaN')
    return result

EB2J=cspice.pxform("ECLIPB1950","J2000",0)
J2EB=cspice.pxform("J2000","ECLIPB1950",0)
CamLat=25
CamLon=-36
CamAngle=45.9
CamDist=3.655
CamUpDenom=3.0
imwidth=backimg.shape[1]
imheight=backimg.shape[0]
CamLook= np.array((0.28,-0.15,0.0))
CamSky= np.array((0.0,0.0,1.0))

#Load stars
maxidx=0
LimitMag=5
for i,this_star in enumerate(BrightStarCatalog):
    star=" "+this_star
    if GetMag(star)>LimitMag:
        maxidx=i
        break

etcal="1977-08-21 00:00:00 TDB"
et=cspice.str2et(etcal)
print(etcal,cspice.etcal(et))
v=np.zeros((4,maxidx+5))

ra=[]
dec=[]
name=[]
for i,this_star in enumerate(BrightStarCatalog):
    star=" "+this_star
    if GetMag(star)>LimitMag:
        break
    this_dec=np.radians(GetDec(star))
    this_ra=np.radians(GetRA(star))
    dec.append(this_dec)
    ra.append(this_ra)
    name.append(GetName(star))
    #name.append("")
    this_vec=np.array((np.cos(this_dec)*np.cos(this_ra),np.cos(this_dec)*np.sin(this_ra),np.sin(this_dec)))
    v[0:3,i]=J2EB @ this_vec

au=150000000
TrajFrame="ECLIPB1950"

def project_stuff(v):
    global CamLoc,CamLook,CamDist,CamLat,CamLon,CamSky,CamAngle,imwidth,imheight,CamUpDenom
    CamLoc = (CamDist * np.array((np.cos(np.radians(CamLat)) * np.cos(np.radians(CamLon)),
                                  np.cos(np.radians(CamLat)) * np.sin(np.radians(CamLon)),
                                  np.sin(np.radians(CamLat)))))
    CamLoc+=CamLook
    print("ImWidth:     ", imwidth)
    print("ImHeight:    ", imheight)
    print("CamAngle:    ", CamAngle)
    print("CamLat(deg): ", CamLat)
    print("CamLon(deg): ", CamLon)
    print("CamLook:     ", CamLook)
    print("CamDist:     ", CamDist)
    print("CamLoc:      ", CamLoc)
    print("CamSky:      ", CamSky)
    print("CamUpDenom:  ", CamUpDenom)
    C=cmatrix(loc=CamLoc,look=CamLook,sky=CamSky)
    Cv=C@v
    prj=project(up=1,right=4/CamUpDenom,angle=CamAngle,width=imwidth,height=imheight,target_c=Cv)
    return prj

def traj(body,et0,et1,dt):
    global TrajFrame,au
    et=et0
    v=np.ones((4,int((et1-et0)/dt)))
    for i in range(v.shape[1]):
        et=et0+i*dt
        v[0:3, i] = cspice.spkezr(str(body), et, TrajFrame, "NONE", "0")[0][0:3] / au
    return v

def circpos(a,theta):
    if theta.size>0:
        result=np.zeros((4,theta.size))
    else:
        result=np.zeros(4,1)
    result[0,...]=a*np.cos(theta)
    result[1,...]=a*np.sin(theta)
    result[3,...]=1
    return result

a   =[0,0.38709927,0.72333566,1,1.52371034,5.20288700]
l0  =[0,252.25032350,181.97909950, -4.55343205, 34.39644051]
ldot=[0,149472.67411175, 58517.81538729, 19140.30268499,  3034.74612775]
def circ(a):
    return circpos(a,np.arange(0,2*np.pi,0.01))

fig=plt.figure()
ax=fig.add_subplot(111)
plt.imshow(backimg)
prj=project_stuff(circ(a[1]))
plot1,=plt.plot(prj[0,...],prj[1,...],'w-')
prj=project_stuff(circ(a[2]))
plot2,=plt.plot(prj[0,...],prj[1,...],'w-')
prj=project_stuff(traj(3,et-365*86400,et,86400))
plot3,=plt.plot(prj[0,...],prj[1,...],'w-')
prj=project_stuff(traj(4,et-686*86400,et,86400))
plot4,=plt.plot(prj[0,...],prj[1,...],'w-')
T=et/(86400*36525)
v[0:3,-4]=cspice.spkezr("4",et,TrajFrame,"NONE","0")[0][0:3]/au
v[3,-4]=1
name.append("Mars")
v[0:3,-4]=cspice.spkezr("3",et,TrajFrame,"NONE","0")[0][0:3]/au
v[3,-4]=1
name.append("Earth")
v[:,-3]=circpos(a[2],np.array(np.radians(l0[2]+ldot[2]*T)))[:,0]
name.append("Venus")
v[:,-2]=circpos(a[1],np.array(np.radians(l0[1]+ldot[1]*T)))[:,0]
name.append("Mercury")
v[:,-1]=circpos(a[0],np.array(np.radians(l0[0]+ldot[0]*T)))[:,0]
name.append("Sun")


prj=project_stuff(v)
starplot,=plt.plot(prj[0,...],prj[1,...],'r+')

if False:
    for i in range(prj.shape[1]):
        if np.isfinite(prj[0,i]):
            if name[i]!="":
                print(i, name[i])
                plt.text(prj[0,i],prj[1,i],name[i],color='blue')

stepsize=10

class Index(object):
    def big(self, event):
        global stepsize
        stepsize*=10
    def sm(self, event):
        global stepsize
        stepsize/=10
    def replot(self):
        global v,starplot,plot1,plot2,plot3,plot4,a
        prj=project_stuff(v)
        starplot.set_xdata(prj[0,...])
        starplot.set_ydata(prj[1,...])
        prj=project_stuff(circ(a[1]))
        plot1.set_xdata(prj[0,...])
        plot1.set_ydata(prj[1,...])
        prj=project_stuff(circ(a[2]))
        plot2.set_xdata(prj[0,...])
        plot2.set_ydata(prj[1,...])
        prj=project_stuff(traj(3, et -365 * 86400, et, 86400))
        plot3.set_xdata(prj[0,...])
        plot3.set_ydata(prj[1,...])
        prj=project_stuff(traj(4, et -365 * 86400, et, 86400))
        plot4.set_xdata(prj[0,...])
        plot4.set_ydata(prj[1,...])
        plt.draw()
    def camlonp(self, event):
        global CamLon, stepsize
        CamLon+=stepsize
        self.replot()
    def camlonm(self, event):
        global CamLon, stepsize
        CamLon-=stepsize
        self.replot()
    def camlatp(self, event):
        global CamLat, stepsize
        CamLat+=stepsize
        self.replot()
    def camlatm(self, event):
        global CamLat, stepsize
        CamLat-=stepsize
        self.replot()
    def lookxp(self, event):
        global CamLook, stepsize
        CamLook[0]+=stepsize
        self.replot()
    def lookxm(self, event):
        global CamLook, stepsize
        CamLook[0]-=stepsize
        self.replot()
    def lookyp(self, event):
        global CamLook, stepsize
        CamLook[1]+=stepsize
        self.replot()
    def lookym(self, event):
        global CamLook, stepsize
        CamLook[1]-=stepsize
        self.replot()
    def anglep(self, event):
        global CamAngle, stepsize
        CamAngle += stepsize
        self.replot()
    def anglem(self, event):
        global CamAngle, stepsize
        CamAngle -= stepsize
        self.replot()
    def upp(self, event):
        global CamUpDenom, stepsize
        CamUpDenom += stepsize
        self.replot()
    def upm(self, event):
        global CamUpDenom, stepsize
        CamUpDenom -= stepsize
        self.replot()
    def distp(self, event):
        global CamDist, stepsize
        CamDist += stepsize
        self.replot()
    def distm(self, event):
        global CamDist, stepsize
        CamDist -= stepsize
        self.replot()

callback=Index()
axs=[]
btns=[]
def makebtn(x,y,name,f):
    ax = plt.axes([x, y, 0.05, 0.05])
    axs.append(ax)
    bx = widgets.Button(ax, name)
    btns.append(bx)
    bx.on_clicked(f)

makebtn(0.80,0.05,'-camlon',callback.camlonm)
makebtn(0.70,0.05,'+camlon',callback.camlonp)
makebtn(0.75,0.00,'-camlat',callback.camlatm)
makebtn(0.75,0.10,'+camlat',callback.camlatp)
makebtn(0.60,0.05,'-lookx',callback.lookxm)
makebtn(0.50,0.05,'+lookx',callback.lookxp)
makebtn(0.55,0.00,'-looky',callback.lookym)
makebtn(0.55,0.10,'+looky',callback.lookyp)
makebtn(0.95,0.00,'-angle',callback.anglem)
makebtn(0.95,0.10,'+angle',callback.anglep)
makebtn(0.90,0.00,'-dist',callback.distm)
makebtn(0.90,0.10,'+dist',callback.distp)
makebtn(0.85,0.00,'-up',callback.upm)
makebtn(0.85,0.10,'+up',callback.upp)
makebtn(0.00,0.00,'/step',callback.sm)
makebtn(0.10,0.00,'*step',callback.big)

plt.show()


