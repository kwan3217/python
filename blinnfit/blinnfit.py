import spiceypy as cspice
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.widgets as widgets
from bsc import *
import os
import sqlite3

fitsky=False

def open_frame_index():
    """
    Open an SQLite database and make sure that the appropriate table(s) are present
    in the database

    :param desc: Used to build the filename of the sqlite database. Different values
      of desc= will result in different sqlite files being used.
    """
    if fitsky:
        dbname = "frame_index.sqlite"
    else:
        dbname="framenfs_index.sqlite"
    #log.info(dbname)
    conn=sqlite3.connect(dbname)
    sql=("create table if not exists frames ("+
      "framenum    integer not null,"+
      "lat_c       real,"+
      "lon_c       real,"+
      "dist        real,"+
      "angle       real,"+
      "right_denom real,"+
      "lookx       real,"+
      "looky       real,"+
      "lookz       real,"+
      "lat_s       real,"+
      "lon_s       real,"+
      "clock       real,"+
      "primary key (framenum))")
    #log.info(sql)
    cur=conn.cursor()
    cur.execute(sql)
    conn.commit()
    return conn

conn=open_frame_index()

os.chdir('/home/jeppesen/workspace/Data/spice/Voyager/')
cspice.furnsh('vgr1.tm')
cspice.furnsh('vgr2.tm')
cspice.furnsh('../generic/lsk/naif0012.tls')
cspice.furnsh('../generic/spk/planets/de430.bsp')
if True:
    framepath='/home/jeppesen/workspace/Data/Prototype/Blinn/SuperTrajectory/'
    framenum=432
    infn=framepath+"frame%04d.png"%framenum
else:
    framepath='/home/jeppesen/workspace/pov/Blinn/'
    infn=framepath+"SuperTrajectory.png"
backimg=mpimg.imread(infn)

goodstars=[1553,1097, 191, 765, 613, 289, 456, 472, 968,1000,
           1292,1443, 411,1400,1243, 979, 801,1013,1492, 791,
            741,1149,1020,1496,1399, 724,1569,1148,1108, 467,
           1475, 887,  47, 666, 133,1414, 241, 909, 845,1109,
           1021, 501,1045, 436,  88,1323, 848,  65, 277,  22,
            980,1611, 379, 495, 755,1585, 931, 300,1585,1371,
           1291, 205, 278,  92, 197, 886, 565, 357, 159,1491,
            494,  37, 391, 515,   1, 175,1290,1320,1035,1199,
            237, 799, 790,1060, 908, 227, 510, 877, 753,1136,
            433,1269, 750,1147, 489, 516, 386,1474,1318, 226,

            580, 397, 661, 850,1146, 432, 220,1058, 308, 567,
            879,1497,1416,1324,1002,1124,1159, 301, 172,1444,
            530, 172,  90, 694,1125, 107, 821, 517,1557,1526,
            122, 458, 682,  94,1189, 811, 744, 263, 737,  98,
           1046,1446, 173,1139, 228,  57,1062, 727, 111,1229,
            216, 128, 318, 491, 552,1007, 563, 454, 252, 346,
            949, 309, 460,1064, 648,1039,1595, 871, 614, 649,
            387,  75, 200, 240, 142, 243, 334, 163,1387, 555,
           1217,1559,  24, 115, 706,  39, 992, 347,  83,1376,

            214, 193, 594, 972, 951, 924, 208, 239,  36,1510,
            584, 505, 705, 595, 194, 983, 784,1421,  13,1390,
           1378, 586, 603,1092,1308,1222,1434,1142, 177, 860,
            618, 633, 862,1360, 192, 566,1589, 591, 738, 100,
            154,1602, 955, 135, 392, 900, 564, 118, 376, 746,
            147, 711, 686,  16,  30,1433,  60,  43, 796, 815,
            785, 718,  52, 285, 320, 230, 976, 835,1513, 960,
             81,1267,1145,1181, 284,1118, 305,1194,1422,1350,
            513, 689, 461, 182, 161, 448, 543, 872, 470, 678,
            479, 371, 660,1183, 275, 557,  93, 700, 375, 732,

            722, 855,1071, 646,1604,1423, 907, 576, 307, 977,
           1184, 677, 805, 953, 306, 659, 843,1317, 106, 143]


def llr2xyz(lat,lon,r=1):
    """Calculate rectangular vector coordinates from spherical coordinates
    :param lat: Latitude in degrees
    :param lon: Longitude in degrees
    :param r:   Radius
    """
    return (r * np.array((np.cos(np.radians(lat)) * np.cos(np.radians(lon)),
                          np.cos(np.radians(lat)) * np.sin(np.radians(lon)),
                          np.sin(np.radians(lat)))))

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
    look_rel = look - loc
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

def project(up=None,right=None,angle=None,width=None,height=None,cx=None,cy=None,target_c=None,clock=None):
    """
    Project the target into the camera field of view
    :param up: Length of Up vector, equivalent to camera{up y*...}
    :param right: Length of Right vector, equivalent to camera{right -x*...}. Note that a right-handed camera should use a positive value for right
    :param angle: Field-of-view angle in degrees, equivalent to camera{angle...}
    :param width: Width of image in pixels, equivalent to image_width
    :param height: Height of image in pixels, equivalent to image_height
    :param cx: center of projection in x direction in pixels. In POV-Ray perspective camera, cx is always implicitly image_width/2
    :param cy: center of projection in y direction in pixels. In POV-Ray perspective camera, cy is always implicitly image_height/2
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
    if cx is None:
        cx=width/2
    if cy is None:
        cy=height/2
    if clock is None:
        clock=0
    rx=linterp(-0.5*right,-width /2,0.5*right,width /2,target_scl[0,...])
    ry=linterp(-0.5*up,   -height/2,0.5*up   ,height/2,target_scl[1,...])
    result[0,...]=rx*np.cos(np.radians(clock))-ry*np.sin(np.radians(clock))+cx
    result[1,...]=rx*np.sin(np.radians(clock))+ry*np.cos(np.radians(clock))+cy
    result[:,result[0,...]<0]=float('NaN')
    result[:,result[1,...]<0]=float('NaN')
    result[:,result[0,...]>width]=float('NaN')
    result[:,result[1,...]>height]=float('NaN')
    return result

EB2J=cspice.pxform("ECLIPB1950","J2000",0)
J2EB=cspice.pxform("J2000","ECLIPB1950",0)

#Load stars
maxidx=0
LimitMag=5
for i,this_star in enumerate(BrightStarCatalog):
    star=" "+this_star
    if GetMag(star)>LimitMag:
        maxidx=i
        break

def frame_to_et(framenum):
    frame0et=cspice.str2et("1970-07-25 20:40:26 TDB")
    h=3.647870235*86400
    return framenum*h+frame0et

et=frame_to_et(framenum)
#et=cspice.str2et(etcal)
print(cspice.etcal(et))
v=np.zeros((4,maxidx+10))

name=[]
for i,this_star in enumerate(BrightStarCatalog):
    star=" "+this_star
    if GetMag(star)>LimitMag:
        break
    dec=np.radians(GetDec(star))
    ra=np.radians(GetRA(star))
    name.append(GetName(star))
    this_vec=np.array((np.cos(dec)*np.cos(ra),
                       np.cos(dec)*np.sin(ra),
                       np.sin(dec)))
    v[0:3,i]=J2EB @ this_vec

au = 150000000
TrajFrame = "ECLIPB1950"

def traj(body,et0,et1,dt):
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

def circ(a):
    return circpos(a,np.arange(0,2*np.pi,0.01))

boxfig=None
boxax=None
boximg=None

def find_star(starm,boxr=10,vis=True):
    import twodgauss
    #Subset the image
    reject=False
    try:
        xo=int(prj[0,starm])
        yo=int(prj[1,starm])
    except ValueError:
        #reject=True
        #print("Reject: No coordinate for star",starm)
        return float('nan'), float('nan')
    if(xo<boxr):
        reject=True
        print("Reject: off left edge %",xo)
    if(xo>backimg.shape[1]-boxr):
        reject=True
        print("Reject: off right edge %",xo)
    if(yo<boxr):
        reject=True
        print("Reject: off top edge %",yo)
    if(yo>backimg.shape[0]-boxr):
        reject=True
        print("Reject: off bottom edge %",yo)
    if reject:
        return float('nan'),float('nan')
    box=backimg[yo-boxr:yo+boxr,xo-boxr:xo+boxr,:]
    box=np.sum(box,2)
    box=box-np.min(box)

    if vis:
        global boxfig,boxax,boximg
        if boximg is None:
            boxfig, boxax = plt.subplots(1, 1,num="box")
            boximg=boxax.imshow(box)
        else:
            boxax.cla()
            boximg=boxax.imshow(box)
        boxax.set_title(name[starm])
        plt.pause(0.001)

    #Try the 2D Gaussian fit on the subset
    try:
        (amp, xo,   yo,   sigx, sigy, theta, ofs, data_fitted)=twodgauss.fit_twoD_Gaussian(box,
         0.5, boxr, boxr, 5,    3,    0,     0
        )
    except RuntimeError:
        print("Reject: fit failed")
        return float('nan'),float('nan')
    xs=box.shape[1]
    ys=box.shape[0]
    x=xo+int(prj[0,starm])-boxr
    y=yo+int(prj[1,starm])-boxr
    print("%4d %7.3f %7.3f %6.3f %6.3f"%(starm,x,y,sigx,sigy))
    if np.abs(sigx)>5:
        print("Reject: bad sigx %f"%sigx)
        reject=True
    if np.abs(sigy)>5:
        print("Reject: bad sigy %f"%sigy)
        reject=True
    if reject:
        x=float('nan')
        y=float('nan')
    if vis:
        color='r' if reject else 'g'
        xg = np.linspace(0, xs - 1, xs)
        yg = np.linspace(0, ys - 1, ys)
        xg, yg = np.meshgrid(xg, yg)
        boxax.contour(xg, yg, data_fitted, 8, colors=color)
        plt.pause(0.001)
    return x,y

def curve_fit_interface(starvec,lat_c,lon_c,angle,right_denom,clock):
    """
    Calculate the pixel positions of the given stars, given these camera parameters
    :param starvec: List of star vectors with homogeneous coordinates of shape (4,M//2)
    :param camlat: Scalar camera latitude in degrees
    :param camlon: Scalar camera longitude in degrees
    :param dist:   Scalar camera distance in AU
    :param angle:  Scalar camera FOV angle in degrees
    :param right_denom: Scalar camera aspect ratio constant
    :param cx:     Scalar image distortion center horizontal coordinate
    :param cy:     Scalar image distortion center vertical coordinate
    :param skylat: Scalar sky vector latitude in degrees
    :param skylon: Scalar sky vector longitude in degrees
    :return: 1D array of shape M, representing a 2D array of pixel coordinates [x_or_y,star] shape (2,M//2)
             raveled so as to work with scipy.optimize.curve_fit. This will be all the x coordinates first,
             then all the y coordinates
    """
    # The curve fitter scipy.optimize.curve_fit takes a function to fit f,
    # independent xdata (can be any object, but f(xdata,*p) must return an
    # array of shape M), dependent ydata (shape M), and an initial guess at
    # a set of parameters p0 (shape N). It returns a set of parameters popt
    # which best fits the data. In our case, the ydata is pixel positions of
    # the stars, and therefore M is twice the number of stars we are trying
    # to fit. The p is camera parameters, and therefore by process of elimination
    # the xdata must be the positions of the stars. In our case, it's easiest to take
    # the vectors of the stars as inputs, so xdata will be an array of shape
    # (4,M//2) and we will return a 1D array of raveled x and y pixel coordinates
    # of each star
    loc = llr2xyz(lat=lat_c,lon=lon_c)
    sky = np.array((0,0,1))
    C = cmatrix(loc=loc, look=np.zeros(3),sky=sky)
    Cv = C @ starvec
    prj = project(up=1, right=4 / right_denom, angle=angle, width=backimg.shape[1], height=backimg.shape[0], target_c=Cv,
              cx=None, cy=None, clock=clock)
    return prj.ravel()

def curve_fitsky_interface(starvec,lat_c,lon_c,angle,right_denom,lat_s,lon_s):
    """
    Calculate the pixel positions of the given stars, given these camera parameters
    :param starvec: List of star vectors with homogeneous coordinates of shape (4,M//2)
    :param camlat: Scalar camera latitude in degrees
    :param camlon: Scalar camera longitude in degrees
    :param dist:   Scalar camera distance in AU
    :param angle:  Scalar camera FOV angle in degrees
    :param right_denom: Scalar camera aspect ratio constant
    :param cx:     Scalar image distortion center horizontal coordinate
    :param cy:     Scalar image distortion center vertical coordinate
    :param skylat: Scalar sky vector latitude in degrees
    :param skylon: Scalar sky vector longitude in degrees
    :return: 1D array of shape M, representing a 2D array of pixel coordinates [x_or_y,star] shape (2,M//2)
             raveled so as to work with scipy.optimize.curve_fit. This will be all the x coordinates first,
             then all the y coordinates
    """
    # The curve fitter scipy.optimize.curve_fit takes a function to fit f,
    # independent xdata (can be any object, but f(xdata,*p) must return an
    # array of shape M), dependent ydata (shape M), and an initial guess at
    # a set of parameters p0 (shape N). It returns a set of parameters popt
    # which best fits the data. In our case, the ydata is pixel positions of
    # the stars, and therefore M is twice the number of stars we are trying
    # to fit. The p is camera parameters, and therefore by process of elimination
    # the xdata must be the positions of the stars. In our case, it's easiest to take
    # the vectors of the stars as inputs, so xdata will be an array of shape
    # (4,M//2) and we will return a 1D array of raveled x and y pixel coordinates
    # of each star
    loc = llr2xyz(lat=lat_c,lon=lon_c)
    sky = llr2xyz(lat=lat_s,lon=lon_s)
    C = cmatrix(loc=loc, look=np.zeros(3),sky=sky)
    Cv = C @ starvec
    prj = project(up=1, right=4 / right_denom, angle=angle, width=backimg.shape[1], height=backimg.shape[0], target_c=Cv,
              cx=None, cy=None)
    return prj.ravel()

class CameraMount(object):
    def __init__(self,conn):
        self.fitplot=None
        self.conn=conn
        self.read(framenum)
        self.width = backimg.shape[1]
        self.height = backimg.shape[0]
        self.cx = self.width / 2
        self.cy = self.height / 2
        self.step_size=10
    def read(self,framenum):
        has_row=False
        #Check if this frame is already recorded
        sql="select lat_c,lon_c,dist,angle,right_denom,lookx,looky,lookz,lat_s,lon_s,framenum,clock from frames order by abs(framenum-?) asc"
        cur = self.conn.cursor()
        for row in cur.execute(sql, (framenum,)):
            if row[10]==framenum:
                self.lat_c        = row[0]
                self.lon_c        = row[1]
                self.dist         = row[2]
                self.angle        = row[3]
                self.right_denom  = row[4]
                self.look=np.array((row[5],
                                    row[6],
                                    row[7]))
                self.lat_s        = row[8]
                self.lon_s        = row[9]
                self.clock        = row[11]
                has_row = True
                break
            elif not has_row:
                fn0=row[10]
                lat_c0       = row[0]
                lon_c0       = row[1]
                dist0        = row[2]
                angle0       = row[3]
                right_denom0 = row[4]
                lookx0       = row[5]
                looky0       = row[6]
                lookz0       = row[7]
                lat_s0       = row[8]
                lon_s0       = row[9]
                clock0        = row[11]
                has_row=True
            else:
                fn1          = row[10]
                lat_c1       = row[0]
                lon_c1       = row[1]
                dist1        = row[2]
                angle1       = row[3]
                right_denom1 = row[4]
                lookx1       = row[5]
                looky1       = row[6]
                lookz1       = row[7]
                lat_s1        = row[8]
                lon_s1        = row[9]
                clock1        = row[11]
                self.lat_c       = linterp(fn0,lat_c0,fn1,lat_c1,framenum)
                self.lon_c       = linterp(fn0,lon_c0,fn1,lon_c1,framenum)
                self.dist        = linterp(fn0,dist0,fn1,dist1,framenum)
                self.angle       = linterp(fn0,angle0,fn1,angle1,framenum)
                self.right_denom = linterp(fn0,right_denom0,fn1,right_denom1,framenum)
                self.look=np.array((linterp(fn0,lookx0,fn1,lookx1,framenum),
                                    linterp(fn0,looky0,fn1,looky1,framenum),
                                    linterp(fn0,lookz0,fn1,lookz1,framenum)))
                self.lat_s       = linterp(fn0,lat_s0,fn1,lat_s1,framenum)
                self.lon_s       = linterp(fn0,lon_s0,fn1,lon_s1,framenum)
                self.clock       = linterp(fn0,clock0,fn1,clock1,framenum)
                break
        if not has_row:
            #initial conditions valid for frame 675
            self.lat_c = 25.186          #Spherical coordinate latitude of camera position relative to its look point, in degrees
            self.lon_c = -35.169         #Spherical coordinate longitude of camera position, in degrees
            self.dist = 3.655            #Spherical coordinate radius of camera position
            self.angle = 45.987          #Angle parameter of camera, equivalent to angle keyword in POV-Ray perspective camera
            self.right_denom = 2.943643  #Camera right vector denominator -- right vector in POV-Ray perspective camera is right -x*4/right_denom
            self.look = np.array(( 0.28,
                                  -0.15,
                                   0.00))
            self.lat_s=90                #Spherical coordinate latitude of sky vector, degrees
            self.lon_s=1.455             #Spherical coordinate longitude of sky vector, degrees. OK but ineffectual if lat_s=+-90deg.
            self.clock=0
    def write(self):
        sql="insert or replace into frames (framenum,lat_c,lon_c,dist,angle,right_denom,lookx,looky,lookz,lat_s,lon_s,clock) values (?,?,?,?,?,?,?,?,?,?,?,?)"
        cur = conn.cursor()
        cur.execute(sql, (framenum,
                          self.lat_c,
                          self.lon_c,
                          self.dist,
                          self.angle,
                          self.right_denom,
                          self.look[0],
                          self.look[1],
                          self.look[2],
                          self.lat_s,
                          self.lon_s,
                          self.clock))
        conn.commit()
    def project_stuff(self,v):
        loc=llr2xyz(lat=self.lat_c,lon=self.lon_c,r=self.dist) + self.look
        sky=llr2xyz(lat=self.lat_s,lon=self.lon_s)
        print("image_width  ", self.width)
        print("image_height ", self.height)
        print("angle        ", self.angle)
        print("lat_c(deg)   ", self.lat_c)
        print("lon_c(deg)   ", self.lon_c)
        print("dist         ", self.dist)
        print("lat_s(deg)   ", self.lat_s)
        print("lon_s(deg)   ", self.lon_s)
        print("right   -x*4/", self.right_denom)
        print("cx           ", self.cx)
        print("cy           ", self.cy)
        print("clock        ", self.clock)
        print("framenum     ", framenum)
        C = cmatrix(loc=loc, look=self.look, sky=sky)
        Cv = C @ v
        prj = project(up=1, right=4 / self.right_denom, angle=self.angle, width=self.width, height=self.height, target_c=Cv, cx=self.cx, cy=self.cy, clock=self.clock)
        return prj
    def big(self, event):
        self.step_size*=10
    def sm(self, event):
        self.step_size/=10
    def replot(self):
        global v,starplot,orbplot,prj,et
        if True:
            for i_planet in range(10):
                if orbcirc[i_planet]:
                    v[:, maxidx + i_planet] = circpos(a[i_planet], np.array(
                        np.radians(l0[i_planet] + ldot[i_planet] * (et / 36525 / 86400))))[:, 0]
                else:
                    v[0:3, maxidx + i_planet] = cspice.spkezr(str(i_planet), et, TrajFrame, "None", "0")[0][0:3] / au
                    v[3, maxidx + i_planet] = 1
                prj = callback.project_stuff(orbtraj[i_planet])
                orbplot[i_planet].set_xdata(prj[0,...])
                orbplot[i_planet].set_ydata(prj[1,...])
        prj=self.project_stuff(v)
        for i,n_o in enumerate(nameobj):
            if n_o is not None:
                n_o.set_visible(np.isfinite(prj[0,i]))
                if np.isfinite(prj[0,i]):
                    n_o.set_position((prj[0,i],prj[1,i]))
        starplot.set_xdata(prj[0,...])
        starplot.set_ydata(prj[1,...])
        plt.pause(0.001)
    def camlonp(self, event):
        self.lon_c+=self.step_size
        self.replot()
    def camlonm(self, event):
        self.lon_c-=self.step_size
        self.replot()
    def camlatp(self, event):
        self.lat_c+=self.step_size
        self.replot()
    def camlatm(self, event):
        self.lat_c-=self.step_size
        self.replot()
    def lookxp(self, event):
        self.look[0]+=self.step_size
        self.replot()
    def lookxm(self, event):
        self.look[0]-=self.step_size
        self.replot()
    def lookyp(self, event):
        self.look[1]+=self.step_size
        self.replot()
    def lookym(self, event):
        self.look[1]-=self.step_size
        self.replot()
    def anglep(self, event):
        self.angle += self.step_size
        self.replot()
    def anglem(self, event):
        self.angle -= self.step_size
        self.replot()
    def upp(self, event):
        self.right_denom += self.step_size
        self.replot()
    def upm(self, event):
        self.right_denom -= self.step_size
        self.replot()
    def clockp(self, event):
        self.clock += self.step_size
        self.replot()
    def clockm(self, event):
        self.clock -= self.step_size
        self.replot()
    def distp(self, event):
        self.dist += self.step_size
        self.replot()
    def distm(self, event):
        self.dist -= self.step_size
        self.replot()
    def framem(self, event):
        global framenum,backimg,et
        framenum-=1
        et = frame_to_et(framenum)
        print(cspice.etcal(et))
        self.read(framenum)
        infn = framepath + "frame%04d.png" % framenum
        backimg = mpimg.imread(infn)
        figimg.set_data(backimg)
        if self.fitplot is not None:
            self.fitplot.set_visible(False)
        self.replot()
    def framep(self, event):
        global framenum,backimg,et
        framenum+=1
        et = frame_to_et(framenum)
        print(cspice.etcal(et))
        self.read(framenum)
        infn = framepath + "frame%04d.png" % framenum
        backimg = mpimg.imread(infn)
        figimg.set_data(backimg)
        if self.fitplot is not None:
            self.fitplot.set_visible(False)
        self.replot()
    def fit(self, event, vis=False):
        """
        Given the current position as an initial guess, find the optimum
        camera parameters and position to fit the stars.
        """
        global ax
        import scipy.optimize
        fitx=np.zeros(len(goodstars))*float('nan')
        fity=np.zeros(len(goodstars))*float('nan')
        fitv=np.zeros((4,len(goodstars)))*float('nan')
        for i,i_goodstar in enumerate(goodstars):
            fitx[i],fity[i]=find_star(i_goodstar,vis=vis)
            fitv[:,i]=v[:,i_goodstar]
        if self.fitplot is None:
            self.fitplot,=ax.plot(fitx,fity,'g*')
        else:
            self.fitplot.set_visible(True)
            self.fitplot.set_xdata(fitx)
            self.fitplot.set_ydata(fity)
        #plt.pause(0.001)
        w=np.where(np.isfinite(fitx))
        fitx=fitx[w]
        fity=fity[w]
        pixdata=np.stack((fitx,fity)).ravel()
        fitv=fitv[:,w[0]] #If we do fitv[:,w] we get a shape (4,1,nstars) instead of the (4,nstars) we want
        #fiti=np.array(goodstars)[w]
        if fitsky:
            p0=np.array((self.lat_c,self.lon_c,self.angle,self.right_denom,self.lat_s,self.lon_s))
            mn=np.array((-90,       -np.inf,   0,         0,               -90,       -np.inf   ))
            mx=np.array(( 90,        np.inf,   90,        np.inf,           90,        np.inf   ))
            if p0[4] > mx[4]:
                p0[4] = 90 - p0[4]
                p0[5] = 180 + p0[5]
        else:
            p0=np.array((self.lat_c,self.lon_c,self.angle,self.right_denom,self.clock))
            mn=np.array((-90,       -np.inf,   0,         0,               -np.inf))
            mx=np.array(( 90,        np.inf,   90,        np.inf,           np.inf))
        try:
            if fitsky:
                (popt,_)=scipy.optimize.curve_fit(curve_fitsky_interface,fitv,pixdata,p0=p0,bounds=(mn,mx))
            else:
                (popt, _) = scipy.optimize.curve_fit(curve_fit_interface, fitv, pixdata, p0=p0, bounds=(mn, mx))
        except ValueError:
            print(p0)
            print(p0>=mn)
            print(p0<=mx)
            raise
        while popt[1]>360:
            popt[1]-=360
        while popt[1]<0:
            popt[1]+=360
        if fitsky:
            while popt[5] > 360:
                popt[5] -= 360
            while popt[5] < 0:
                popt[5] += 360
            (self.lat_c,self.lon_c,self.angle,self.right_denom,self.lat_s,self.lon_s)=popt
        else:
            (self.lat_c,self.lon_c,self.angle,self.right_denom,self.clock)=popt
        self.write()
        self.replot()
    def autop(self, event):
        for i in range(2171-94):
            self.framep(event)
            self.fit(event,vis=False)
            self.fit(event,vis=False)
    def autom(self, event):
        for i in range(100):
            self.framem(event)
            self.fit(event,vis=False)
            self.fit(event,vis=False)



callback=CameraMount(conn)
fig=plt.figure()
ax=fig.add_subplot(111)
figimg=plt.imshow(backimg)
#                 0      1             2           3            4            5            6            7            8            9
plname=          ["Sun","Mercury"     ,"Venus"    ,"Earth"     ,"Mars"      ,"Jupiter"   ,"Saturn"    ,"Uranus"    ,"Neptune"   ,"Pluto"]
orbcirc=         [True,  True         ,True       ,True        ,True        ,False       ,False       ,True        ,True        ,False]
a      =np.array((0,       0.38709927,  0.72333566,  1.0       ,  1.52371034,  5.20288700,  9.53667594, 19.18916464, 30.06992276, 39.48211675))
l0     =np.array((0,     252.25032350,181.97909950,100.46457166, -4.55343205, 34.39644051, 49.95424423,313.23810451,-55.12002969,238.92903833))
ldot   =a*0
ldot[1:]=100*360/np.sqrt(a[1:]**3) #Don't try to calculate orbital period of Sun
T=360*365*100*86400/ldot #Orbit period in seconds
#ldot   =np.array((0, 149472.67411175,58517.81538729, 19140.30268499,  3034.74612775]

orbplot=[None]*10 #Sun is planet 0 in this case
orbtraj=[None]*10
name+=plname
for i_planet in range(10):
    if orbcirc[i_planet]:
        v[:, maxidx+i_planet] = circpos(a[i_planet], np.array(np.radians(l0[i_planet] + ldot[i_planet] * (et/36525/86400))))[:, 0]
        orbtraj[i_planet]=circ(a[i_planet])
    else:
        v[0:3,maxidx+i_planet]=cspice.spkezr(str(i_planet),et,TrajFrame,"None","0")[0][0:3]/au
        v[3,maxidx+i_planet]=1
        orbtraj[i_planet]=traj(i_planet,et-T[i_planet]-86400,et,T[i_planet]/360)
    prj=callback.project_stuff(orbtraj[i_planet])
    orbplot[i_planet],=plt.plot(prj[0,...],prj[1,...],'w-')

prj=callback.project_stuff(v)
starplot,=plt.plot(prj[0,...],prj[1,...],'r+')

nameobj=[None]*prj.shape[1]
if True:
    for i in range(prj.shape[1]):
        if name[i]!="":
            print(i, name[i])
            nameobj[i]=plt.text(prj[0,i] if np.isfinite(prj[0,i]) else 0,
                                prj[1,i] if np.isfinite(prj[0,i]) else 0," %d %s"%(i,name[i]),
                                color='yellow' if i in goodstars else 'blue',
                                visible=np.isfinite(prj[0,i]))
axs=[]
btns=[]
def makebtn(x,y,name,f):
    ax = plt.axes([x, y, 0.05, 0.05])
    axs.append(ax)
    bx = widgets.Button(ax, name)
    btns.append(bx)
    bx.on_clicked(f)

makebtn(0.80,0.05,'-camlon',callback.camlonm)
makebtn(0.80,0.00,'-clock',callback.clockm)
makebtn(0.80,0.10,'+clock',callback.clockp)
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
makebtn(0.00,0.05,'<frame',callback.framem)
makebtn(0.10,0.05,'>frame',callback.framep)
makebtn(0.05,0.00,'fit',callback.fit)
makebtn(0.00,0.10,'<auto',callback.autom)
makebtn(0.10,0.10,'auto>',callback.autop)

plt.show()




