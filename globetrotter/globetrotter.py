from datasystem.python_lib.log import log
log("import numpy as np")
import numpy as np
import matplotlib.image as mpimg
log("import matplotlib.pyplot as plt")
import matplotlib.pyplot as plt
from datasystem import ANC
import csv

def dms2rad(d,m,s):
    return np.radians(d+m/60+s/3600)

def ll2xyz(lat,lon):
    clat=np.cos(lat)
    clon=np.cos(lon)
    slat=np.sin(lat)
    slon=np.sin(lon)
    return np.array((clat*clon,clat*slon,slat))

def xyz2ll(xyz):
    nxyz=xyz/np.linalg.norm(xyz)
    return (np.arcsin(nxyz[2]),np.arctan2(nxyz[1],nxyz[0]))

def rq2ll(lat0,lon0,r,q):
    lat=np.arcsin(np.sin(lat0)*np.cos(r)+np.cos(lat0)*np.sin(r)*np.cos(q))
    lon=np.arctan2(np.sin(q)*np.sin(r)*np.cos(lat0),np.cos(r)-np.sin(lat0)*np.sin(lat))+lon0
    lon[lon>np.pi]-=2*np.pi
    lon[lon<-np.pi]+=2*np.pi
    return lat,lon

def ll2rq(lat0,lon0,lat,lon):
    r=np.arccos(np.sin(lat0)*np.sin(lat)+np.cos(lat0)*np.cos(lat)*np.cos(lon-lon0))
    q=np.arctan2(np.sin(lon-lon0)*np.cos(lat),np.cos(lat0)*np.sin(lat)-np.sin(lat0)*np.cos(lat)*np.cos(lon-lon0))
    try:
        q[q<0]+=2*np.pi
        q[q>2*np.pi]-=2*np.pi
    except:
        if(q<0):
            q+=2*np.pi
        elif(q>2*np.pi):
            q+=2*np.pi
    return (r,q)

def rq2xy(r,q,xsize,ysize,scl=3/np.pi,rot=0):
    x=r*np.sin(q-rot)
    x*=scl*xsize/2
    x=xsize//2+x
    y=r*np.cos(q-rot)
    y*=scl*xsize/2 #Yes, use xsize, because this makes it the same scale even if the image is not a square
    y=ysize//2-y
    try:
        x[x<0]=0
        x[x>=xsize]=xsize-1
    except:
        if(x<0):
            x=0
        elif(x>=xsize):
            x=xsize-1
    try:
        y[y<0]=0
        y[y>=ysize]=ysize-1
    except:
        if(y<0):
            y=0
        elif(y>=ysize):
            y=ysize-1
    return (x,y)

def xy2rq(x,y,xsize,ysize,scl=3/np.pi,rot=0):
    r=np.sqrt((x-xsize/2)**2+(ysize/2-y)**2)/(scl*xsize/2)
    q=np.arctan2(x-xsize/2,(ysize/2-y))+rot
    q[q<0]+=2*np.pi
    q[q>2*np.pi]-=2*np.pi
    return (r,q)
    
def main(lat0=dms2rad( 42, 21, 47),
    #Logan Airport from Wikipedia
    lon0=dms2rad(-71,- 0,-23),
    #Dubai Airport from Wikipedia
    lat1=dms2rad( 25, 15, 10),
    lon1=dms2rad( 55, 21, 52),
    xsize=2000,
    ysize=1000
             ):
    """
    Conventions:
    
    x - Horizontal pixel coordinate in final image. Ranges from 0 (left) to xsize-1 (top)
    y - Vertical   pixel coordinate in final image. ranges from 0 (top) to ysize-1 (bottom)
    img - final image. A 2D image has indices [y,x] and a 3D color image has indices [y,x,color]
    r - distance from center of projection in radians
    q - bearing from center of projection in radians, ranges from [0,2*pi). 0=north, pi/2(90deg)=east, pi(180deg)=south, 3pi/2(270deg)=west
    lat - latitude in radians, ranges from [-pi/2,pi/2]
    lon - longitude in radians, ranges from [-pi,pi)
    rot - angle that the map is rotated, such that the destination is directly to the right of the center. OR, with rot=0, no rotation is performed, which
          results in the north pole being directly above the center. Positive rot rotate the map counterclockwise
          relative to rot=0
          
    Which means that each pixel has three coordinates:
       x,y - pixel coordinate from the top
       r,q - distance and bearing from center. Bearing takes rot into account.
       lat,lon -   
    """
    xyz0=ll2xyz(lat0,lon0)
    xyz1=ll2xyz(lat1,lon1)
    (latm,lonm)=xyz2ll((xyz0+xyz1)/2)
    (_,rot)=ll2rq(latm,lonm,lat1,lon1)
    rot-=np.pi/2
    #rot=1
    #Now, use latm and lonm as the center of an azimutha equidistant map
    log("Loading Earth map")
    Map=np.flipud(mpimg.imread(ANC+"EarthMap.png").astype(np.float32))
    #Normalize map such that brightest point is 1.0 . This normalizes across
    #all channels 
    #Map/=255
    log("Loading track")
    tracklat=[]
    tracklon=[]
    with open(ANC+'track.csv') as inf:
        reader=csv.reader(inf,delimiter=',')
        for row in reader:
            tracklat.append(float(row[2]))
            tracklon.append(float(row[1]))
    #Set up the destination image
    log("Calculating projection")
    x=np.arange(0,xsize)
    y=np.arange(0,ysize)
    x.shape=(1,xsize)
    y.shape=(ysize,1)
    (r,q)=xy2rq(x,y,xsize=xsize,ysize=ysize,rot=rot)
    (lat,lon)=rq2ll(latm,lonm,r,q)
    xpix=((np.degrees(lon)+180)*(Map.shape[1]-1)/360).astype(np.int16)
    ypix=((np.degrees(lat)+90)*(Map.shape[0]-1)/180).astype(np.int16)
    xpix[xpix>Map.shape[1]-1]=Map.shape[1]-1
    xpix[xpix<0]=0
    ypix[ypix>Map.shape[0]-1]=Map.shape[0]-1
    ypix[ypix<0]=0
    img=Map[ypix,xpix,:]
    plt.imshow(img)
    (r0,q0)=ll2rq(latm,lonm,lat0,lon0)
    (r1,q1)=ll2rq(latm,lonm,lat1,lon1)
    (x,y)=rq2xy(np.array([r0,r1]),np.array([q0,q1]),xsize=xsize,ysize=ysize,rot=rot)
    plt.plot(x,y,'k')
    #for i in range(24):
    #    (r,q)=ll2rq(latm,lonm,np.radians(np.arange(-90,91)),np.radians(i*15))
    #    (x,y)=rq2xy(r,q,xsize=xsize,ysize=ysize,rot=rot)
    #    plt.plot(x,y,'b')
    #for i in range(12):
    #    (r,q)=ll2rq(latm,lonm,np.radians((i-6)*15),np.radians(np.arange(0,361)))
    #    (x,y)=rq2xy(r,q,xsize=xsize,ysize=ysize,rot=rot)
    #    plt.plot(x,y,'b')
    (r,q)=ll2rq(latm,lonm,np.radians(np.array(tracklat)),np.radians(np.array(tracklon)))
    (x,y)=rq2xy(r,q,xsize=xsize,ysize=ysize,rot=rot)
    plt.plot(x,y,'r')
    plt.axis([0,xsize,ysize,0])
    plt.axis('off')
    plt.axis('equal')    
    plt.show()
    pass

if __name__ == "__main__":
    main()
