"""
Generate and draw path objects
"""
from bezier import Bezier
import numpy as np

class Path:
    def __init__(self,x0,y0):
        self.xdata=[x0]
        self.ydata=[y0]
        self.xm2=None
        self.ym2=None
    def lineto(self,x1,y1):
        self.xm2=self.xdata[-1]
        self.ym2=self.ydata[-1]
        self.xdata.append(x1)
        self.ydata.append(y1)
    def curveto(self,x1,y1,x2,y2,x3,y3):
        self.xm2=x2
        self.ym2=y2
        x0=self.xdata[-1]
        y0=self.ydata[-1]
        B = Bezier(((x0,y0), (x1, y1), (x2, y2), (x3, y3)))
        for i in range(1, 101):
            t = i / 100
            p = B.eval(t)
            self.xdata.append(p[0])
            self.ydata.append(p[1])
    def smoothcurveto(self,x2,y2,x3,y3):
        x1=self.xdata[-1]+(self.xdata[-1]-self.xm2)
        y1=self.ydata[-1]+(self.ydata[-1]-self.ym2)
        self.curveto(x1,y1,x2,y2,x3,y3)
    def pathto(self,path):
        """
        :param path: Iterable of tuples
          * Element 0 must be a two-element tuple of numbers, indicating the start point x0,y0
          * Successive elements must be either:
             - two-element (x1,y1) line segment. Initial point is final point of previous segment.
             - four-element (x2,y2,x3,y3) smooth Bezier curve. Initial point is final point of previous segment,
               and initial node is mirror reflection of initial point of previous line segment or final node of previous
               curve segment.
             - six-element (x1,y1,x2,y2,x3,y3) Bezier curve. Initial point is final point of previous segment
        :return: xdata and ydata suitable for Line2D constructor
        """
        self.xdata = [path[0][0]]
        self.ydata = [path[0][1]]
        # Collect segments -- break each Bezier curve into 100 segments.
        for pathseg in path[1:]:
            if len(pathseg) == 2:
                self.lineto(*pathseg)
            elif len(pathseg) == 4:
                self.smoothcurveto(*pathseg)
            elif len(pathseg) == 6:
                self.curveto(*pathseg)
            else:
                raise ValueError("Don't know what to do with path data of length %d -- %s"%(len(pathseg),str(pathseg)))
    def getpathdata(self,sx=1,sy=1,dx=0,dy=0):
        xdata=np.array(self.xdata)
        ydata=np.array(self.ydata)
        #Scale data around initial point (more useful for fonts)
        xdata=(xdata-xdata[0])*sx+xdata[0]+dx
        ydata=(ydata-ydata[0])*sy+ydata[0]+dy
        return xdata,ydata
    def stroke(self,pb,sx=1,sy=1,dx=0,dy=0,**kwags):
        xdata,ydata=self.getpathdata(sx=sx,sy=sy,dx=dx,dy=dy)

    def fill(self,pb,dx=0,dy=0,**kwargs):
        pass