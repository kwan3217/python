import matplotlib.pyplot as plt
import matplotlib.lines as lines
import matplotlib.patches as patches
from bezier import Bezier
import pytest
import numpy as np

class PictureBox():
    def __init__(self,figure,w,h,dpi=100,**kwargs):
        self.figure=figure
        self.w=w
        self.h=h
        plt.figure(figure,figsize=(w/dpi,h/dpi),dpi=dpi,**kwargs)
    def stroke(self,xdata,ydata,**kwargs):
        """
        :param xdata: numpy array of x coordinates
        :param ydata: numpy array of y coordinates
        :param kwargs:
          Passed to the artist constructor. Consider adding things like "color" etc.
        :return: None
        """
        fig=plt.figure(self.figure)
        l = lines.Line2D(xdata/self.w, 1-ydata/self.h, transform=fig.transFigure, figure=fig, **kwargs)
        fig.lines.extend([l])
        plt.pause(0.001)
    def fill(self,xdata,ydata,**kwargs):
        """
        :param path: Iterable of tuples, passed to translate_path
        :param kwargs:
          Passed to the artist constructor. Consider adding things like "color" etc.
        :return: None
        """
        fig=plt.figure(self.figure)
        l = patches.Polygon(np.array([xdata/self.w, 1-ydata/self.h]).T, transform=fig.transFigure, figure=fig, **kwargs)
        fig.lines.extend([l])
        plt.pause(0.001)
    def line(self,x0,y0,x1,y1,**kwargs):
        self.stroke(np.array([x0,x1]),np.array([y0,y1]),**kwargs)
    def savepng(self,oufn):
        plt.figure(self.figure)
        plt.savefig(oufn)
    def clear(self):
        plt.figure(self.figure)
        plt.clf()

def test_PictureBox():
    pb=PictureBox("PictureBox",1280,720)
    pb.stroke(np.array([0,200]),np.array([0,90]),color="#ff0000")
    pb.fill(np.array([100,200,200,100]),np.array([100,100,200,200]),color="#ff8000")
    pb.savepng("test_PictureBox1.png")
    pb.clear()
    pb.fill(np.array([150,250,250,150]),np.array([150,150,250,250]),color="#ffff00")
    pb.savepng("test_PictureBox2.png")
    plt.show()

if __name__=="__main__":
    test_PictureBox()


