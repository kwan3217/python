import numpy as np

class Bezier:
    def __init__(self,r):
        """

        :param r: Iterable of numpy arrays representing the vectors.
          Elements 0 and 3 are the endpoints which the curve will touch at
          t=0 and t=1, while elements 1 and 2 are the nodes which control
          the shape.
        """
        self.r=[np.array(rn) for rn in r]
        self.c=3*(self.r[1]-self.r[0])
        self.b=3*self.r[2]-3*self.r[1]-self.c
        self.a=self.r[3]-self.r[0]-self.b-self.c
    def b3(self,t):
        return t**3
    def b2(self,t):
        return -3*t**3+3*t**2
    def b1(self,t):
        return self.b2(1-t)
    def b0(self,t):
        return self.b0(1-t)
    def eval(self,t):
        return tuple((((self.a)*t+self.b)*t+self.c)*t+self.r[0])
