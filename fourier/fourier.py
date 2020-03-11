"""
Follow the leader, 3b1b. Draw a Fourier transform of a given path
"""

from cmath import exp

class Fourier:
    def __init__(self,coef=None,path=None):
        if coef is None:
            #Do something to trace a path into Fourier coefficients
        else:
            self.coef=coef
    def n(self,i):
        """
        Given a coefficient index, return the order n
        :param i:
        :return:
        """
        return i//2 if i%2==0 else -i//2
    def eval(self,t):
        """

        :param t: Runs from 0 to 2*pi
        :return:
        """
        result=0
        for i,c in coef:
            result=result+coef*exp(1j*t*n(i))
        return result

from cairo import *


