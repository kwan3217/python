'''
Created on Dec 30, 2016

@author: chrisj

Functions to calculate the moments of inertia of various primitive objects.
Larger objects can be built up from these by applying the parallel axis theorem
to transform each primitive from its own center of mass to that of the
composite object, then simply summing the resulting moments.

All of these functions take mass as their last parameter, which always has a
default value of 1. This means that you can ignore the last argument, and get
a *normalized* moment of inertia. The vessel parts take advantage of this
by multiplying this normalized inertia by the mass they are already passed, so
that the mass doesn't have to be passed twice (thereby keeping the One
Definitive Source rule).
'''
import numpy as np
from math import pi

def CylinderV(r,l):
    """
    Calculate volume of right circular cylinder

    :param r: radius of cylinder
    :param l: length of cylinder
    :return: volume of cylinder
    """
    V=r**2*np.pi*l
    return V

def CylinderA(r,l):
    """
    Calculate surface area of cylinder, not including end caps

    :param r: radius of cylinder
    :param l: length of cylinder
    :return: Area of side of cylinder
    """
    A=r*2*np.pi*l
    return A

def HollowCylinderV(rOuter,rInner,l):
    V=CylinderV(rOuter,l)-CylinderV(rInner,l)
    return V

def SolidCylinderM(r,l,rho=1):
    M=CylinderV(r,l)*rho
    return M

def ThinCylinderM(r,l,rho=1):
    M=CylinderA(r,l)*rho
    return M

def HollowCylinderM(rOuter,rInner,l,rho=1):
    M=HollowCylinderV(rOuter,rInner,l)*rho
    return M

def SolidCylinderI(r, l, m=1):
    """
    Calculate inertia tensor for solid cylinder.

    Inertia tensor is relative to cylinder center of mass. Cylinder axis is
    along X axis.

    Parameters
    ----------
    r : real
        radius of cylinder
    l : real
        length of cylinder
    m : real
        Total mass of cylinder

    Returns
    -------
    I : numpy 3x3 matrix
        Inertia tensor with respect to the center of mass and these axes.

    Notes
    -----
    from https://en.wikipedia.org/wiki/List_of_moments_of_inertia
    """
    IAxis=m*r*r/2
    IDiameter=m*(3*r*r+l*l)/12
    return np.matrix([[IAxis,0,0],[0,IDiameter,0],[0,0,IDiameter]])

def ThinCylinderI(r, l, m=1):
    """
    Calculate inertia tensor for a thin-walled cylinder with no end caps.

    Inertia tensor is relative to cylinder center of mass. Cylinder axis is
    along X axis.

    Parameters
    ----------
    r : real
        radius of cylinder
    l : real
        length of cylinder
    m : real
        Total mass of cylinder

    Returns
    -------
    I : numpy 3x3 matrix
        Inertia tensor with respect to the center of mass and these axes.

    Notes
    -----
    from https://en.wikipedia.org/wiki/List_of_moments_of_inertia
    """
    IAxis=r*r
    IDiameter=(3*r*r+2*l*l)/6
    return np.matrix([[IAxis,0,0],[0,IDiameter,0],[0,0,IDiameter]])*m

def HollowCylinderI(rOuter, rInner, l, m=1):
    """
    Calculate inertia tensor for a thick-walled hollow cylinder.

    Inertia tensor is relative to cylinder center of mass. Cylinder axis is
    along X axis.

    Parameters
    ----------
    rOuter : real
        outer radius of cylinder
    rInner : real
        inner radius of cylinder (must be >=rOuter, silently gives
        non-physical results otherwise, but =rOuter is OK)
    l : real
        length of cylinder
    m : real
        Total mass of cylinder

    Returns
    -------
    I : numpy 3x3 matrix
        Inertia tensor with respect to the center of mass and these axes.

    Notes
    -----
    calculated by subtracting the inertia of a cylinder of given inner radius
    from the inertia of a cylinder of given outer radius. In the thin-shell
    limit, this calls ThinCylinderI
    """
    if rOuter==rInner:
        return ThinCylinderI(rOuter,l,m)
    OuterVolume=pi*rOuter*rOuter*l
    InnerVolume=pi*rInner*rInner*l
    TotalVolume=OuterVolume-InnerVolume
    Density=m/TotalVolume
    OuterMass=OuterVolume*Density
    InnerMass=InnerVolume*Density
    return SolidCylinderI(rOuter, l, OuterMass)-SolidCylinderI(rInner,l,InnerMass)

def SphereV(r):
    """
    Calculate volume of cylinder
    """
    return r**3*np.pi*4.0/3.0

def SphereA(r,l):
    """
    Calculate surface area of cylinder, not including end caps
    """
    A=r*2*np.pi*4
    return A

def HollowSphereV(rOuter,rInner):
    V=SphereV(rOuter)-SphereV(rInner)
    return V

def SolidSphereM(r,rho=1):
    M=SphereV(r)*rho
    return M

def ThinSphereM(r,rho=1):
    return SphereA(r)*rho

def HollowSphereM(rOuter,rInner,rho=1):
    return HollowSphereV(rOuter,rInner)*rho

def SolidSphereI(r,m=1):
    """
    Calculate inertia tensor for a solid sphere.

    Inertia tensor is relative to sphere center of mass.

    Parameters
    ----------
    r : real
        radius of sphere
    m : real
        Total mass of cylinder

    Returns
    -------
    I : numpy 3x3 matrix
        Inertia tensor with respect to the center of mass and these axes.
    """
    I=2*m*r*r/5;
    I=np.identity(3)*I
    return I

def ThinSphereI(r,m=1):
    """
    Calculate inertia tensor for a thin-walled sphere.

    Inertia tensor is relative to sphere center of mass.

    Parameters
    ----------
    m : real
        Total mass of sphere
    r : real
        radius of sphere

    Returns
    -------
    I : numpy 3x3 matrix
        Inertia tensor with respect to the center of mass and these axes.
    """
    I=2*m*r*r/3;
    I=np.identity(3)*I
    return I

def HollowSphereI(rOuter, rInner,m=1):
    """
    Calculate inertia tensor for a solid sphere.

    Inertia tensor is relative to sphere center of mass.

    Parameters
    ----------
    m : real
        Total mass of cylinder
    rOuter : real
        Outer radius of sphere
    rInner : real
        inner radius of sphere (must be >=rOuter, silently gives
        non-physical results otherwise)

    Returns
    -------
    I : numpy 3x3 matrix
        Inertia tensor with respect to the center of mass and these axes.
    """
    if rOuter==rInner:
        return ThinSphereI(rOuter,m)
    OuterVolume=4.0/3.0*pi*rOuter*rOuter*rOuter
    InnerVolume=4.0/3.0*pi*rInner*rInner*rInner
    TotalVolume=OuterVolume-InnerVolume
    Density=m/TotalVolume
    OuterMass=OuterVolume*Density
    InnerMass=InnerVolume*Density
    I=SolidSphereI(rOuter,OuterMass)-SolidSphereI(rInner,InnerMass)
    return I

def ConeA(r,h):
    """
    Calculate the area of "side" of a right circular cone, not counting the base cap
    :param r: Radius of cone base
    :param h: height of the cone
    :return: Area of the side of the cone
    """
    s=np.sqrt(h**2+r**2)
    A=pi*r*s
    return A

def ConeV(r,h):
    """
    Calculate the volume of a right circular cone
    :param r: Radius of the cone base
    :param h: Height of the cone
    :return: Volume of the cone
    """
    V=pi*r*r*h/3
    return V

def SolidConeI(r,h,m=1):
    """
    :param r:
    :param h:
    :param m:
    :return:
    """
    #From Wolfram Alpha: "moment of inertia tensor of solid cone"
    Iaxis=3*r**2/10
    Iperp=(3*a**2+2*h**2)/20


def TriangularPrismI(a, b, c,m=1):
    """
    Calculate inertia tensor for a right triangular prism.

    Parameters
    ----------
    m : real
        Total mass of prism
    a : real
        length of prism base along +x axis
    b : real
        length of prism base along +y axis
    c : real
        height of prism along +z axis

    Notes
    -----
    I is with respect to the center of mass and these axes.
    To describe a prism out -x or -y, use a negative a or b
    """
    I=np.matrix([[2*b*b+3*c*c,        a*b,           0],
                 [a*b,        2*a*a+3*c*c,           0],
                 [0,                    0, 2*a*a+2*b*b]])*m/36
    return I

def RectangularPrismI(a, b, c,m=1):
    """
    Calculate inertia tensor for a right rectangular prism

    Parameters
    ----------
    m : real
        Total mass of prism
    a : real
        length of prism along +x axis
    b : real
        length of prism along +y axis
    c : real
        length of prism along +z axis

    Notes
    -----
    I is with respect to the center of mass and these axes.
    """
    I=np.matrix([[b*b+c*c,       0,       0],
                 [      0, a*a+c*c,       0],
                 [      0,       0, a*a+b*b]])*m/12;
    return I

def MonteCarloI(rho,x0,x1,y0,y1,z0,z1,n=1000000,n_inter=1000,*arg,**kwarg):
    """
    Calculate the moment of inertia of an arbitrary object by the Monte Carlo method
    :param rho: Function that takes rho(x,y,z) and returns density at that point
    :param x0: lower x bound. Will silently give wrong answers if x0>x1
    :param x1: upper x bound
    :param y0: lower y bound. Will silently give wrong answers if y0>y1
    :param y1: upper y bound
    :param z0: lower z bound. Will silently give wrong answers if z0>z1
    :param z1: upper z bound
    :param n: Number of Monte Carlo samples to use
    :return: A tuple.
      First element scalar estimate of total mass.
      Second element is 3x3 matrix of estimate of inertia tensor.
      Third element is estimate of center of mass
    """
    import random
    m_acc=0
    I_acc=np.zeros((3,3))
    CoM_acc=np.zeros(3)
    box_v=(x1-x0)*(y1-y0)*(z1-z0)
    for i in range(int(n)):
        x=random.uniform(x0,x1)
        y=random.uniform(y0,y1)
        z=random.uniform(z0,z1)
        r=np.array((x,y,z))
        this_rho=rho(x,y,z,*arg,**kwarg)
        m_acc+=this_rho
        CoM_acc+=r*this_rho
        I_acc+=this_rho*np.array(((y**2+z**2,-x*y,-x*z),
                                  (-y*x,x**2+z**2,-y*z),
                                  (-z*x,-z*y,x**2+y**2)))
        dV=box_v/(i+1)
        m=m_acc*dV
        I=I_acc*dV
        CoM=CoM_acc*dV/m
        if(i % n_inter==0):
            print("i: %d dv: %f m: %f"%(i,dV,m))
            print("CoM: ",CoM)
            print("I: ", I)
    print("i: %d dv: %f m: %f" % (n, dV, m))
    print("CoM: ", CoM)
    print("I: ", I)
    return (m,I,CoM)

def cubeRho(x,y,z):
    return 1

def sphereRho(x,y,z,r=1):
    if(x**2+y**2+z**2<r**2):
        return 1
    return 0

def spheroidRho(x,y,z,a=1,b=1,zz0=None,zz1=None):
    if zz0 is None:
        zz0=-b
    if zz1 is None:
        zz1=b
    if z<zz0:
        return 0
    if z>zz1:
        return 0
    if(x**2/a**2+y**2/a**2+z**2/b**2<1):
        return 1
    return 0

if __name__=="__main__":
    (m_est,I_est,CoM_est)=MonteCarloI(spheroidRho,-1,1,-1,1,-1,1,n=10000000,n_inter=100000,zz0=0)
    m=8
    CoM=np.zeros(3)
    I=m*np.array(((4/6,0,0),
                  (0,4/6,0),
                  (0,0,4/6)))
    print("calc m: %f" %(m))
    print("calc I: ",I)
    print("calc CoM: ",CoM)
    print("diff m: %f" %(m-m_est))
    print("diff I: ",I-I_est)
    print("diff CoM: ",CoM-CoM_est)
