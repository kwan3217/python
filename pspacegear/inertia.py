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

def SphereA(r,z0=None,z1=None):
    """
    Calculate surface area of segment of sphere. Calculates for a whole sphere by default

    :param r: Radius of sphere
    :param z0: Lower bound of segment, -r by default
    :param z1: Upper bound of segment, r by default
    :return: Area of segment, not including circular end caps.
    """
    if z0 is None:
        z0=-r
    if z1 is None:
        z1=r
    A=2*np.pi*r*(z1-z0)
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

def ThinSphereI(r,m=1,rho=None,z0=None,z1=None):
    """
    Calculate inertia tensor for a thin-walled spherical segment (whole sphere by default). Segment is cut by planes
    of constant z, so the z axis remains an axis of wheel symmetry for the segment.

    Inertia tensor is relative to the center of the sphere this is a segment of. This might
    not be the center of mass, if the segment isn't symmetrical about the equator of the sphere.

    :param r: Radius of sphere
    :param m: Mass of sphere
    :param rho: Areal density of sphere. Optional, but if passed, m is ignored.
    :param z0: Lower bound of segment, bottom of sphere by default
    :param z1: Upper bound of segment, top of sphere by default
    :return: Inertia tensor with respect to the center of the sphere.
    """
    if rho is None:
        A=SphereA(r,z0,z1)
        rho=m/A
    z1int=(z1-z0)
    z3int=(z1**3-z0**3)/3
    Izz=2*np.pi*rho*r**3*(z1int-z3int)
    Ixx=Izz/2+2*np.pi*rho*r*z3int
    I=np.array((Ixx,0,0),
               (0,Ixx,0),
               (0,0,Izz))
    return I

def HollowSphereI(rOuter, rInner,m=1):
    """
    Calculate inertia tensor for a hollow thick-walled sphere.

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

def SpheroidV(a,b,z0=None,z1=None):
    if z0 is None:
        z0=-b
    if z1 is None:
        z1=b
    z1int=z1-z0
    z3int=(z1**3-z0**3)/3
    V=np.pi*a**2*(z1int-z3int/b**2)
    return V

def SpheroidI(a,b,z0=None,z1=None,m=1,rho=None):
    """
    Calculate the moment of inertia of a segment of a solid spheroid.

    This segment is a solid of revolution, bounded by an ellipse spun
    around the z axis, which is also one of the axes of the ellipse. It
    is optionally also bound by planes z=z0 and z=z1. This is useful for
    doing things like the fuel in a spheroid fuel tank.
    :param a: Equatorial radius of spheroid, may be greater than (oblate),
              equal to (sphere) or less than (prolate) the polar radius
    :param b: Polar radius of spheroid
    :param z0: Optional lower bound of integration. If included, the
               spheroid is assumed to be cut off below this level. May
               be any value between -b and b inclusive. Do not exceed
               this bound, or this function will silently give a wrong
               answer.
    :param z1: Optional upper bound of integration. If included, the
               spheroid is assumed to be cut off above this level. May
               be any value between -b and b, but must be greater than
               z0.
    :param m: mass of the spheroid segment, defaults to 1, ignored if
              rho is set
    :param rho: density of the spheroid segment. If not set, then the
                mass m and volume are used to calculate it. If it is
                set, m is ignored (its value isn't used).
    :return:
    """
    V=SpheroidV(a,b,z0=z0,z1=z1)
    #If rho has a value, then the value of m is ignored
    if rho is None:
        rho=m/V
    if z0 is None:
        z0=-b
    if z1 is None:
        z1=b
    z1int=z1-z0
    z3int=(z1**3-z0**3)/3
    z5int=(z1**5-z0**5)/5
    Izz=rho*np.pi*a**4/2*(z1int-2*z3int/b**2+z5int/b**4)
    Ixx=Izz/2+rho*np.pi*a**2*(z3int-z5int/b**2)
    return np.array(((Ixx,0,0),
                     (0,Ixx,0),
                     (0,0,Izz)))

def SpheroidCoM(a,b,z0=None,z1=None):
    V=SpheroidV(a,b,z0=z0,z1=z1)
    if z0 is None:
        z0 = -b
    if z1 is None:
        z1 = b
    z2int = (z1 ** 2 - z0 ** 2) / 2
    z4int = (z1 ** 4 - z0 ** 4) / 4
    zcom=np.pi*a**2/V*(z2int-z4int/b**2)
    return np.array((0,0,zcom))

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
    if x**2/a**2+y**2/a**2+z**2/b**2<1:
        return 1
    return 0

if __name__=="__main__":
    a=2
    b=1
    z0=-b
    z1=0
    V=SpheroidV(a,b,z0=z0,z1=z1)
    m=V
    CoM=SpheroidCoM(a,b,z0=z0,z1=z1)
    I=SpheroidI(a,b,z0=0,rho=1)
    print("calc m: %f" %(m))
    print("calc I: ",I)
    print("calc CoM: ",CoM)
    (m_est,I_est,CoM_est)=MonteCarloI(spheroidRho,-a,a,-a,a,z0,z1,n=10000000,n_inter=1000000,a=a,b=b,zz0=z0,zz1=z1)
    print("diff m: %f" %(m-m_est))
    print("diff I: ",I-I_est)
    print("diff CoM: ",CoM-CoM_est)
