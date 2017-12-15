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
    Calculate volume of cylinder
    """
    return r**2*np.pi*l

def CylinderA(r,l):
    """
    Calculate surface area of cylinder, not including end caps
    """
    return r*2*np.pi*l

def HollowCylinderV(rOuter,rInner,l):
    return CylinderV(rOuter,l)-CylinderV(rInner,l)

def SolidCylinderM(r,l,rho=1):
    return CylinderV(r,l)*rho

def ThinCylinderM(r,l,rho=1):
    return CylinderA(r,l)*rho

def HollowCylinderM(rOuter,rInner,l,rho=1):
    return HollowCylinderV(rOuter,rInner,l)*rho

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
    return r*2*np.pi*4

def HollowSphereV(rOuter,rInner):
    return SphereV(rOuter)-SphereV(rInner)

def SolidSphereM(r,rho=1):
    return SphereV(r)*rho

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
    return np.identity(3)*I

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
    return np.identity(3)*I

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
    return SolidSphereI(rOuter,OuterMass)-SolidSphereI(rInner,InnerMass)

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
    return np.matrix([[2*b*b+3*c*c,        a*b,           0],
                      [a*b,        2*a*a+3*c*c,           0],
                      [0,                    0, 2*a*a+2*b*b]])*m/36

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
    return np.matrix([[b*b+c*c,       0,       0],
                      [      0, a*a+c*c,       0],
                      [      0,       0, a*a+b*b]])*m/12;
