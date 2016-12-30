#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 29 14:48:24 2016

@author: chrisj
"""

import numpy as np
from inertia import *

class mass:
    """
    Describes a mass element of a vessel.
    
    This class is usable by itself for structural elements and other
    things such as fuel tank walls. The fuel itself should be based on a 
    subclass which knows how to calculate its own mass properties and be
    updated.
    
    Parameters
    ----------
    m : real
      Mass of this element in global mass units, usually kg
    I : 3x3 numpy matrix
      Normaized inertia tensor about center of mass of mass element, in global 
      length units, usually m. This is the actual inertia tensor divided by the
      mass of the object
    CoM : 3-element numpy vector
      Center of mass of this element in station coordinates
    
    Notes
    -----
    Use the fields of this class directly. Subclasses may recalculate these
    fields when timestep is called.
    """
    def __init__(self,m,I,CoM):
        self.CoM=CoM
        self.m=m
        self.I=I*m
        
class vessel:
    """
    Describes an independent vessel, a single rigid body which responds
    to dynamics, as opposed to purely kinematics like would come from 
    a Spice kernel.
    
    A vessel consists of:
      *A state vector which describes the position and orientation of the
       vessel
      *A list of actuators, each of which implements a non-gravitational
       force on the vessel
      *A list of masses, each of which is used to calculate the center of mass
       and moment of inertia of the vessel
       
    """
    E3=np.identity(3)
    
    def __init__(self,actuators,masses):
        self.actuators=actuators
        self.masses=masses
        
    def timestep(self,t,delta_t):
        """
        Update the properties of the vessel. This is used to do such things as:
            *drain fuel from fuel tanks
            *calculate actuator control values, say from an autopilot
            
        Mostly this is delegated to the lists of actuators and masses that
        this vessel holds.
        
        Parameters
        ----------
        t : real
            range time in seconds. Vessel properties should be updated to be 
            valid for this point in time, in other words, this is the end of
            the time step
        delta_t : real
            length of time step in seconds. This is used to allow such things
            as subtracting fuel from propellant tanks given their engines'
            current mass flow rate. This function must be prepared for delta_t
            to be zero or negative, a
            
        Returns
        -------
        None
        
        Notes
        -----
        This function must be prepared to be called during a Runge-Kutta
        substep, so it must be able to handle a zero or negative delta_t. In
        most cases, a zero delta_t will just be a no-op.
        """
        pass
        
    def inertia(self):
        """
        Get the mass and moment of inertia matrix about the center of mass
        
        Returns
        -------
        m : real
            mass of vessel in global mass units, usually kg
        I : 3x3 numpy matrix
            Inertia tensor in global mass and length units, usually kg and m
        CoM : 3-element numpy vector
            Center of mass in station coordinates (useful for drawing the 
            vessel but not needed in dynamics)
        """
        #vessel mass.
        m=0
        CoM=np.array([0,0,0])
        for mass in self.masses:
            m=m+mass.m
            CoM=CoM+np.array((mass.m,))*mass.CoM
        #Divide the weighted sum of center of masses by the total mass to get
        #the actual center of mass
        CoM=CoM/m
        #Now that we know where the total center of mass is, use the parallel
        #axis theorem to figure the moment of inertia of each element
        #relative to that center of mass
        #https://en.wikipedia.org/wiki/Parallel_axis_theorem#Tensor_generalization
        I=np.matrix([[0,0,0],[0,0,0],[0,0,0]])
        for mass in self.masses:
            R=mass.CoM-CoM
            I=I+mass.I+m*(np.inner(R,R)*vessel.E3-np.outer(R,R))
        return (m,I,CoM)
    def actuate(self):
        """
        Get the total non-gravitational force and moment relative to the center
        of mass of the vessel
        
        Returns
        -------
        f : 3-element numpy vector
            Force on center of mass in global mass, length, and time units,
            usually N
        M : 3-element numpy vector
            Moment about center of mass in global mass, length, and time units,
            usually N-m
        """
        pass
    
v=vessel([],[mass(100,SolidSphereI(100,1),np.array([-1,0,0])),
             mass(100,ThinSphereI (100,1),np.array([ 1,0,0]))])
print(v.inertia())

