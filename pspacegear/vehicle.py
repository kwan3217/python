#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 29 14:48:24 2016

@author: chrisj
"""

import numpy as np
from inertia import SolidSphereI, ThinSphereI
import math

class Mass:
    """
    Describes a mass element of a vehicle.
    
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
        """
        Initialize a mass with constant mass properties.

        :param m:  mass of element in kg
        :param I:   Normalized inertia tensor about center of mass in m^2.
                    This is actual inertia tensor divided by mass
        :param CoM: Center of mass in station coordinates
        """
        self.CoM=CoM
        self.m=m
        self.I=I*m
    def timestep(self,t,delta_t):
        """
        Placeholder virtual method
        :param t: Range time at end of step
        :param delta_t: Time step size
        """
        pass
    def inertia(self):
        return (self.m,self.I,self.CoM)

class Propellant(Mass):
    """
    Usable by itself as a point propellant resource or one with constant normalized
    MoI and CoM (say a gaseous tank). Subclasses will define time-variable
    moment of inertia and/or center of mass.
    """
    def __init__(self,m_full,I,CoM,propType,initialLoad=1.0):
        """

        :param M: Fully-loaded mass of propellant resource
        :param I: Normalized moment of inertia of full propellant resource
        :param CoM: Center of mass of full propellant resource in station coordinates
        :param propType: Object describing what kind of propellant this stores
        :param initialLoad: Initial fraction of propellant loaded. Default is fully loaded.
        """
        Mass.__init__(m_full*initialLoad,I,CoM)
        self.m_full=m_full
        self.propType=propType
    def demand(self,mdot,t,delta_t):
        """
        Demand a certain amount of propellant from the propellant resource. Each
        actuator that uses a resource will call this during its own timestep.
        Note that calling this method represents actual consumption/unloading/loading
        of the propellant resource, and therefore directly and immediately
        changes the mass properties.

        :param mdot: Mass flow rate demanded. Negative flow rate is the usual case
                     for an engine consuming propellant. Positive flow rate is such
                     as (re)filling the resource.
        :param t:    Range time at end of time step
        :param delta_t: Time step size
        :return: Fraction of demand that can be met. Will usually be 1.0 or 0.0,
                 but may be less. Might be less for a pressure-fed engine where
                 available flow rate drops with quantity, or might be fractional
                 during the last time step that a resource has fuel. For instance,
                 if there is 1kg of fuel left, and 2kg/s over 1 second is demanded,
                 only 0.5 of that demand can be met.
        """
        m_full=self.M/self.load #Calculate a full load
        delta_m=mdot*delta_t
        new_m_requested=self.M+delta_m
        if new_m_requested<0:
            result=self.M/delta_m
            self.M=0
        if new_m_requested>m_full:
            result=(m_full-self.M)/delta_m
            self.M=m_full
        self.load=self.M/m_full
        #A non-pointmass subclass would call this method, recalculate I and CoM based on
        #the new mass, then return the result of this method.
        return result
    def timestep(self,t,delta_t):
        """
        Adjust mass properties caused independent of propellant consumption. Most
        propellant resources are inert, but maybe this is venting, leakage, etc.
        :param t: Range time at end of time step
        :param delta_t: Time step size
        """
        Mass.timestep(t,delta_t)

class Substance:
    """
    Describe a substance which occupies space and has mass.
    """
    def __init__(self,name,phase=2,density=math.nan):
        """
        :param name: Name of substance, used to match up substance
                     consumption and supply. IE if a propellant resource
                     supplies a substance called "RP1" and an engine
                     demands a substance called "RP1" then life is good.
                     Conversely, if the engine demands "LH2" then we can
                     detect this and throw the appropriate error.
        :param density: Density of the substance in mass units per cubic length unit (kg/m^3 for SI)
                        Gas doesn't have a definite density, so it's density must be ignored.
        :param phase:   Phase of the substance. 0 represents solid, 1 represents liquid, 2 represents gas.
                        At the moment there is no practical difference between solid and liquid. The
                        critical property shared by solid and liquid is that density is a meaningful
                        concept for the substance itself, relatively independent of its storage
                        conditions.
        """
        self.name=name
        self.density=density
        self.phase=phase
Al  =Substance("Al" ,phase=0,density=2700) #Solid Aluminum, intended to be used for structures
RP1 =Substance("RP1",phase=1,density= 810) #Room-temperature Rocket Propellant 1
dRP1=Substance("RP1",phase=1,density= 810  *1.04) #Chilled RP1, used in Falcon 9. Reported 2.5-4% increase in density, but compatible with any consumer of RP1
#LOX data from https://www.gpo.gov/fdsys/pkg/GOVPUB-C13-26d428ad4ca587866a90da5f71b4a727/pdf/GOVPUB-C13-26d428ad4ca587866a90da5f71b4a727.pdf
LOX =Substance("LOX",phase=1,density=1141.6)      #"Normal" temperature liquid oxygen, 1 bar at boiling point, 90K
dLOX=Substance("LOX",phase=1,density=1254.8)      #Chilled LOX. SpaceX reports using LOX at 66K
H2O =Substance("H2O",phase=1,density=1000) #Water, according to the intended definition of the kilogram
#LH2 data from http://webbook.nist.gov/cgi/fluid.cgi?Action=Load&ID=C1333740&Type=SatT&Digits=5&PLow=.5&PHigh=1.5&PInc=.1&RefState=DEF&TUnit=K&PUnit=atm&DUnit=kg/m3&HUnit=kJ/mol&WUnit=m/s&VisUnit=uPa*s&STUnit=N/m
LH2 =Substance("LH2",phase=1,density=  70.85)  #Liquid hydrogen at boiling at 1 bar, 20K
#Solid propellant
AP  =Substance("AP",phase=0,density=1950) #Ammonium Perchlorate, according to Wikipedia
PBAN=Substance("PBAN",phase=1,density= 930)  #PBAN in liquid form, before casting
#The following is intended to calculate the density of AP-PBAN-Al composite solid propellant,
#as used in the Space Shuttle. Unfortunately, you can't do this, since components can change
#density as things are cast. That makes this only a reasonable justifiable estimate.
APwt=69.06
Alwt=16
PBwt=12.04
#Neglect the epoxy and iron oxide
SRF =Substance("SRF",phase=0,density=(AP.density*APwt+PBAN.density*PBwt+Al.density*Alwt)/(APwt+PBwt+Alwt))
print(SRF.name,SRF.phase,SRF.density)

class Actuator:
    def __init__(self,CoA):
        """

        :param CoA: Initial center of action in station coordinates
        """
        self.CoA=CoA
    def control(self,u):
        """
        Set the control input vector u
        :param u: Control input vector. The number and meanings of the components
                  depend on the type of actuator. It might be a single element for
                  a throttleable but non-steerable engine, elements for aileron,
                  elevator, and rudder deflection, throttle and pointing vector
                  for a steerable engine, etc.
        """
        self.u=u
    def actuate(self,x,t,delta_t):
        """
        Calculate the force and moment generated by this actuator
        :param x: State vector, included so that conversions from world to body can be done
        :param t: range time
        :param delta_t: time step size
        :return: A tuple with the following items:
          f: Force vector on the center of mass in the reference frame
          M: Moment about the center of mass in the body frame.

        If the actuator applies a linear force off-center, it is up to this function
        to calculate the moment and add it to any pure moment applied by the actuator directly
        """
        f=np.zeros(3)
        M=np.zeros(3)
        return (f,M)

class Rocket(Actuator):
    def __init__(self,name,CoA,fhat,throttle=1.0,mdot=None,f0=None,f1=None,ve0=None,ve1=None,prop=None,mix=None):
        """

        :param CoA:
        :param fhat:
        :param throttle:
        :param mdot: Total mass flow rate of engine in kg/s. Should be a positive number.
        :param f0:   Vacuum thrust in N. This is the thrust when the ambient pressure is 0.
        :param f1:   Sea-level thrust in N. This is the thrust when the ambient pressure is 101325Pa.
        :param ve0:  Vacuum effective exhaust velocity in m/s
        :param ve1:  Sea-level effective exhaust velocity in m/s

        ve=f/mdot
        """
        self.name=name
        #Primary engine performance parameters are mdot, ve0, and ve1.
        if mdot is None:
            #Don't have mdot, so must have either two thrusts or two ves
            if ve0 is None:
                #Given both thrusts and sea-level ve
                #Figure mdot and vacuum ve
                mdot=f1/ve1
                f0=mdot*ve0
            elif ve1 is None:
                #Given both thrusts and vacuum ve
                #Figure mdot and sea-level ve
                mdot=f0/ve0
                ve1=f1/mdot
            elif f0 is None:
                #Given both ve and sea-level thrust
                #Figure mdot and vacuum thrust
                mdot=f1/ve1
                f0=mdot*ve0
            elif f1 is None:
                #Given both ve and vacuum thrust
                #Figure mdot and sea-level thrust
                mdot=f0/ve0
                f1=mdot*ve1
            else:
                mdot1=f1/ve1
                mdot0=f0/ve0
                print("Engine is overdetermined, mdot0=",mdot0," mdot1=",mdot1)
                mdot=mdot0
        else:
            if f0 is None:
                #Have mdot but not f0, must have ve0
                f0=mdot*ve0
            if f1 is None:
                #Have mdot but not f1, must have ve1
                f1=mdot*ve1
            if ve0 is None:
                #Have mdot but not ve0, must have f0
                ve0=f0/mdot
            if ve1 is None:
                #Have mdot but not ve1, must have f1
                ve1=mdot*f1
        #We've finally figured out the right parameters, so keep them
        #with this class:
        self.mdot=mdot
        self.ve0=ve0
        self.ve1=ve1
        self.f0=self.mdot*self.ve0
        self.f1=self.mdot*self.ve1

#Merlin 1D engine. Stats from http://www.spacex.com/falcon9
M1D=Rocket("M1D",np.array([0,0,0]),np.array([1,0,0]),f0=8227000/9,f1=7607000/9,ve0=3050,ve1=2770)
print(M1D.mdot,M1D.ve0,M1D.ve1,M1D.f0,M1D.f1)
#Merlin 1D Vacuum engine.
M1Dvac=Rocket("M1Dvac",np.array([0,0,0]),np.array([1,0,0]),f0=934000,f1=0,ve0=348*9.80665)
print(M1Dvac.mdot,M1Dvac.ve0,M1Dvac.ve1,M1Dvac.f0,M1Dvac.f1)

class vehicle(Mass):
    """
    Describes an independent vehicle, a single rigid body which responds
    to dynamics, as opposed to purely kinematics like would come from 
    a Spice kernel.
    
    A vehicle consists of:
      *A state vector which describes the position and orientation of the
       vehicle
      *A list of actuators, each of which implements a non-gravitational
       force on the vehicle
      *A list of masses, each of which is used to calculate the center of mass
       and moment of inertia of the vehicle
       
    """
    E3=np.identity(3) #Identity 3x3 matrix
    
    def __init__(self,actuators,masses):
        self.actuators=actuators
        self.masses=masses
        
    def timestep(self,t,delta_t):
        """
        Update the properties of the vehicle. This is used to do such things as:
            *drain fuel from fuel tanks
            *calculate actuator control values, say from an autopilot
            
        Mostly this is delegated to the lists of actuators and masses that
        this vehicle holds.
        
        Parameters
        ----------
        t : real
            range time in seconds. vehicle properties should be updated to be 
            valid for this point in time, in other words, this is the end of
            the time step
        delta_t : real
            length of time step in seconds. This is used to allow such things
            as subtracting fuel from propellant tanks given their engines'
            current mass flow rate.
            
        Returns
        -------
        None
        
        Notes
        -----
        This function must be prepared to be called during a Runge-Kutta
        substep, so it must be able to handle a zero or negative delta_t. In
        most cases, a zero delta_t will just be a no-op.

        Example use during an RK4 timestep. Range time is 10s at the beginning of the step. Step size
        is 1s, so 0.5s for substeps.

        1) Run timestep(10.0,0.0) to set properties for first RK4 step
        2) Run timestep(10.5,0.5) to set properties for second step
        3) Run timestep(10.5,0.0) to set properties for third step
        4) Run timestep(11.0,0.5) to set properties for fourth step
        """
        (F,M)=actuate(t,delta_t)
        for mass in self.masses:
            mass.timestep(t,delta_t)

    def inertia(self):
        """
        Get the mass and moment of inertia matrix about the center of mass
        
        Returns
        -------
        m : real
            mass of vehicle in global mass units, usually kg
        I : 3x3 numpy matrix
            Inertia tensor about center of mass in global mass and length
            units, usually kg and m
        CoM : 3-element numpy vector
            Center of mass in station coordinates (useful for drawing the 
            vehicle but not needed in dynamics)
        """
        #vehicle mass in station coordinates
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
            I=I+mass.I*mass.m+m*(np.inner(R,R)*vehicle.E3-np.outer(R,R))
        return (m,I,CoM)
    def actuate(self,t,delta_t):
        """
        Get the total non-gravitational force and moment relative to the center
        of mass of the vehicle
        
        Returns
        -------
        f : 3-element numpy vector
            Force on center of mass in body frame in global mass, length, and 
            time units, usually N
        M : 3-element numpy vector
            Moment about center of mass in body frame in global mass, length, 
            and time units, usually N-m
        """
        f=np.zeros(3)
        M=np.zeros(3)
        for actuator in self.actuators:
            (this_f,this_M)=actuator.actuate(t,delta_t)
            f+=this_f
            m+=this_M
        return (f,M)
    
v=vehicle([],[Mass(1,SolidSphereI(1),np.array([-1,0,0])),
             Mass(1,ThinSphereI (1),np.array([ 1,0,0]))])
print(v.inertia())

