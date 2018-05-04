"""
Powered Explicit Guidance - see the wiki
"""

import numpy as np #Use numpy math and trig functions, in case we get passed an array
import numpy.log as ln

def b(n,T,tau,ve):
    """
    nth moment of integrated rocket acceleration. Zeroth moment is effectively Delta-V
    bn(T)=\int_0^T t^n a(t) dt
    :param T: Time to burn in time units (SI - seconds)
    :param tau: Time needed to burn the entire vessel as fuel, given the current mass flow rate (SI - seconds)
    :param ve: Effective exhaust velocity (SI - m/s)
    :param n: moment order to return
    :return: nth moment of acceleration
    """
    if n==0:
        return -ve*ln(1-T/tau)
    else:
        return b(n-1,T,tau,ve)*tau-(ve*T**n)/n

def c(n,T,tau,ve):
    """
    nth moment of double-integrated rocket acceleration. Zeroth moment is ideal distance moved during acceleration.
    bn(T)=\int_0^T \int_0^t s^n a(s) ds dt
    :param T: Time to burn in time units (SI - seconds)
    :param tau: Time needed to burn the entire vessel as fuel, given the current mass flow rate (SI - seconds)
    :param ve: Effective exhaust velocity (SI - m/s)
    :param n: moment order to return
    :return: nth moment of acceleration
    """
    if n==0:
        return b(0,T,tau,ve)*T-b(1,T,tau,ve)
    else:
        return c(n-1,T,tau,ve)*tau-(ve*T**(n+1))/(n*(n+1))

def peg(A,B,T,dt,a,ve,r,rdot,vq,rT,rdotT,vqT,mu,n=1):
    """
    Execute Powered Explicit Guidance

    :param A: A steering constant, unitless
    :param B: B steering constant, in inverse time units (SI - 1/s)
    :param T: Time to go in time units (SI - s)
    :param dt: Time since last execution of PEG in time units (SI - s)
    :param a: Current rocket acceleration
    :param ve: Current effective exhaust velocity
    :param r: Current distance from center of planet
    :param rdot: Current vertical speed
    :param vq: Current horizontal speed
    :param rT: Target distance from center of planet
    :param rdotT: Target vertical speed
    :param vqT: Target horizontal speed
    :param mu: Gravitational parameter for planet
    :param n: Number of times to run to convergence
    :return: A tuple
      New A constant
      New B constant
      new T constant
      fdotr - vertical component of direction of thrust
      fdotq - horizontal component of direction of thrust
    """
    #Update steering constants with time
    A+=B*dt
    T-=dt

    #Intermediate values which are determined from given parameters
    tau=ve/a #Time to burn the entire rocket as if it was fuel
    omega=vq/r
    omegaT=vqT/rT
    aT=a/(1-T/tau) #Acceleration at end of burn

    #Iterate as requested
    for i in range(n):
        #Calculate new T
        rbar=(r+rT)/2 #Eqn 26
        hT=vqT*rT
        h=vq*r
        deltah=hT-h #Eqn 33
        fh=0 #No yaw steering for now
        fhdot=0
        #Approximate value and derivative of vertical component of fhat
        fr=A+(mu/r**2-omega**2*r)/a          #Eqn 22b
        frT=A+B*T+(mu/rT**2-omegaT**2*rT)/aT
        frdot=(frT-fr)/T                     #Eqn 22c corrected
        #Approximate value and derivatives of downrange component of fhat
        fq=1-fr**2/2-fh**2/2                 #Eqn 25a
        fqdot=-fr*frdot-fh*fhdot             #Eqn 25b
        fqdotdot=-frdot**2/2-fhdot**2/2      #Eqn 25c
        deltaV=(deltah/rbar+ve*T*(fqdot+fqdotdot*tau)+fqdotdot*ve*T**2/2)(fq+fqdot*tau+fqdotdot*tau**2) #Eqn 36
        T=tau*(1-np.exp(-deltaV/ve))         #Eqn 37b

        #Calculate new A and B
        kb=rdotT-rdot
        kc=rT-r-rdot*T
        b0=b(0,T,tau,ve)
        b1=b(1,T,tau,ve)
        c0=c(0,T,tau,ve)
        c1=c(1,T,tau,ve)
        B=(kc*b0-c0*kb)/(c1*b0-c0*b1)
        A=kb/b0-b1/b0*B

    #Calculate vector components of fhat
    fdotr=A+(mu/r**2-omega**2*r)/a
    fdotq=np.sqrt(1-fdotr**2)

    #Return the results
    return A,B,T,fdotr,fdotq

