"""
Powered Explicit Guidance - see the wiki
"""

import numpy as np #Use numpy math and trig functions, in case we get passed an array
import collections

def b(n,T,tau,ve):
    """
    nth moment of integrated rocket acceleration. Zeroth moment is ideal Delta-V
    bn(T)=\int_0^T t^n a(t) dt
    :param T: Time to burn in time units (SI - seconds)
    :param tau: Time needed to burn the entire vessel as fuel, given the current mass flow rate (SI - seconds)
    :param ve: Effective exhaust velocity (SI - m/s)
    :param n: moment order to return
    :return: nth moment of acceleration
    """
    if n==0:
        return -ve*np.log(1-T/tau)
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


PegExtra = collections.namedtuple("PegExtra",
                                  ["A","B","T","rdot","vq","aT", "rbar", "hT", "h", "deltah", "fr", "frT", "frdot", "fq", "fqdot", "fqdotdot",
                                   "N1", "N2", "N3", "N", "D0", "D1", "D2", "D", "deltaV","fdotr","fdotq"])
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
    omegaT = vqT / rT
    pegextra=None
    if T>10:
        #If T<10, then just coast on the time-updated constants from before
        #Iterate as requested
        for i in range(n):
            #Calculate new T
            aT = a / (1 - T / tau)  # Acceleration at end of burn
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
            # Eqn 36
            N1=deltah/rbar
            N2=ve*T*(fqdot+fqdotdot*tau)
            N3=fqdotdot*ve*T**2/2
            N=N1+N2+N3
            D0=fq
            D1=fqdot*tau
            D2=fqdotdot*tau**2
            D=D0+D1+D2
            #if N/D>0:
            #deltaV=N/D
            #else:
            deltaV=N1 #Equivalent to just thrusting horizontal
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
        # Calculate vector components of fhat
        fdotr = A + (mu / r ** 2 - omega ** 2 * r) / a
        fdotq = np.sqrt(1 - fdotr ** 2)
        if not np.isfinite(fdotq):
            print("Something happened!")
        pegextra=PegExtra(A=A,B=B,T=T,rdot=rdot,vq=vq,aT=aT,rbar=rbar,hT=hT,h=h,deltah=deltah,fr=fr,frT=frT,frdot=frdot,fq=fq,
                          fqdot=fqdot,fqdotdot=fqdotdot,N1=N1,N2=N2,N3=N3,N=N,D0=D0,D1=D1,D2=D2,D=D,deltaV=deltaV,fdotr=fdotr,fdotq=fdotq)
    else:
        # Calculate vector components of fhat
        fdotr = A + (mu / r ** 2 - omega ** 2 * r) / a
        fdotq = np.sqrt(1 - fdotr ** 2)
        if not np.isfinite(fdotq):
            print("Something happened!")
        pegextra=PegExtra(A=A,B=B,T=T,rdot=rdot,vq=vq,aT=None,rbar=None,hT=None,h=None,deltah=None,fr=None,frT=None,frdot=None,fq=None,
                          fqdot=None,fqdotdot=None,N1=None,N2=None,N3=None,N=None,D0=None,D1=None,D2=None,D=None,deltaV=None,fdotr=None,fdotq=None)
    #Return the results
    return A,B,T,fdotr,fdotq,pegextra

