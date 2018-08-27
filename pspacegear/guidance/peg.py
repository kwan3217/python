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
                                  ["A","B","T","rdot","vq","aT", "rbar", "hT", "h", "deltah", "fr", "frT", "fq", "fqT", "N1", "deltaV","fdotr","fdotq"])
def peg(A,B,T,dt,a,ve,r,rdot,vq,rT,rdotT,vqT,mu,n=1,warn=True):
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
    rbar = (r + rT) / 2  # Eqn 26
    hT = vqT * rT
    h = vq * r
    deltah = hT - h  # Eqn 33
    if T>10:
        #If T<10, then just coast on the time-updated constants from before
        #Iterate as requested
        for i in range(n):
            #Calculate new T
            aT = a / (1 - T / tau)  # Acceleration at end of burn
            #Approximate value and derivative of vertical component of fhat
            fr=A+(mu/r**2-omega**2*r)/a          #Eqn 22b
            frT=A+B*T+(mu/rT**2-omegaT**2*rT)/aT
            # Eqn 36
            N1=deltah/rbar
            #Do this on our own. Look at fq (newly calculated, not eqn25a) and fqT
            #to get representative fraction of thrust in downrange direction. If
            #fr and frT are opposite signs, then we went through horizontal, and some of our thrust was more effective
            #than that on either end. The representative value is then the weighted mean of fq, fqT and 1 weighted
            #twice as heavy. If they are the same sign, then we just use the mean of fq and fqT
            fq =np.sqrt(1-fr **2)
            fqT=np.sqrt(1-frT**2)
            if fr*frT>0:
                #Same sign
                fqR=(fq+fqT)/2
            else:
                #Opposite sign
                w1=2
                fqR=(fq+fqT+1*w1)/(2+w1)
            deltaV=N1/fqR
            nextT=tau*(1-np.exp(-deltaV/ve))         #Eqn 37b

            #Calculate new A and B
            kb=rdotT-rdot
            kc=rT-r-rdot*nextT
            b0=b(0,nextT,tau,ve)
            b1=b(1,nextT,tau,ve)
            c0=c(0,nextT,tau,ve)
            c1=c(1,nextT,tau,ve)
            nextB=(kc*b0-c0*kb)/(c1*b0-c0*b1)
            nextA=kb/b0-b1/b0*nextB
            aT = a / (1 - nextT / tau)  # Acceleration at end of burn
            nextfr=nextA+(mu/r**2-omega**2*r)/a          #Eqn 22b
            nextfrT=nextA+nextB*nextT+(mu/rT**2-omegaT**2*rT)/aT
            if np.abs(nextfr)>1 or np.abs(nextfrT)>1:
                if warn:
                    print("Being asked to fly an impossible trajectory: fr=%f,frT=%f"%(nextfr,nextfrT))
                    warn=False
                break
            else:
                A=nextA
                B=nextB
                T=nextT
        # Calculate vector components of fhat
        fdotr = A + (mu / r ** 2 - omega ** 2 * r) / a
        fdotq = np.sqrt(1 - fdotr ** 2)
        if not np.isfinite(fdotq):
            print("Something happened!")
        pegextra=PegExtra(A=A,B=B,T=T,rdot=rdot,vq=vq,aT=aT,rbar=rbar,hT=hT,h=h,deltah=deltah,fr=fr,frT=frT,N1=N1,fq=fq,fqT=fqT,deltaV=deltaV,fdotr=fdotr,fdotq=fdotq)
    else:
        # Calculate vector components of fhat
        fdotr = A + (mu / r ** 2 - omega ** 2 * r) / a
        fdotq = np.sqrt(1 - fdotr ** 2)
        if not np.isfinite(fdotq):
            print("Something happened!")
        pegextra=PegExtra(A=A,B=B,T=T,rdot=rdot,vq=vq,N1=None,aT=None,rbar=None,hT=None,h=None,deltah=None,fr=None,frT=None,fq=None,fqT=None,deltaV=None,fdotr=fdotr,fdotq=fdotq)
    #Return the results
    return A,B,T,fdotr,fdotq,pegextra,warn

def peg_startup(a,ve,r,rdot,vq,rT,rdotT,vqT,mu,n=10):
    """
    Estimate starting values for steering constants A, B, T

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
      Initial guess for A constant
      Initial guess for B constant
      Initial guess for T constant
    """

    #Intermediate values which are determined from given parameters
    tau=ve/a #Time to burn the entire rocket as if it was fuel
    omega=vq/r
    omegaT = vqT / rT
    #If T<10, then just coast on the time-updated constants from before
    #Iterate as requested

    #Assume that the initial thrust vector is parallel to the current motion
    fr=rdot/np.sqrt(rdot**2+vq**2)
    #Calculate initial A from that
    A=fr-(mu/r**2-omega**2*r)/a
    #Assume that the final thrust vector is parallel to the target motion
    frT=rdotT/np.sqrt(rdotT**2+vqT**2)
    #Calculate T from that
    rbar = (r + rT) / 2  # Eqn 26
    hT = vqT * rT
    h = vq * r
    deltah = hT - h  # Eqn 33
    N1 = deltah / rbar
    # Do this on our own. Look at fq (newly calculated, not eqn25a) and fqT
    # to get representative fraction of thrust in downrange direction. If
    # fr and frT are opposite signs, then we went through horizontal, and some of our thrust was more effective
    # than that on either end. The representative value is then the weighted mean of fq, fqT and 1 weighted
    # twice as heavy. If they are the same sign, then we just use the mean of fq and fqT
    fq = np.sqrt(1 - fr ** 2)
    fqT = np.sqrt(1 - frT ** 2)
    if fr * frT > 0:
        # Same sign
        fqR = (fq + fqT) / 2
    else:
        # Opposite sign
        w1 = 2
        fqR = (fq + fqT + 1 * w1) / (2 + w1)
    deltaV = N1 / fqR
    T = tau * (1 - np.exp(-deltaV / ve))  # Eqn 37b
    #Now that we have A and T, calculate B to hit assumed frT
    aT = a / (1 - T / tau)  # Acceleration at end of burn
    B=(frT-A- (mu / rT ** 2 - omegaT ** 2 * rT) / aT)/T

    A,B,T,_,_,_,_=peg(A=A, B=B, T=T, dt=0, a=a, ve=ve, r=r, rdot=rdot, vq=vq, rT=rT, rdotT=rdotT, vqT=vqT, mu=mu, n=n, warn=False)

    return A,B,T

if __name__=="__main__":
    a    =      3.843015
    ve   =   4417.895825
    Re   =6378137
    mu   = 398600.4415e9
    rdot0=    750
    rdot1=   1500
    drdot=     10
    vq0  =   4000
    vq1  =   5000
    dvq  =    100
    h0   = 130000
    h1   = 200000
    dh   =   1000
    rT   = 200000+Re
    rdotT=   -100
    vqT  =   7783.691397
    nh   =len(range(   h0,   h1,dh   ))
    nrdot=len(range(rdot0,rdot1,drdot))
    nvq  =len(range(  vq0,  vq1,dvq  ))
    print(nh,nrdot,nvq,nh*nrdot*nvq)
    A=np.zeros((nh,nrdot,nvq))
    B=np.zeros((nh,nrdot,nvq))
    T=np.zeros((nh,nrdot,nvq))
    for ih,h in enumerate(range(h0,h1,dh)):
        print(ih,h)
        for irdot,rdot in enumerate(range(rdot0,rdot1,drdot)):
            for ivq,vq in enumerate(range(vq0,vq1,dvq)):
                A[ih, irdot, ivq],B[ih,irdot,ivq],T[ih,irdot,ivq]=peg_startup(a=a,ve=ve,mu=mu,r=h+Re,rdot=rdot,vq=vq,rT=rT,rdotT=rdotT,vqT=vqT)
    pass