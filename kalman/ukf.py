"""
Implementation of the K2 project in Python.
https://omoikane.kwansystems.org/wiki/index.php/Kalman_Filter/Implementation

This implements the Unscented Kalman filter, pretty much as described on the library page, except with the following
differences:

1) Since we are using python, we use numpy instead of Eigen
2) Since passing a function is easy, we use it instead of object orientation.
"""

import numpy as np
import scipy

def update(t0=None, t1=None, x0=None, P0=None, Q=None, f=None, z=None, R=None, h=None,**kwargs):
    """
    Use the unscented Kalman filter to do one time and measurement update for one measurement.
    :param t0: Time that a priori state estimate is valid
    :param t1: Time of measurement  
    :param x0: a priori state estimate, n element 1D array
    :param P0: a priori state covariance, nxn element 2D array (matrix)
    :param Q:  process noise covariance, must match the size of P
    :param f:  State propagation function. Must be of the form:
        def f(t0=None,t1=None,x0=None,...other named parameters from **kwargs...,**kwargs (to soak up the ones we don't care about))
            #:param t0: Time that initial state is valid
            #:param t1: Time to propagate to
            #:param x0: State at t0
            #:return: State at t1
            ... #Perform whatever calculations are necessary to do the time update...
            return x1
    :param z: Measurement vector at t1, m element 1D array
    :param R: Measurement noise covariance, mxm element 2D array (matrix)
    :param h: Measurement function. Must be of the form:
        def h(t=None,x=None,...other named parameters from update **kwargs...,**kwargs)
            #:param t: Time of measurement
            #:param x: State vector
            #:return: Measurement vector
            ... #Perform whatever calculations are necessary to construct the measurement...
            return z
    :param kwargs: 
    :return: A tuple with the following elements:
      xp: a posteriori state estimate, updated both in time and with the measurement
      Pp: a posteriori state covariance, likewise updated
    """
    W0=1.0/3.0             #Weighting of center sigma point, Time Update step 1
    n=x0.size              #Number of elements in state vector
    Wj=(1-W0)/(2*n)        #Weighting of other sigma points, Time Update step 1
    Sigma=np.linalg.cholesky(P0*n/(1-W0)) #Sigma matrix, used to create sigma points, Time Update step 1
    xxm=[]                #<Sigma points of time-updated estimate
    zzh=[]                #<Sigma points of predicted measurement
    xxm.append(f(t0=t0,t1=t1,x0=x0,**kwargs)) #Propagate center sigma point, Time Update step 2
    zzh.append(h(x=xxm[0],t=t1,**kwargs))     #Predict measurement for center sigma point, Measurement Update step 1
    xm=W0*xxm[0]          #<Time-updated estimate, initial accumulator
    zh=W0*zzh[0]          #<Predicted measurement, initial accumulator
    #VectorXd xx;             #<Sigma points of un-updated estimate, needed because of flaw in Eigen
    for j in range(n):
        xx=x0+Sigma[:,j]      #Create positive sigma point, Time Update step 1
        xxm.append(f(t0=t0,t1=t1,x0=xx,**kwargs)) #Propagate positive sigma point, Time Update step 2
        zzh.append(h(x=xxm[j+1],t=t1,**kwargs))            #Predict measurement for positive sigma point, Measurement Update step 1
        xm+=Wj*xxm[j+1] #Accumulate time-updated measurement, Time Update step 3
        zh+=Wj*zzh[j+1]
    for j in range(n):
        xx=x0-Sigma[:,j]      #Create negative sigma point, Time Update step 1
        xxm.append(f(t0=t0,t1=t1,x0=xx,**kwargs)) #Propagate negative sigma point, Time Update step 2
        zzh.append(h(x=xxm[j+n+1],t=t1,**kwargs)) #Predict measurement for negative sigma point, Measurement Update step 1
        xm+=Wj*xxm[j+n+1] #Accumulate time-updated measurement, Time Update step 3
        zh+=Wj*zzh[j+n+1]
    #Note - for any column vector V, np.outer(V,V) is effectively V@V.T. The problem is that
    #we don't really have column vectors, we have rank-1 vectors which are neither row nor column.
    #To do V@V.T, we have to write V[:,np.newaxis]@V[:,np.newaxis].T
    Pm   =Q+W0*np.outer(xxm[0]-xm,xxm[0]-xm) #<Time-updated covariance, initial accumulator
    Gamma=R+W0*np.outer(zzh[0]-zh,zzh[0]-zh) #<Residual covariance, initial accumulator
    S    =  W0*np.outer(xxm[0]-xm,zzh[0]-zh) #<Cross-covariance, initial accumulator
    for j in range(n):
        Pm   +=Wj*np.outer(xxm[j  +1]-xm,xxm[j  +1]-xm) #Update the state covariance with the covariance of the positive sigma point
        Pm   +=Wj*np.outer(xxm[j+n+1]-xm,xxm[j+n+1]-xm) #Update the state covariance with the covariance of the negative sigma point
        Gamma+=Wj*np.outer(zzh[j  +1]-zh,zzh[j  +1]-zh) #Update the residual covariance with the covariance of the positive sigma point
        Gamma+=Wj*np.outer(zzh[j+n+1]-zh,zzh[j+n+1]-zh) #Update the residual covariance with the covariance of the negative sigma point
        S    +=Wj*np.outer(xxm[j  +1]-xm,zzh[j  +1]-zh) #Update the cross-covariance with the covariance of the positive sigma point
        S    +=Wj*np.outer(xxm[j+n+1]-xm,zzh[j+n+1]-zh) #Update the cross-covariance with the covariance of the negative sigma point
    K=S@np.linalg.inverse(Gamma)  #<Kalman gain
    Pp=Pm-K@Gamma@(K.T)   #<Measurement-updated covariance
    #Note that the actual measurement isn't involved until down here, well past the point
    #where the covariances and gain are calculated.
    y=z-zh               #<Measurement innovation
    xp=xm+K@y            #<Measurement-updated estimate
    return (xp,Pp)

def slow_chol(a):
    """
    Compute the Cholesky decomposition of a matrix. This computes
    the matrix L such that a=L@L.T . As such, it is often used
    as a matrix "square root". Called "slow" because it uses pure
    Python, rather than a call to an external library. Direct implementation
    of equations (2.9.4) and (2.9.5) in Numerical Recipes.
    :param a: A square positive definite matrix
    :return: Lower triangular factor.
    """
    if a.shape[0]!=a.shape[1]:
        raise ValueError("Not a square matrix")
    L=np.zeros(a.shape)#*float('NaN') #For debugging purposes, initially set all values to NaN so there is a train
                                      #wreck if we read a vaule before it is written.
    n=a.shape[0]
    for i in range(n):
        Lii=a[i,i]
        for k in range(0,i):
            Lii-=L[i,k]**2
        if Lii<0:
            raise ValueError("Not positive semidefinite, about to sqrt a negative")
        L[i,i]=np.sqrt(Lii)
        for j in range(i+1,n):
            Lji=a[i,j]
            for k in range(0,i):
                Lji-=L[i,k]*L[j,k]
            if L[i,i]==0:
                raise ValueError("Not positive definite, about to divide by zero")
            L[j,i]=Lji/L[i,i]
    return L

def updatea(t0=None, t1=None, x0=None, P0=None, Q=None, f=None, z=None, R=None, h=None,**kwargs):
    """
    Use the augmented unscented Kalman filter to do one time and measurement update for one measurement.
    :param t0: Time that a priori state estimate is valid
    :param t1: Time of measurement
    :param x0: a priori state estimate
    :param P0: a priori state covariance
    :param Q:  process noise covariance. For the augmented UKF, this doesn't have to match the size of P.
    :param f:  State propagation function. Must be of the form:
        def f(t0=None,t1=None,x0=None,...other named parameters from **kwargs...,**kwargs (to soak up the ones we don't care about))
            #:param t0: Time that initial state is valid
            #:param t1: Time to propagate to
            #:param x0: State at t0
            #:return: State at t1
            ... #Perform whatever calculations are necessary to do the time update...
            return x1
    :param z: Measurement vector at t1
    :param R: Measurement noise covariance
    :param h: Measurement function. Must be of the form:
        def h(t=None,x=None,...other named parameters from **kwargs...,**kwargs)
            #:param t: Time of measurement
            #:param x: State vector
            #:return: Measurement vector
            ... #Perform whatever calculations are necessary to construct the measurement...
            return z
    :param kwargs:
    :return: A tuple with the following elements:
      xp: a posteriori state estimate, updated both in time and with the measurement
      Pp: a posteriori state covariance, likewise updated
    """
    W0=1.0/3.0             #Weighting of center sigma point, Time Update step 1
    n=x0.size              #Number of elements in state vector
    m=z.size               #Number of elements in measurement vector
    q=Q.shape[0]           #Number of elements in process noise vector
    na=n+q
    Wj=(1-W0)/(2*na)        #Weighting of other sigma points, Time Update step 1
    #Construct the augmented covariance
    #Pa0=[P0  ,Pxv0
    #     Pxv0,Q   ]
    #Pxv0 is the cross correlation between the state and process noise (assumed to be zero here)
    Pa0=np.zeros((na,na))
    Pa0[0:n,0:n]=P0
    Pa0[n:,n:]=Q
    Sigma=slow_chol(Pa0*na/(1-W0)) #Sigma matrix, used to create sigma points, Time Update step 1
    xx0a=np.zeros((na,2*na+1)) *float('nan')             #<Sigma points of a priori estimate
    xxm=np.zeros((n,2*na+1))   *float('nan')             #<Sigma points of time-updated estimate
    zzh=np.zeros((m,2*na+1))   *float('nan')             #<Sigma points of predicted measurement
    q0=np.zeros(q)
    x0a=np.concatenate((x0,q0))
    #Create center sigma point
    xx0a[:,0]=x0a
    xxm[:,0]=f(t0=t0,t1=t1,x0=x0a[0:n],v=x0a[n:],**kwargs) #Propagate center sigma point, Time Update step 2
    zzh[:,0]=h(x=xxm[:,0],t=t1,**kwargs)     #Predict measurement for center sigma point, Measurement Update step 1
    xm=W0*xxm[:,0]          #<Time-updated estimate, initial accumulator
    zh=W0*zzh[:,0]          #<Predicted measurement, initial accumulator
    for j in range(na):
        xx0a[:,j+1]=x0a+Sigma[:,j]      #Create positive sigma point, including both state and process noise vector, Time Update step 1
        #Note that xxm is not augmented
        xxm[:,j+1]=f(t0=t0,t1=t1,x0=xx0a[0:n,j+1],v=xx0a[n:,j+1],**kwargs) #Propagate positive sigma point, Time Update step 2
        zzh[:,j+1]=h(x=xxm[:,j+1],t=t1,**kwargs)            #Predict measurement for positive sigma point, Measurement Update step 1
        xm+=Wj*xxm[:,j+1] #Accumulate time-updated measurement, Time Update step 3
        zh+=Wj*zzh[:,j+1]
        xx0a[:,j+na+1]=x0a-Sigma[:,j]      #Create negative sigma point, Time Update step 1
        xxm[:,j+na+1]=f(t0=t0,t1=t1,x0=xx0a[0:n,j+na+1],v=xx0a[n:,j+na+1],**kwargs) #Propagate negative sigma point, Time Update step 2
        zzh[:,j+na+1]=h(x=xxm[:,j+na+1],t=t1,**kwargs) #Predict measurement for negative sigma point, Measurement Update step 1
        xm+=Wj*xxm[:,j+na+1] #Accumulate time-updated measurement, Time Update step 3
        zh+=Wj*zzh[:,j+na+1]
    #Note - for any column vector V, np.outer(V,V) is effectively V@V.T. The problem is that
    #we don't really have column vectors, we have rank-1 vectors which are neither row nor column.
    #To do V@V.T, we have to write V[:,np.newaxis]@V[:,np.newaxis].T
    #Note that even though we have na sigma points, Pm is not augmented, only n x n, not na x na.
    #Likewise, Gamma is m x m.
    # Don't add Q, since we take care of that with the augmented sigma points
    Pm  =   W0*np.outer(xxm[:,0]-xm,xxm[:,0]-xm) #<Time-updated covariance, initial accumulator.
    Gamma=R+W0*np.outer(zzh[:,0]-zh,zzh[:,0]-zh) #<Residual covariance, initial accumulator
    S    =  W0*np.outer(xxm[:,0]-xm,zzh[:,0]-zh) #<Cross-covariance, initial accumulator
    for j in range(na*2):
        Pm   +=Wj*np.outer(xxm[:,j+1]-xm,xxm[:,j+1]-xm) #Update the state covariance with the covariance of the sigma point
        Gamma+=Wj*np.outer(zzh[:,j+1]-zh,zzh[:,j+1]-zh) #Update the residual covariance with the covariance of the sigma point
        S    +=Wj*np.outer(xxm[:,j+1]-xm,zzh[:,j+1]-zh) #Update the cross-covariance with the covariance of the sigma point
    K=S@np.linalg.inv(Gamma)  #<Kalman gain
    Pp=Pm-K@Gamma@(K.T)   #<Measurement-updated covariance
    #Note that the actual measurement isn't involved until down here, well past the point
    #where the covariances and gain are calculated.
    y=z-zh               #<Measurement innovation
    xp=xm+K@y            #<Measurement-updated estimate
    return (xp,Pp)
