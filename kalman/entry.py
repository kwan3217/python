import ukf
import numpy as np
import matplotlib.pyplot as plt

def entryxdot(x,R0=6374,Gm0=3.9860e5,H0=13.406,beta0=0.59783,v=None,**kwargs):
    """
    Evaluate the equations of motion (IE time derivative of the state vector)
    for the reentry problem.
    :param x: State vector. First two elements are position, second two are velocity, last is ballistic term
    :param R0: Radius of Earth. Used to establish the reference altitude in the model atmosphere
    :param Gm0: Gravitational parameter of the Earth
    :param H0:  Scale height of the model atmosphere
    :param beta0: Nominal ballistic coefficient
    :param v: Optional random walk vector
    :param kwargs: Soak up extra named parameters intended for other functions
    :return:
    """
    R=np.linalg.norm(x[0:2])
    V=np.linalg.norm(x[2:4])
    beta=beta0*np.exp(x[4])
    G=-Gm0/R**3
    D=-beta*V*np.exp((R0-R)/H0)
    result=np.array((x[2],x[3],D*x[2]+G*x[0],D*x[3]+G*x[1],0))
    if v is not None:
        result[2:5]+=v
    return result

def entryF(t0=None, t1=None, x0=None,v=None,steps=1,**kwargs):
    """
    Integrate the equations of motion. We will use a simple Euler method to begin with
    :param t0: Time that initial state is valid
    :param t1: Time to propagate to
    :param x0: State at t0
    :param v: Random walk vector, presumed constant from t0 to t1
    :param steps: Number of Euler steps to take
    :return: State at t1
    """
    x1=x0.copy()
    dt=(t1-t0)/steps
    for i in range(steps):
        x1+=dt*entryxdot(x0,v=v,**kwargs)
    return x1

def entryH(t=None,x=None,xr=6374,yr=0,w=None,**kwargs):
    """
    Calculate observation (range and azimuth) of target
    :param t: Time that the state is valid (not used)
    :param x: State at t
    :param xr: Radar x coordinate
    :param yr: Radar y coordinate
    :param w: Optional noise vector. If not passed, produce perfect measurements
    :param kwargs:
    :return: Two-element measurement vector. First element is range, second is azimuth in radians
    """
    result=np.array((np.sqrt((x[0]-xr)**2+(x[1]-yr)**2),np.arctan2(x[1]-yr,x[0]-xr)))
    if w is not None:
        result+=w
    return result

#Test Cholesky decomposition
a=np.array(((  4, 12,-16),
            ( 12, 37,-43),
            (-16,-43, 98)))
b=np.linalg.cholesky(a)
print(b)
print(b@b.T)
c=ukf.slow_chol(a)
print(c)
print(c@c.T)
print(c-b)

t0=0
t1=200
dt=0.1
n=int(t1/dt)
x0t=np.array((6500.4,349.14,-1.8093,-6.7967,0))
P0mc=1e-6*np.eye(5)
P0mc[4,4]=0
Q=np.zeros((3,3))
Q[0,0]=2.4064e-5
Q[1,1]=2.4064e-5
AQ=np.zeros((3,3))
AQ[0:2,0:2]=np.linalg.cholesky(Q[0:2,0:2])
x0=np.array((6500.4,349.14,-1.8093,-6.7967,0.6932))
P0=1e-6*np.eye(5)
P0[4,4]=1
R=np.zeros((2,2))
R[0,0]=0.001**2
R[1,1]=0.00017**2
AR=np.linalg.cholesky(R)
x=x0.copy()
P=P0.copy()
t=t0
xt=x0t.copy()
ts=np.zeros(n)
xts=np.zeros((5,n))
zts=np.zeros((2,n))
zs=np.zeros((2,n))
xs=np.zeros((5,n))
for i in range(n):
    #Generate the actual trajectory
    v=AQ@np.random.randn(3) #Generate the random walk
    #v=None
    xt=entryF(t0=t,t1=t+dt,x0=xt,v=v,steps=2) #Integrate the trajectory
    ts[i]=t+dt
    xts[:,i]=xt
    #Generate the measurement
    w=AR@np.random.randn(2) #Generate the measurement noise
    z=entryH(t=t,x=xt,w=w,steps=2)  #Generate the measurement
    zts[:,i]=z
    #Run the UKF
    (x,P)=ukf.updatea(t0=t,t1=t+dt,x0=x,P0=P,Q=Q,f=entryF,z=z,R=R,h=entryH,steps=2)
    xs[:,i]=x
    z=entryH(t=t,x=x,steps=2)
    zs[:,i]=z
    t=t+dt

for i in range(5):
    plt.figure(i)
    plt.plot(ts, xts[i, :]-xs[i, :])

plt.figure(6)
plt.plot(xts[0,:],xts[1,:],'b-',np.sqrt(6374**2-xts[1,:]**2),xts[1,:],'g-')
plt.plot(xts[0,::100],xts[1,::100],'b+')
plt.plot(xs[0,:],xs[1,:],'r-')
plt.plot(xs[0,::100],xs[1,::100],'r+')
plt.axis('equal')
plt.show()