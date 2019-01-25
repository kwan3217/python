import numpy as np
import matplotlib.pyplot as plt

def q(a,b,yi,yip1,xi,xip1,t):
    return (1-t)*yi+t*yip1+t*(1-t)*(a*(1-t)+b*t)

def qp(a,b,yi,yip1,xi,xip1,t):
    return (yip1-yi)/(xip1-xi)+(1-2*t)*(a*(1-t)+b*t)/(xip1-xi)+t*(1-t)*(b-a)/(xip1-xi)

def qpp(a,b,yi,yip1,xi,xip1,t):
    return 2*(b-2*a+(a-b)*3*t)/((xip1-xi)**2)

def clamped_spline(x,y,K0=0,Kn=0,test=True):
    """
    Fit n cubic polynomials through n+1 points. Across each point, the polynomials will be smooth (equal first and
    second derivative). At the first and last points, the
    :param x: array of n+1 scalar independent variable values
    :param y: array of n+1,... dependent values at the given independent values
    :return: A tuple of coefficients for n cubic polynomials
    """
    #Find the first derivative at each point by solving a tridiagonal matrix equation
    n=len(x)-1
    if test:
        A=np.zeros((n+1,n+1),dtype=float)
        Ap=np.zeros((n+1,n+1),dtype=float)
        B=np.zeros(n+1,dtype=float)
        Bp=np.zeros(n+1,dtype=float)
    fi=1
    gi=0
    hi=K0
    gp=np.zeros(n)*float('nan')
    gp[0]=gi/fi
    hp=np.zeros(n+1)*float('nan')
    hp[0]=hi/fi
    if test:
        A[0,0]=fi
        A[0,1]=gi
        Ap[0,0]=1
        Ap[0,1]=gp[0]
        B[0]=hi
        Bp[0]=hp[0]
    for i in range(1,n+1):
        gi = float('nan')
        if i<n:
            ei = 1 / (x[i] - x[i - 1])
            gi = 1 / (x[i + 1] - x[i])
            fi = 2 * (ei + gi)
            gp[i]=gi/(fi-ei*gp[i-1])
            hi=3*((y[i+1]-y[i])/((x[i+1]-x[i])**2)+(y[i]-y[i-1])/((x[i]-x[i-1])**2))
        else:
            ei=0
            fi=1
            hi=Kn
        hp[i] = (hi - ei * hp[i - 1]) / (fi - ei * gp[i - 1])
        if test:
            A[i, i - 1] = ei
            A[i, i] = fi
            Ap[i,i]=1
            if i<n-1:
                A[i, i + 1] = gi
                Ap[i, i + 1] = gp[i]
            B[i] = hi
            Bp[i] = hp[i]
    k=np.zeros(n+1)*float('nan')
    k[n]=hp[n]
    for i in range(n-1,-1,-1):
        k[i]=hp[i]-gp[i]*k[i+1]
    #Now that we have the slopes, calculate the a and b for each segment
    a= k[0:n  ]*(x[1:n+1]-x[0:n])-(y[1:n+1]-y[0:n])
    b=-k[1:n+1]*(x[1:n+1]-x[0:n])+(y[1:n+1]-y[0:n])
    return (a,b)

x=np.array([1,2,3,4])
y=np.array([0,1,1,0])
(a,b)=clamped_spline(x,y)
plt.plot(x,y,'k*')
for i in range(0,len(a)):
    t=np.arange(0,1.01,0.01)
    print(i,x[i],y[i],q(a[i],b[i],y[i],y[i+1],x[i],x[i+1],0))
    plt.plot(x[i]+t*(x[i+1]-x[i]),q  (a[i],b[i],y[i],y[i+1],x[i],x[i+1],t),'b-')
    plt.plot(x[i]+t*(x[i+1]-x[i]),qp (a[i],b[i],y[i],y[i+1],x[i],x[i+1],t),'g-')
    plt.plot(x[i]+t*(x[i+1]-x[i]),qpp(a[i],b[i],y[i],y[i+1],x[i],x[i+1],t),'r-')
plt.axis('equal')
plt.show()