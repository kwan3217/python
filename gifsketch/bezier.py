import matplotlib.pyplot as plt
import numpy as np

def r_conic(a,e,theta):
    return

#Theta is intended to be 0 at apoapse, pi at periapse, and back to tau at next apoapse
theta=np.arange(0,6.28,0.01)
ra=2.0
i_e1=range(0,99)
i_e2=range(990,1001)
es=np.hstack((np.array(i_e1)/100,np.array(i_e2)/1000))
for i_e,e in enumerate(es):
    rp=ra*(1-e)/(1+e)
    a=(ra+rp)/2
    #We use 1-ecos instead of 1+ecos as in the front cover equations, since
    #the front cover wants \theta=0 to be zero at periapse instead of apoapse.
    #To convert this, we want 1+ecos(\theta+pi), and cos(\theta+pi)=-cos\theta
    r=(a*(1-e**2))/(1-e*np.cos(theta))
    r[0]=ra
    v=np.sqrt(2.0/r-1.0/a)
    T=2*np.pi*np.sqrt(a**3)
    Tfall=2*np.pi*np.sqrt((ra/2)**3)/2
    Tcirc=2*np.pi*np.sqrt(ra**3)
    vcirca=np.sqrt(1/ra)
    try:
        vcircp=np.sqrt(1/rp)
    except ZeroDivisionError:
        vcircp=float('Inf')
    va=np.sqrt(2.0/ra-1.0/a)
    try:
        vp=np.sqrt(2.0/rp-1.0/a)
    except ZeroDivisionError:
        vp=float('Inf')
    x=r*np.cos(theta)
    y=r*np.sin(theta)
    rf2=ra-rp
    print(ra,rp,a,e,T,Tcirc,Tfall/Tcirc)
    plt.figure("r,v")
    plt.cla()
    plt.plot(theta,r,theta,v)
    plt.figure("orbit")
    plt.cla()
    plt.plot(np.cos(theta),np.sin(theta),'k-')
    plt.plot([0,rf2],[0,0],'r+')
    draw_v=False
    draw_c=False
    if draw_v:
        plt.plot([-rp,-rp],[0,-np.min((vp,2.5))],'r-')
        if draw_c:
            plt.plot([-rp-0.05,-rp-0.05],[0,-np.min((vcircp,2.5))],'g-')
        plt.plot([ra,ra],[0,va],'r-')
        if draw_c:
            plt.plot([ra+0.05,ra+0.05],[0,vcirca],'g-')
    plt.plot(x,y,'b-')
    plt.axis('equal')
    plt.axis([-3,3,-2.2,2.2])
    ax=plt.gca()
    ax.set_autoscale_on(False)
    plt.text( ra, 0.2,"ra=%5.3f"%ra)
    plt.text( ra, 0.0,"e=%5.3f"%e)
    plt.text( ra,-0.2,"T=%5.3f"%T)
    plt.text( ra,-0.4,"va=%5.3f"%va,color='r' if draw_v else 'k')
    if draw_c:
        plt.text( ra,-0.6,"vc=%5.3f"%vcirca,color='g' if draw_v else 'k')
    plt.text(-rp, 0.2,"rp=%5.3f"%rp)


    plt.text(-rp,-0.4,"vp=%5.3f"%vp,color='r' if draw_v else 'k')
    if draw_c:
        plt.text(-rp,-0.6,"vc=%5.3f"%vcircp,color='g' if draw_v else 'k')
    plt.savefig("Frames/eccentricity/ecc%03d"%i_e)
    plt.pause(0.001)
#Repeat the last frame a few times
for i in range(len(es),len(es)+30):
    plt.savefig("Frames/eccentricity/ecc%03d"%i)
