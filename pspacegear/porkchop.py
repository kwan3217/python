import spiceypy as cspice
import bmw
import numpy as np
import collections

PorkchopResult=collections.namedtuple("PorkchopResult",["dr","dv_planet","dv_sc","ar","av_planet","av_sc"])
def porkchop(etd0, etd1, nd, eta0, eta1, na, bodyd, bodya,frame="J2000"):
    mu=cspice.gdpool("BODY10_GM",0,1)[0]
    au=cspice.gdpool("AU",0,1)[0]
    dr       =np.zeros(( 1,nd,3))
    dv_planet=np.zeros(( 1,nd,3))
    dv_sc    =np.zeros((na,nd,3))
    ar       =np.zeros((na, 1,3))
    av_planet=np.zeros((na, 1,3))
    av_sc    =np.zeros((na,nd,3))
    for id in range(nd):
        print(id)
        if nd==1:
            etd=etd0
        else:
            etd = etd0 + (etd1 - etd0) * id / (nd-1)
        (dx,_)=cspice.spkezr(str(bodyd), etd, frame,"NONE","SUN")
        dr       [0,id,:]=dx[0:3]
        dv_planet[0,id,:]=dx[3:6]
        for ia in range(na):
            if na==1:
                eta=eta0
            else:
                eta = eta0 + (eta1 - eta0) * ia / (na-1)
            if id==0:
                (ax,_) = cspice.spkezr(str(bodya), eta, frame, "NONE", "SUN")
                ar       [ia, 0, :] = ax[0:3]
                av_planet[ia, 0, :] = ax[3:6]
            (dv_sc[ia,id,:],av_sc[ia,id,:])=bmw.gauss(dx[0:3],ar[ia,0,:],eta-etd,l_DU=au,mu=mu,Type=-1)
    return PorkchopResult(dr=dr,dv_planet=dv_planet,dv_sc=dv_sc,
                          ar=ar,av_planet=av_planet,av_sc=av_sc)

def exercise_porkchop():
    import os
    import matplotlib.pyplot as plt
    os.chdir("/Users/jeppesen/workspace/Data/spice/generic/")
    cspice.furnsh("generic.tm")
    etd0=cspice.str2et("2018-05-05 00:00:00 TDB")
    etd1=cspice.str2et("2018-06-09 00:00:00 TDB")
    id=0
    eta0=cspice.str2et("2018-11-26 00:00:00 TDB")-30*86400
    ia=30
    eta1=eta0+60*86400
    nd=int((etd1-etd0)//86400+1)
    na=int((eta1-eta0)//86400+1)
    result=porkchop(etd0,etd1,nd,eta0,eta1,na,"3","4",frame="J2000")
    dvd=result.dv_sc-result.dv_planet
    avd=result.av_sc-result.av_planet
    dc3=np.sum(dvd**2,axis=2)
    ac3=np.sum(avd**2,axis=2)
    cs=plt.contour(dc3,range(21))
    plt.clabel(cs)
    #plt.plot(range(nd),c3.transpose())

    rpark=cspice.gdpool("BODY399_RADII",0,1)[0]+185
    dmu=cspice.gdpool("BODY399_GM",0,1)[0]
    amu=cspice.gdpool("BODY499_GM",0,1)[0]
    rentry=cspice.gdpool("BODY499_RADII",0,1)[0]+125
    vpark=np.sqrt(dmu/rpark)
    dvinf=np.sqrt(dc3)
    avinf=np.sqrt(ac3)
    da=-dmu/dc3
    aa=-amu/ac3
    dvp=np.sqrt(2*dmu/rpark-dmu/da)
    departv=dvp-vpark
    dla=np.degrees(np.arctan2(dvd[:,:,2],np.sqrt(dvd[:,:,0]**2+dvd[:,:,1]**2)))
    plt.show()

if __name__=="__main__":
    exercise_porkchop()