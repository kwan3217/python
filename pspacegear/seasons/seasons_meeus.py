"""
Python version of algorithm from Astronomical Algorithms, chapter 26. Only uses
the tables that are valid from 1000 AD to 3000 AD (Table 26.B)
"""

import numpy as np

def JDE0(year,season):
    y=(year-2000)/1000.0
    return season[0]+season[1]*y+season[2]*y**2+season[3]*y**3+season[4]*y**4

#polynomials describing the rough JDE of equinox as a polynomial function
#of millenia from AD 2000 (Table 26.B)
Mar=[2451623.80984,365242.37404,+0.05169,-0.00411,-0.00057]
Jun=[2451716.56767,365241.62603,+0.00325,+0.00888,-0.00030]
Sep=[2451810.21715,365242.01767,-0.11575,+0.00337,+0.00078]
Dec=[2451900.05952,365242.74049,+0.06223,-0.00823,-0.00032]

dtYear=list(range(1974,2019,1))
         #xxx0    xxx1    xxx2    xxx3    xxx4    xxx5    xxx6    xxx7    xxx8    xxx9
dtTable=[                                44.4841,45.4761,46.4567,47.5214,48.5344,49.5861, #1974-1979
         50.5387,51.3808,52.1668,52.9565,53.7882,54.3427,54.8712,55.3222,55.8197,56.3000, #1980-1989
         56.8553,57.5653,58.3092,59.1218,59.9845,60.7853,61.6287,62.2950,62.9659,63.4673, #1990-1999
         63.8285,64.0908,64.2998,64.4734,64.5736,64.6876,64.8452,65.1464,65.4573,65.7768, #2000-2009
         66.0699,66.3246,66.6030,66.9069,67.2810,67.6439,68.1024,68.5927,68.9676          #2010-2018
]

def T(year,season):
    return (JDE0(year,season)-2451545.0)/36525.0

def W(T):
    return 35999.373*T-2.47

def Deltalam(W):
    return 1+0.0334*np.cos(np.radians(W))+0.0007*np.cos(np.radians(2*W))

def DeltaT(year):
    """
    Calculate difference between ephemeris time and universal time in days
    :param year:
    :return:
    """
    #From https://eclipse.gsfc.nasa.gov/SEhelp/deltatpoly2004.html
    return -(26.92+0.32217*(year-2000)+0.005589*(year-2000)**2)/86400.0
    #return -32.184/86400.0
    #return -(66.0+(year-2000)*1.0)/86400.0

As=[485,203,199,182,156,136, 77, 74, 70, 58, 52, 50, 45, 44, 29, 18, 17, 16, 14, 12, 12, 12,  9,  8]
Bs=[324.96,337.23,342.08, 27.85, 73.14,171.52,222.54,296.72,243.58,119.81,297.17, 21.02,247.54,325.15, 60.93,155.12,288.79,198.04,199.76, 95.39,287.11,320.81,227.73, 15.45]
Cs=[  1934.136, 32964.467,    20.186,445267.112, 45036.886, 22518.443, 65928.934,  3034.906,  9037.513, 33718.147,   150.678,  2281.226, 29929.562, 31555.956,  4443.417, 67555.328,  4562.452, 62894.029, 31436.921, 14577.848, 31931.756, 34777.259,  1222.114, 16859.074]

def S(T):
    result=0
    for A,B,C in zip(As,Bs,Cs):
        result+=A*np.cos(np.radians(B+C*T))
    return result

def JDE(year,season):
    jde0=JDE0(year,season)
    t=T(year,season)
    s=S(t)
    w=W(t)
    dl=Deltalam(w)
    return jde0+s*0.00001/dl+DeltaT(year)

def caldat(JD):
    Z=int(JD+0.5)
    F=(JD+0.5)-Z
    if Z<2299161:
        A=Z
    else:
        alpha=int((Z-1867216.25)/36524.25)
        A=Z+1+alpha-int(alpha/4)
    B=A+1524
    C=int((B-122.1)/365.25)
    D=int(365.25*C)
    E=int((B-D)/30.6001)
    day=B-D-int(30.6001*E)
    if E<14:
        month=E-1
    else:
        month=E-13
    if month>2:
        year=C-4716
    else:
        year=C-4715
    sod=F*86400
    min=int(sod/60)
    sec=sod-min*60
    hour=int(min/60)
    min=min-hour*60
    return (year,month,day,hour,min,sec)

if __name__=="__main__":
    with open("/home/jeppesen/Python seasons table.txt","w") as ouf:
        print("                  d  h                      d  h  m           d  h  m",file=ouf)
        for year in range(1992,2021):
            mar_equ=caldat(JDE(year,Mar))
            jun_sol=caldat(JDE(year,Jun))
            sep_equ=caldat(JDE(year,Sep))
            dec_sol=caldat(JDE(year,Dec))
            print("%04d                        %04d"%(year,year),file=ouf)
            print("Perihelion  Jan   X XX    Equinoxes  Mar   %02d %02d %02d %02d Sept  %02d %02d %02d %02d"%(mar_equ[2],mar_equ[3],mar_equ[4],mar_equ[5],sep_equ[2],sep_equ[3],sep_equ[4],sep_equ[5]),file=ouf)
            print("Aphelion    July  X XX    Solstices  June  %02d %02d %02d %02d Dec   %02d %02d %02d %02d"%(jun_sol[2],jun_sol[3],jun_sol[4],jun_sol[5],dec_sol[2],dec_sol[3],dec_sol[4],dec_sol[5]),file=ouf)
            print(file=ouf)
            #print(caldat(JDE(year,Mar)),caldat(JDE(year,Jun)),caldat(JDE(year,Sep)),caldat(JDE(year,Dec)))