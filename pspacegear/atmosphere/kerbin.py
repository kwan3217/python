import numpy as np
from collections import namedtuple

"""
 Implementation of an atmosphere made up of earth-composition air, from a
 table of altitude, pressure, and density (like that given in the Kerbin atmosphere page on the KSP wiki).
 We take those tables, interpolate to the given altitude, then try to calculate as many of the other
 parameters as possible, as physically realistically as possible.
"""

# Physical constants and functions relating to the Earth atmosphere at Z=0 (sea level).
# These are declared private since they are from the USSA1976 model, *NOT* the 2002 CODATA 
# recommended values.

# Boltzmann constant, J/K
_k = 1.380622e-23
# Avogadro constant, 1/kmol 
_Na =  6.022169e26
#Ideal gas constant, J/(K kmol). This is k*Na in theory, but this uses the USSA1976 value.
_Rs =      8314.32
# Mean molecular weight of air at sea level, kg/kmol
_M0 =      28.9644
# Standard latitude, deg
_Lat0 = 45.5425
# acceleration of gravity at standard lat, m/s^2
_g0 =      9.82
# Kerbin radius at standard lat, m */
_r0 =       600000
# Celsius scale zero point, K */
_Td =       273.15
# Air pressure at g0, Pa */
_P0 =      101325.
# Temperature at g0, K */
_T0 =       288.15
# Effective collision diameter, m */
_sigma = 3.65e-10
# Ratio of Specific heats for calorically perfect diatomic gas
# http://www.engapplets.vt.edu/fluids/shockcal/shockCalHelp.html#Theory
_gamma =     1.40
# Viscosity coefficient beta, kg/(s m K^1/2)
_beta  = 1.458e-6
# Sutherland's constant, K
_S     =    110.4



def _geopotential(Z):
    """
    Calculate geopotential at a given altitude above the standard latitude
    :param Z: Geometric altitude, m
    :return: Geopotential, m'
    """
    return _r0*Z/(_r0+Z)

def _gravity(Z):
    """
    Calculate magnitude of acceleration of gravity at a given
    altitude above the standard latitude

    :param Z: geometric height, m
    :return: acceleration of gravity, m/s^2
    """
    return _g0*(_r0/(_r0+Z))**2

table_Z  =np.array((     0, 2500, 5000, 7500,10000,15000,20000,25000  ,30000  ,40000   ,50000   ,60000    ,70000))
table_P  =np.array((101325,69015,45625,29126,17934, 6726, 2549,  993.6,  404.1,   79.77,   15.56,    2.387,0    ))
table_rho=np.array(())
table_H=geopotential(table_Z)

# Number of geopotential reference levels
_bmax = 7;
# Maximum geometric altitude of lower atmosphere model, m
_Zmax=86000;
# Geopotential of reference levels, m'
_H =  (      0, 11000,  20000,  32000, 47000,   51000,   71000,_geopotential(_Zmax))
# Lapse rates of reference layeres, K/m'
_Lm = (-0.0065,   0.0, 0.0010, 0.0028,   0.0, -0.0028, -0.0020,                 0.0)
# Hydrostatic constant. Surface gravity * molecular weight / Gas Constant. Constant across lower atmosphere */
_As = _g0*_M0/_Rs

# molecular weight, kg/kmol
_Mgas = {"N2" : 28.0134,
         "O2" : 31.9988,
         "Ar" : 39.948 ,
         "CO2": 44.00995,
         "Ne" : 20.183  ,
         "He" :  4.0026,
         "Kr" : 83.80,
         "Xe" :131.30,
         "CH4": 16.04303,
         "H2" :  2.01594,
         "O3" :47.998}
# Fractional volume composition, unitless. Assumed constant across entire lower atmosphere (except for ozone)
_Fgas = {"N2" :  0.78084,
         "O2" :  0.209476,
         "Ar" :  9.34e-3,
         "CO2":314.00e-6,
         "Ne" : 18.18e-6,
         "He" :  5.24e-6,
         "Kr" :  1.14e-6,
         "Xe" : 87.00e-9,
         "CH4":  2.00e-6,
         "H2" :500.00e-9,
         "O3" :float('NaN')}  #Ozone is calculated separately
  
# Table of Ozone number density (molecule/m^3) vs geometric altitude (m)
_OzoneNumberDensityData=(                                    6.80e+17, #0
                         6.80e+17,5.80e+17,5.70e+17,6.50e+17,1.13e+18, #10000
	                     2.02e+18,2.35e+18,2.95e+18,4.04e+18,4.77e+18, #20000
	                     4.86e+18,4.54e+18,4.03e+18,3.24e+18,2.52e+18, #30000
                         2.03e+18,1.58e+18,1.22e+18,8.73e+17,6.07e+17, #40000
	                     3.98e+17,2.74e+17,1.69e+17,1.03e+17,6.64e+16, #50000
                         3.84e+16,2.55e+16,1.61e+16,1.12e+16,7.33e+15, #60000
	                     4.81e+15,3.17e+15,1.72e+15,7.50e+14,5.40e+14, #70000
	                     2.20e+14,1.70e+14,0.00e+14,0.00e+14,0.00e+14, #80000
	                     0.00e+14,0.00e+14,0.00e+14)

def _TLapse(T0, H0, L, H):
    """
    Calculates temperature at given geopotential in a layer with a
    nonzero temperature lapse rate.

    :param T0: Temperature at base of layer, K
    :param H0: Geopotential of base of layer, m'
    :param L:  Lapse rate, K/m'
    :param H:  Geopotential at altitude at which to calculate temperature, m'
    :return:   Temperature at H, K
    """
    return T0+L*(H-H0)

def _PNoLapse(P0, H0, T0, H):
    """
    Calculates total pressure at given geopotential in a layer with a
    zero temperature lapse rate
    :param P0: Pressure at base of layer, Pa
    :param H0: Geopotential of base of layer, m'
    :param T0: Temperature across layer, K
    :param H:  Geopotential at altitude at which to calculate pressure, m'
    :return:   Pressure at H, Pa
    """
    return P0*np.exp(-_As*(H-H0)/T0)
  
def _PLapse(P0, T0, L, T):
    """
    Calculates pressure at given temperature in a layer with a
    nonzero temperature lapse rate. How bizarre! Pressure is a function
    of altitude, but in this formula it is a function of temp, which
    makes it an indirect function of alt.
    :param P0: Pressure at base of layer, Pa
    :param T0: Temperature of base of layer, K
    :param L:  Lapse rate, K/m'
    :param T:  Temperature at altitude at which to calculate pressure, K
    :return:   Pressure at T, Pa
    """
    return P0*np.exp(_As*np.log(T0/T)/L)

# Pressure at reference levels, Pa */
_P = [_P0]*(_bmax+1)
# Temperature at reference levels, K */
_T = [_T0]*(_bmax+1)
# Calculate base temperatures and pressures, given lapse rates and
# temperature and pressure at Z=0
for b in range(_bmax):
    #Calculate temperature at base of next layer, by linear extrapolation
    _T[b+1] = _TLapse(_T[b],_H[b],_Lm[b],_H[b+1])
    if _Lm[b]== 0.0: #If there is no temp gradient in this layer
        _P[b+1]=_PNoLapse(_P[b],_H[b],_T[b],_H[b+1])
    else:           #If there is
        _P[b+1] = _PLapse(_P[b],_T[b],_Lm[b],_T[b+1])

#Zone of varying molecular weight, from 80-86km
_f = (1.000000, 0.999996, 0.999989, 0.999971, 0.999941, 0.999909,
      0.999870, 0.999829, 0.999786, 0.999741, 0.999694, 0.999641, 0.999579)
_Zinc = 500 # m, incremental height
_Zm = 80000;# m, initial altitude
def _mol(Z):
    """
    Calculate mean molecular weight M from 80 to 86 km

    :param Z: Altitude to calculate at, m
    :return:  Molecular weight at Z, kg/kmol
    """
    if (Z < _Zm):
        return _M0
    else:
      Zi = (Z-_Zm)/_Zinc
      I = int(Zi)
      return _M0*((_f[I+1]-_f[I])*(Zi-I) + _f[I])

AirProperties=namedtuple("AirProperties",["Z",    #Geometric (if it matters) altitude, m. USSA1976 Table 1 property.
                                          "H",    #Geopotential, m'.                               Table 1
                                                  #One "meter" (m') of geopotential is a change in height at any altitude
                                                  #which results in the same change in potential energy as a 1m change in
                                                  #height at sea level Z=0. So, going from H=10000m' to H=20000m' has the
                                                  #same change in potential energy as going from H=0m' to H=10000m', and
                                                  #almost exactly 10000 times as much as going from Z=0m to Z=1m (ignoring
                                                  #the change in gravity over 1m).
                                          "T",    #Temperature, K.                                 Table 1
                                          "Tm",   #Molecular-scale Temperature, K'.                Table 1
                                                  #This is the temperature scaled by the molecular weight, and is the important
                                                  #temperature for some purposes
                                          "P",    #Pressure, Pa.                                   Table 1
                                          "rho",  #Density, kg/m^3 Density, kg/m^3.                Table 1
                                          "M",    #Molecular weight, kg/kmol                       Table 1
                                          "g",    #Acceleration of gravity, m/s^2.                 Table 2
                                          "Hp",   #Pressure Scale Height, m.                       Table 2
                                          "n",    #Molecular number density, 1/m^3.                Table 2
                                          "V",    #Average molecular velocity, m/s.                Table 2
                                          "L",    #Mean free path length, m.                       Table 2
                                          "f",    #Mean Colision frequency, Hz.                    Table 2
                                          "gamma",#Ratio of specific heats. This is usually assumed constant for given
                                                  #gas compostion, so not USSA1976, but could be   Table 2
                                          "Cs",   #Speed of sound, m/s.                            Table 3 (lower atmosphere)
                                          "mu",   #Dynamic Viscosity, Ns/(m^2).                    Table 3
                                                  #Other than being resistance to flow, I don't know what it means
                                          "eta",  #Kinematic Viscosity, m^2/s                      Table 3
                                          "kt",   #Thermal conductivity, W/(m*K).                  Table 3
                                          "ngas"  #Number density of component gases, 1/m^3.       Table 4 (upper atmosphere)
                                                  #This is a dictionary. Key is chemical formula of gas, value is
                                                  #number density of that gas. It is only tabulated for the upper atmosphere
                                                  #because it is constant for all gases except O3 in the lower atmosphere.
])

def lower_atmosphere(Z):
    """
    Calculate the properties of the lower atmosphere (<86000m' geopotential)
    :param Z:
    :return:
    """
    H = _geopotential(Z)
    g = _gravity(Z)
    for bi in range(_bmax):
        if H>=_H[bi]:
            b=bi
    ngas = {}
    if Z<=_Zmax:
        Tm = _TLapse(_T[b],_H[b],_Lm[b],H)
        if(_Lm[b]== 0.0):
          P = _PNoLapse(_P[b],_H[b],_T[b],H)
        else:
          P = _PLapse(_P[b],_T[b],_Lm[b],Tm)
        M = _mol(Z)
        T = Tm*M/_M0
        rho = (P/Tm)*(_M0/_Rs)
        Hp = (_Rs/_M0)*(Tm/g)                     #pressure scale height, m
        n = (_Na/_Rs)*(P/T)                       #number density, 1/m^3
        V = np.sqrt((8*_Rs/(np.pi*_M0))*Tm)       #mean air-particle speed, m/s
        L = 1.0/((np.sqrt(2)*np.pi*_sigma**2)*n)  #mean free path, m
        f=V/L                                     #collision frequency, 1/s
        Cs = np.sqrt((_gamma*_Rs/_M0)*Tm)         #speed of sound, m/s
        mu = _beta*np.sqrt(T)/(1+_S/T)            #dynamic viscosity, Ns/m^2
        eta = mu/rho;                             #kinematic viscosity, m^2/s
        kt = 2.64638e-3*np.sqrt(T)/(1+245.4*np.exp(-12*np.log(10)/T)/T);   #coef. of thermal conductivity, W/(m K)
        for gasname in _Fgas:
            ngas[gasname]=_Fgas[gasname]*n
        # Calculate for ozone
        # ngas["O3"]=OzoneNumberDensity.Interp(Z,0);
        # if ngas["O3"]<0:
        #     ngas["O3"]=0;
        gamma=_gamma
    else:
        # These values describe a vacuum, with zero where appropriate (like pressure)
        # and NaN where the quantity is not defined (like temperature)
        Tm = float('nan')
        P = 0.0
        T = float('nan')
        M = float('nan')
        rho = 0.0
        Hp= float('nan')
        n = float('nan')
        V = float('nan')
        L = float('nan')
        f = float('nan')
        Cs= float('nan')
        mu= float('nan')
        eta=float('nan')
        kt =float('nan')
        gamma =float('nan')
        ngas={}
    A=AirProperties(Z=Z,H=H,Tm=Tm,P=P,M=M,T=T,rho=rho,g=g,Hp=Hp,n=n,V=V,L=L,f=f,Cs=Cs,mu=mu,eta=eta,kt=kt,gamma=gamma,ngas=ngas)
    return A

if __name__=="__main__":
    AP=lower_atmosphere(0)
    print(AP)
    AP = lower_atmosphere(5334.0*0.3048) #Boulder (Folsom field level)
    print(AP)
    AP = lower_atmosphere(7220.0*0.3048) #Top of the Third Flatiron
    print(AP)
    AP = lower_atmosphere(14200.0*0.3048) #Top of Mt Evans
    print(AP)
    AP = lower_atmosphere(29035.0*0.3048) #Top of Mt Everest
    print(AP)
    AP = lower_atmosphere(49112.48)       #"50km" height in calc_albedo model (2.4e23 /cm^2 coldens)
    print(AP)
    AP = lower_atmosphere(83000)          #PMC Cloud Deck
    print(AP)
    AP = lower_atmosphere(100000)         #Above the top of the model - should return vacuum
    print(AP)
