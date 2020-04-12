"""
Make the equivalent to DE430 for the Kerbal system.
"""
import subprocess
import numpy as np
import os

#Spice number follows convention of real solar system, except with 1000 added. This allows real and Kerbal
#solar systems to be loaded at the same time.
#If no orbital elements are provided, we will create a kernel with a constant zero vector state relative to parent
#We use the name Kerbol instead of Sun so that it doesn't conflict with the real Sun. Likewise Kerbol system barycenter
#rather than solar system barycenter (should it be Kerbolar system? Or is that pushing consistency too far?)
elements={
#   Sp#   Name                         GM/(km,s)            r_e/km  parent a/km         e       i/deg ap/deg lan/deg M/rad
    1000: ("Kerbol system barycenter", 1.17233279483249e+9, 261600,    0,  None,                                                     None),
    1199: ("Moho",                     1.68609378654509e+2,    250, 1001,  None,                                                     1210000),
    1299: ("Eve",                      8.17173022921085e+3,    700, 1002,  None,                                                       80500),
    1201: ("Gilly ",                   8.28944981471635e-3,     13, 1002, (   31500    ,0.55   ,12    , 10, 80, 0.899999976158143),    28255),
    1399: ("Kerbin ",                  3.53160000000000e+3,    600, 1003,  None,                                                     -6*3600),
    1301: ("Mun",                      6.51383975207806e+1,    200, 1003, (   12000    ,0.00   , 0    , 0, 0, 1.7),                  "sync" ),
    1302: ("Minmus ",                  1.76580002631247e+0,     60, 1003, (   47000    ,0.00   , 6    , 38, 78, 0.9),                  40400),
    1499: ("Duna ",                    3.01363211975098e+2,    320, 1004,  None,                                                     "moonsync/1401"),
    1401: ("Ike",                      1.85683685731441e+1,    130, 1004, (    3200    ,0.03   , 0.2  , 0, 0, 1.70000004768372),     "sync" ),
    1599: ("Dres",                     2.14844886000000e+1,    138, 1005,  None,                                                       34800),
    1699: ("Jool",                     2.82528004209995e+5,   6000, 1006,  None,                                                       36000),
    1601: ("Laythe ",                  1.96200002923608e+3,    500, 1006, (   27184    ,0.00   , 0    , 0, 0, 3.14000010490418),     "sync" ),
    1602: ("Vall",                     2.07481499473751e+2,    300, 1006, (   43152    ,0.00   , 0    , 0, 0, 0.899999976158143),    "sync" ),
    1603: ("Tylo",                     2.82528004209995e+3,    600, 1006, (   68500    ,0.00   , 0.025, 0, 0, 3.14000010490418),     "sync" ),
    1604: ("Bop",                      2.48683494441491e+0,     65, 1006, (  128500    ,0.235  ,15    ,25, 10, 0.899999976158143),   "sync" ),
    1605: ("Pol",                      7.21702080000000e-1,     44, 1006, (  179890    ,0.17085, 4.25 ,15, 2, 0.9),                  "sync" ),
    1799: ("Eeloo",                    7.44108145270496e+1,    210, 1007,  None,                                                       19460),
    1010: ("Kerbol",                   1.17233279483249e+9, 261600, 1000,  None,                                                      432000),
    1001: ("Moho Barycenter",          1.68609378654509e+2,    250, 1000, ( 5263138.304,0.2    , 7, 15, 70, 3.14000010490418),       None),
    1002: ("Eve Barycenter",           8.17173022921085e+3,    700, 1000, ( 9832684.544,0.01   , 2.1, 0, 15, 3.14000010490418),      None),
    1003: ("Kerbin Barycenter",        3.53160000000000e+3,    600, 1000, (13599840.256,0      , 0, 0, 0, 3.14),                     None),
    1004: ("Duna Barycenter",          3.01363211975098e+2,    320, 1000, (20726155.264,0.051  , 0.06, 0, 135.5, 3.14000010490418),  None),
    1005: ("Dres Barycenter",          2.14844886000000e+1,    138, 1000, (40839348.203,0.145  , 5, 90, 280, 3.14),                  None),
    1006: ("Jool Barycenter",          2.82528004209995e+5,   6000, 1000, (68773560.320,0.05   , 1.304, 0, 52, 0.100000001490116),   None),
    1007: ("Eeloo Barycenter",         7.44108145270496e+1,    210, 1000, (90118820.000,0.26   , 6.15, 260, 50, 3.14),               None)
}

ksp_version="0.18.2"
ke_version=18 #Intended to be xyz for ksp version x.y.z, which will be a problem if any of those exceed one digit

elements_header="""
Spice kernel covering the natural objects in the Kerbol system. This segment
is for the position of %s (Spice ID %d) relative to %s (Spice ID %d).

Using data from spreadsheet, extracted from KSP version %s

The K0001 epoch (equivalent to J2000) is used, as this matches the in-game
time scale. The clock ticks at 1 SI second per second, and starts at 
Year 1, Day 1, 00:00:00 UT which is the exact instant that the game starts
and matches all UT in seconds in the persistence file.

\\begindata
INPUT_DATA_TYPE = 'ELEMENTS'
OUTPUT_SPK_TYPE = 15
OBJECT_ID=%d
CENTER_ID=%d
REF_FRAME_NAME='J2000'
PRODUCER_ID='C. Jeppesen, Kwan Systems'
DATA_ORDER='EPOCH A E INC PER NOD MEAN'
TIME_WRAPPER='# ETSECONDS'
INPUT_DATA_UNITS = ('ANGLES=DEGREES' 'DISTANCES=km')
DATA_DELIMITER=','
LINES_PER_RECORD=1
CENTER_GM=%.14e
CENTER_J2=0
CENTER_EQ_RADIUS=%.0f
PRECESSION_TYPE='NO PRECESSION'
CENTER_POLE_RA=0
CENTER_POLE_DEC=90
START_TIME='2000-Jan-01 12:00:00 TDB'
STOP_TIME='3000-Jan-01 12:00:00 TDB'
LEAPSECONDS_FILE='kerb0000.tls'
"""

state_header="""
Spice kernel covering the natural objects in the Kerbol system. This segment
is for the position of %s (Spice ID %d) relative to %s (Spice ID %d).

Using data from spreadsheet, extracted from KSP version %s

The K0001 epoch (equivalent to J2000) is used, as this matches the in-game
time scale. The clock ticks at 1 SI second per second, and starts at 
Year 1, Day 1, 00:00:00 UT which is the exact instant that the game starts
and matches all UT in seconds in the persistence file.

\\begindata
INPUT_DATA_TYPE = 'STATES'
OUTPUT_SPK_TYPE = 13
OBJECT_ID=%d
CENTER_ID=%d
REF_FRAME_NAME='J2000'
PRODUCER_ID='C. Jeppesen, Kwan Systems'
DATA_ORDER='EPOCH X Y Z VX VY VZ'
TIME_WRAPPER='# ETSECONDS'
INPUT_DATA_UNITS = ('ANGLES=DEGREES' 'DISTANCES=km')
DATA_DELIMITER=','
LINES_PER_RECORD=1
START_TIME='2000-Jan-01 12:00:00 TDB'
STOP_TIME='3000-Jan-01 12:00:00 TDB'
LEAPSECONDS_FILE='kerb0000.tls'
POLYNOM_DEGREE=1
"""
K0001=0.0
K1001=31556995200.0 #1000 Earth years after K0001, in seconds

try:
    os.remove("ke%03d.bsp"%ke_version)
except FileNotFoundError:
    pass
with open("kerb018.tpc","wt") as tpc:
    print("\\begindata",file=tpc)
    for i_spice in elements:
        name=elements[i_spice][0]
        gm=elements[i_spice][1]
        re=elements[i_spice][2]
        i_parent=elements[i_spice][3]
        if i_parent!=0:
            parent_name=elements[i_parent][0]
            parent_gm=elements[i_parent][1]
            parent_re=elements[i_parent][2]
        else:
            parent_name="Inertial Origin"
            parent_gm=0
            parent_re=0
        this_elements = elements[i_spice][4]
        rot_period=elements[i_spice][5]
        print("NAIF_BODY_CODE+=(%d)" %i_spice,file=tpc)
        print("NAIF_BODY_NAME+=('%s')"%name.strip().upper(),file=tpc)
        if(i_spice>=1010):
            print("BODY%04d_RADII=(%d,%d,%d)"%(i_spice,re,re,re),file=tpc)
        print("BODY%04d_GM=(%.15e)" %(i_spice,gm), file=tpc)
        if rot_period is None:
            pass
        else:
            print("BODY%04d_POLE_RA=(0,0,0)"%i_spice, file=tpc)
            print("BODY%04d_POLE_DEC=(0,0,0)"%i_spice, file=tpc)
            if rot_period=="sync":
                #Synchronous with planet -- calculate orbital period and use it
                orb_period=2*np.pi*np.sqrt(this_elements[0]**3/parent_gm)
                rot_period=orb_period
            elif type(rot_period)==str and rot_period.startswith("moonsync"):
                #Synchronous with a moon -- calculate orbital period of moon and use it
                i_child=int(rot_period.split("/")[1])
                child_elements = elements[i_child][4]
                orb_period=2*np.pi*np.sqrt(child_elements[0]**3/gm)
                rot_period=orb_period
            elif rot_period<0:
                #Solar day - calculate orbital period and use formula to calculate sidereal rot_period
                parent_element=elements[i_parent][4]
                grandparent_gm=elements[elements[i_parent][3]][1]
                orb_period=2*np.pi*np.sqrt(parent_element[0]**3/grandparent_gm)
                rot_period=-rot_period
                n=orb_period/rot_period
                d=n+1
                rot_period=rot_period*n/d
            else:
                #Already sidereal day, so we're good.
                pass
            w=360*86400.0/rot_period
            print("BODY%04d_POLE_W=(0,%.15f,0)"%(i_spice,w), file=tpc)
        mkspkfn = "mkspk%03d_%04d.txt" % (ke_version,i_spice)
        datafn = "data%03d_%04d.txt" %(ke_version, i_spice)
        with open(mkspkfn, "wt") as mkspk:
            with open(datafn, "wt") as data:
                if this_elements is None:
                    print(state_header%(name,i_spice,parent_name,i_parent,ksp_version,i_spice,i_parent),file=mkspk)
                    for t in (K0001,K1001):
                        print("%f,0,0,0,0,0,0"%t,file=data)
                else:
                    print(elements_header % (name, i_spice, parent_name, i_parent, ksp_version,i_spice, i_parent,parent_gm,parent_re), file=mkspk)
                    print("0,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e" % (
                        this_elements[0],this_elements[1],this_elements[2],
                        this_elements[3],this_elements[4],np.degrees(this_elements[5])
                    ), file=data)
        subprocess.call(["mkspk", "-append", "-setup", mkspkfn, "-input", datafn, "-output", "ke%03d.bsp" % ke_version])


