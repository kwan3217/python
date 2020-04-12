"""
Make the equivalent to DE430 for the Kerbal system.
"""
import subprocess
import numpy as np
import os
import krpc
from collections import namedtuple

#Spice number follows convention of real solar system, except with 1000 added. This allows real and Kerbal
#solar systems to be loaded at the same time.
#If no orbital elements are provided, we will create a kernel with a constant zero vector state relative to parent
#We use the name Kerbol instead of Sun so that it doesn't conflict with the real Sun. Likewise Kerbol system barycenter
#rather than solar system barycenter (should it be Kerbolar system? Or is that pushing consistency too far?)
ElementTuple=namedtuple("ElementTuple",["a","e","i","lan","ap","M"])
BodyTuple=namedtuple("BodyTuple",["spice_name","gm","re","i_parent","elements","rot_period"])
SystemTuple=namedtuple("SystemTuple",["spice_name","ksp_name","i_parent"])

sys_structure={
#   Sp#   Spice Name                              KSP name  Parent
    1000: SystemTuple("Kerbol system barycenter", "Sun",       0),
    1199: SystemTuple("Moho",                     "Moho",   1001),
    1299: SystemTuple("Eve",                      "Eve",    1002),
    1201: SystemTuple("Gilly",                    "Gilly",  1002),
    1399: SystemTuple("Kerbin",                   "Kerbin", 1003),
    1301: SystemTuple("Mun",                      "Mun"   , 1003),
    1302: SystemTuple("Minmus",                   "Minmus", 1003),
    1499: SystemTuple("Duna",                     "Duna"  , 1004),
    1401: SystemTuple("Ike",                      "Ike"   , 1004),
    1599: SystemTuple("Dres",                     "Dres"  , 1005),
    1699: SystemTuple("Jool",                     "Jool"  , 1006),
    1601: SystemTuple("Laythe",                   "Laythe", 1006),
    1602: SystemTuple("Vall",                     "Vall"  , 1006),
    1603: SystemTuple("Tylo",                     "Tylo"  , 1006),
    1604: SystemTuple("Bop",                      "Bop"   , 1006),
    1605: SystemTuple("Pol",                      "Pol"   , 1006),
    1799: SystemTuple("Eeloo",                    "Eeloo" , 1007),
    1010: SystemTuple("Kerbol",                   "Sun"   , 1000),
    1001: SystemTuple("Moho Barycenter",          "Moho"  , 1000),
    1002: SystemTuple("Eve Barycenter",           "Eve"   , 1000),
    1003: SystemTuple("Kerbin Barycenter",        "Kerbin", 1000),
    1004: SystemTuple("Duna Barycenter",          "Duna"  , 1000),
    1005: SystemTuple("Dres Barycenter",          "Dres"  , 1000),
    1006: SystemTuple("Jool Barycenter",          "Jool"  , 1000),
    1007: SystemTuple("Eeloo Barycenter",         "Eeloo" , 1000)
}

ksp_version="1.0.5"
ke_version=105 #Intended to be xyz for ksp version x.y.z, which will be a problem if any of those exceed one digit

elements_header="""
Spice kernel covering the natural objects in the Kerbol system. This segment
is for the position of {spice_name} (Spice ID {i_spice}) relative to {parent_name} (Spice ID {i_parent}).

Using data extracted via kRPC, extracted from KSP version {ksp_version}

The K0001 epoch (equivalent to J2000) is used, as this matches the in-game
time scale. The clock ticks at 1 SI second per second, and starts at 
Year 1, Day 1, 00:00:00 UT which is the exact instant that the game starts
and matches all UT in seconds in the persistence file.

The actual data is as follows:
Parent GM/(km**3/s**2):          {parent_gm:.17e}
Parent radius/km:                {parent_re:.17e}
Semimajor axis/km:               {a:.17e}
Eccentricity:                    {e:.17e}
Inclination/rad:                 {i:.17e}
Argument of Periapse/rad:        {ap:.17e}
Longitude of Ascending Node/rad: {lan:.17e}
Mean Anomaly at epoch UT=0s/rad: {M:.17e}

Note that KSP uses meters and radians. Spice uses km and
no explicit angle unit (uses unit vectors instead).

\\begindata
INPUT_DATA_TYPE = 'ELEMENTS'
OUTPUT_SPK_TYPE = 15
OBJECT_ID={i_spice}
CENTER_ID={i_parent}
REF_FRAME_NAME='J2000'
PRODUCER_ID='C. Jeppesen, Kwan Systems'
DATA_ORDER='EPOCH A E INC PER NOD MEAN'
TIME_WRAPPER='# ETSECONDS'
INPUT_DATA_UNITS = ('ANGLES=radians' 'DISTANCES=km')
DATA_DELIMITER=','
LINES_PER_RECORD=1
CENTER_GM={parent_gm:.17e}
CENTER_J2=0
CENTER_EQ_RADIUS={parent_re:.17e}
PRECESSION_TYPE='NO PRECESSION'
CENTER_POLE_RA=0
CENTER_POLE_DEC=90
START_TIME='2000-Jan-01 12:00:00 TDB'
STOP_TIME='3000-Jan-01 12:00:00 TDB'
LEAPSECONDS_FILE='kerb0000.tls'
"""

state_header="""
Spice kernel covering the natural objects in the Kerbol system. This segment
is for the position of {spice_name} (Spice ID {i_spice}) relative to {parent_name} (Spice ID {i_parent})
and uses a constant zero vector offset.

Using data from spreadsheet, extracted from KSP version {ksp_version}

The K0001 epoch (equivalent to J2000) is used, as this matches the in-game
time scale. The clock ticks at 1 SI second per second, and starts at 
Year 1, Day 1, 00:00:00 UT which is the exact instant that the game starts
and matches all UT in seconds in the persistence file.

\\begindata
INPUT_DATA_TYPE = 'STATES'
OUTPUT_SPK_TYPE = 13
OBJECT_ID={i_spice}
CENTER_ID={i_parent}
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

def mk_ke(ke_version,ksp_version,elements):
    try:
        os.remove("ke%03d.bsp"%ke_version)
    except FileNotFoundError:
        pass
    with open("kerb%03d.tpc"%ke_version,"wt") as tpc:
        print("\\begindata",file=tpc)
        for i_spice in elements:
            name=elements[i_spice].spice_name
            gm=elements[i_spice].gm
            re=elements[i_spice].re
            i_parent=elements[i_spice].i_parent
            if i_parent!=0:
                parent_name=elements[i_parent].spice_name
                parent_gm=elements[i_parent].gm
                parent_re=elements[i_parent].re
            else:
                parent_name="Inertial Origin"
                parent_gm=0
                parent_re=0
            this_elements = elements[i_spice].elements
            rot_period=elements[i_spice].rot_period
            print("NAIF_BODY_CODE+=(%d)" %i_spice,file=tpc)
            print("NAIF_BODY_NAME+=('%s')"%name.strip().upper(),file=tpc)
            if(i_spice>=1010):
                print("BODY%04d_RADII=(%d,%d,%d)"%(i_spice,re/1e3,re/1e3,re/1e3),file=tpc)
            print("BODY%04d_GM=(%.17e)" %(i_spice,gm/1e9), file=tpc)
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
            mkspkfn = "mkspk%03d_%04d.txt" % (ke_version, i_spice)
            datafn = "data%03d_%04d.txt" % (ke_version, i_spice)
            with open(mkspkfn,"wt") as mkspk:
                with open(datafn,"wt") as data:
                    if this_elements is None:
                        print(state_header.format(spice_name=name,i_spice=i_spice,parent_name=parent_name,i_parent=i_parent,ksp_version=ksp_version),file=mkspk)
                        for t in (K0001,K1001):
                            print("%f,0,0,0,0,0,0"%t,file=data)
                    else:
                        print(elements_header.format(spice_name=name,i_spice=i_spice,parent_name=parent_name,i_parent=i_parent,ksp_version=ksp_version,
                                                     parent_gm=parent_gm/1e9,parent_re=parent_re/1e3,
                                                     a=this_elements.a/1e3,e=this_elements.e,i=this_elements.i,
                                                     lan=this_elements.lan,ap=this_elements.ap,M=this_elements.M), file=mkspk)
                        print("0,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e" % (
                            this_elements[0]/1e3,this_elements[1],this_elements[2],
                            this_elements[3],this_elements[4],np.degrees(this_elements[5])
                        ), file=data)
            subprocess.call(["mkspk","-append","-setup","kerbol_mkspk.txt","-input","kerbol_data.txt","-output","ke%03d.bsp"%ke_version])

def extract_elements(sys_structure):
    conn = krpc.connect(address="127.0.0.1", name="ephemeris")
    ksp_bodies = conn.space_center.bodies
    elements={}
    for i_spice in sys_structure:
        s=sys_structure[i_spice]
        gm = ksp_bodies[s.ksp_name].gravitational_parameter
        re = ksp_bodies[s.ksp_name].equatorial_radius
        rot_period=ksp_bodies[s.ksp_name].rotational_period
        if i_spice==1000 or i_spice==1010:
            #Sun
            el=None
        elif i_spice>1100:
            #planet/moon
            i_planet=(i_spice-1000)//100
            i_moon=i_spice%100
            if i_moon==99:
                #planet state relative to barycenter
                el=None
            else:
                el=ElementTuple(a=ksp_bodies[s.ksp_name].orbit.semi_major_axis,
                                e=ksp_bodies[s.ksp_name].orbit.eccentricity,
                                i=ksp_bodies[s.ksp_name].orbit.inclination,
                                lan=ksp_bodies[s.ksp_name].orbit.longitude_of_ascending_node,
                                ap=ksp_bodies[s.ksp_name].orbit.argument_of_periapsis,
                                M=ksp_bodies[s.ksp_name].orbit.mean_anomaly_at_epoch)
        else:
            #planet barycenter index
            el = ElementTuple(a=ksp_bodies[s.ksp_name].orbit.semi_major_axis,
                              e=ksp_bodies[s.ksp_name].orbit.eccentricity,
                              i=ksp_bodies[s.ksp_name].orbit.inclination,
                              lan=ksp_bodies[s.ksp_name].orbit.longitude_of_ascending_node,
                              ap=ksp_bodies[s.ksp_name].orbit.argument_of_periapsis,
                              M=ksp_bodies[s.ksp_name].orbit.mean_anomaly_at_epoch)
        elements[i_spice]=BodyTuple(spice_name=s.spice_name,
                                    gm=gm,
                                    re=re,
                                    i_parent=s.i_parent,
                                    elements=el,
                                    rot_period=rot_period)
    return elements

if __name__=="__main__":
    elements=extract_elements(sys_structure)
    for i_spice in elements:
        print(i_spice,":",elements[i_spice])
    mk_ke(ke_version,ksp_version,elements)
