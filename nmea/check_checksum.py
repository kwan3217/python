import re
from math import cos,sin,radians,acos
import matplotlib.pyplot as plt
import os
import os.path

re_nmea=re.compile(".*\\$([^*]*)\*([0-9A-F][0-9A-F]).*")
re_gga=re.compile(r"""
(?P<gpstype>..)GGA,  #Type of GPS system
(?P<time>[0-9]{6}(?:\.[0-9]+)?),               #Time, with optional fraction of second
(?P<lat>[0-9]{4}\.[0-9]+),(?P<NS>[NS]),         #Latitude and hemisphere
(?P<lon>[0-9]{5}\.[0-9]+),(?P<EW>[EW]),         #Longitude and hemisphere
(?P<fixq>[0-8])?,                             #(Optional) fix quality
(?P<nsat>[0-9]+)?,                            #(Optional) number of satellites in the fix
(?P<HDOP>[0-9]+.[0-9]+)?,                      #(Optional) Horizontal dillution of precision
(?P<alt>-?[0-9]+.[0-9]+)?,(?P<altUnit>[M])?,      #(Optional) altitude and units
(?P<geoid>-?[0-9]+.[0-9]+)?,(?P<geoidUnit>[M])?,  #(Optional) geoid altitude above ellipsoid and units
(?P<DGPStime>[^,]*),(?P<DGPSsta>.*) #Slots for DGPS update time and station number
""",re.VERBOSE)
re_rmb=re.compile(r"""
(?P<gpstype>..)RMB,  #Type of GPS system
(?P<valid>[AV])?,                             #Data valid
(?P<xtk>[0-9]+\.[0-9]+),(?P<xtkdir>[LR]),     #Crosstrack error and direction to steer to correct
(?P<origin>\S+),                              #Origin waypoint name
(?P<dest>\S+),                                #Destination waypoint name
(?P<lat>[0-9]{4,4}\.[0-9]+),(?P<NS>[NS]),     #Destination Latitude and hemisphere
(?P<lon>[0-9]{5,5}\.[0-9]+),(?P<EW>[EW]),     #Destination Longitude and hemisphere
(?P<range_nmi>[0-9]+.[0-9]+),                 #Range to destination, nautical miles
(?P<abs_bear>[0-9]+.[0-9]+),                  #absolute bearing, degrees east of true north
(?P<vmg>[0-9]+.[0-9]+),                       #speed towards the destination, knots
(?P<alt>-?[0-9]+.[0-9]+)?,(?P<altUnit>[M])?,  #(Optional) altitude and units
(?P<arrive>[AV])?                             #Arrival alarm
""",re.VERBOSE)
re_wpl=re.compile(r"""
(?P<gpstype>..)WPL,                           #Type of GPS system
(?P<lat>[0-9]{4}\.[0-9]+),(?P<NS>[NS]),       #Latitude and hemisphere
(?P<lon>[0-9]{5}\.[0-9]+),(?P<EW>[EW]),       #Longitude and hemisphere
(?P<name>\S+)                                 #waypoint name
""",re.VERBOSE)
re_rmc=re.compile(r"""
(?P<gpstype>..)RMC,  #Type of GPS system
(?P<time>[0-9]{6}(?:\.[0-9]+)?),               #Time, with optional fraction of second
(?P<valid>[AV])?,                            #fix valid
(?P<lat>[0-9]{4}\.[0-9]+),(?P<NS>[NS]),      #Latitude and hemisphere
(?P<lon>[0-9]{5}\.[0-9]+),(?P<EW>[EW]),      #Longitude and hemisphere
(?P<spd>[0-9]+.[0-9]+)?,                     #(Optional) speed
(?P<hdg>[0-9]+.[0-9]+)?,                     #(Optional) true heading
(?P<date>[0-9]{6}),                          #Date in DDMMYY
(?P<magvar>[0-9]+\.[0-9]+)?,(?P<magvarDir>[EW])?      #(optional) magnetic variation
""",re.VERBOSE)

old_lat = None
old_lon = None
old_time = None
old_date = None
old_spd = None
high_alt=None
high_lineno=None
wpl_dict = {}
lineno = 0
bad_alt=True

def check_checksum(infn):
    oufn = infn + "fix.nmea"
    global old_lat,old_lon,old_time,old_date,old_spd,wpl_dict,lineno,bad_alt,high_alt,high_lineno
    old_lat = None
    old_lon = None
    old_time = None
    old_date = None
    old_spd = None
    wpl_dict={}
    lineno=0
    def get_lat(latin,hemi):
        lat=float(latin)
        latdeg=int(lat)//100
        min=lat-latdeg*100
        lat=latdeg+min/60
        if hemi=="S" or hemi=="W":
            lat=-lat
        return lat

    def sod(timein):
        h=int(timein[0:2])
        m=int(timein[2:4])
        s=float(timein[4:])
        return h*3600+m*60+s

    def lla2xyz(latdeg,londeg):
        return(cos(radians(latdeg))*cos(radians(londeg)),
               cos(radians(latdeg))*sin(radians(londeg)),
               sin(radians(latdeg)))

    def dotp(a,b):
        return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]

    def dist(a,b):
        if(a==b):
            return 0
        d=dotp(a,b)
        if d>1:
            #print("dot product out of range")
            return 0
        if d<-1:
            print("dot product out of range (negative)")
            return 6378137*3.1415926535897932
        return acos(d)*6378137

    speed=[]
    bad_alt=True

    def handle_gga(gga_match):
        """

        :param gga_match:
        :return:
        """
        global old_lat,old_lon,old_time,old_spd,bad_alt,lineno,high_alt,high_lineno
        lat = get_lat(gga_match.group('lat'), gga_match.group('NS'))
        lon = get_lat(gga_match.group('lon'), gga_match.group('EW'))
        time = sod(gga_match.group('time'))
        alt = gga_match.group('alt')
        geoid = gga_match.group('geoid')
        bad_alt = (alt == "-" + geoid) or ("-" + alt == geoid)
        if bad_alt:
            print("Bad altitude on line ", lineno, data)
            return False
        if old_lat is not None:
            dt = time - old_time
            dd = dist(lla2xyz(lat, lon), lla2xyz(old_lat, old_lon))
            if dt == 0:
                if dd > 10:
                    print("Position step on line ", lineno)
                    return False
                else:
                    spd = 0
                    acc = 0
            else:
                spd = dd / dt
                acc = (spd - old_spd) / dt
            speed.append(spd)
            if abs(acc) > 99:
                print("Position glitch on line ", lineno, data)
                return False
            else:
                old_lat = lat
                old_lon = lon
                old_time = time
                old_spd = spd
                if high_alt is not None and alt>high_alt:
                    high_alt=alt
                    high_lineno=lineno
        else:
            old_lat = lat
            old_lon = lon
            old_time = time
            old_spd = 0
            if alt is not None:
                high_alt = alt
                high_lineno = lineno
        return True

    def handle_rmc(rmc_match):
        global old_lat,old_lon,old_time,old_spd,lineno,high_alt,high_lineno
        lat = get_lat(rmc_match.group('lat'), rmc_match.group('NS'))
        lon = get_lat(rmc_match.group('lon'), rmc_match.group('EW'))
        time = sod(rmc_match.group('time'))
        if old_lat is not None:
            dt = time - old_time
            dd = dist(lla2xyz(lat, lon), lla2xyz(old_lat, old_lon))
            if dt == 0:
                if dd > 0:
#                    print("Position step on line ", lineno)
                    return False
                else:
                    spd = 0
                    acc = 0
            else:
                spd = dd / dt
                acc = (spd - old_spd) / dt
            speed.append(spd)
            if abs(acc) > 99:
#                print("Position glitch on line ", lineno, data)
                return False
            else:
                old_lat = lat
                old_lon = lon
                old_time = time
                old_spd = spd
        else:
            old_lat = lat
            old_lon = lon
            old_time = time
            old_spd = 0
        return True

    def handle_wpl(match):
        global wpl_dict
        lat=match.group('lat')
        ns=match.group('NS')
        lon = match.group('lon')
        ew=match.group('EW')
        name=match.group('name')
        if match.group('name') not in wpl_dict:
            print("Found waypoint %s at (%s%s,%s%s)"%(name,lat,ns,lon,ew))
            wpl_dict[name]=(lat,ns,lon,ew)
            return True
        else:
            if (lat,ns,lon,ew)!=wpl_dict[name]:
                print_tuple=(name,)+wpl_dict[name]+(lat,ns,lon,ew)
                print("Waypoint %s changed from (%s%s,%s%s) to (%s%s,%s%s)"%print_tuple)
                wpl_dict[name]=(lat,lon)
                return True
            else:
                return False

    with open(oufn,"w",encoding="cp437") as ouf:
        with open(infn,"r",encoding="cp437") as inf:
            npos=0
            lineno=1
            for line in inf:
                result=re_nmea.match(line)
                if result is not None:
                    data=result.group(1)
                    cksum_stored=result.group(2)
                    cksum_calc=0x00
                    for char in data:
                        cksum_calc=cksum_calc^ord(char)
                    if int("0x"+cksum_stored,16)==cksum_calc:
                        write_line=True
                        gga_match = re_gga.match(data)
                        if gga_match is not None:
                            if not handle_gga(gga_match):
                                write_line=False
                            else:
                                npos+=1
                        elif data[2:5]=="GGA":
                            #print("Bad GGA at line ",lineno,data)
                            write_line=False
                        rmc_match = re_rmc.match(data)
                        if rmc_match is not None:
                            write_line=not bad_alt
                            if not handle_rmc(rmc_match):
                                write_line=False
                            else:
                                npos+=1
                        elif data[2:5]=="RMC":
                            #print("Bad RMC at line ",lineno,data)
                            write_line=False
                        wpl_match = re_wpl.match(data)
                        if wpl_match is not None:
                            if not handle_wpl(wpl_match):
                                write_line=False
                        elif data[2:5]=="WPL":
                            print("Bad WPL at line ",lineno,data)
                            write_line=False
                        if data[0:4]=="PKWN":
                            write_line=False #PKWN data is probably valid, but Google Earth doesn't care
                        if write_line:
                            print("$"+data+"*"+cksum_stored,file=ouf)
                    else:
                        print("Problem with checksum in line %d"%lineno)
                        pass
                lineno+=1
    print("Number of positions found in %s: %d"%(os.path.basename(infn),npos))
    if npos==0:
        pass
        #os.remove(infn)
        #os.remove(oufn)

if __name__=="__main__":
    import glob
    #infns=glob.glob("/home/jeppesen/workspace/Data/recover_gps/SensorLogs/NMEA*.txt")
    infns=glob.glob("/home/jeppesen/Desktop/MNSensorLogs/*NMEA*.txt")
    for infn in infns:
        if ".fix." not in infn:
            check_checksum(infn)
            print(high_alt,high_lineno)