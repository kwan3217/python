import re
from math import cos,sin,radians,acos
import matplotlib.pyplot as plt

re_nmea=re.compile(".*\\$([^*]*)\*([0-9A-F][0-9A-F]).*")
re_gga=re.compile(r"""
(?P<gpstype>..)GGA,  #Type of GPS system
(?P<time>[0-9]{6,6}(?:\.[0-9]+)),               #Time, with optional fraction of second
(?P<lat>[0-9]{4,4}\.[0-9]+),(?P<NS>[NS]),         #Latitude and hemisphere
(?P<lon>[0-9]{5,5}\.[0-9]+),(?P<EW>[EW]),         #Longitude and hemisphere
(?P<fixq>[0-8])?,                             #(Optional) fix quality
(?P<nsat>[0-9]+)?,                            #(Optional) number of satellites in the fix
(?P<HDOP>[0-9]+.[0-9]+),                      #(Optional) Horizontal dillution of precision
(?P<alt>-?[0-9]+.[0-9]+)?,(?P<altUnit>[M])?,      #(Optional) altitude and units
(?P<geoid>-?[0-9]+.[0-9]+)?,(?P<geoidUnit>[M])?,  #(Optional) geoid altitude above ellipsoid and units
, #Slots for DGPS update time and station number
""",re.VERBOSE)
re_rmc=re.compile(r"""
(?P<gpstype>..)RMC,  #Type of GPS system
(?P<time>[0-9]{6}(?:\.[0-9]+)),               #Time, with optional fraction of second
(?P<valid>[AV])?,                            #fix valid
(?P<lat>[0-9]{4}\.[0-9]+),(?P<NS>[NS]),      #Latitude and hemisphere
(?P<lon>[0-9]{5}\.[0-9]+),(?P<EW>[EW]),      #Longitude and hemisphere
(?P<spd>[0-9]+.[0-9]+)?,                     #(Optional) speed
(?P<hdg>[0-9]+.[0-9]+)?,                     #(Optional) true heading
(?P<date>[0-9]{6}),                          #Date in DDMMYY
(?P<magvar>[0-9]{5}\.[0-9]+)?,(?P<magvarDir>[EW])?,      #(optional) magnetic variation
""",re.VERBOSE)

infn="/media/jeppesen/PassportData/Florida5/Southwest/20150303/SOUTHW05.nmea"
oufn=infn[:-4]+"fix.nmea"
old_lat=None
old_lon=None
old_time=None
old_date=None
old_spd=None

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
    if abs(d)>1:
        print("dot product out of range")
    return acos(dotp(a,b))*6378137

speed=[]
bad_alt=True

with open(oufn,"w",encoding="cp437") as ouf:
    with open(infn,"r",encoding="cp437") as inf:
        lineno=1
        for line in inf:
            result=re_nmea.match(line)
            if result is not None:
                data=result.group(1)
                cksum_stored=int("0x"+result.group(2),16)
                cksum_calc=0x00
                for char in data:
                    cksum_calc=cksum_calc^ord(char)
                if cksum_stored==cksum_calc:
                    write_line=True
                    gga_match = re_gga.match(data)
                    if gga_match is not None:
                        lat=get_lat(gga_match.group('lat'),gga_match.group('NS'))
                        lon=get_lat(gga_match.group('lon'),gga_match.group('EW'))
                        time=sod(gga_match.group('time'))
                        alt=gga_match.group('alt')
                        geoid=gga_match.group('geoid')
                        bad_alt=(alt=="-"+geoid) or ("-"+alt==geoid)
                        if bad_alt:
                            write_line = False
                            print("Bad altitude on line ", lineno, data)
                        if old_lat is not None:
                            dt=time-old_time
                            dd=dist(lla2xyz(lat,lon),lla2xyz(old_lat,old_lon))
                            if dt==0:
                                if dd>0:
                                    print("Position step on line ",lineno)
                                    write_line=False
                                    spd=-1
                                    acc=-1
                                else:
                                    spd=0
                                    acc=0
                            else:
                                spd=dd/dt
                                acc=(spd-old_spd)/dt
                            speed.append(spd)
                            if abs(acc)>99:
                                print("Position glitch on line ",lineno,data)
                                write_line=False
                            else:
                                old_lat=lat
                                old_lon=lon
                                old_time=time
                                old_spd=spd
                        else:
                            old_lat=lat
                            old_lon=lon
                            old_time=time
                            old_spd=0
                    elif data[2:5]=="GGA":
                        print("Bad GGA at line ",lineno,data)
                        write_line=False
                    rmc_match = re_rmc.match(data)
                    if rmc_match is not None:
                        write_line=not bad_alt
                        lat=get_lat(rmc_match.group('lat'),rmc_match.group('NS'))
                        lon=get_lat(rmc_match.group('lon'),rmc_match.group('EW'))
                        time=sod(rmc_match.group('time'))
                        if old_lat is not None:
                            dt=time-old_time
                            dd=dist(lla2xyz(lat,lon),lla2xyz(old_lat,old_lon))
                            if dt==0:
                                if dd>0:
                                    print("Position step on line ",lineno)
                                    write_line=False
                                    spd=-1
                                    acc = -1
                                else:
                                    spd = 0
                                    acc = 0
                            else:
                                spd = dd / dt
                                acc = (spd - old_spd) / dt
                            speed.append(spd)
                            if abs(acc) > 99:
                                print("Position glitch on line ", lineno, data)
                                write_line = False
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
                    elif data[2:5]=="RMC":
                        print("Bad RMC at line ",lineno,data)
                        write_line=False
                    if write_line:
                        ouf.write(line)
                else:
#                    print("Problem with checksum in line %d"%lineno)
                    pass
            lineno+=1
            if lineno==25173:
                print("On line ",lineno)

plt.plot(speed)
plt.show()
