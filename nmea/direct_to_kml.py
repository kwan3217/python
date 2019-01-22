import re
from math import cos,sin,radians,acos
import matplotlib.pyplot as plt
import os
import os.path

re_nmea=re.compile(".*\\$([^*]*)\*([0-9A-F][0-9A-F])?.*")
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
re_pkwne=re.compile(r"""
PKWNE,
(?P<type>[0-9]+),                            #Not sure what this is, perhaps button number?
(?P<hour>[0-9]2):                            #Hour
(?P<minute>[0-9]2):                          #Minute
(?P<second>[0-9]2(?:\.[0-9]+)?),             #Second with optional fractional part
(?P<tag>\S+)                                 #waypoint tag
""",re.VERBOSE)


kml_header="""
<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2" xmlns:gx="http://www.google.com/kml/ext/2.2" xmlns:kml="http://www.opengis.net/kml/2.2" xmlns:atom="http://www.w3.org/2005/Atom">
<Document>
	<name>2018-12-22.kml</name>
	<Style id="s_ylw-pushpin">
		<IconStyle>
			<scale>1.1</scale>
			<Icon>
				<href>http://maps.google.com/mapfiles/kml/pushpin/ylw-pushpin.png</href>
			</Icon>
			<hotSpot x="20" y="2" xunits="pixels" yunits="pixels"/>
		</IconStyle>
	</Style>
	<StyleMap id="multiTrack">
		<Pair>
			<key>normal</key>
			<styleUrl>#multiTrack_n</styleUrl>
		</Pair>
		<Pair>
			<key>highlight</key>
			<styleUrl>#multiTrack_h</styleUrl>
		</Pair>
	</StyleMap>
	<Style id="multiTrack_h">
		<IconStyle>
			<scale>1.2</scale>
			<Icon>
				<href>http://earth.google.com/images/kml-icons/track-directional/track-0.png</href>
			</Icon>
		</IconStyle>
		<LineStyle>
			<color>99ffac59</color>
			<width>8</width>
		</LineStyle>
	</Style>
	<StyleMap id="m_ylw-pushpin">
		<Pair>
			<key>normal</key>
			<styleUrl>#s_ylw-pushpin</styleUrl>
		</Pair>
		<Pair>
			<key>highlight</key>
			<styleUrl>#s_ylw-pushpin_hl</styleUrl>
		</Pair>
	</StyleMap>
	<Style id="s_ylw-pushpin_hl">
		<IconStyle>
			<scale>1.3</scale>
			<Icon>
				<href>http://maps.google.com/mapfiles/kml/pushpin/ylw-pushpin.png</href>
			</Icon>
			<hotSpot x="20" y="2" xunits="pixels" yunits="pixels"/>
		</IconStyle>
	</Style>
	<Style id="multiTrack_n">
		<IconStyle>
			<Icon>
				<href>http://earth.google.com/images/kml-icons/track-directional/track-0.png</href>
			</Icon>
		</IconStyle>
		<LineStyle>
			<color>99ffac59</color>
			<width>6</width>
		</LineStyle>
	</Style>
	<Folder>
		<name>%s</name>
		<open>1</open>
"""
kml_footer="""
	</Folder>
</Document>
</kml>
"""

#Parameter: Name of track (from filename of event log)
track_header="""
		<Placemark>
			<name>%s</name>
			<styleUrl>#multiTrack</styleUrl>
			<gx:Track>
"""
track_footer="""
			</gx:Track>
		</Placemark>
"""
#parameters: Time of track point in ISO8601 2018-12-25T12:34:56Z format
when="""
				<when>%s</when>
"""
# 4. Signed longitude in decimal degrees, positive east (so all locations in lower 48 are negative)
# 5. Signed latitude in decimal degrees, positive north
# 6. Altitude above geoid in meters
coord="""
				<gx:coord>%f %f %f</gx:coord>
"""
#parameters:
# 1. Name of placemark
# 2. Description (write the time of the placemark in original local time format)
# 3. Time of placemark in ISO8601 2018-12-25T12:34:56Z format
# 4. Signed longitude in decimal degrees, positive east (so all locations in lower 48 are negative)
# 5. Signed latitude in decimal degrees, positive north
# 6. Altitude above geoid in meters
pkwne="""
		<Placemark>
			<name>%s</name>
			<description>%s</description>
			<LookAt>
				<gx:TimeStamp><when>%s</when>
</gx:TimeStamp>
			</LookAt>
			<styleUrl>#m_ylw-pushpin</styleUrl>
			<Point>
				<altitudeMode>absolute</altitudeMode>
				<gx:drawOrder>1</gx:drawOrder>
				<coordinates>%f,%f,%f</coordinates>
			</Point>
		</Placemark>
"""

lat = []
lon = []
time = []
date=[]
alt=[]
spd=[]
lineno = 0
bad_alt=True

def check_checksum(infn):
    oufn = infn[:-4] + ".fix.nmea"
    global lat,lon,time,alt,spd,lineno,bad_alt
    old_lat = None
    old_lon = None
    old_time = None
    old_date = None
    old_spd = None
    printlat=""
    printlon=""
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

    def handle_pos(match):
        global lat,lon,alt,time,date,spd
        this_lat = get_lat(match.group('lat'), match.group('NS'))
        this_lon = get_lat(match.group('lon'), match.group('EW'))
        this_alt = match.group('alt')
        geoid = match.group('geoid')
        if this_alt is not None:
            bad_alt = (this_alt == "-" + geoid) or ("-" + this_alt == geoid)
            if bad_alt:
                print("Bad altitude on line ", lineno, data)
                return
        this_time = sod(match.group('time'))
        if len(time)>0:
            dt = this_time-time[-1]
            dd = dist(lla2xyz(this_lat, this_lon), lla2xyz(lat[-1], lon[-1]))
            if dt == 0:
                if dd > 10:
                    print("Position step on line ", lineno)
                return
            else:
                this_spd = dd / dt
                this_acc = (this_spd - spd[-1]) / dt
            speed.append(this_spd)
            if abs(this_acc) > 99:
                print("Position glitch on line ", lineno, data)
                return
        lat.append(this_lat)
        lon.append(this_lon)
        time.append(this_time)
        alt.append(None)
        spd.append(this_spd)

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

    def handle_pkwne(match):
        global printlat,printlon
        tag=match.group('name')
        return "GPWPL,%s,%s,%s"%(printlat,printlon,tag)

    def calc_checksum(data):
        cksum_calc = 0x00
        for char in data:
            cksum_calc = cksum_calc ^ ord(char)
        return cksum_calc

    with open(oufn,"w",encoding="cp437") as ouf:
        with open(infn,"r",encoding="cp437") as inf:
            npos=0
            lineno=1
            for line in inf:
                if lineno == 123116:
                    print("break!")
                result=re_nmea.match(line)
                if result is not None:
                    data=result.group(1)
                    cksum_stored=result.group(2)
                    cksum_calc=calc_checksum(data)
                    if cksum_stored is None or int("0x"+cksum_stored,16)==cksum_calc:
                        write_line=True
                        gga_match = re_gga.match(data)
                        if gga_match is not None:
                            if not handle_pos(gga_match):
                                write_line=False
                            else:
                                npos+=1
                        elif data[2:5]=="GGA":
                            #print("Bad GGA at line ",lineno,data)
                            write_line=False
                        rmc_match = re_rmc.match(data)
                        if rmc_match is not None:
                            write_line=not bad_alt
                            if not handle_pos(rmc_match):
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
                        pkwne_match = re_pkwne.match(data)
                        if pkwne_match is not None:
                            data=handle_pkwne(pkwne_match)
                        elif data[0:5]=="PKWNE":
                            print("Bad PKWNE at line ",lineno,data)
                            write_line=False
                        if data[0:4]=="PKWN" and data[4]!="E":
                            write_line=False #PKWN data is probably valid, but Google Earth doesn't care
                        if write_line:
                            data="GP"+data[2:]
                            cksum_new=calc_checksum(data)
                            print("$"+data+"*%02X"%cksum_new,file=ouf)
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
    #infns=glob.glob("/home/chrisj/Florida7/*.nmea")
    infns=glob.glob("C:\\Users\\Chris Jeppesen\\Desktop\\Florida7SensorNMEA\\Saturday*.nmea")
    for infn in infns:
        print(infn)
        if "fix" not in infn:
            check_checksum(infn)
            print(high_alt,high_lineno)
