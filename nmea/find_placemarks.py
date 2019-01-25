kml_header="""<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2" xmlns:gx="http://www.google.com/kml/ext/2.2" xmlns:kml="http://www.opengis.net/kml/2.2" xmlns:atom="http://www.w3.org/2005/Atom">
<Document>
	<name>%s</name>
	<Style id="s_ylw-pushpin">
		<IconStyle>
			<scale>1.1</scale>
			<Icon>
				<href>http://maps.google.com/mapfiles/kml/pushpin/ylw-pushpin.png</href>
			</Icon>
			<hotSpot x="20" y="2" xunits="pixels" yunits="pixels"/>
		</IconStyle>
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
"""
kml_footer="""</Document>
</kml>
"""
placemark="""	<Placemark>
		<name>%s</name>
		<LookAt>
			<gx:TimeSpan>
				<begin>%s</begin>
				<end>%s</end>
			</gx:TimeSpan>
			<longitude>%18.13f</longitude>
			<latitude>%18.13f</latitude>
			<altitude>%18.13f</altitude>
			<heading>7.898471913806163e-07</heading>
			<tilt>44.99871477414683</tilt>
			<range>999.9406734086464</range>
			<gx:altitudeMode>absolute</gx:altitudeMode>
		</LookAt>
		<styleUrl>#m_ylw-pushpin</styleUrl>
		<Point>
			<gx:drawOrder>1</gx:drawOrder>
			<altitudeMode>absolute</altitudeMode>
			<coordinates>%18.13f,%18.13f,%18.13f</coordinates>
		</Point>
	</Placemark>
"""

import os.path
import sys
import re
gga=re.compile("\$..GGA,([^,]*),([^,]*),([^,]*)")
if len(sys.argv)<2:
    infn="android_gps/ToCaliforniaNMEA20180411T052713.txt"
else:
    infn=sys.argv[1]
oufn = os.path.basename(infn)
oufn = oufn.split(".")
basebase = ".".join(oufn[:-1])
oufn = os.path.dirname(infn) + "/" + basebase + ".kml"
print(infn)
print(oufn)
print(basebase)
with open(infn,"r") as inf:
    with open(oufn,"w") as ouf:
        ouf.write(kml_header % basebase)
        #Search for the placemarks here
        for line in inf:
            #Get rid of timestamp, it is local time anyway
            line=line[12:]
            #Does it match the GGA pattern?
            if line[3:6]=="GGA":
                parts=line.split(",")
                this_q=int(parts[6])
                if this_q>0:
                    this_time=float(parts[1])
                    this_lat=float(parts[2])
#                    if parts[3]!="N":
                    this_lon=float(parts[4])
#                if parts[5]!="E":
                    this_alt=float(parts[9])
                    #convert to ISO time format
                    h=int(this_time/10000)
                    this_time-=h*10000
                    m=int(this_time/100)
                    this_time-=m*100
                    s=this_time
                    time="%02d:%02d:%05.2fZ"%(h,m,s)
                    d=int(this_lat/100)
                    this_lat-=d*100
                    m=this_lat
                    lat=d+m/60.0
                    if parts[3]!="N":
                        lat=-lat
                    d=int(this_lon/100)
                    this_lon-=d*100
                    m=this_lon
                    lon=d+m/60.0
                    if parts[5]!="E":
                        lon=-lon
                    alt=this_alt
            if line[3:6]=="RMC":
                parts = line.split(",")
                this_q = parts[2]
                if this_q =="A":
                    this_time = float(parts[1])
                    this_lat = float(parts[3])
                    this_lon = float(parts[5])
                    this_date = int(parts[9])
                    #convert to ISO time format
                    h=int(this_time/10000)
                    this_time-=h*10000
                    m=int(this_time/100)
                    this_time-=m*100
                    s=this_time
                    time="%02d:%02d:%05.2fZ"%(h,m,s)
                    # convert to ISO date format
                    d = int(this_date // 10000)
                    this_date -= d * 10000
                    m = int(this_date // 100)
                    this_date -= m * 100
                    y = this_date
                    date = "20%02d-%02d-%02dT" % (y,m,d)
                    #convert lat to decimal degrees
                    d = int(this_lat / 100)
                    this_lat -= d * 100
                    m = this_lat
                    lat = d + m / 60.0
                    if parts[4] != "N":
                        lat = -lat
                    #convert lon to decimal degrees
                    d = int(this_lon / 100)
                    this_lon -= d * 100
                    m = this_lon
                    lon = d + m / 60.0
                    if parts[6] != "E":
                        lon = -lon
            if line[1:6]=="PKWNE":
                tag=line.split(",")[1].split("*")[0]
                ouf.write(placemark%(tag,
                                     date+time,
                                     date+time,
                                     lon,lat,alt,
                                     lon,lat,alt))
        ouf.write(kml_footer)

