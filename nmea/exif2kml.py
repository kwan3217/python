'''
Created on Aug 3, 2017

@author: chrisj
'''
import exifread
import glob
import time
import datetime
import pytz
import calendar
from timezonefinder import TimezoneFinder

def getgps(tags):
    latdeg = tags["GPS GPSLatitude"].values[0].num / tags["GPS GPSLatitude"].values[0].den
    latmin = tags["GPS GPSLatitude"].values[1].num / tags["GPS GPSLatitude"].values[1].den
    latsec = tags["GPS GPSLatitude"].values[2].num / tags["GPS GPSLatitude"].values[2].den
    latNS = tags["GPS GPSLatitudeRef"].values
    lat = (latdeg + latmin / 60 + latsec / 3600) * (1 if latNS == "N" else -1)
    londeg = tags["GPS GPSLongitude"].values[0].num / tags["GPS GPSLongitude"].values[0].den
    lonmin = tags["GPS GPSLongitude"].values[1].num / tags["GPS GPSLongitude"].values[1].den
    lonsec = tags["GPS GPSLongitude"].values[2].num / tags["GPS GPSLongitude"].values[2].den
    lonEW = tags["GPS GPSLongitudeRef"].values
    lon = (londeg + lonmin / 60 + lonsec / 3600) * (1 if lonEW == "E" else -1)
    GPSDate = tags["GPS GPSDate"].values
    GPSDate = GPSDate.split(":")
    GPSYear = int(GPSDate[0])
    GPSMonth = int(GPSDate[1])
    GPSDay = int(GPSDate[2])
    GPSHour = int(tags["GPS GPSTimeStamp"].values[0].num / tags["GPS GPSTimeStamp"].values[0].den)
    GPSMin = int(tags["GPS GPSTimeStamp"].values[1].num / tags["GPS GPSTimeStamp"].values[1].den)
    GPSSec = int(tags["GPS GPSTimeStamp"].values[2].num / tags["GPS GPSTimeStamp"].values[2].den)
    GPSssec=tags["GPS GPSTimeStamp"].values[2].num / tags["GPS GPSTimeStamp"].values[2].den-GPSSec
    GPSValid=True
    if "GPS GPSStatus" in tags and tags["GPS GPSStatus"].values=="V":
        GPSValid=False
    return (lat,lon,(GPSYear,GPSMonth,GPSDay,GPSHour,GPSMin,GPSSec),GPSValid)

def getimgtime(tags):
    ImageTimestamp = tags["Image DateTime"].values.split(" ")
    ImageDate = ImageTimestamp[0].split(":")
    ImageYear = int(ImageDate[0])
    ImageMonth = int(ImageDate[1])
    ImageDay = int(ImageDate[2])
    ImageTime = ImageTimestamp[1].split(":")
    ImageHour = int(ImageTime[0])
    ImageMin = int(ImageTime[1])
    ImageSec = int(ImageTime[2])
    return (ImageYear,ImageMonth,ImageDay,ImageHour,ImageMin,ImageSec)

def is_img_in_utc(gpstime,imgtime):
    pass

def make_aware(naive,tzname='UTC'):
    naive_datetime=datetime.datetime(naive[0],naive[1],naive[2],naive[3],naive[4],naive[5])
    tz=pytz.timezone(tzname)
    return tz.localize(naive_datetime)

if __name__ == '__main__':
    tf=TimezoneFinder()
    path="/home/jeppesen/Desktop/DCA/"
    with open(path+"Geotag.kml","w") as ouf:
        ouf.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        ouf.write('<kml xmlns="http://www.opengis.net/kml/2.2" xmlns:gx="http://www.google.com/kml/ext/2.2" xmlns:kml="http://www.opengis.net/kml/2.2" xmlns:atom="http://www.w3.org/2005/Atom">\n')
        ouf.write('  <Document>\n')
        ouf.write('    <name>%s</name>\n'%path)
        ouf.write('    <Style id="camera">\n')
        ouf.write('      <IconStyle>\n')
        ouf.write('        <Icon>\n')
        ouf.write('          <href>:/camera_mode.png</href>\n')
        ouf.write('        </Icon>\n')
        ouf.write('      </IconStyle>\n')
        ouf.write('    </Style>\n')
        for infn in glob.iglob(path+"/**/*.[Jj][Pp]*",recursive=True):
            ouf.write("    <!-- %s -->\n"%infn)
            with open(infn,"rb") as inf:
                tags=exifread.process_file(inf)
                #for k,v in tags.items():
                #    if k!="JPEGThumbnail":
                #        print(k,v)
                if "GPS GPSLatitude" in tags:
                    print(infn)
                    (lat,lon,gpstime,gpsvalid)=getgps(tags)
                    tzname=tf.timezone_at(lat=lat,lng=lon)
                    imgtime=getimgtime(tags)
                    g=make_aware(gpstime)
                    l=make_aware(imgtime,tzname=tzname)
                    l_utc=l.astimezone(pytz.timezone("UTC"))

                    diff=(g-l).total_seconds()
                    warning=""
                    if not gpsvalid:
                        warning="Warning! GPS data is marked invalid.<br>"
                    if abs(diff)>8*3600:
                        warning=warning+("Warning! GPS and image timestamps are too different. Old GPS fix?"
                                 "Timestamp difference: %s%02d:%02d:%02d (%ds)<br />"%("+" if diff>0 else "-",abs(diff)//3600,(abs(diff)//60)%60,abs(diff)%60,diff))
                    elif abs(diff)>5:
                        warning=warning+("Warning! GPS and image timestamps are too different. Camera reported image time in UTC?"
                                 "Timestamp difference: %s%02d:%02d:%02d (%ds)<br />"%("+" if diff>0 else "-",abs(diff)//3600,(abs(diff)//60)%60,abs(diff)%60,diff))
                    if warning=="":
                        warning="GPS looks good<br />"
                    width = tags["EXIF ExifImageWidth"].values[0]
                    height = tags["EXIF ExifImageLength"].values[0]
                    if width>height:
                        height=height*640//width
                        width=640
                    else:
                        width=width*640//height
                        height=640
                    warning=("GPS time: %04d-%02d-%02dT%02d:%02d:%02dZ<br />"%gpstime+
                             "Image time: %s<br />"%l.strftime("%Y-%m-%dT%H:%M:%S%z")+warning)
                    print(warning)
                    ouf.write("    <Placemark>\n")
                    ouf.write("      <name>%s</name>\n"%infn.split("/")[-1])
                    ouf.write("      <description><![CDATA[<img src=\"%s\" width=\"%d\" height=\"%d\" />%s]]></description>\n" % (infn,width,height,warning))
                    ouf.write("      <styleUrl>camera</styleUrl>\n")
                    ouf.write("      <TimeStamp><when>%s</when></TimeStamp>\n"%l.strftime("%Y-%m-%dT%H:%M:%S%z"))
                    ouf.write("      <Point>\n")
                    ouf.write("        <coordinates>%f,%f</coordinates>\n"%(lon,lat))
                    ouf.write("      </Point>\n")
                    ouf.write("    </Placemark>\n")
        ouf.write('  </Document>\n')
        ouf.write('</kml>\n')
        