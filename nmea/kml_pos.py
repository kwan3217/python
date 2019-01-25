"""
Read a KML file and return a position given a time
"""
from xml.dom import minidom
from datetime import datetime
import pytz
import glob
import exifread
from timezonefinder import TimezoneFinder
import re
tf = TimezoneFinder()
import os.path

def parseKML(infn):
    """
    Find all tracks in a KML file

    :param infn: Filename of KML file to parse
    :return: A list of tracks, in the order that they are found.
    A track is represented by a tuple:
      # name - name of track
      # firstTime - datetime object representing earliest time in the track
      # lastTime - datetime object representing latest time in the track
      # altMode - string representing altitude mode of height coordinate
      # track - a dictionary of tuples, indexed by time. Each tuple is:
         # East longitude in degrees
         # North latitude in degrees
         # Altitude above reference in meters

    All times have tzinfo set to pytz.utc. Currently the date parser only handles times with time zone specified as Z
    (not even +00:00 will work for now).
    """
    xmldoc=minidom.parse(infn)
    tracks=[]
    for tracknode in xmldoc.getElementsByTagName('gx:Track'):
        namenode = tracknode.parentNode.getElementsByTagName("name")[0]
        name=namenode.firstChild.data
        times = []
        coords=[]
        firstTime=None
        lastTime=None
        altModenode=tracknode.getElementsByTagName('altitudeMode')
        if altModenode is not None and len(altModenode)>0:
            altMode=altModenode[0].firstChild.data
        else:
            altMode=None
        for whennode in tracknode.getElementsByTagName('when'):
            try:
                time=pytz.utc.localize(datetime.strptime(whennode.firstChild.data,'%Y-%m-%dT%H:%M:%SZ'))
            except ValueError:
                try:
                    time=pytz.utc.localize(datetime.strptime(whennode.firstChild.data,'%Y-%m-%dT%H:%M:%S.%fZ'))
                except ValueError:
                    #Google Earth will write a bare date if the UTC time is exactly midnight
                    time = pytz.utc.localize(datetime.strptime(whennode.firstChild.data, '%Y-%m-%d'))
            times.append(time)
            if firstTime is None:
                firstTime=time
                lastTime=time
            else:
                if firstTime>time:
                    firstTime=time
                if lastTime<time:
                    lastTime=time
        for coordnode in tracknode.getElementsByTagName('gx:coord'):
            coordstr=coordnode.firstChild.data.split()
            coords.append((float(coordstr[0]),float(coordstr[1]),float(coordstr[2])))
        if len(coords)!=len(times):
            raise ValueError("Coords and times don't match up: %d times, %d coords"%(len(times),len(coords)))
        track={}
        for time,coord in zip(times,coords):
            track[time]=coord
        tracks.append((name,firstTime,lastTime,altMode,track))
    return tracks

def getgps(exif):
    """
    Get GPS information from EXIF data
    :param tags:
    :return:
    """
    latdeg = exif["GPS GPSLatitude"].values[0].num / exif["GPS GPSLatitude"].values[0].den
    latmin = exif["GPS GPSLatitude"].values[1].num / exif["GPS GPSLatitude"].values[1].den
    latsec = exif["GPS GPSLatitude"].values[2].num / exif["GPS GPSLatitude"].values[2].den
    latNS = exif["GPS GPSLatitudeRef"].values
    lat = (latdeg + latmin / 60 + latsec / 3600) * (1 if latNS == "N" else -1)
    londeg = exif["GPS GPSLongitude"].values[0].num / exif["GPS GPSLongitude"].values[0].den
    lonmin = exif["GPS GPSLongitude"].values[1].num / exif["GPS GPSLongitude"].values[1].den
    lonsec = exif["GPS GPSLongitude"].values[2].num / exif["GPS GPSLongitude"].values[2].den
    lonEW = exif["GPS GPSLongitudeRef"].values
    lon = (londeg + lonmin / 60 + lonsec / 3600) * (1 if lonEW == "E" else -1)
    GPSDate = exif["GPS GPSDate"].values
    GPSDate = GPSDate.split(":")
    GPSYear = int(GPSDate[0])
    GPSMonth = int(GPSDate[1])
    GPSDay = int(GPSDate[2])
    GPSHour = int(exif["GPS GPSTimeStamp"].values[0].num / exif["GPS GPSTimeStamp"].values[0].den)
    GPSMin = int(exif["GPS GPSTimeStamp"].values[1].num / exif["GPS GPSTimeStamp"].values[1].den)
    GPSSec = int(exif["GPS GPSTimeStamp"].values[2].num / exif["GPS GPSTimeStamp"].values[2].den)
    GPSssec=exif["GPS GPSTimeStamp"].values[2].num / exif["GPS GPSTimeStamp"].values[2].den-GPSSec
    GPSValid=True
    if "GPS GPSStatus" in exif and exif["GPS GPSStatus"].values=="V":
        GPSValid=False
    alt=1
    altMode="relativeToGround"
    return (lat,lon,alt,(GPSYear,GPSMonth,GPSDay,GPSHour,GPSMin,GPSSec),GPSValid,altMode)

def getimgtime(exif):
    ImageTimestamp = exif["Image DateTime"].values.split(" ")
    ImageDate = ImageTimestamp[0].split(":")
    ImageYear = int(ImageDate[0])
    ImageMonth = int(ImageDate[1])
    ImageDay = int(ImageDate[2])
    ImageTime = ImageTimestamp[1].split(":")
    ImageHour = int(ImageTime[0])
    ImageMin = int(ImageTime[1])
    ImageSec = int(ImageTime[2])
    return (ImageYear,ImageMonth,ImageDay,ImageHour,ImageMin,ImageSec)

def make_aware(naive,tzname='UTC'):
    naive_datetime=datetime(naive[0],naive[1],naive[2],naive[3],naive[4],naive[5])
    tz=pytz.timezone(tzname)
    return tz.localize(naive_datetime)

def searchKML(tracks,utc):
    bestdt=None
    besttime=None
    bestcoord=None
    bestaltMode=None
    besttrackname=None
    for track in tracks:
        #0 name - name of track
        #1 firstTime - datetime object representing earliest time in the track
        #2 lastTime - datetime object representing latest time in the track
        #3 altMode - string representing altitude mode of height coordinate
        #4 track - a dictionary of tuples, indexed by time. Each tuple is:
            #0 East longitude in degrees
            #1 North latitude in degrees
            #2 Altitude above reference in meters
        if track[1]<=utc and track[2]>=utc:
            if utc in track[4]:
                #We have an exact match
                return (track[0],utc,track[4][utc],track[3])
            else:
                t0=None
                coord0=None
                for t1,coord1 in track[4].items():
                    if t0 is None:
                        t0=t1
                        coord0=coord1
                        continue
                    if t1>utc:
                        #return first time after given utc
                        dt0=(t0 - utc).total_seconds()
                        t0=t1
                        dt1=(t1-utc).total_seconds()
                        if abs(dt0)<abs(dt1):
                            t=t0
                            dt=dt0
                            coord=coord0
                        else:
                            t=t1
                            dt=dt1
                            coord=coord1
                        if bestdt is None or abs(dt)<abs(bestdt):
                            bestdt=dt
                            besttime=t
                            bestcoord=coord
                            bestaltMode=track[3]
                            besttrackname=track[0]
                        break #Stop searching this track
    if bestcoord is not None:
        return(besttrackname,besttime,bestcoord,bestaltMode)
    else:
        return None

def guess_from_filename(infn):
    """
    Guess the image timestamp from the filename by matching it against certain regular expressions
    :param infn:
    :return:
    """
    match=re.match("20([0-9]{2})([0-9]{2})([0-9]{2})_([0-9]{2})([0-9]{2})([0-9]{2}).jpg",infn)
    if match is None:
        return None

def printTimeDiff(imgtime,tracktime):
    """

    :param imgtime: Image timestamp
    :param tracktime: GPS best estimate timestamp
    :return: String describing how much behind the GPS is from the image. Negative indicates GPS data is from after the image
    """
    diff = (imgtime - tracktime).total_seconds()
    return "%s%02d:%02d:%02d (%ds)"%("+" if diff>0 else "-", abs(diff)//3600, (abs(diff)//60)%60, abs(diff)%60, diff)

def find_pic(tracks,infn,tz_guess="America/New_York"):
    """
    Given tracks and an image file name, try to place the image in space and time
    :param tracks: Set of tracks, as returned by parseKML
    :param infn: name of image file to place
    :return:
    *If the image cannot be placed, None
    *If the image can be placed, a tuple:
       # datetime object representing time of image.
           Timezone will be set to UTC, appropriate local time adjustment already made
       # tuple of width and height of image in pixels
       # string representing altMode of altitude coordinate
       # tuple of east longitude in degrees, north latitude in degrees, and altitude above reference surface
       # Warning string describing any problems found in the placing routine

    This routine parses the EXIF image data for image time. It also looks for GPS data in the image, but
    only for limited purposes:
    * If the image has no GPS data, the timestamp in the EXIF (if present) or filename is used. The time zone
      guess is used to correct to UTC
    * If the image has GPS data, the position is used to search the time zone geometry database and adjust the
      image time
    Once this is done, it searches for the GPS position closest in time to the image time with guesstimated time zone.
    It then uses that position to get the time zone again, and complains if the time zone changed.
    """
    global tf
    with open(infn, "rb") as inf:
        exif = exifread.process_file(inf)
    # for k,v in tags.items():
    #    if k!="JPEGThumbnail":
    #        print(k,v)
    if "Image DateTime" in exif:
        imgtime = getimgtime(exif)
    else:
        imgtime=guess_from_filename(infn)
        if imgtime is None:
            return None
    gpswarning=""
    if "GPS GPSLatitude" in exif:
        print(infn)
        (lat, lon, alt, gpstime, gpsvalid, altMode) = getgps(exif)
        tzname = tf.timezone_at(lat=lat, lng=lon)
        g = make_aware(gpstime)
        l = make_aware(imgtime, tzname=tzname)

        gpswarning=""
        if not gpsvalid:
            gpswarning = "Warning! GPS data is marked invalid.<br />\n"
        gpswarning+="GPS and image time difference: %s<br />\n" % printTimeDiff(l,g)
    else:
        gpstime=None
        tzname=tz_guess
        l = make_aware(imgtime, tzname=tzname)
    l_utc = l.astimezone(pytz.utc)
    if "EXIF ExifImageWidth" in exif:
        width = exif["EXIF ExifImageWidth"].values[0]
        height = exif["EXIF ExifImageLength"].values[0]
    else:
        return None
    warning="Image time: %s<br />\n"%l.strftime("%Y-%m-%dT%H:%M:%S%z")
    if gpstime is not None:
        warning+="GPS time: %04d-%02d-%02dT%02d:%02d:%02dZ<br />\n" % gpstime + gpswarning
    #Now that we have our best guess for image time in UTC, check the tracks
    sk=searchKML(tracks,l_utc)
    if sk is None:
        #No position found in the tracks, fall back to GPS position
        if gpstime is None:
            #If no GPS position, we can't find it
            return None
        warning+="No position found, using GPS data recorded in image<br />\n"
        trackName=infn
        trackTime=g
    else:
        (trackName,trackTime,(lon,lat,alt),altMode)=sk
    warning+="Track: %s<br />\n"%trackName
    warning+="Track time: %s (%s difference)<br />\n"%(trackTime.strftime("%Y-%m-%dT%H:%M:%S%z"),printTimeDiff(l_utc,trackTime))

    #Check if the time zone changed. If so, image is ambiguous (need human help)
    newtzname = tf.timezone_at(lat=lat, lng=lon)
    if newtzname!=tzname:
        warning+="Time zone uncertain, found %s and %s\n"%(tzname,newtzname)
    print(warning)
    return (l_utc,(width,height),altMode,(lon,lat,alt),warning)

def kml_header(oufn):
    ouf=open(oufn, "wt")
    ouf.write('<?xml version="1.0" encoding="UTF-8"?>\n')
    ouf.write('<kml xmlns="http://www.opengis.net/kml/2.2" xmlns:gx="http://www.google.com/kml/ext/2.2" xmlns:kml="http://www.opengis.net/kml/2.2" xmlns:atom="http://www.w3.org/2005/Atom">\n')
    ouf.write('  <Document>\n')
    ouf.write('    <name>%s</name>\n' % oufn)
    ouf.write('    <Style id="camera">\n')
    ouf.write('      <IconStyle>\n')
    ouf.write('        <Icon>\n')
    ouf.write('          <href>:/camera_mode.png</href>\n')
    ouf.write('        </Icon>\n')
    ouf.write('      </IconStyle>\n')
    ouf.write('    </Style>\n')
    return ouf

def kml_picture(tracks,infn,tz_guess=None):
    """
    Write one placemark for one image into the target KML
    :param ouf: open file to write to
    :param infn: Filename of image to write
    :param tz_guess: Name of the guesstimated time zone of the image. For instance, if you
      know that the image was taken in California, pass "US/Pacific".

    """
    result="    <!-- %s -->\n" % infn
    fp=find_pic(tracks,infn)
    if fp is not None:
        (imgtime, (width, height), altMode, (lon, lat, alt), warning)=fp

        if width > height:
            height = height * 640 // width
            width = 640
        else:
            width = width * 640 // height
            height = 640
        result+="    <Placemark>\n"
        result +="      <name>%s %s</name>\n" % (infn.split("/")[-1],imgtime.strftime("%Y-%m-%dT%H:%M:%SZ"))
        result +="      <description><![CDATA[<img src=\"%s\" width=\"%d\" height=\"%d\" /><br />%s]]></description>\n" % (
                infn, width, height, warning)
        result +="      <styleUrl>camera</styleUrl>\n"
        result +="      <TimeStamp><when>%s</when></TimeStamp>\n" % imgtime.strftime("%Y-%m-%dT%H:%M:%SZ")
        result +="      <Point>\n"
        if altMode is not None:
            result +="        <altitudeMode>%s</altitudeMode>\n"%altMode
        result +="        <coordinates>%f,%f,%f</coordinates>\n" % (lon, lat, alt)
        result +="      </Point>\n"
        result +="    </Placemark>\n"
        return (imgtime,result)
    return None

def kml_footer(ouf):
    ouf.write('  </Document>\n')
    ouf.write('</kml>\n')

print("Parsing KML")
tracks=parseKML("/home/jeppesen/workspace/Data/kml/Florida5.kml")
print("Tracks found: %d"%len(tracks))

path = "/home/jeppesen/Florida5/"
print("Attempting to find place for each image in "+path)
with kml_header(path + "Geotag.kml") as ouf:
    kml={}
    for infn in glob.iglob(path + "/**/*.[Jj][Pp]*", recursive=True):
        print(infn)
        kp=kml_picture(tracks,infn,tz_guess="America/New_York")
        if kp is not None:
            (time,this_kml)=kp
            if time in kml:
                kml[time][os.path.basename(infn)]=this_kml
            else:
                kml[time]={os.path.basename(infn):this_kml}
    for k in sorted(kml):
        for base_infn in kml[k]:
            this_kml=kml[k][base_infn]
            ouf.write(this_kml)
    kml_footer(ouf)
