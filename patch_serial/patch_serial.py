import itertools

data={}
time={}
count={}
line={}

with open("scratch.txt","r") as inf:
    linenum=1
    for line1, line2 in itertools.zip_longest(*[inf] * 2):
        if line1 is not None:
            line1=line1.strip(" \t\n\r")
        else:
            print("None line found at line %d"%linenum)
        if line2 is not None:
            line2=line2.strip(" \t\n\r")
        else:
            print("None line found at line %d" % linenum)
        if line1=="":
            print("Blank line found on line %d"%linenum)
        linenum+=1
        if line2=="":
            print("Blank line found on line %d"%linenum)
        linenum+=1
        parts=line1.split()
        this_pkt=int(parts[0])
        this_time=parts[1]
        pkt_key=this_time
        if this_pkt in count:
            count[pkt_key]+=1
            if data[pkt_key]!=line2:
                print("Different data with the same packet number %s"%pkt_key)
                print("Line %5d Was: %s"%(line[pkt_key],data[pkt_key]))
                print("Line %5d Is:  %s"%(linenum,line2))

            if time[pkt_key]!=this_time:
                print("Different time with the same line number")
                print("Time %5d Was: %s"%(line[pkt_key],time[pkt_key]))
                print("Time %5d Is:  %s"%(linenum,this_time))
        else:
            count[pkt_key]=1
        line[pkt_key]=linenum
        time[pkt_key]=this_time
        data[pkt_key]=line2

#Second pass - sort by time in the GPRMC or GPGGA lines. Drop the other lines
sorted_lines={}


with open("scratch.nmea","w") as ouf:
    for k in sorted(data):
        ouf.write(data[k]+"\n")
