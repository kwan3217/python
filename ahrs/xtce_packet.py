import xml.etree.ElementTree as ET
import collections
import struct
import glob

def readPktDesc(xml):
    tree=ET.parse(xml)
    namespaces={'xtce':'http://www.omg.org/space/xtce'}
    #Dig out the TelemetryMetaData tag
    tlm=tree.getroot().findall('xtce:TelemetryMetaData',namespaces)[0]
    #From this we need:
    #* ParameterTypeSet to get names of types, which include bit sizes
    pts=tlm.find('xtce:ParameterTypeSet',namespaces)
    #* ParameterSet, which defines fields and references types defined in 
    #  ParameterTypeSet
    ps=tlm.find('xtce:ParameterSet',namespaces)
    #* ContainerSet, which defines headers and packets. Headers have their abstract
    #  attribute set to "true"
    cs=tlm.find('xtce:ContainerSet',namespaces)
    #So: Make three dictionaries:
    #* A types dictionary made from ParameterTypeSet which is indexed by the name
    #  attribute of the ParameterTypeSet children
    types={}
    for ptype in list(pts):
        #Dig out the FloatDataEncoding or IntegerDataEncoding
        fde=ptype.find("xtce:FloatDataEncoding",namespaces)
        ide=ptype.find("xtce:IntegerDataEncoding",namespaces)
        sde=ptype.find("xtce:StringDataEncoding",namespaces)
        bde=ptype.find("xtce:BinaryDataEncoding",namespaces)
        if fde is not None:
            #it's a float scalar  
            bits=int(fde.attrib["sizeInBits"])
            if bits==32:
                typestr="<f"
            else:
                typestr="<d"
        elif ide is not None:
            #It's an int scalar
            bits=int(ide.attrib["sizeInBits"])
            signed=(ide.attrib["encoding"]!="unsigned")
            if bits<=8:
                typestr="B"
            elif bits<=16:
                typestr=">H"
            else:
                typestr=">L"
            #We can do this because the unsigned version of all of the above are
            #the same as the upper-case version of the signed
            if signed:
                typestr=typestr.lower()
        elif sde is not None:
            typestr="s"
            bits=int(sde.find("xtce:SizeInBits",namespaces).find("xtce:Fixed",namespaces).find("xtce:FixedValue",namespaces).text)
        elif bde is not None:
            bits=int(bde.find("xtce:SizeInBits",namespaces).find("xtce:FixedValue",namespaces).text)
            typestr="c"
        else:
            raise ValueError("Not an int or a float: "+ptype.tag)
        types[ptype.attrib["name"]]={'bits':bits,'typestr':typestr}    
    #* A fields dictionary made from ParameterSet which is indexed by the name 
    #  attribute of the ParameterSet children
    fields={}
    for p in list(ps):
        fields[p.attrib["name"]]=types[p.attrib["parameterTypeRef"]]
    #* A header dictionary made from ContainerSet children which is indexed by the
    #  name attribute. This includes both headers and packets. We will also make 
    #  a dictionary indexed by apid.
    headers={}
    apids={}
    for h in list(cs):
        headerfields=collections.OrderedDict()
        b=h.find("xtce:BaseContainer",namespaces)
        if b is not None:
            headerfields.update(headers[b.attrib["containerRef"]])
        e=h.find("xtce:EntryList",namespaces)
        for f in list(e):
            if f.tag.endswith("ParameterRefEntry"):
                name=f.attrib["parameterRef"]
                name=name.split("__")[-1]
                headerfields[name]={}
                headerfields[name].update(fields[f.attrib["parameterRef"]])
                r=f.find("xtce:RepeatEntry",namespaces)
                if r is not None:
                    headerfields[name]["repeat"]=int(r.find("xtce:Count",namespaces).find("xtce:FixedValue",namespaces).text)
                else:
                    headerfields[name]["repeat"]=1
            else:
                raise ValueError("Unhandled tag "+f.tag)
        headers[h.attrib["name"]]=headerfields
        if "abstract" not in h.attrib:
            apid=int(b.find("xtce:RestrictionCriteria",namespaces).find("xtce:ComparisonList",namespaces).find("xtce:Comparison",namespaces).attrib["value"])
            apids[apid]=headerfields
        elif len(apids)==0:
            apids[0]=headerfields
    packets={}
    for apid,packet in apids.items():
        bitofs=0
        fieldnames=[]
        fieldtype=[] #Will accumulate field type strings
        bitstart=[] #will accumulate field start points in bits
        bitlen=[] #will accumualte field lengths in bits
        repeat=[] #will accululate repeat counts
        for name,field in packet.items():
            fieldnames.append(name)
            fieldtype.append(field["typestr"])
            bitstart.append(bitofs)
            bitofs+=field["bits"]
            bitlen.append(field["bits"])
            repeat.append(field["repeat"])
    #    for i in range(len(fieldnames)):
    #        print(fieldnames[i],bitstart[i],bitlen[i],fieldtype[i],repeat[i])
        packets[apid]=(fieldnames,bitstart,bitlen,fieldtype,repeat)    
    return packets

def getPacket(buf,pktDef):
    packet=collections.OrderedDict()
    for i in range(len(pktDef[0])):
        fieldname=pktDef[0][i]
        bitstart=pktDef[1][i]
        bitlen=pktDef[2][i]
        fieldtype=pktDef[3][i]
        repeat=pktDef[4][i]
        if fieldtype in ['s','c']:
            if repeat<0:
                val=buf[bitstart//8:]
            else:
                val=buf[bitstart//8:bitstart//8+bitlen//8]
            if fieldtype=='s':
                val=val.decode("utf-8","ignore")
        else:
            val=struct.unpack_from(fieldtype,buf,bitstart//8)[0]
            #In CCSDS, the MSB is numbered 0, since it is transmitted first.
            #Boo to that, but what can we do? So, a general bit field will have
            #some padding before, some number of data bits, and some pad after.
            #The value needs to be shifted right the number of bits of the
            #padding after, then masked to the number of bits in the value.
            bufbitlen=struct.calcsize(fieldtype)*8
            bitsbefore=bitstart%8
            bitsafter=bufbitlen-bitsbefore-bitlen
            if bitsafter>0:
                val=val>>bitsafter
            if bufbitlen>bitlen:
                val=val & ((1<<bitlen)-1)
        packet[fieldname]=val 
    return packet

def getNextPacket(inf,pktDefs,apid_filter=None):
    CCSDSheader=pktDefs[0]
    nbitheader=0
    for bitlen in CCSDSheader[2]:
        nbitheader+=bitlen
    nbytes=nbitheader//8
    if (nbitheader%8)!=0:
        nbytes+=1
    while True:
        b=inf.read(nbytes)
        header=getPacket(b,CCSDSheader)
        pktlen=header["PKT_LEN"]+1
        b=b+inf.read(pktlen)
        if apid_filter is None or header["PKT_APID"] in apid_filter:
            return getPacket(b,pktDefs[header["PKT_APID"]])

