import xml.etree.ElementTree as ET

tree=ET.parse('xtce_pocketometer20150317.xml')
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
            typestr=">f"
        else:
            typestr=">d"
    elif ide is not None:
        #It's an int scalar
        bits=int(ide.attrib["sizeInBits"])
        signed=(ide.attrib["encoding"]!="unsigned")
        if bits<=8:
            typestr="B"
        elif bits<=16:
            typestr=">I"
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
    types[ptype.attrib["name"]]={'bits':bits,'typestr':typestr,'ptype':ptype}    
#* A fields dictionary made from ParameterSet which is indexed by the name 
#  attribute of the ParameterSet children
fields={}
for p in list(ps):
    fields[p.attrib["name"]]={'parameter':p,'type':types[p.attrib["parameterTypeRef"]]}
#* A header dictionary made from ContainerSet children which is indexed by the
#  name attribute and only includes abstract containers
headers={}
for h in list(cs):
    if "abstract" in h.attrib:
        headerfields=[]
        
        fields[h.attrib["name"]]={header:}
#* A packet dictionary made from ContainerSet children which is indexed by the
#  name attribute and doesn't include abstract containers
print(fields["CONFIG__DIRENTRY"])
for name,field in fields.items():
    print(name,field)
pass