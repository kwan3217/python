import struct
import os.path
from collections import OrderedDict
if __name__ == '__main__':
    t_u8=1
    t_i16=2
    t_i32=3
    t_float=4
    t_double=5
    t_string=7
    t_bin=0x0A
    t_u16=0x0C
    t_u32=0x0D
                #0  1    2    3    4    5   6  7   8  9  A   B   C    D
                #   u8  i16  i32  f32  f64    str       bin    u16   u32
    lookupType=['','B','>h','>i','>f','>d','','s','','','c','','>H','>I']
    lookupSize=[ 0, 1,   2,   4,   4,   8,  0, 0,  0, 0, 0,  0,  2,   4 ]
    infn="/Users/jeppesen/QuaternionCompass/testOutputData/20171004T230307/packets.sds"
    base=os.path.dirname(infn)
    with open(infn,"rb") as inf:
        done=False
        pktdesc=OrderedDict() #Dictionary of packet descriptions keyed on apid (numeric)
        #Preload 11 (bootstrap) and 12 (not complete)
        apidfirst=OrderedDict()
        def csv(pkt):
            oufn=base+'/'+pktdesc[pkt['apid']]['name']+".csv"
            if pkt['apid'] not in apidfirst:
                with open(oufn,"w") as ouf:
                    apidfirst[pkt['apid']]=False
                    for k,v in pkt.items():
                        if k!='apid':
                            print(k+",",file=ouf,end='')
                    print("",file=ouf)
            with open(oufn,"a") as ouf:
                for k,v in pkt.items():
                    if k!='apid':
                        print(str(v)+",",file=ouf,end='')
                print("",file=ouf)                    
        def dump(pkt):
            oufn=base+'/'+pktdesc[pkt['apid']]['name']+".dump"
            with open(oufn,"ab" if pkt['apid'] in apidfirst else "wb") as ouf:
                ouf.write(pkt['DumpData'])
            apidfirst[pkt['apid']]=False
        def addField(pkt):
            """
            Add a field description to a packet. If the packet describes a whole packet, check if the packet description
            exists yet. If it does, just change the name, otherwise create it fresh.
            """
            print(pkt)
            if pkt['descType']==0:
                if pkt['descApid'] in pktdesc:
                    pktdesc[pkt['descApid']]['name']=pkt['descName']
                else:
                    pktdesc[pkt['descApid']]={'name':pkt['descName'],'handler':csv if pkt['descName']!='Dump' else dump,'fields':OrderedDict()}
            else:
                pktdesc[pkt['descApid']]['fields'][pkt['descName']]={'pos':pkt['descPos'],'type':pkt['descType']}
        pktdesc[11]={'name':  'PacketDescription',
                     'handler':addField,
                     'fields':OrderedDict({'descApid':{'pos': 6,'type':t_u16},
                               'descPos': {'pos': 8,'type':t_u16},
                               'descType':{'pos':10,'type':t_u8},
                               'descName':{'pos':11,'type':t_string}})}
        pktdesc[12]={'name':'CCSDS self-documentation','handler':print,'fields':{'Documentation':{'pos': 6,'type':t_string}}}
        while not done:
            header=inf.read(6)
            apid=struct.unpack(">H",header[0:2])[0] & ((1<<11)-1)
            pktlen=struct.unpack(">H",header[4:6])[0]+1
            body=inf.read(pktlen)
            pkt={'apid':apid}
            for k,v in pktdesc[apid]['fields'].items():
                typ =lookupType[pktdesc[apid]['fields'][k]['type']]
                pos =pktdesc[apid]['fields'][k]['pos' ]
                size=lookupSize[pktdesc[apid]['fields'][k]['type']]
                if typ=='s':
                    pkt[k]=body[pos-6:].decode("utf-8")
                elif typ=='c':
                    pkt[k]=body[pos-6:]
                else:
                    pkt[k]=struct.unpack(typ,body[pos-6:pos-6+size])[0]
            if pktdesc[apid]['handler'] is not None:
                pktdesc[apid]['handler'](pkt)
                