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
    t_i64=0x0E
    t_u64=0x0F
                #0  1    2    3    4    5   6  7  8   9  A   B   C    D    E    F
                #   u8  i16  i32  f32  f64    str       bin    u16   u32  i64  u64
    lookupType=['','B','>h','>i','>f','>d','','s','','','c','','>H','>I','>q','>Q']
    lookupSize=[ 0, 1,   2,   4,   4,   8,  0, 0,  0,  0, 0,  0,  2,   4,  8 ,  8 ]
    infn="/home/jeppesen/workspace/python/parseYukari/20150312_POCKET03.SDS"
    base=os.path.dirname(infn)+"/"+os.path.basename(infn)[:-4]
    try:
        os.mkdir(base)
    except:
        pass
    with open(infn,"rb") as inf:
        inf.seek(8)
        done=False
        pktdesc=OrderedDict() #Dictionary of packet descriptions keyed on apid (numeric)
        #Preload 11 (bootstrap) and 12 (not complete)
        apidfirst=OrderedDict()
        #These routines are inner because they have access to variables defined above
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
        yukaripktdesc=OrderedDict()
        yukaripktdesc[11]={'name':  'PacketDescription',
                     'handler':addField,
                     'fields':OrderedDict({'descApid':{'pos': 6,'type':t_u16},
                               'descPos': {'pos': 8,'type':t_u16},
                               'descType':{'pos':10,'type':t_u8},
                               'descName':{'pos':11,'type':t_string}})}
        yukaripktdesc[12]={'name':'CCSDS self-documentation','handler':print,'fields':{'Documentation':{'pos': 6,'type':t_string}}}
        southwestpktdesc=OrderedDict()
        southwestpktdesc[0x00C]={'name':'Version',
                                 'handler':csv,
                                 'fields':OrderedDict({'hw_type'  :{'pos': 6,'type':t_i32},
                                                       'hw_serial':{'pos':10,'type':t_i32},
                                                       'mamcr'    :{'pos':14,'type':t_u8},
                                                       'mamtim'   :{'pos':15,'type':t_u8},
                                                       'pll0stat' :{'pos':16,'type':t_u16},
                                                       'vpbdiv'   :{'pos':18,'type':t_u8},
                                                       'fosc'     :{'pos':19,'type':t_u32},
                                                       'cclk'     :{'pos':23,'type':t_u32},
                                                       'pclk'     :{'pos':27,'type':t_u32},
                                                       'preint'   :{'pos':31,'type':t_u32},
                                                       'prefrac'  :{'pos':35,'type':t_u32},
                                                       'ccr'      :{'pos':39,'type':t_u8},
                                                       'version'  :{'pos':40,'type':t_string} })}
        southwestpktdesc[0x022]={'name':'ConfigFile',
                                 'handler':csv,
                                 'fields':OrderedDict({'base'     :{'pos':6+0x00,'type':t_string,'len':8},
                                                       'ext'      :{'pos':6+0x08,'type':t_string,'len':3},
                                                       'attr'     :{'pos':6+0x0B,'type':t_u8},
                                                       'resv1'    :{'pos':6+0x0C,'type':t_u8},
                                                       'ctenths'  :{'pos':6+0x0D,'type':t_u8},
                                                       'ctime'    :{'pos':6+0x0E,'type':t_u16,'end':'<'},
                                                       'cdate'    :{'pos':6+0x10,'type':t_u16,'end':'<'},
                                                       'adate'    :{'pos':6+0x12,'type':t_u16,'end':'<'},
                                                       'clusterM' :{'pos':6+0x14,'type':t_u16,'end':'<'},
                                                       'wtime'    :{'pos':6+0x16,'type':t_u16,'end':'<'},
                                                       'wdate'    :{'pos':6+0x18,'type':t_u16,'end':'<'},
                                                       'clusterL' :{'pos':6+0x1A,'type':t_u16,'end':'<'},
                                                       'size'     :{'pos':6+0x1C,'type':t_u32,'end':'<'},
                                                       'configContents':{'pos':6+0x20,'type':t_string}})}
        southwestpktdesc[0x028]={'name':'ConfigParse',
                                 'handler':csv,
                                 'fields':OrderedDict({'gyroSens'  :{'pos': 6,'type':t_i32},
                                                       'gyroODR'   :{'pos':10,'type':t_i32},
                                                       'gyroBW'    :{'pos':14,'type':t_i32},
                                                       'P'         :{'pos':18,'type':t_i32},
                                                       'Ps'        :{'pos':22,'type':t_i32},
                                                       'I'         :{'pos':26,'type':t_i32},
                                                       'Is'        :{'pos':30,'type':t_i32},
                                                       'D'         :{'pos':34,'type':t_i32},
                                                       'Ds'        :{'pos':38,'type':t_i32},
                                                       'nWaypoints':{'pos':42,'type':t_i32}})}
        southwestpktdesc[0x003]={'name':'source',
                                 'handler':dump,
                                 'fields':OrderedDict({'ofs'       :{'pos': 6,'type':t_u16},
                                                       'DumpData'  :{'pos': 8,'type':t_bin}})}
        southwestpktdesc[0x012]={'name':'sdinfo',
                                 'handler':csv,
                                 'fields':OrderedDict({'manufacturer':{'pos': 6,'type':t_u8},
                                                       'oem'         :{'pos': 7,'type':t_string,'len':3+1},
                                                       'product'     :{'pos':11,'type':t_string,'len':6+1},
                                                       'revision'    :{'pos':18,'type':t_u8},
                                                       'serial'      :{'pos':19,'type':t_u32,'end':'<'},
                                                       'mfg_year'    :{'pos':23,'type':t_u8},
                                                       'mfg_month'   :{'pos':24,'type':t_u8},
                                                       'capacity'    :{'pos':25,'type':t_u64,'end':'<'},
                                                       'flag_copy'   :{'pos':33,'type':t_u8},
                                                       'flag_wp'     :{'pos':34,'type':t_u8},
                                                       'flag_wp_temp':{'pos':34,'type':t_u8},
                                                       'format'      :{'pos':34,'type':t_u8}})}
        southwestpktdesc[0x00E]={'name':'HMC5883config',
                                 'handler':csv,
                                 'fields':OrderedDict({'reg0':{'pos': 6,'type':t_u8},
                                                       'reg1':{'pos': 7,'type':t_u8},
                                                       'reg2':{'pos': 8,'type':t_u8},
                                                       'reg9':{'pos': 9,'type':t_u8},
                                                       'regA':{'pos':10,'type':t_u8},
                                                       'regB':{'pos':11,'type':t_u8},
                                                       'regC':{'pos':12,'type':t_u8}})}
        southwestpktdesc[0x002]={'name':'BMP180config',
                                 'handler':csv,
                                 'fields':OrderedDict({'ac1':{'pos': 6,'type':t_i16},
                                                       'ac2':{'pos': 8,'type':t_i16},
                                                       'ac3':{'pos':10,'type':t_i16},
                                                       'ac4':{'pos':12,'type':t_u16},
                                                       'ac5':{'pos':14,'type':t_u16},
                                                       'ac6':{'pos':16,'type':t_u16},
                                                       'b1' :{'pos':18,'type':t_i16},
                                                       'b2' :{'pos':20,'type':t_i16},
                                                       'mb' :{'pos':22,'type':t_i16},
                                                       'mc' :{'pos':24,'type':t_i16},
                                                       'md' :{'pos':26,'type':t_i16}})}
        southwestpktdesc[0x020]={'name':'L3G4200Dconfig',
                                 'handler':csv,
                                 'fields':OrderedDict({'reg20':{'pos': 6,'type':t_u8},
                                                       'reg21':{'pos': 7,'type':t_u8},
                                                       'reg22':{'pos': 8,'type':t_u8},
                                                       'reg23':{'pos': 9,'type':t_u8},
                                                       'reg24':{'pos':10,'type':t_u8}})}
        southwestpktdesc[0x023]={'name':'ADXL345config',
                                 'handler':csv,
                                 'fields':OrderedDict({'reg31':{'pos': 6,'type':t_u8}})}

        southwestpktdesc[0x024]={'name':'AccGyro',
                                 'handler':csv,
                                 'fields':OrderedDict({'TC0'   :{'pos': 6,'type':t_u32},
                                                       'gx'    :{'pos':10,'type':t_i16},
                                                       'gy'    :{'pos':12,'type':t_i16},
                                                       'gz'    :{'pos':14,'type':t_i16},
                                                       'ax'    :{'pos':16,'type':t_i16},
                                                       'ay'    :{'pos':18,'type':t_i16},
                                                       'az'    :{'pos':20,'type':t_i16},
                                                       't'     :{'pos':22,'type':t_u8},
                                                       'status':{'pos':23,'type':t_u8},
                                                       'TC1'   :{'pos':24,'type':t_u32}})}
        southwestpktdesc[0x004]={'name':'Mag',
                                 'handler':csv,
                                 'fields':OrderedDict({'TC0'   :{'pos': 6,'type':t_u32},
                                                       'bx'    :{'pos':10,'type':t_i16},
                                                       'bz'    :{'pos':12,'type':t_i16},
                                                       'by'    :{'pos':14,'type':t_i16}})}
        southwestpktdesc[0x00A]={'name':'BMP',
                                 'handler':csv,
                                 'fields':OrderedDict({'TC0'   :{'pos': 6,'type':t_u32},
                                                       'traw'  :{'pos':10,'type':t_i16},
                                                       'praw'  :{'pos':12,'type':t_u32},
                                                       'tcal'  :{'pos':16,'type':t_i16},
                                                       'pcal'  :{'pos':18,'type':t_u32},
                                                       'TC1'   :{'pos':22,'type':t_u32}})}
        southwestpktdesc[0x01A]={'name':'gps',
                                 'handler':dump,
                                 'fields':OrderedDict({'TC0'   :{'pos': 6,'type':t_u32},
                                                       'DumpData':{'pos':10,'type':t_bin}})}
        southwestpktdesc[0x039]={'name':'rmc',
                                 'handler':csv,
                                 'fields':OrderedDict({'TC'   :{'pos': 6,'type':t_u32},
                                                       'HMS'  :{'pos':10,'type':t_float,'end':'<'},
                                                       'lat'  :{'pos':14,'type':t_float,'end':'<'},
                                                       'lon'  :{'pos':18,'type':t_float,'end':'<'},
                                                       'spd'  :{'pos':22,'type':t_float,'end':'<'},
                                                       'hdg'  :{'pos':26,'type':t_float,'end':'<'},
                                                       'DMY'  :{'pos':30,'type':t_u32}})}
        southwestpktdesc[0x038]={'name':'gga',
                                 'handler':csv,
                                 'fields':OrderedDict({'TC'   :{'pos': 6,'type':t_u32},
                                                       'HMS'  :{'pos':10,'type':t_float,'end':'<'},
                                                       'lat'  :{'pos':14,'type':t_float,'end':'<'},
                                                       'lon'  :{'pos':18,'type':t_float,'end':'<'},
                                                       'alt'  :{'pos':22,'type':t_float,'end':'<'}})}
        pocketpktdesc=OrderedDict()
        pocketpktdesc[0x001]={'name':'Version',
                                 'handler':csv,
                                 'fields':OrderedDict({'hw_type'  :{'pos': 6,'type':t_i32},
                                                       'hw_serial':{'pos':10,'type':t_i32},
                                                       'mamcr'    :{'pos':14,'type':t_u8},
                                                       'mamtim'   :{'pos':15,'type':t_u8},
                                                       'pll0stat' :{'pos':16,'type':t_u16},
                                                       'vpbdiv'   :{'pos':18,'type':t_u8},
                                                       'fosc'     :{'pos':19,'type':t_u32},
                                                       'cclk'     :{'pos':23,'type':t_u32},
                                                       'pclk'     :{'pos':27,'type':t_u32},
                                                       'preint'   :{'pos':31,'type':t_u32},
                                                       'prefrac'  :{'pos':35,'type':t_u32},
                                                       'ccr'      :{'pos':39,'type':t_u8},
                                                       'version'  :{'pos':40,'type':t_string} })}
        pocketpktdesc[0x002]={'name':'source',
                                 'handler':dump,
                                 'fields':OrderedDict({'start'     :{'pos': 6,'type':t_u32},
                                                       'end'       :{'pos':10,'type':t_u32},
                                                       'base'      :{'pos':14,'type':t_u32},
                                                       'DumpData'  :{'pos':18,'type':t_bin}})}
        pocketpktdesc[0x003]={'name':'image',
                                 'handler':dump,
                                 'fields':OrderedDict({'start'     :{'pos': 6,'type':t_u32},
                                                       'end'       :{'pos':10,'type':t_u32},
                                                       'base'      :{'pos':14,'type':t_u32},
                                                       'DumpData'  :{'pos':18,'type':t_bin}})}
        pocketpktdesc[0x004]={'name':'ConfigFile',
                                 'handler':csv,
                                 'fields':OrderedDict({'base'     :{'pos':6+0x00,'type':t_string,'len':8},
                                                       'ext'      :{'pos':6+0x08,'type':t_string,'len':3},
                                                       'attr'     :{'pos':6+0x0B,'type':t_u8},
                                                       'resv1'    :{'pos':6+0x0C,'type':t_u8},
                                                       'ctenths'  :{'pos':6+0x0D,'type':t_u8},
                                                       'ctime'    :{'pos':6+0x0E,'type':t_u16,'end':'<'},
                                                       'cdate'    :{'pos':6+0x10,'type':t_u16,'end':'<'},
                                                       'adate'    :{'pos':6+0x12,'type':t_u16,'end':'<'},
                                                       'clusterM' :{'pos':6+0x14,'type':t_u16,'end':'<'},
                                                       'wtime'    :{'pos':6+0x16,'type':t_u16,'end':'<'},
                                                       'wdate'    :{'pos':6+0x18,'type':t_u16,'end':'<'},
                                                       'clusterL' :{'pos':6+0x1A,'type':t_u16,'end':'<'},
                                                       'size'     :{'pos':6+0x1C,'type':t_u32,'end':'<'},
                                                       'configContents':{'pos':6+0x20,'type':t_string}})}
        pocketpktdesc[0x005]={'name':'MPU6050config',
                                 'handler':csv,
                                 'fields':OrderedDict({'self_test_x':{'pos': 6,'type':t_u8},
                                                       'self_test_y':{'pos': 7,'type':t_u8},
                                                       'self_test_z':{'pos': 8,'type':t_u8},
                                                       'self_test_a':{'pos': 9,'type':t_u8},
                                                       'smplrt_div' :{'pos':10,'type':t_u8},
                                                       'config'     :{'pos':11,'type':t_u8},
                                                       'gyro_config':{'pos':12,'type':t_u8},
                                                       'acc_config' :{'pos':13,'type':t_u8},
                                                       'mot_thr'    :{'pos':14,'type':t_u8},
                                                       'int_pin_cfg':{'pos':15,'type':t_u8},
                                                       'int_enable' :{'pos':16,'type':t_u8},
                                                       'user_ctrl'  :{'pos':17,'type':t_u8},
                                                       'pwr_mgmt_1' :{'pos':18,'type':t_u8},
                                                       'whoami'     :{'pos':19,'type':t_u8}})}
        pocketpktdesc[0x006]={'name':'AK8975config',
                                 'handler':csv,
                                 'fields':OrderedDict({'whoami':{'pos': 6,'type':t_u8},
                                                       'reg01':{'pos': 7,'type':t_u8},
                                                       'reg10':{'pos': 8,'type':t_u8},
                                                       'reg11':{'pos': 9,'type':t_u8},
                                                       'reg12' :{'pos':10,'type':t_u8}})}
        pocketpktdesc[0x007]={'name':'BMP180config',
                                 'handler':csv,
                                 'fields':OrderedDict({'ac1':{'pos': 6,'type':t_i16},
                                                       'ac2':{'pos': 8,'type':t_i16},
                                                       'ac3':{'pos':10,'type':t_i16},
                                                       'ac4':{'pos':12,'type':t_u16},
                                                       'ac5':{'pos':14,'type':t_u16},
                                                       'ac6':{'pos':16,'type':t_u16},
                                                       'b1' :{'pos':18,'type':t_i16},
                                                       'b2' :{'pos':20,'type':t_i16},
                                                       'mb' :{'pos':22,'type':t_i16},
                                                       'mc' :{'pos':24,'type':t_i16},
                                                       'md' :{'pos':26,'type':t_i16}})}
        pocketpktdesc[0x008]={'name':'sdinfo',
                                 'handler':csv,
                                 'fields':OrderedDict({'manufacturer':{'pos': 6,'type':t_u8},
                                                       'oem'         :{'pos': 7,'type':t_string,'len':3+1},
                                                       'product'     :{'pos':11,'type':t_string,'len':6+1},
                                                       'revision'    :{'pos':18,'type':t_u8},
                                                       'serial'      :{'pos':19,'type':t_u32,'end':'<'},
                                                       'mfg_year'    :{'pos':23,'type':t_u8},
                                                       'mfg_month'   :{'pos':24,'type':t_u8},
                                                       'capacity'    :{'pos':25,'type':t_u64,'end':'<'},
                                                       'flag_copy'   :{'pos':33,'type':t_u8},
                                                       'flag_wp'     :{'pos':34,'type':t_u8},
                                                       'flag_wp_temp':{'pos':34,'type':t_u8},
                                                       'format'      :{'pos':34,'type':t_u8}})}
        pocketpktdesc[0x009]={'name':'MPU',
                                 'handler':csv,
                                 'fields': OrderedDict({'TC0': {'pos': 6, 'type': t_u32},
                                                        'ax': {'pos': 10, 'type': t_i16},
                                                        'ay': {'pos': 12, 'type': t_i16},
                                                        'az': {'pos': 14, 'type': t_i16},
                                                        'gx': {'pos': 16, 'type': t_i16},
                                                        'gy': {'pos': 18, 'type': t_i16},
                                                        'gz': {'pos': 20, 'type': t_i16},
                                                        'mt': {'pos': 22, 'type': t_i16},
                                                        'ti': {'pos': 24, 'type': t_u8},
                                                        'tt': {'pos': 25, 'type': t_u8}})}
        pocketpktdesc[0x00A]={'name':'AK',
                                 'handler':csv,
                                 'fields': OrderedDict({'TC0': {'pos': 6, 'type': t_u32},
                                                        'bx_mpu': {'pos': 10, 'type': t_i16},
                                                        'by_mpu': {'pos': 12, 'type': t_i16},
                                                        'bz_mpu': {'pos': 14, 'type': t_i16},
                                                        'st1': {'pos': 16, 'type': t_u8},
                                                        'st2': {'pos': 17, 'type': t_u8}})}
        pocketpktdesc[0x00B]={'name':'BMP',
                                 'handler':csv,
                                 'fields':OrderedDict({'TC0'   :{'pos': 6,'type':t_u32},
                                                       'traw'  :{'pos':10,'type':t_i16},
                                                       'praw'  :{'pos':12,'type':t_u32},
                                                       'tcal'  :{'pos':16,'type':t_float,'end':"<"},
                                                       'pcal'  :{'pos':20,'type':t_u32}})}
        pocketpktdesc[0x00D]={'name':'gps',
                              'handler':dump,
                              'fields':OrderedDict({'TC0'     :{'pos': 6,'type':t_u32},
                                                    'DumpData':{'pos':10,'type':t_bin}})}
        pocketpktdesc[0x00F]={'name':'rmc',
                                 'handler':csv,
                                 'fields':OrderedDict({'TC'   :{'pos': 6,'type':t_u32},
                                                       'HMS'  :{'pos':10,'type':t_float,'end':'<'},
                                                       'lat'  :{'pos':14,'type':t_float,'end':'<'},
                                                       'lon'  :{'pos':18,'type':t_float,'end':'<'},
                                                       'spd'  :{'pos':22,'type':t_float,'end':'<'},
                                                       'hdg'  :{'pos':26,'type':t_float,'end':'<'},
                                                       'DMY'  :{'pos':30,'type':t_u32}})}
        pocketpktdesc[0x00E]={'name':'gga',
                                 'handler':csv,
                                 'fields':OrderedDict({'TC'   :{'pos': 6,'type':t_u32},
                                                       'HMS'  :{'pos':10,'type':t_float,'end':'<'},
                                                       'lat'  :{'pos':14,'type':t_float,'end':'<'},
                                                       'lon'  :{'pos':18,'type':t_float,'end':'<'},
                                                       'alt'  :{'pos':22,'type':t_float,'end':'<'}})}
        pktdesc=pocketpktdesc
        while not done:
            header=inf.read(6)
            apid=struct.unpack(">H",header[0:2])[0] & ((1<<11)-1)
            pktlen=struct.unpack(">H",header[4:6])[0]+1
            body=inf.read(pktlen)
            pkt={'apid':apid}
            fields=pktdesc[apid]['fields']
            for name,field in fields.items():
                typ =lookupType[field['type']]
                pos =           field['pos']
                size=lookupSize[field['type']]
                if typ=='s':
                    if 'len' in field:
                        end=pos-6+field['len']
                    else:
                        end=-1
                    try:
                        pkt[name]=body[pos-6:end].decode("utf-8")
                    except:
                        pkt[name]=body[pos-6:end].decode("cp437")
                elif typ=='c':
                    if 'len' in field:
                        end=pos-6+field['len']
                        pkt[name] = body[pos - 6:end]
                    else:
                        pkt[name] = body[pos - 6:]
                else:
                    if 'end' in field:
                        typ=field['end']+typ[1:]
                    pkt[name]=struct.unpack(typ,body[pos-6:pos-6+size])[0]
            if pktdesc[apid]['handler'] is not None:
                pktdesc[apid]['handler'](pkt)
                