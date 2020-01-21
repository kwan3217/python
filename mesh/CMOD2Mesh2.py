"""
Code to import and export Celestia Models (CMOD). Initial target will be
a Wavefront obj, which can be imported by Blender.
"""

from enum import Enum,auto
import numpy as np
import struct

#These were not present in Java, where there was instead a bunch of read[Float|Int](),
#SwapEndian(), and bit masking to handle Java's choice of endianness and lack of
#unsigned values. The following just make things easier to read.
def read_and_dump(Inf,len):
    b=bytearray(len)
    for i in range(len):
#        if(Inf.tell()%32==0):
#            print("\n%08X"%Inf.tell(),end='')
        b[i]=Inf.read(1)[0]
#        print("%02X "%b[i],end='')
    return b
def read_int16(Inf):
    return struct.unpack("<h", read_and_dump(Inf,2))[0]
def read_int32(Inf):
    return struct.unpack("<i", read_and_dump(Inf,4))[0]
def read_uint16(Inf):
    return struct.unpack("<H", read_and_dump(Inf,2))[0]
def read_uint32(Inf):
    return struct.unpack("<I", read_and_dump(Inf,4))[0]
def read_float32(Inf):
    return struct.unpack("<f", read_and_dump(Inf,4))[0]
def read_str(Inf,l):
    return str(read_and_dump(Inf,l),encoding='cp437')

#This is a common idiom in the readers, that makes them easier to read.
def check_type(Inf,expected_type):
    type=read_uint16(Inf)
    if type!=expected_type:
        raise ValueError("Expected type %d"%expected_type)


class CMOD2Mesh2:
    # These were all static final ints in Java
    CMOD_Material_ID = 1001  # 0x03E9
    CMOD_EndMaterial_ID = 1002  # 0x03EA
    CMOD_Diffuse_ID = 1003  # 0x03EB
    CMOD_Specular_ID = 1004  # 0x03EC
    CMOD_SpecularPower_ID = 1005  # 0x03ED
    CMOD_Opacity_ID = 1006  # 0x03EE
    CMOD_Texture_ID = 1007  # 0x03EF
    CMOD_Mesh_ID = 1009  # 0x03F1
    CMOD_EndMesh_ID = 1010  # 0x03F2
    CMOD_VertexDesc_ID = 1011  # 0x03F3
    CMOD_EndVertexDesc_ID = 1012  # 0x03F4
    CMOD_Vertices_ID = 1013  # 0x03F5
    CMOD_Emissive_ID = 1014  # 0x03F6
    CMOD_Blend_ID = 1015  # 0x03F7
    CMOD_Float1_ID = 1
    CMOD_Float2_ID = 2
    CMOD_Float3_ID = 3
    CMOD_Float4_ID = 4
    CMOD_String_ID = 5
    CMOD_Uint32_ID = 6
    CMOD_Color_ID = 7
    #Occasionally in Java, we used for some Enum x, x[i] to
    #get the i'th member of the enum. In Python, use x.__members__[i].
    class VertexSemantic(Enum):
        def __init__(self, LpovrayStart=None):
            self.povrayStart = LpovrayStart
        Position= "vertex_vectors"
        Color0=   "Color0"
        Color1=   "Color1"
        Normal=   "normal_vectors"
        Tangent=  "Tangent"
        Texture0= "uv_vectors"
        Texture1= "Texture1"
        Texture2= "Texture2"
        Texture3= "Texture3"
        PointSize="Color0"

    class VertexReader:
        def __init__(self, stride):
            self.stride = stride
        def read(self, inf):
            result = np.zeros(self.stride)
            for i in range(self.stride):
                result[i] = read_float32(inf)
            return result
    class UByte4Reader:
        def __init__(self,stride):
            self.stride=stride
        def read(self, inf):
            result = np.zeros(self.stride)
            for i in range(self.stride):
                result[i] = struct.unpack("B", inf.read(1))[0]
            return result

    VertexFormatList=[
        VertexReader(1), #Float1
        VertexReader(2), #Float2
        VertexReader(3), #Float3
        VertexReader(4), #Float4
        UByte4Reader(4)  #UByte4
    ]
    VertexFormatNameList=["Float1",
                          "Float2",
                          "Float3",
                          "Float4",
                          "UByte4"]

    class PrimitiveType(Enum):
        TriList = 0
        TriStrip = 1
        TriFan = 2
        LineList = 3
        LineStrip = 4
        PointList = 5
        SpriteList = 6

    class CMOD_Float1:
        def __init__(self, Ldata):
            try:  # try using Ldata as a file-like object
                check_type(Ldata,CMOD2Mesh2.CMOD_Float1_ID)
                self.data = read_float32(Ldata)
            except AttributeError:
                # Not a file-like, assume that it contains the data itself
                self.data = Ldata

    class CMOD_Float2:
        def __init__(self, Inf):
            check_type(Inf, CMOD2Mesh2.CMOD_Float2_ID)
            self.data1 = read_float32(Inf)
            self.data2 = read_float32(Inf)

    class CMOD_Float3:
        def __init__(self, Inf):
            check_type(Inf, CMOD2Mesh2.CMOD_Float3_ID)
            self.data1 = read_float32(Inf)
            self.data2 = read_float32(Inf)
            self.data3 = read_float32(Inf)

    class CMOD_Float4:
        def __init__(self, Inf):
            check_type(Inf, CMOD2Mesh2.CMOD_Float4_ID)
            self.data1 = read_float32(Inf)
            self.data2 = read_float32(Inf)
            self.data3 = read_float32(Inf)
            self.data4 = read_float32(Inf)

    class CMOD_String:
        def __init__(self,Inf):
            check_type(Inf, CMOD2Mesh2.CMOD_String_ID)
            l=struct.unpack("<h",Inf.read(2))[0]
            self.s=read_str(Inf,l)

    class CMOD_Uint32:
        def __init__(self,Inf):
            check_type(Inf, CMOD2Mesh2.CMOD_Uint32_ID)
            self.data=struct.unpack("<I",Inf.read(4))[0]

    class CMOD_Color:
        def __str__(self):
            return "rgb <%f,%f,%f>"%(self.red,self.green,self.blue)
        def __init__(self,Lred,Lgreen=None,Lblue=None):
            if Lgreen is not None:
                self.red  =Lred
                self.green=Lgreen
                self.blue =Lblue
            else:
                check_type(Lred, CMOD2Mesh2.CMOD_Color_ID)
                self.red  =read_float32(Lred)
                self.green=read_float32(Lred)
                self.blue =read_float32(Lred)
                print(str(self),end='')
    class CMOD_Material:
        def str_header(self):
            return "texture {\n"
        def str_diffuse_header(self):
            return "  pigment {color "
        def str_diffuse_footer(self):
            return "} //diffuse\n"
        def str_diffuse(self):
            return self.str_diffuse_header()+"%s"%self.diffuse+self.str_diffuse_footer()
        def str_specular_header(self):
            return "  finish {reflection "
        def str_specular_footer(self):
            return "} //specular\n"
        def str_specular(self):
            return self.str_specular_header()+"%s"%self.specular+self.str_specular_footer()
        def str_emissive_header(self):
            return "  finish {ambient "
        def str_emissive_footer(self):
            return "} //emissive\n"
        def str_emissive(self):
            return (self.str_emissive_header()+"%s"%self.emissive+self.str_emissive_footer()) if self.emissive is not None else "  //no emissive\n"
        def str_phong(self):
            return "  finish {phong %f phong_size %f} //specularPower\n"%(self.specular.red, self.specularPower.data) if self.specular.red>0 else "  //no specularPower\n"
        def str_opacity(self):
            return "  //opacity %f\n"%self.opacity.data
        def str_blend(self):
            return "  //blend %d\n"%self.blend
        def str_footer(self):
            return "}\n"
        def __str__(self):
            result =self.str_header()
            result+=self.str_diffuse()
            result+=self.str_specular()
            result+=self.str_emissive()
            result+=self.str_phong()
            result+=self.str_blend()
            result+=self.str_opacity()
            result+=self.str_footer()
            return result
        def __init__(self,Inf):
            self.diffuse=CMOD2Mesh2.CMOD_Color(0,0,0)
            self.specular=CMOD2Mesh2.CMOD_Color(0,0,0)
            self.emissive=None #CMOD_Color
            self.specularPower=CMOD2Mesh2.CMOD_Float1(1)
            self.opacity=CMOD2Mesh2.CMOD_Float1(1)
            self.blend=0 #int16
            self.texture=[None]*4 #CMOD_String[4]
            print(self.str_header(),end='')
            type=read_uint16(Inf)
            while type!=CMOD2Mesh2.CMOD_EndMaterial_ID:
                if type==CMOD2Mesh2.CMOD_Diffuse_ID:
                    print(self.str_diffuse_header(),end='')
                    self.diffuse=CMOD2Mesh2.CMOD_Color(Inf)
                    print(self.str_diffuse_footer(),end='')
                elif type==CMOD2Mesh2.CMOD_Specular_ID:
                    print(self.str_specular_header(),end='')
                    self.specular=CMOD2Mesh2.CMOD_Color(Inf)
                    print(self.str_specular_footer(),end='')
                elif type==CMOD2Mesh2.CMOD_Emissive_ID:
                    print(self.str_emissive_header(),end='')
                    self.emissive=CMOD2Mesh2.CMOD_Color(Inf)
                    print(self.str_emissive_footer(),end='')
                elif type==CMOD2Mesh2.CMOD_Opacity_ID:
                    self.opacity=CMOD2Mesh2.CMOD_Float1(Inf)
                    print(self.str_opacity(),end='')
                elif type==CMOD2Mesh2.CMOD_SpecularPower_ID:
                    self.specularPower=CMOD2Mesh2.CMOD_Float1(Inf)
                    print(self.str_phong(),end='')
                elif type==CMOD2Mesh2.CMOD_Blend_ID:
                    self.blend=read_int16(Inf)
                    print(self.str_blend(),end='')
                elif type==CMOD2Mesh2.CMOD_Texture_ID:
                    texType=read_int16(Inf)
                    self.texture[texType]=CMOD2Mesh2.CMOD_String(Inf)
                    print("  #declare Texture[%d]=\"%s\"\n"%(texType,self.texture[texType].s));
                type=read_uint16(Inf)
            print(self.str_footer())

    class CMOD_VertexDescription:
        def str_header(self):
            return "  /*VertexDesc {\n"
        def str_footer(self):
            return "  }*/\n"
        def __str__(self):
            result =self.str_header()
            for i in range(len(self.semantic)):
                result+="    %s (%s) //%d\n"%(str(self.semantic[i]),str(self.format[i]),i)
            result+=self.str_footer()
            return result
        def __init__(self,Inf):
            print(self.str_header(),end='')
            self.semantic=[] #List<VertexSemantic>
            self.format=[]   #List<VertexFormat>
            check_type(Inf,CMOD2Mesh2.CMOD_VertexDesc_ID)
            type=read_int16(Inf)
            i=0
            while type!=CMOD2Mesh2.CMOD_EndVertexDesc_ID:
                fmt=read_int16(Inf)
                self.semantic.append(CMOD2Mesh2.VertexSemantic[list(CMOD2Mesh2.VertexSemantic.__members__)[type]])
                self.format.append(CMOD2Mesh2.VertexFormatList[fmt])
                print("    type[%d]%s fmt[%d](%s) //%d\n"%(type,str(self.semantic[i]),fmt,CMOD2Mesh2.VertexFormatNameList[fmt],i),end='')
                i+=1
                type=read_int16(Inf)
            print(self.str_footer())
    class CMOD_Vertices:
        def __str__(self):
            result=""
            for s in self.desc.semantic:
                d=self.data[s]
                result+="  %s {\n"%s.povrayStart
                result+="    %d"%len(d)
                for i in range(len(d)):
                    result+=",\n    <%f"%d[i][0]
                    for j in range(1,len(d[i])):
                        result+=",%f"%d[i][j]
                    result+=">"
            result+="\n  }"
            return result
        def __init__(self,Inf,Ldesc):
            self.desc=Ldesc #CMOD_VertexDescription
            self.length=None #long
            self.data={} #Map<VertexSemantic,float[][]>
            check_type(Inf,CMOD2Mesh2.CMOD_Vertices_ID)
            self.length=read_uint32(Inf)
            for i in range(len(self.desc.format)):
                self.data[self.desc.semantic[i]]=[[]]*self.length
            for i in range(self.length):
                for j in range(len(self.desc.format)):
                    self.data[self.desc.semantic[j]][i]=self.desc.format[j].read(Inf)
            print(str(self))

    class PrimitiveGroup:
        def __str__(self):
            result ="  face_indices {\n"
            result+="    %d"%(len(self.data)//3)
            for i in range(len(self.data)//3):
                result+=",\n    <%d"%self.data[i*3+0]
                for j in range(1,3):
                    result+=",%d"%self.data[i*3+j]
                result+=">"
            result+="\n  }"
            return result
        def __init__(self,Ltype,Inf):
            self.type=list(CMOD2Mesh2.PrimitiveType.__members__)[Ltype]    #PrimitiveType
            self.materialIdx=read_int32(Inf) #int
            indexCount=read_int32(Inf)
            self.data = [0]*indexCount  # int[]
            for i in range(indexCount):
                self.data[i]=read_int32(Inf)
            print(str(self))
    class CMOD_Mesh:
        def str_header(self):
            return "mesh2 {\n"
        def str_footer(self):
            return "} //mesh2\n"
        def __str__(self):
            result =self.str_header()
            result+=str(self.V)
            for G in self.groups:
                result+=str(G)
            result+=self.str_footer()
            return result
        def __init__(self,Inf):
            print(self.str_header())
            self.desc=CMOD2Mesh2.CMOD_VertexDescription(Inf) #CMOD_VertexDescription
            self.V=CMOD2Mesh2.CMOD_Vertices(Inf,self.desc)   #CMOD_Vertices
            self.groups=[] #List<PrimitiveGroup>
            type=read_int16(Inf)
            while type!=CMOD2Mesh2.CMOD_EndMesh_ID:
                self.groups.append(CMOD2Mesh2.PrimitiveGroup(type,Inf))
                type=read_int16(Inf)
            print(self.str_footer())
    def __str__(self):
        i=0
        result="#declare materialArray=array[%d];"%len(self.materials)
        for M in self.materials:
            result+="#declare material%03d="%i
            if i==769:
                M.specular.red=0
                M.specular.green=0
                M.specular.blue=0
                M.diffuse.red=1
                M.diffuse.green=1
                M.diffuse.blue=1
                M.specularPower.data=4.0
            result+=str(M)
            result+="#declare materialArray[%d]=material%03d;"%(i,i)
            result+="\n"
            i+=1
        i=0;
        result+="#declare partArray=array[%d];"%len(self.meshes)
        for M in self.meshes:
            result+="#declare partArray[%d]=%s\n"%(i,str(M))
            i+=1
        result+="#declare Voyager=union {\n"
        result+="  #local I=0;\n"
        result+="  #while(I<dimension_size(partArray,1))\n"
        result+="    object{partArray[I] texture {materialArray[I]}}\n"
        result+="    #local I=I+1;\n"
        result+="  #end\n"
        result+="}\n"
        return result
    def __init__(self,Inf):
        self.header="" #String
        self.materials=[] #LinkedList<CMOD_Material>
        self.meshes=[] #LinkedList<CMOD_Mesh>
        b=read_str(Inf,16)
        type=read_int16(Inf)
        seenMesh=False
        materials=0
        while type>0:
            if type==CMOD2Mesh2.CMOD_Material_ID:
                if seenMesh:
                    raise ValueError("All materials must be before any meshes")
                print("#declare material[%3d]="%materials,end='')
                materials+=1
                self.materials.append(CMOD2Mesh2.CMOD_Material(Inf))
            elif type==CMOD2Mesh2.CMOD_Mesh_ID:
                seenMesh=True
                self.meshes.append(CMOD2Mesh2.CMOD_Mesh(Inf))
            try:
                type=read_int16(Inf)
            except:
                type=-1

if __name__=="__main__":
    with open("Data/Mesh/voyager.cmod", "rb") as Inf:
        try:
            Imp = CMOD2Mesh2(Inf)
        except:
            print("0x%08X"%Inf.tell())
            raise
        else:
            with open("Data/Mesh/voyager.mesh2") as ouf:
                print(str(Imp))
