"""
Registers
"""
r=[0]*256  # 256 registers. Is that too many?
PC=0xFF
LR=0xFE
SP=0xFD
FLAGS=0xFC
mem=[0]*65536  # Initially we only deal with the lower 64kiw of memory. At some point
               # we should probably make this a list of lists, and only populate them
               # as they are used. To do otherwise would claim 16GiB of space.
               
def MOV(D,A):
    """Copy data from rA to rD"""
    r[D]=r[A]
    
def readMem(addr):
    """Read main memory, returning the value at the given address."""                       
    return mem[addr]

def writeMem(addr,value):
    """Write main memory at the given address with the given value"""
    mem[addr]=value
    
def LDRimm(D,addr):
    """Load Register. Use the given address, read memory at that point and set rD to that value"""
    r[D]=readMem(addr)
    
def LDRReg(D,A,B):
    """Load Register. Use rA+B as address, read memory and set rD to that value"""
    r[D]=readMem(r[A]+B)
    
def LDRRegReg(D,A,B):
    """Load Register. Use rA+rB as address, read memory and set rD to that value"""
    r[D]=readMem(r[A]+r[B])
    
def POP(D,A):
    """Pop stack. Use rA as address, read memory and set rD to that value, and then increment rA. Register rA will most commonly be SP."""
    r[D]=readMem(r[A])
    r[A]=r[A]+1

def STRimm(D,addr):
    """Store Register. Use the given address, write value in rD to that point in memory"""
    writeMem(addr,r[D])
    
def STRReg(D,A,B):
    """Store Register. Use rA+B as address, write value in rD to that point in memory"""
    writeMem(r[A]+B,r[D])
    
def STRRegReg(D,A,B):
    """Store Register. Use rA+rB as address, write value in rD to that point in memory"""
    writeMem(r[A]+r[B],r[D])
    
def PUSH(D,A):
    """Push stack. Decrement rA, use rA as address, write value in rD to memory at that address. Register rA will most commonly be SP."""
    r[A]=r[A]-1
    writeMem(r[A],r[D])

writeMem(0x300,0x1234)
LDRimm(0x12,0x300) #Should have 0x1234 in rx12    
print(r[0x12])        