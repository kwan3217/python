import spiceypy as cspice
import os

old=os.getcwd()
os.chdir('../../../Data/spice/Ranger/')
cspice.furnsh('Ranger7Background.tm')
os.chdir(old)

cspice.furnsh("Ranger7.bsp")
cspice.furnsh("Ranger7.bc")
cspice.furnsh("Ranger7.tf")
cspice.furnsh("Ranger7.tsc")
timpact=     -1117751615.876467
print(cspice.etcal(timpact))
print(cspice.etcal(timpact-35.32527123))
(state,ltime)=cspice.spkezr("-1007",timpact,"J2000","LT","399012")
print(state,ltime)
print(cspice.etcal(timpact+ltime))
print(cspice.etcal(timpact-35.32527123+ltime))
tdt=cspice.str2et("1964-Jul-31 13:25:00 TDT")
tdb=cspice.str2et("1964-Jul-31 13:25:00 TDB")
print(tdt,tdb,tdb-tdt)
