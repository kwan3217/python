EESchema Schematic File Version 2
LIBS:power
LIBS:device
LIBS:transistors
LIBS:conn
LIBS:linear
LIBS:regul
LIBS:74xx
LIBS:cmos4000
LIBS:adc-dac
LIBS:memory
LIBS:xilinx
LIBS:microcontrollers
LIBS:dsp
LIBS:microchip
LIBS:analog_switches
LIBS:motorola
LIBS:texas
LIBS:intel
LIBS:audio
LIBS:interface
LIBS:digital-audio
LIBS:philips
LIBS:display
LIBS:cypress
LIBS:siliconi
LIBS:opto
LIBS:atmel
LIBS:contrib
LIBS:valves
LIBS:kwire
LIBS:switches
LIBS:KwanSystems
LIBS:xordriver-cache
EELAYER 25 0
EELAYER END
$Descr A4 11693 8268
encoding utf-8
Sheet 2 2
Title ""
Date ""
Rev ""
Comp ""
Comment1 ""
Comment2 ""
Comment3 ""
Comment4 ""
$EndDescr
$Comp
L VCC #PWR05
U 1 1 597217F5
P 4300 2700
F 0 "#PWR05" H 4300 2550 50  0001 C CNN
F 1 "VCC" H 4300 2850 50  0000 C CNN
F 2 "" H 4300 2700 50  0000 C CNN
F 3 "" H 4300 2700 50  0000 C CNN
	1    4300 2700
	1    0    0    -1  
$EndComp
$Comp
L GND #PWR06
U 1 1 59721809
P 4300 4250
F 0 "#PWR06" H 4300 4000 50  0001 C CNN
F 1 "GND" H 4300 4100 50  0000 C CNN
F 2 "" H 4300 4250 50  0000 C CNN
F 3 "" H 4300 4250 50  0000 C CNN
	1    4300 4250
	1    0    0    -1  
$EndComp
$Comp
L NMOS Q202
U 1 1 59722202
P 4200 3600
F 0 "Q202" H 4240 3600 50  0000 L CNN
F 1 "NMOS" H 4335 3550 50  0001 L CNN
F 2 "KwanSystems:SOT-23" H 4400 3700 50  0001 C CNN
F 3 "" H 4200 3600 50  0000 C CNN
	1    4200 3600
	1    0    0    -1  
$EndComp
$Comp
L PMOS Q201
U 1 1 59722235
P 4200 3000
F 0 "Q201" H 4240 3000 50  0000 L CNN
F 1 "PMOS" H 4335 2950 50  0001 L CNN
F 2 "KwanSystems:SOT-23" H 4400 3100 50  0001 C CNN
F 3 "" H 4200 3000 50  0000 C CNN
	1    4200 3000
	1    0    0    -1  
$EndComp
Wire Wire Line
	4000 3000 4000 3850
Wire Wire Line
	4300 3200 4300 3400
Connection ~ 4300 3300
Wire Wire Line
	3650 3300 4000 3300
Connection ~ 4000 3300
Text HLabel 3650 3300 0    60   Input ~ 0
A
Text HLabel 5850 3300 2    60   Output ~ 0
Y
Wire Wire Line
	4300 2700 4300 2800
Wire Wire Line
	4300 3800 4300 4250
Text HLabel 3650 3950 0    60   Input ~ 0
B
$Comp
L PMOS Q205
U 1 1 59722D7D
P 5400 3000
F 0 "Q205" H 5440 3000 50  0000 L CNN
F 1 "PMOS" H 5535 2950 50  0001 L CNN
F 2 "KwanSystems:SOT-23" H 5600 3100 50  0001 C CNN
F 3 "" H 5400 3000 50  0000 C CNN
	1    5400 3000
	1    0    0    -1  
$EndComp
$Comp
L NMOS Q206
U 1 1 59722DE2
P 5400 3600
F 0 "Q206" H 5440 3600 50  0000 L CNN
F 1 "NMOS" H 5535 3550 50  0001 L CNN
F 2 "KwanSystems:SOT-23" H 5600 3700 50  0001 C CNN
F 3 "" H 5400 3600 50  0000 C CNN
	1    5400 3600
	1    0    0    -1  
$EndComp
Wire Wire Line
	4300 3300 4500 3300
Wire Wire Line
	4500 3300 4500 2800
Wire Wire Line
	4500 2800 4950 2800
Wire Wire Line
	4650 3950 3650 3950
Wire Wire Line
	4650 3000 4650 3950
Connection ~ 4650 3600
Wire Wire Line
	4000 3850 4950 3850
Wire Wire Line
	4950 3850 4950 3800
Connection ~ 4000 3600
Wire Wire Line
	4950 3200 4950 3400
Wire Wire Line
	4950 3300 5200 3300
Wire Wire Line
	5200 3000 5200 3600
Connection ~ 4950 3300
Connection ~ 5200 3300
Wire Wire Line
	5500 3800 5500 4250
Wire Wire Line
	4300 4250 5700 4250
Wire Wire Line
	4300 2700 5500 2700
Wire Wire Line
	5500 2700 5500 2800
Wire Wire Line
	5500 3200 5500 3400
$Comp
L LED_Small D201
U 1 1 59725495
P 5700 3400
F 0 "D201" H 5650 3525 50  0000 L CNN
F 1 "Y" H 5700 3300 50  0000 L CNN
F 2 "KwanSystems:D_0603" V 5700 3400 50  0001 C CNN
F 3 "" V 5700 3400 50  0000 C CNN
	1    5700 3400
	0    -1   -1   0   
$EndComp
Wire Wire Line
	5500 3300 5850 3300
Connection ~ 5500 3300
Connection ~ 5700 3300
$Comp
L RESISTOR R201
U 1 1 59725531
P 5700 3800
F 0 "R201" H 5650 3850 45  0000 L BNN
F 1 "RESISTOR" H 5650 3700 45  0001 L BNN
F 2 "KwanSystems:SMD_0402" H 5755 3950 20  0001 C CNN
F 3 "" H 6125 3500 60  0001 C CNN
	1    5700 3800
	0    1    1    0   
$EndComp
Wire Wire Line
	5700 3500 5700 3700
Wire Wire Line
	5700 4250 5700 3900
Connection ~ 5500 4250
Text Label 4550 2800 0    60   ~ 0
#A
Text Label 5000 3300 0    60   ~ 0
#Y
Text Label 5700 3650 0    60   ~ 0
Y-
Wire Wire Line
	4950 3000 5050 3000
Wire Wire Line
	5050 3000 5050 2700
Connection ~ 5050 2700
Wire Wire Line
	4950 3600 5050 3600
Wire Wire Line
	5050 3600 5050 4250
Connection ~ 5050 4250
$Comp
L PMOS_BODY Q203
U 1 1 59766780
P 4850 3000
F 0 "Q203" H 4980 3040 50  0000 L CNN
F 1 "PMOS_BODY" H 4985 2950 50  0001 L CNN
F 2 "" H 5050 3100 50  0000 C CNN
F 3 "" H 4850 3000 50  0000 C CNN
	1    4850 3000
	1    0    0    -1  
$EndComp
$Comp
L NMOS_BODY Q204
U 1 1 597667A9
P 4850 3600
F 0 "Q204" H 4980 3640 50  0000 L CNN
F 1 "NMOS_BODY" H 4985 3550 50  0001 L CNN
F 2 "" H 5050 3700 50  0000 C CNN
F 3 "" H 4850 3600 50  0000 C CNN
	1    4850 3600
	1    0    0    -1  
$EndComp
Text Notes 6100 1500 0    60   ~ 0
Theory of operation: This is a chain of 3 inverters. \nThe first one (Qx01 and Qx02)is used to control \nthe power rails of the second one (Qx03 and Qx04).
Text Notes 6100 2950 0    60   ~ 0
Since the power rails of the second inverter aren't\nthe normal VCC, there are two issues:
Text Notes 6100 3900 0    60   ~ 0
1) A normal discrete FET has its body and source tied\ntogether. This sets the reference point for the gate. If\nthe power rails are switched, then the NMOS will need\nits gate higher than high and the PMOS lower than low\nin order to turn on the switch. We overcome this by\nusing transistors where the body is tied directly to VCC\nor GND independent of the source. This is easy on an\nintegrated circuit, but discrete transistors with 4\nterminals are rare. 
Text Notes 6100 4550 0    60   ~ 0
2) Since the middle inverter has its power from the\noutput of another inverter rather than from the rails,\nthe voltage swing of the middle gate won't be as\nwide. This will cause trouble with fanout. We\novercome this by adding the third inverter (Qx05\nand Qx06), which is powered directly from the rails.
Text Notes 6100 5400 0    60   ~ 0
\nDetails in Lee and Sobelman, \n"New Low-Voltage Circuits for XOR and XNOR", \nDOI 10.1109/SECON.1997.598676 . This circuit is \ndrawn in figure 2, and is NOT the recommended\ncircuit in that paper. Their recommendation uses\n7 transistors (odd to have an odd number in\nCMOS) and is optimized for speed, while I\noptimize to minimize transistor count.
Text Notes 8800 3500 0    60   ~ 0
An alternate solution is to use a CD4007.\nThis is used in the XOR4007.sch alternate\npage. I prefer the discrete solution, since\nit is easier to hide stuff in a 14-pin IC,\ntherefore harder to trust.
Text Notes 6100 2400 0    60   ~ 0
When A is false, Qx01 is on and Qx02 is off.\nThis puts high voltage on the top rail of\nthe second inverter. The low voltage on\nthe bottom comes straight from the input.\nThe inverter acts as normal, so the output\nof the middle inverter is the opposite of\nthe B input. This is fed through the last\ninverter as a buffer, so we get a copy of\nthe B input as the final output.
Text Notes 8800 2600 0    60   ~ 0
When A is true, Qx01 is off and Qx02 is on.\nThis puts low voltage ob the top rail of the\nsecond inverter. The bottom voltage \ncomes from the A input directly, and is\ntherefore high. Since the rails are reversed,\nthe second inverter is no longer an inverter,\nsince a true B signal will turn on the NMOS\nand let through the high voltage on the\nbottom rail, and turn off the PMOS and\nblock the low voltage on the top rail, and\nvice versa for false B.
$EndSCHEMATC
