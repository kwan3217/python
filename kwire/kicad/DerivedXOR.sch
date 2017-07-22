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
EELAYER 26 0
EELAYER END
$Descr A4 11693 8268
encoding utf-8
Sheet 10 18
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
L VCC #PWR801
U 1 1 597217F5
P 4300 2700
AR Path="/596E3F12/59691506/596FACCF/59722BD8/597217F5" Ref="#PWR801"  Part="1" 
AR Path="/596E3F12/59691506/596FACCF/597233AA/597217F5" Ref="#PWR901"  Part="1" 
AR Path="/596E3F12/59691506/59723831/597217F5" Ref="#PWR1001"  Part="1" 
AR Path="/596E3F12/596919C0/596FACCF/59722BD8/597217F5" Ref="#PWR1601"  Part="1" 
AR Path="/596E3F12/596919C0/596FACCF/597233AA/597217F5" Ref="#PWR1701"  Part="1" 
AR Path="/596E3F12/596919C0/59723831/597217F5" Ref="#PWR1801"  Part="1" 
F 0 "#PWR1001" H 4300 2550 50  0001 C CNN
F 1 "VCC" H 4300 2850 50  0000 C CNN
F 2 "" H 4300 2700 50  0000 C CNN
F 3 "" H 4300 2700 50  0000 C CNN
	1    4300 2700
	1    0    0    -1  
$EndComp
$Comp
L GND #PWR802
U 1 1 59721809
P 4300 4250
AR Path="/596E3F12/59691506/596FACCF/59722BD8/59721809" Ref="#PWR802"  Part="1" 
AR Path="/596E3F12/59691506/596FACCF/597233AA/59721809" Ref="#PWR902"  Part="1" 
AR Path="/596E3F12/59691506/59723831/59721809" Ref="#PWR1002"  Part="1" 
AR Path="/596E3F12/596919C0/596FACCF/59722BD8/59721809" Ref="#PWR1602"  Part="1" 
AR Path="/596E3F12/596919C0/596FACCF/597233AA/59721809" Ref="#PWR1702"  Part="1" 
AR Path="/596E3F12/596919C0/59723831/59721809" Ref="#PWR1802"  Part="1" 
F 0 "#PWR1002" H 4300 4000 50  0001 C CNN
F 1 "GND" H 4300 4100 50  0000 C CNN
F 2 "" H 4300 4250 50  0000 C CNN
F 3 "" H 4300 4250 50  0000 C CNN
	1    4300 4250
	1    0    0    -1  
$EndComp
$Comp
L NMOS Q802
U 1 1 59722202
P 4200 3600
AR Path="/596E3F12/59691506/596FACCF/59722BD8/59722202" Ref="Q802"  Part="1" 
AR Path="/596E3F12/59691506/596FACCF/597233AA/59722202" Ref="Q902"  Part="1" 
AR Path="/596E3F12/59691506/59723831/59722202" Ref="Q1002"  Part="1" 
AR Path="/596E3F12/596919C0/596FACCF/59722BD8/59722202" Ref="Q1602"  Part="1" 
AR Path="/596E3F12/596919C0/596FACCF/597233AA/59722202" Ref="Q1702"  Part="1" 
AR Path="/596E3F12/596919C0/59723831/59722202" Ref="Q1802"  Part="1" 
F 0 "Q1002" H 4240 3600 50  0000 L CNN
F 1 "NMOS" H 4335 3550 50  0001 L CNN
F 2 "KwanSystems:SOT-23" H 4400 3700 50  0001 C CNN
F 3 "" H 4200 3600 50  0000 C CNN
	1    4200 3600
	1    0    0    -1  
$EndComp
$Comp
L PMOS Q801
U 1 1 59722235
P 4200 3000
AR Path="/596E3F12/59691506/596FACCF/59722BD8/59722235" Ref="Q801"  Part="1" 
AR Path="/596E3F12/59691506/596FACCF/597233AA/59722235" Ref="Q901"  Part="1" 
AR Path="/596E3F12/59691506/59723831/59722235" Ref="Q1001"  Part="1" 
AR Path="/596E3F12/596919C0/596FACCF/59722BD8/59722235" Ref="Q1601"  Part="1" 
AR Path="/596E3F12/596919C0/596FACCF/597233AA/59722235" Ref="Q1701"  Part="1" 
AR Path="/596E3F12/596919C0/59723831/59722235" Ref="Q1801"  Part="1" 
F 0 "Q1001" H 4240 3000 50  0000 L CNN
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
L PMOS Q803
U 1 1 59722D45
P 4850 3000
AR Path="/596E3F12/59691506/596FACCF/59722BD8/59722D45" Ref="Q803"  Part="1" 
AR Path="/596E3F12/59691506/596FACCF/597233AA/59722D45" Ref="Q903"  Part="1" 
AR Path="/596E3F12/59691506/59723831/59722D45" Ref="Q1003"  Part="1" 
AR Path="/596E3F12/596919C0/596FACCF/59722BD8/59722D45" Ref="Q1603"  Part="1" 
AR Path="/596E3F12/596919C0/596FACCF/597233AA/59722D45" Ref="Q1703"  Part="1" 
AR Path="/596E3F12/596919C0/59723831/59722D45" Ref="Q1803"  Part="1" 
F 0 "Q1003" H 4890 3000 50  0000 L CNN
F 1 "PMOS" H 4985 2950 50  0001 L CNN
F 2 "KwanSystems:SOT-23" H 5050 3100 50  0001 C CNN
F 3 "" H 4850 3000 50  0000 C CNN
	1    4850 3000
	1    0    0    -1  
$EndComp
$Comp
L PMOS Q805
U 1 1 59722D7D
P 5400 3000
AR Path="/596E3F12/59691506/596FACCF/59722BD8/59722D7D" Ref="Q805"  Part="1" 
AR Path="/596E3F12/59691506/596FACCF/597233AA/59722D7D" Ref="Q905"  Part="1" 
AR Path="/596E3F12/59691506/59723831/59722D7D" Ref="Q1005"  Part="1" 
AR Path="/596E3F12/596919C0/596FACCF/59722BD8/59722D7D" Ref="Q1605"  Part="1" 
AR Path="/596E3F12/596919C0/596FACCF/597233AA/59722D7D" Ref="Q1705"  Part="1" 
AR Path="/596E3F12/596919C0/59723831/59722D7D" Ref="Q1805"  Part="1" 
F 0 "Q1005" H 5440 3000 50  0000 L CNN
F 1 "PMOS" H 5535 2950 50  0001 L CNN
F 2 "KwanSystems:SOT-23" H 5600 3100 50  0001 C CNN
F 3 "" H 5400 3000 50  0000 C CNN
	1    5400 3000
	1    0    0    -1  
$EndComp
$Comp
L NMOS Q804
U 1 1 59722DAB
P 4850 3600
AR Path="/596E3F12/59691506/596FACCF/59722BD8/59722DAB" Ref="Q804"  Part="1" 
AR Path="/596E3F12/59691506/596FACCF/597233AA/59722DAB" Ref="Q904"  Part="1" 
AR Path="/596E3F12/59691506/59723831/59722DAB" Ref="Q1004"  Part="1" 
AR Path="/596E3F12/596919C0/596FACCF/59722BD8/59722DAB" Ref="Q1604"  Part="1" 
AR Path="/596E3F12/596919C0/596FACCF/597233AA/59722DAB" Ref="Q1704"  Part="1" 
AR Path="/596E3F12/596919C0/59723831/59722DAB" Ref="Q1804"  Part="1" 
F 0 "Q1004" H 4890 3600 50  0000 L CNN
F 1 "NMOS" H 4985 3550 50  0001 L CNN
F 2 "KwanSystems:SOT-23" H 5050 3700 50  0001 C CNN
F 3 "" H 4850 3600 50  0000 C CNN
	1    4850 3600
	1    0    0    -1  
$EndComp
$Comp
L NMOS Q806
U 1 1 59722DE2
P 5400 3600
AR Path="/596E3F12/59691506/596FACCF/59722BD8/59722DE2" Ref="Q806"  Part="1" 
AR Path="/596E3F12/59691506/596FACCF/597233AA/59722DE2" Ref="Q906"  Part="1" 
AR Path="/596E3F12/59691506/59723831/59722DE2" Ref="Q1006"  Part="1" 
AR Path="/596E3F12/596919C0/596FACCF/59722BD8/59722DE2" Ref="Q1606"  Part="1" 
AR Path="/596E3F12/596919C0/596FACCF/597233AA/59722DE2" Ref="Q1706"  Part="1" 
AR Path="/596E3F12/596919C0/59723831/59722DE2" Ref="Q1806"  Part="1" 
F 0 "Q1006" H 5440 3600 50  0000 L CNN
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
L LED_Small D801
U 1 1 59725495
P 5700 3400
AR Path="/596E3F12/59691506/596FACCF/59722BD8/59725495" Ref="D801"  Part="1" 
AR Path="/596E3F12/59691506/596FACCF/597233AA/59725495" Ref="D901"  Part="1" 
AR Path="/596E3F12/59691506/59723831/59725495" Ref="D1001"  Part="1" 
AR Path="/596E3F12/596919C0/596FACCF/59722BD8/59725495" Ref="D1601"  Part="1" 
AR Path="/596E3F12/596919C0/596FACCF/597233AA/59725495" Ref="D1701"  Part="1" 
AR Path="/596E3F12/596919C0/59723831/59725495" Ref="D1801"  Part="1" 
F 0 "D1001" H 5650 3525 50  0000 L CNN
F 1 "LED_Small" H 5525 3300 50  0000 L CNN
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
L RESISTOR R801
U 1 1 59725531
P 5700 3800
AR Path="/596E3F12/59691506/596FACCF/59722BD8/59725531" Ref="R801"  Part="1" 
AR Path="/596E3F12/59691506/596FACCF/597233AA/59725531" Ref="R901"  Part="1" 
AR Path="/596E3F12/59691506/59723831/59725531" Ref="R1001"  Part="1" 
AR Path="/596E3F12/596919C0/596FACCF/59722BD8/59725531" Ref="R1601"  Part="1" 
AR Path="/596E3F12/596919C0/596FACCF/597233AA/59725531" Ref="R1701"  Part="1" 
AR Path="/596E3F12/596919C0/59723831/59725531" Ref="R1801"  Part="1" 
F 0 "R1001" H 5650 3850 45  0000 L BNN
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
$EndSCHEMATC