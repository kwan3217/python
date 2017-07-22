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
Sheet 1 18
Title ""
Date ""
Rev ""
Comp ""
Comment1 ""
Comment2 ""
Comment3 ""
Comment4 ""
$EndDescr
$Sheet
S 2550 1300 750  700 
U 596E3F12
F0 "bit01" 60
F1 "2bitaddsub.sch" 60
F2 "SUB" I L 2550 1950 60 
F3 "A0" I L 2550 1350 60 
F4 "A1" I L 2550 1450 60 
F5 "B0" I L 2550 1550 60 
F6 "B1" I L 2550 1650 60 
F7 "C0" I L 2550 1850 60 
F8 "S0" O R 3300 1350 60 
F9 "S1" O R 3300 1450 60 
F10 "C2" O R 3300 1950 60 
$EndSheet
Wire Wire Line
	1850 1350 2550 1350
Wire Wire Line
	2050 1450 2550 1450
Wire Wire Line
	3300 1350 3950 1350
Wire Wire Line
	3300 1950 3550 1950
Text Label 2400 1350 0    60   ~ 0
A0
Text Label 2400 1450 0    60   ~ 0
A1
Text Label 3500 1350 0    60   ~ 0
S0
Text Label 3450 1950 0    60   ~ 0
K
$Comp
L GND #PWR102
U 1 1 5970F080
P 1450 3250
F 0 "#PWR102" H 1450 3000 50  0001 C CNN
F 1 "GND" H 1450 3100 50  0000 C CNN
F 2 "" H 1450 3250 50  0000 C CNN
F 3 "" H 1450 3250 50  0000 C CNN
	1    1450 3250
	1    0    0    -1  
$EndComp
$Comp
L VCC #PWR101
U 1 1 5970F0C4
P 1450 1150
F 0 "#PWR101" H 1450 1000 50  0001 C CNN
F 1 "VCC" H 1450 1300 50  0000 C CNN
F 2 "" H 1450 1150 50  0000 C CNN
F 3 "" H 1450 1150 50  0000 C CNN
	1    1450 1150
	1    0    0    -1  
$EndComp
$Comp
L PWR_FLAG #FLG101
U 1 1 5970F13A
P 1450 1150
F 0 "#FLG101" H 1450 1245 50  0001 C CNN
F 1 "PWR_FLAG" H 1450 1330 50  0000 C CNN
F 2 "" H 1450 1150 50  0000 C CNN
F 3 "" H 1450 1150 50  0000 C CNN
	1    1450 1150
	0    -1   -1   0   
$EndComp
$Comp
L PWR_FLAG #FLG102
U 1 1 5970F1C6
P 1450 3250
F 0 "#FLG102" H 1450 3345 50  0001 C CNN
F 1 "PWR_FLAG" H 1450 3430 50  0000 C CNN
F 2 "" H 1450 3250 50  0000 C CNN
F 3 "" H 1450 3250 50  0000 C CNN
	1    1450 3250
	0    -1   -1   0   
$EndComp
$Comp
L LED_Small D101
U 1 1 5970F37E
P 3550 2250
F 0 "D101" H 3500 2375 50  0000 L CNN
F 1 "K" H 3600 2300 50  0000 L CNN
F 2 "KwanSystems:D_0603" V 3550 2250 50  0001 C CNN
F 3 "" V 3550 2250 50  0000 C CNN
	1    3550 2250
	0    -1   -1   0   
$EndComp
Wire Wire Line
	1450 3250 3950 3250
Connection ~ 3750 3250
$Comp
L 8CONTROL S101
U 1 1 5970FD98
P 1650 1350
F 0 "S101" H 1600 1350 50  0000 C CNN
F 1 "A" H 1800 1400 50  0000 C CNN
F 2 "Housings_DIP:DIP-16_W7.62mm" H 1650 1350 50  0001 C CNN
F 3 "" H 1650 1350 50  0000 C CNN
	1    1650 1350
	1    0    0    -1  
$EndComp
$Comp
L 8CONTROL S101
U 2 1 597100DB
P 1850 1450
F 0 "S101" H 1800 1450 50  0000 C CNN
F 1 "A" H 2000 1500 50  0000 C CNN
F 2 "Housings_DIP:DIP-16_W7.62mm" H 1850 1450 50  0001 C CNN
F 3 "" H 1850 1450 50  0000 C CNN
	2    1850 1450
	1    0    0    -1  
$EndComp
$Comp
L 8CONTROL S102
U 1 1 59710157
P 1650 1550
F 0 "S102" H 1600 1550 50  0000 C CNN
F 1 "B" H 1800 1600 50  0000 C CNN
F 2 "Housings_DIP:DIP-16_W7.62mm" H 1650 1550 50  0001 C CNN
F 3 "" H 1650 1550 50  0000 C CNN
	1    1650 1550
	1    0    0    -1  
$EndComp
$Comp
L 8CONTROL S102
U 2 1 5971022B
P 1850 1650
F 0 "S102" H 1800 1650 50  0000 C CNN
F 1 "B" H 2000 1700 50  0000 C CNN
F 2 "Housings_DIP:DIP-16_W7.62mm" H 1850 1650 50  0001 C CNN
F 3 "" H 1850 1650 50  0000 C CNN
	2    1850 1650
	1    0    0    -1  
$EndComp
Wire Wire Line
	1950 1550 2550 1550
Wire Wire Line
	2050 1650 2550 1650
Text Label 2400 1550 0    60   ~ 0
B0
Text Label 2400 1650 0    60   ~ 0
B1
$Comp
L 8CONTROL S103
U 1 1 59710E50
P 1650 1850
F 0 "S103" H 1600 1850 50  0000 C CNN
F 1 "C" H 1800 1900 50  0000 C CNN
F 2 "Housings_DIP:DIP-16_W7.62mm" H 1650 1850 50  0001 C CNN
F 3 "" H 1650 1850 50  0000 C CNN
	1    1650 1850
	1    0    0    -1  
$EndComp
$Comp
L 8CONTROL S103
U 2 1 59710EBD
P 1850 1950
F 0 "S103" H 1800 1950 50  0000 C CNN
F 1 "C" H 2000 2000 50  0000 C CNN
F 2 "Housings_DIP:DIP-16_W7.62mm" H 1850 1950 50  0001 C CNN
F 3 "" H 1850 1950 50  0000 C CNN
	2    1850 1950
	1    0    0    -1  
$EndComp
Wire Wire Line
	1850 1850 2550 1850
Wire Wire Line
	2050 1950 2550 1950
Text Label 2400 1850 0    60   ~ 0
C0
Text Label 2350 1950 0    60   ~ 0
SUB
$Comp
L 8INDICATOR D102
U 2 1 59711E88
P 3950 2250
F 0 "D102" H 4050 2300 50  0000 L CNN
F 1 "S" H 4050 2200 50  0000 L CNN
F 2 "" V 3950 2250 50  0000 C CNN
F 3 "" V 3950 2250 50  0000 C CNN
	2    3950 2250
	0    1    1    0   
$EndComp
$Comp
L 8INDICATOR D102
U 1 1 5971200C
P 3750 2250
F 0 "D102" H 3850 2300 50  0000 L CNN
F 1 "S" H 3850 2200 50  0000 L CNN
F 2 "" V 3750 2250 50  0000 C CNN
F 3 "" V 3750 2250 50  0000 C CNN
	1    3750 2250
	0    1    1    0   
$EndComp
Wire Wire Line
	3300 1450 3750 1450
Text Label 3500 1450 0    60   ~ 0
S1
Connection ~ 3550 3250
$Comp
L RP4 R101
U 1 1 59731188
P 3550 2800
F 0 "R101" V 3508 2858 45  0000 L CNN
F 1 "RP4" V 3592 2858 45  0000 L CNN
F 2 "0402-RES" H 3605 2950 20  0001 C CNN
F 3 "" H 3975 2500 60  0001 C CNN
	1    3550 2800
	0    1    1    0   
$EndComp
$Comp
L RP4 R101
U 2 1 59731246
P 3750 2800
F 0 "R101" V 3708 2858 45  0000 L CNN
F 1 "RP4" V 3792 2858 45  0000 L CNN
F 2 "0402-RES" H 3805 2950 20  0001 C CNN
F 3 "" H 4175 2500 60  0001 C CNN
	2    3750 2800
	0    1    1    0   
$EndComp
Wire Wire Line
	3550 1950 3550 2150
Wire Wire Line
	3550 2350 3550 2700
Wire Wire Line
	3550 2900 3550 3250
Wire Wire Line
	3750 2900 3750 3250
Wire Wire Line
	3750 2350 3750 2700
$Comp
L RP4 R101
U 3 1 59731490
P 3950 2800
F 0 "R101" V 3908 2858 45  0000 L CNN
F 1 "RP4" V 3992 2858 45  0000 L CNN
F 2 "0402-RES" H 4005 2950 20  0001 C CNN
F 3 "" H 4375 2500 60  0001 C CNN
	3    3950 2800
	0    1    1    0   
$EndComp
Wire Wire Line
	3950 2350 3950 2700
Wire Wire Line
	3950 3250 3950 2900
Wire Wire Line
	3950 1350 3950 2150
Wire Wire Line
	3750 1450 3750 2150
$EndSCHEMATC
