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
Sheet 1 9
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
S 2550 1300 600  400 
U 596E3F12
F0 "bit0" 60
F1 "1bitaddsub.sch" 60
F2 "A" I L 2550 1350 60 
F3 "B" I L 2550 1450 60 
F4 "C_IN" I L 2550 1550 60 
F5 "S" O R 3150 1350 60 
F6 "C_OUT" O R 3150 1450 60 
F7 "SUB" I L 2550 1650 60 
$EndSheet
Wire Wire Line
	2000 1350 2550 1350
Wire Wire Line
	2250 1450 2550 1450
Wire Wire Line
	2000 1550 2550 1550
Wire Wire Line
	3150 1350 3950 1350
Wire Wire Line
	3150 1450 3550 1450
Text Label 2300 1350 0    60   ~ 0
A
Text Label 2300 1450 0    60   ~ 0
B
Text Label 2300 1550 0    60   ~ 0
C
Text Label 3250 1350 0    60   ~ 0
S
Text Label 3250 1450 0    60   ~ 0
K
Wire Wire Line
	2250 1650 2550 1650
Text Label 2300 1650 0    60   ~ 0
SUB
$Comp
L GND #PWR01
U 1 1 5970F080
P 1450 2400
F 0 "#PWR01" H 1450 2150 50  0001 C CNN
F 1 "GND" H 1450 2250 50  0000 C CNN
F 2 "" H 1450 2400 50  0000 C CNN
F 3 "" H 1450 2400 50  0000 C CNN
	1    1450 2400
	1    0    0    -1  
$EndComp
$Comp
L VCC #PWR02
U 1 1 5970F0C4
P 1450 1150
F 0 "#PWR02" H 1450 1000 50  0001 C CNN
F 1 "VCC" H 1450 1300 50  0000 C CNN
F 2 "" H 1450 1150 50  0000 C CNN
F 3 "" H 1450 1150 50  0000 C CNN
	1    1450 1150
	1    0    0    -1  
$EndComp
$Comp
L PWR_FLAG #FLG03
U 1 1 5970F13A
P 1450 1150
F 0 "#FLG03" H 1450 1245 50  0001 C CNN
F 1 "PWR_FLAG" H 1450 1330 50  0000 C CNN
F 2 "" H 1450 1150 50  0000 C CNN
F 3 "" H 1450 1150 50  0000 C CNN
	1    1450 1150
	0    -1   -1   0   
$EndComp
$Comp
L PWR_FLAG #FLG04
U 1 1 5970F1C6
P 1450 2400
F 0 "#FLG04" H 1450 2495 50  0001 C CNN
F 1 "PWR_FLAG" H 1450 2580 50  0000 C CNN
F 2 "" H 1450 2400 50  0000 C CNN
F 3 "" H 1450 2400 50  0000 C CNN
	1    1450 2400
	0    -1   -1   0   
$EndComp
$Comp
L LED_Small D101
U 1 1 5970F37E
P 3650 1450
F 0 "D101" H 3600 1575 50  0000 L CNN
F 1 "K" H 3700 1500 50  0000 L CNN
F 2 "KwanSystems:D_0603" V 3650 1450 50  0001 C CNN
F 3 "" V 3650 1450 50  0000 C CNN
	1    3650 1450
	-1   0    0    1   
$EndComp
Wire Wire Line
	1450 2400 3750 2400
Wire Wire Line
	3750 2400 4150 2400
Connection ~ 3750 2400
$Comp
L 8CONTROL S101
U 1 1 5970FD98
P 1800 1350
F 0 "S101" H 1700 1350 50  0000 C CNN
F 1 "A" H 1900 1400 50  0000 C CNN
F 2 "Housings_DIP:DIP-16_W7.62mm" H 1800 1350 50  0001 C CNN
F 3 "" H 1800 1350 50  0000 C CNN
	1    1800 1350
	1    0    0    -1  
$EndComp
$Comp
L 8CONTROL S102
U 1 1 597100DB
P 2050 1450
F 0 "S102" H 1950 1450 50  0000 C CNN
F 1 "B" H 2150 1500 50  0000 C CNN
F 2 "Housings_DIP:DIP-16_W7.62mm" H 2050 1450 50  0001 C CNN
F 3 "" H 2050 1450 50  0000 C CNN
	1    2050 1450
	1    0    0    -1  
$EndComp
$Comp
L 8CONTROL S103
U 1 1 59710157
P 1800 1550
F 0 "S103" H 1700 1550 50  0000 C CNN
F 1 "C" H 1900 1600 50  0000 C CNN
F 2 "Housings_DIP:DIP-16_W7.62mm" H 1800 1550 50  0001 C CNN
F 3 "" H 1800 1550 50  0000 C CNN
	1    1800 1550
	1    0    0    -1  
$EndComp
$Comp
L 8CONTROL S103
U 2 1 5971022B
P 2050 1650
F 0 "S103" H 2000 1650 50  0000 C CNN
F 1 "C" H 2150 1700 50  0000 C CNN
F 2 "Housings_DIP:DIP-16_W7.62mm" H 2050 1650 50  0001 C CNN
F 3 "" H 2050 1650 50  0000 C CNN
	2    2050 1650
	1    0    0    -1  
$EndComp
$Comp
L 8INDICATOR D102
U 1 1 59711DE8
P 4050 1350
F 0 "D102" H 4150 1400 50  0000 L CNN
F 1 "S" H 4150 1300 50  0000 L CNN
F 2 "Housings_DIP:DIP-16_W7.62mm" V 4050 1350 50  0001 C CNN
F 3 "" V 4050 1350 50  0000 C CNN
	1    4050 1350
	1    0    0    -1  
$EndComp
$Comp
L CONN_01X02 J101
U 1 1 59724BC8
P 1050 1800
F 0 "J101" H 1050 1950 50  0000 C CNN
F 1 "CONN_01X02" V 1150 1800 50  0000 C CNN
F 2 "Pin_Headers:Pin_Header_Straight_1x02_Pitch2.54mm" H 1050 1800 50  0001 C CNN
F 3 "" H 1050 1800 50  0000 C CNN
	1    1050 1800
	-1   0    0    1   
$EndComp
Wire Wire Line
	1250 1750 1450 1750
Wire Wire Line
	1450 1750 1450 1150
Wire Wire Line
	1250 1850 1450 1850
Wire Wire Line
	1450 1850 1450 2400
$Comp
L RESISTOR R102
U 1 1 59724FB3
P 4150 1650
F 0 "R102" H 4100 1700 45  0000 L BNN
F 1 "RESISTOR" H 4100 1550 45  0001 L BNN
F 2 "KwanSystems:SMD_0402" H 4205 1800 20  0001 C CNN
F 3 "" H 4575 1350 60  0001 C CNN
	1    4150 1650
	0    1    1    0   
$EndComp
Wire Wire Line
	4150 1350 4150 1550
Wire Wire Line
	4150 2400 4150 1750
$Comp
L RESISTOR R101
U 1 1 59725163
P 3750 1800
F 0 "R101" H 3700 1850 45  0000 L BNN
F 1 "RESISTOR" H 3700 1700 45  0001 L BNN
F 2 "KwanSystems:SMD_0402" H 3805 1950 20  0001 C CNN
F 3 "" H 4175 1500 60  0001 C CNN
	1    3750 1800
	0    1    1    0   
$EndComp
Wire Wire Line
	3750 1450 3750 1700
Wire Wire Line
	3750 1900 3750 2400
$EndSCHEMATC
