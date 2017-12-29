KPL/FK


Ranger 7 Frames Kernel
===============================================================================

   This frame kernel contains complete set of frame definitions for the
   Ranger 7 spacecraft, its structures and science instruments. This frame
   kernel also contains name - to - NAIF ID mappings for Ranger 7 science
   instruments and s/c structures (see the last section of the file.)


Version and Date
-------------------------------------------------------------------------------

   Version 0.0 -- May 22, 2009 -- Boris Semenov, NAIF

      Initial Release with spacecraft, HGA and UHF frames based on
      temporary spacecraft ID -1007.


References
-------------------------------------------------------------------------------

   1. ``Frames Required Reading''

   2. ``Kernel Pool Required Reading''

   3. ``C-Kernel Required Reading''

   4. E-mail from Gina Signori, LMCO re. s/c, HGA and UHF frames; 05/22/09


Implementation Notes
-------------------------------------------------------------------------------

   This file is used by the SPICE system as follows: programs that make
   use of this frame kernel must ``load'' the kernel, normally during
   program initialization using the SPICELIB routine FURNSH. This file
   was created and may be updated with a text editor or word processor.


MAVEN Frames
-------------------------------------------------------------------------------

   The following MAVEN frames are defined in this kernel file:

           Name                  Relative to           Type       NAIF ID
      ======================  ===================  ============   =======

   Spacecraft frame:
   -----------------
      RANGER7_SPACECRAFT        rel.to ECI_TOD      CK             -1007000


MAVEN Frames Hierarchy
-------------------------------------------------------------------------------

   The diagram below shows MAVEN frames hierarchy:


                               "ECI_TOD" INERTIAL
        +------------------------------------------------------------+
        |                              |                             |
        | <--pck                       |<-fixed                      | <--pck
        |                              |                             |
                                       |
                                       |
                                       |      
                                       |      
                                       |      
                                       |      
                                       |<--ck 
                                       |           
                               "RANGER7_SPACECRAFT"   
                               --------------------

   (*) BFR -- body-fixed rotating frame

Spacecraft Bus Frame
-------------------------------------------------------------------------------
 
   The spacecraft frame is defined by the s/c design as follows [from
   TBD]:

      -  Z axis is parallel to the nominal HGA boresight;
 
      -  X axis is parallel to the nominal UHF antenna boresight;

      -  Y axis completes the right hand frame;

      -  the origin of the frame is centered on the launch vehicle
         separation plane.

   These diagrams illustrate the s/c frame:

      TBD
 
   Since the S/C bus attitude is provided by a C kernel (see [3] for
   more information), this frame is defined as a CK-based frame.

   \begindata

      FRAME_RANGER7_SPACECRAFT       = -1007000
      FRAME_-1007000_NAME            = 'RANGER7_SPACECRAFT'
      FRAME_-1007000_CLASS           = 3
      FRAME_-1007000_CLASS_ID        = -1007000
      FRAME_-1007000_CENTER          = -1007
      CK_-1007000_SCLK               = -1007
      CK_-1007000_SPK                = -1007

