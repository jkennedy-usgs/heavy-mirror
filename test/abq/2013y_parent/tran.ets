# Evapotranspiration Segments package for MF2000
# Mult arrays are the proportion of each cell area with a vegitation zone where ET from the water
#   table can occur.  The Parval is then the annual or seasonal max ET rate.
#   Annual rate = annual max ET/365.25 days; Irr season rate = annual max ET/230 days.
#   ET surface elev is land surface elev => max rate if water table is simulated to be at
#   land surface (LS), or higher.  Extinction depth is 30 feet (~10 meters) below land surface,
#   Divided into three segments:  where, at 9 ft below LS (30% of ETSX) PETM is 40%, and at
#   16 feet below LS (PETM (14 ft below ET surface)
         1        50        11         3    NETSOP IETSCB NPETS NETSEG
ETrg35           ETS 0.137E-01         1    PARNAM PARTYP Parval NCLU #1935 Rio Grande max annual ET rate
ETrg35M          ALL         1              Mltarr Zonarr IZ #ALL => no zone array
ETrg55           ETS 0.137E-01         1    PARNAM PARTYP Parval NCLU #1955 Rio Grande max annual ET rate
ETrg55M          ALL         1              Mltarr Zonarr IZ #ALL => no zone array
ETrg75           ETS 0.137E-01         1    PARNAM PARTYP Parval NCLU #1975 Rio Grande max annual ET rate
ETrg75M          ALL         1              Mltarr Zonarr IZ #ALL => no zone array
ETrg92           ETS 0.137E-01         1    PARNAM PARTYP Parval NCLU #1992 Rio Grande max annual ET rate
ETrg92M          ALL         1              Mltarr Zonarr IZ #ALL => no zone array
ETrg92i          ETS 0.217E-01         1    PARNAM PARTYP Parval NCLU #1992 Rio Grande max irr season ET rate
ETrg92M          ALL         1              Mltarr Zonarr IZ #ALL => no zone array
ETjem55          ETS 0.137E-01         1    PARNAM PARTYP Parval NCLU #1955 Jemez R max annual ET rate
ETjem55M         ALL         1              Mltarr Zonarr IZ #ALL => no zone array
ETjem75          ETS 0.137E-01         1    PARNAM PARTYP Parval NCLU #1975 Jemez R max annual ET rate
ETjem75M         ALL         1              Mltarr Zonarr IZ #ALL => no zone array
ETjem75i         ETS 0.217E-01         1    PARNAM PARTYP Parval NCLU #1975 Jemez R max irr season ET rate
ETjem75M         ALL         1              Mltarr Zonarr IZ #ALL => no zone array
ETjemres         ETS 0.137E-01         1    PARNAM PARTYP Parval NCLU #1975 Jemez R ET rate masked by Jemez Res.
ETjem75M         JEMETzn     1 2            Mltarr Zonarr IZ #zones 1,2 et; zones 3,4 generall covered by Reservoir
ETjemresi        ETS 0.217E-01         1    PARNAM PARTYP Parval NCLU #1975 Jemez R ET rate masked by Jemez Res.
ETjem75M         JEMETzn     1 2            Mltarr Zonarr IZ #same as above except for irr season
ETnone           ETS      0.00         1    PARNAM PARTYP Parval NCLU #used to turn off ET in winter
NONE             ALL         1              #all values null to turn of ET in winter 
         1         2         1         0         1   INETSS INETSR INETSX INIETS INSGDF
OPEN/CLOSE ../externalfiles/arrays/SURF_elev.arr 1. (10f8.1) -1 # ET surface elev array; 0.00 outside active model
ETrg35  -1                                  # ET along Rio Grande for steady state
ETjem55 -1                                  # ET along Jemez R for steady state (1935 info not available)
CONSTANT         30.                        ETSX
CONSTANT        0.30                        PXDP   1st segment from 0 to 9 ft below ET surface 
CONSTANT        0.40                        PETM        40% max et at 9 ft below ET surface
CONSTANT        0.53                        PXDP   2nd segment from 9 to 16 ft below ET surface,
CONSTANT        0.15                        PETM        15% max et at 16 ft below ET surface (last segment 16-30 ft below LS)
        -1        -1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- use previous for 1900-04
        -1        -1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- use previous for 1905-09
        -1        -1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- use previous for 1910-14
        -1        -1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- use previous for 1915-19
        -1        -1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- use previous for 1920-24
        -1        -1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- use previous for 1925-29
        -1        -1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- use previous for 1930-34
        -1        -1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- use previous for 1935-39
        -1        -1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- use previous for 1940-44
        -1         2        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new params but use previous for all others
ETrg55  -1                                  # ET along Rio Grande for 1945-49 (switch to 1955 data)
ETjem55 -1                                  # ET along Jemez R for 1945-49
        -1        -1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- use previous for 1950-54
        -1        -1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- use previous for 1955-59
        -1        -1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- use previous for 1960-64
        -1         2        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new params but use previous for all others
ETrg75  -1                                  # ET along Rio Grande for 1965-69 (switch to 1975 data)
ETjem75 -1                                  # ET along Jemez R for 1965-69 (switch to 1975 data)
        -1        -1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- use previous for 1970-74
        -1        -1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- use previous for 1975
        -1        -1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- use previous for 1976
        -1        -1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- use previous for 1977
        -1        -1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- use previous for 1978
        -1         2        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new params but use previous for all others
ETrg75  -1                                  # ET along Rio Grande for 1979 
ETjemres -1                                   # ET along Jemez R for 1979 masked by Jemez Reservoir
        -1        -1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- use previous for 1980
        -1        -1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- use previous for 1981
        -1        -1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- use previous for 1982
        -1        -1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- use previous for 1983
        -1         2        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new params but use previous for all others
ETrg92   -1                                 # ET along Rio Grande for 1984 (switch to 1992 data)
ETjemres -1                                 # ET along Jemez R for 1984 masked by Jemez Reservoir (1992 info not available)
        -1        -1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- use previous for 1985
        -1        -1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- use previous for 1986
        -1        -1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- use previous for 1987
        -1        -1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- use previous for 1988
        -1        -1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- use previous for 1989
        -1         1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new param but use previous for all others
ETnone   -1                                 # no ET in winter 01/90 - 03/90
        -1         2        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new params but use previous for all others
ETrg92i  -1                                 # ET along Rio Grande for 1990 irrigation season
ETjemresi -1                                 # ET along Jemez R for 1990 irrigation season (03/90 - 10/90)
        -1         1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new param but use previous for all others
ETnone    -1                                # no ET in winter 11/90 - 03/91
        -1         2        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new params but use previous for all others
ETrg92i   -1                                # ET along Rio Grande for 1991 irrigation season
ETjemresi -1                                # ET along Jemez R for 1991 irrigation season (03/91 - 10/91)
        -1         1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new param but use previous for all others
ETnone    -1                                # no ET in winter 11/91 - 03/92
        -1         2        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new params but use previous for all others
ETrg92i   -1                                # ET along Rio Grande for 1992 irrigation season
ETjemresi -1                                # ET along Jemez R for 1992 irrigation season (03/92 - 10/92)
        -1         1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new param but use previous for all others
ETnone    -1                                # no ET in winter 11/92 - 03/93
        -1         2        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new params but use previous for all others
ETrg92i   -1                                # ET along Rio Grande for 1993 irrigation season
ETjemresi -1                                # ET along Jemez R for 1993 irrigation season (03/93 - 10/93)
        -1         1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new param but use previous for all others
ETnone    -1                                # no ET in winter 11/93 - 03/94
        -1         2        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new params but use previous for all others
ETrg92i   -1                                # ET along Rio Grande for 1994 irrigation season
ETjemresi -1                                # ET along Jemez R for 1994 irrigation season (03/94 - 10/94)
        -1         1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new param but use previous for all others
ETnone    -1                                # no ET in winter 11/94 - 03/95
        -1         2        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new params but use previous for all others
ETrg92i   -1                                # ET along Rio Grande for 1995 irrigation season
ETjemresi -1                                # ET along Jemez R for 1995 irrigation season (03/95 - 10/95)
        -1         1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new param but use previous for all others
ETnone    -1                                # no ET in winter 11/95 - 03/96
        -1         2        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new params but use previous for all others
ETrg92i   -1                                # ET along Rio Grande for 1996 irrigation season
ETjemresi -1                                # ET along Jemez R for 1996 irrigation season (03/96 - 10/96)
        -1         1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new param but use previous for all others
ETnone    -1                                # no ET in winter 11/96 - 03/97
        -1         2        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new params but use previous for all others
ETrg92i   -1                                # ET along Rio Grande for 1997 irrigation season
ETjemresi -1                                # ET along Jemez R for 1997 irrigation season (03/97 - 10/97)
        -1         1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new param but use previous for all others
ETnone    -1                                # no ET in winter 11/97 - 03/98
        -1         2        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new params but use previous for all others
ETrg92i   -1                                # ET along Rio Grande for 1998 irrigation season
ETjemresi -1                                # ET along Jemez R for 1998 irrigation season (03/98 - 10/98)
        -1         1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new param but use previous for all others
ETnone    -1                                # no ET in winter 11/98 - 03/99
        -1         2        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new params but use previous for all others
ETrg92i   -1                                # ET along Rio Grande for 1999 irrigation season
ETjemresi -1                                # ET along Jemez R for 1999 irrigation season (03/99 - 10/99)
        -1         1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new param but use previous for all others
ETnone    -1                                # no ET in winter 11/99 - 03/00   
        -1         2        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new params but use previous for all others
ETrg92i   -1                                # ET along Rio Grande for 2000 irrigation season
ETjemresi -1                                # ET along Jemez R for 2000 irrigation season (03/00 - 10/00)
        -1         1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new param but use previous for all others
ETnone    -1                                # no ET in winter 11/00 - 03/01    
        -1         2        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new params but use previous for all others
ETrg92i   -1                                # ET along Rio Grande for 2001 irrigation season
ETjemresi -1                                # ET along Jemez R for 2001 irrigation season (03/01 - 10/01)
        -1         1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new param but use previous for all others
ETnone    -1                                # no ET in winter 11/01 - 03/02    
        -1         2        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new params but use previous for all others
ETrg92i   -1                                # ET along Rio Grande for 2002 irrigation season
ETjem75i  -1                               # ET along Jemez R for 2002 irrigation season (03/02 - 10/02) jem res drained by oct 27, 2001
        -1         1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new param but use previous for all others
ETnone    -1                                # no ET in winter 11/02 - 03/03    
        -1         2        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new params but use previous for all others
ETrg92i   -1                                # ET along Rio Grande for 2003 irrigation season
ETjem75i  -1                               # ET along Jemez R for 2003 irrigation season (03/03 - 10/03) jem res drained by oct 27, 2001
        -1         1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new param but use previous for all others
ETnone    -1                                # no ET in winter 11/03 - 03/04 
        -1         2        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new params but use previous for all others
ETrg92i   -1                             
ETjem75i  -1                             
        -1         1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new param but use previous for all others
ETnone    -1                             # no ET in winter 11/04 - 03/05
        -1         2        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new params but use previous for all others
ETrg92i   -1                             
ETjem75i  -1                             
        -1         1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new param but use previous for all others
ETnone    -1                             # no ET in winter 11/05 - 03/06  
        -1         2        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new params but use previous for all others
ETrg92i   -1                             
ETjem75i  -1                             
        -1         1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new param but use previous for all others
ETnone    -1                             # no ET in winter 11/06 - 03/07
        -1         2        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new params but use previous for all others
ETrg92i   -1                             
ETjem75i  -1                             
        -1         1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new param but use previous for all others
ETnone    -1                             # no ET in winter 11/07 - 03/08
        -1         2        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new params but use previous for all others
ETrg92i   -1                             
ETjem75i  -1                             
        -1         1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new param but use previous for all others
ETnone    -1                             # no ET in winter 11/08 - 03/09
        -1         2        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new params but use previous for all others StressPeriod 71
ETrg92i   -1                             
ETjem75i  -1                             
        -1         1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new param but use previous for all others  StressPeriod 72
ETnone    -1                             # no ET in winter 11/09 - 03/10
        -1         2        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new params but use previous for all others  StressPeriod 73
ETrg92i   -1                             
ETjem75i  -1                             
        -1         1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new param but use previous for all others  StressPeriod 74
ETnone    -1                             # no ET in winter 11/10 - 03/11
        -1         2        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new params but use previous for all others  StressPeriod 75
ETrg92i   -1                             
ETjem75i  -1                             
        -1         1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new param but use previous for all others  StressPeriod 76
ETnone    -1                             # no ET in winter 11/11 - 03/12
        -1         2        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new params but use previous for all others  StressPeriod 77
ETrg92i   -1                             
ETjem75i  -1                             
        -1         1        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new param but use previous for all others  StressPeriod 78
ETnone    -1                             # no ET in winter 11/12 - 03/13
        -1         2        -1        -1        -1    # INETSS INETSR INETSX INIETS INSGDF -- read new params but use previous for all others  StressPeriod 79
ETrg92i   -1                             
ETjem75i  -1                             
