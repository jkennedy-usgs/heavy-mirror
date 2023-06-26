program HEAVY
! ********************************************************************
! Program for calculating the gravitational attraction of mass change
!   defined by MODFLOW head change and model geometry
!
! Josh Larsen, USGS California Water Science Center
! Jeff Kennedy, USGS Arizona Water Science Center
!   jkennedy@usgs.gov
!
! Model read routines are taken from MODFLOW
!
! Forward gravity calculation from Leiraio et al., 2009,
! Calculation of the temporal gravity variation from spatially
! variable water storage change in soils and aquifers, Journal of
! Hydrology
! ********************************************************************
    use GLOBAL
    use HVYGLOBAL
    use HVYDATA
    implicit none
!    include 'openspec.inc'

    character*20 ACCESS, FORM, ACTION(2)
    data ACCESS/'STREAM'/
    data FORM/'UNFORMATTED'/
    data(ACTION(I), I=1, 2)/'READ', 'READWRITE'/

!-------ASSIGN VERSION NUMBER AND DATE
    character*40 VERSION
    character*10 HEAVYNAM
    parameter(VERSION="1.0.0 2023-02-25")
    parameter(HEAVYNAM="-HEAVY-")
    character*80 HEADNG(2)
    character*200 FNAME
    character*200 MESSAGE
    integer IBDT(8)
    character*32 arg

    integer argi
    integer ARGNO, GIDX, LAYER_TO_CALC, TIMESTEP_TO_START

    character*4 CUNIT(NIUNIT)
    data CUNIT/'BAS6', 'DIS ', 'LPF ', 'UPW ', 'PVAL', 'ZONE', 'MULT', &
    &             'HEAD', 'OC  ', 'HVY ', 90*'    '/

    integer :: I, NC, INUNIT, MAXUNIT, KP
    integer :: IGRID, NTSTEPS, GGRID
    real     :: START, FINISH
    character*32 GRIDSTR
    allocate (GARR(4))
    GARR = 0.0
    GGRID = 0
    GRIDSTR = ''
    ARGNO = 1
    LAYER_TO_CALC = 0
    TIMESTEP_TO_START = 0

    ! Parse command-line arguments - fortran 2003 code
    ! These arguments are all before the name file argument,
    ! and are optional. If none exist, name file is read
    ! with call to GETNAMFIL
    do argi = 1, command_argument_count()
        call get_command_argument(argi, arg)
        select case (arg)
            case ('-v', '--version')
                print *, VERSION
                stop

            case ('-h', '--help')
                print *, 'Heavy, a program for MODFLOW gravity modeling'
                print *
                print *, 'https://doi.org/10.5066/P9IIHXN3'
                print *, '------------------------------------------------------------------------'
                print *, 'Command-line flags:'
                print *
                print *, '    -q    Carry out QA calculations'
                print *, '    -g    Calculate gravity on a grid of evenly spaced stations.'
                print *, '            For example, "g -20" creates a 20x20 grid of stations'
                print *, '            in the x and y directions. See the  User Guide for more'
                print *, '            information.'
                print *, '    -t1   Calculate gravity change relative to the end of the first'
                print *, '            timestep instead of initial conditions.'
                print *, '    -l    Calculate gravity change only for storage change in the '
                print *, '            specified layer.'
                print *
                stop

            case ('-q', '--qc')
                print *, 'Heavy QA routine'
                CALL DEBUGFORSBERG()
                CALL DEBUGFORSBERGGRID()
                CALL DEBUGPOINTMASS()
                stop

            case ('-g', '--grid')
                call get_command_argument(argi + 1, arg)

                GRIDSTR = arg
                ! Look for '+', denotes a subset of the model grid
                GIDX = index(string=arg, substring="+")

                ! Parse arguments for subset grid
                ! Values must be single digits, i.e., 0.2, not 0.25
                if (GIDX .gt. 0) then
                    read (arg(GIDX + 1:GIDX + 3), '(f4.0)') GARR(1)
                    read (arg(GIDX + 5:GIDX + 7), '(f4.0)') GARR(2)
                    read (arg(GIDX + 9:GIDX + 11), '(f4.0)') GARR(3)
                    read (arg(GIDX + 13:GIDX + 15), '(f4.0)') GARR(4)
                    read (arg(1:GIDX - 1), '(I10)') GGRID
                    if ((GARR(2) .lt. GARR(1)) .or. (GARR(4) .lt. GARR(3))) then
                        print *, "When using the sub-grid option, the order is left-right-bottom-top. &
                        & The second element must be larger than the first, and the fourth must be larger than the third."
                        stop
                    end if
                ! GGRID = number of stations in x and y dirs
                else
                    read (arg, '(I10)') GGRID
                end if

                ARGNO = ARGNO + 2

            case ('-t1', '--timestep')
                TIMESTEP_TO_START = 1
                ARGNO = ARGNO + 1

            case ('-l', '--layer')
                call get_command_argument(argi + 1, arg)
                read (arg, '(I2)') LAYER_TO_CALC
                ARGNO = ARGNO + 2

            case default
    !            No recognized arguments, proceed to getting name file
        end select
    end do

    call CPU_TIME(START)
!
!2------WRITE BANNER TO SCREEN AND DEFINE CONSTANTS.
    write (*, 1) VERSION
1   format(/, 34x, 'HEAVY-', A, /, &
    & 25X, 'U.S. GEOLOGICAL SURVEY MODELING TOOL FOR', /, &
    & 20X, 'THE FORWARD MODELING OF MICROGRAVITY WITH MODFLOW',/)
    INUNIT = 99
!
!3------GET THE NAME OF THE NAME FILE
    call GETNAMFIL(FNAME, ARGNO)
    MAXUNIT = INUNIT
!
!4------OPEN NAME FILE
    open (unit=INUNIT, file=FNAME, status='OLD', action=ACTION(1))
    NC = index(FNAME, ' ')
    write (*, 2) ' Using NAME file: ', FNAME(1:NC)
2   format(A, A)
!
!5------Get current date and time, assign to IBDT, and write to screen
    call DATE_AND_TIME(VALUES=IBDT)
    write (*, 3) (IBDT(I), I=1, 3), (IBDT(I), I=5, 7)
3   format(1X, 'Run start date and time (yyyy/mm/dd hh:mm:ss): ', &
     & I4, '/', I2.2, '/', I2.2, 1X, I2, ':', I2.2, ':', I2.2,/)
!
!6------Allocate variables
    allocate (NCOL, NROW, NLAY, NPER, NBOTM, NCNFBD)
    allocate (ITMUNI, LENUNI, ITRSS, INBAS) !, INDIS)
    allocate (IOUT)
    allocate (IUNIT(NIUNIT))
    IUHEAD = 8
!
!7------Attempt 2 use the AR method of MODFLOW BAS
    IGRID = 1
    CALL GWF2BAS7AR(INUNIT, CUNIT, VERSION, 2, 6, 7, MAXUNIT, &
                    &    IGRID, 9, HEADNG, 5, HEAVYNAM)

    allocate (HEAD(NCOL, NROW, NLAY))
    allocate (NTOFF(NPER))

    NTOFF(1) = 0
    do KP = 1, NPER
        if (KP .GT. 1) NTOFF(KP) = NTOFF(KP - 1) + NSTP(KP - 1)
    end do

    NTSTEPS = 0
    do KP = 1, NPER
        NTSTEPS = NTSTEPS + NSTP(KP)
    end do
    allocate (NTSTP(NTSTEPS))
    allocate (SIMTIME(0:NTSTEPS))
!8------Allocate and read flow packages
!     IGRID = 1
    if (IUNIT(3) .GT. 0) call GWF2LPF7AR(IUNIT(3), IGRID)
!-------TODO: test UPW!
    if (IUNIT(4) .GT. 0) call GWF2UPW1AR(IUNIT(4), IGRID)
!9------Check the head file
!     IF(IUNIT(8).GT.0) CALL CHECKHEADFILE(IGRID,1)
!10-----Allocate arrays and read HVY package
    call GWF2HVYAR(IUNIT(10), GGRID, GRIDSTR)

!
!11----Calculate gravity
    if (IUNIT(8) .GT. 0) then
        call GWF2HVYHDLP(LAYER_TO_CALC)
    else
        MESSAGE = "BINARY HEAD FILE NOT SUPPLIED, exiting..."
        call USTOP(MESSAGE)
    end if
!
!12----GET Elapsed time
    call CPU_TIME(FINISH)
    print '("Elapsed Run Time: ", f10.3, " seconds")', FINISH - START

    write (IOUT, *) "  "
    write (IOUT, *) 'NORMAL TERMINATION OF HEAVY SIMULATION'


end program


SUBROUTINE GETNAMFIL(FNAME, ARGNO)
!******************************************************************
!    GET THE NAME OF THE NAME FILE: borrowed from mfnwt
!******************************************************************
!    SPECIFICATIONS
!------------------------------------------------------------------
    character*(*) FNAME
    character*200 COMLIN
    integer ARGNO
    logical EXISTS
!------------------------------------------------------------------
!
! Get name file from command line or user interaction.
    FNAME = ' '
    COMLIN = ' '
! *** Subroutines GETARG and GETCL are extensions to Fortran 90/95 that
! *** allow a program to retrieve command-line arguments.  To enable
! *** Modflow-2000 to read the name of a Name file from the command
! *** line, either GETARG or GETCL must be called, but not both.  As
! *** distributed, the call to GETARG is uncommented.  For compilers
! *** that support GETCL but not GETARG, comment out the call to GETARG
! *** and uncomment the call to GETCL.  The calls to both GETARG and
! *** GETCL may be commented out for compilers that do not support
! *** either extension.
!   call GETARG(ARGNO, COMLIN)
!   CALL GETCL(COMLIN)
    call get_command_argument(ARGNO, COMLIN)
    ICOL = 1
    if (COMLIN .NE. ' ') then
        FNAME = COMLIN
    else
15      write (*, *) ' Enter the name of the NAME FILE: '
        read (*, '(A)') FNAME
        call URWORD(FNAME, ICOL, ISTART, ISTOP, 0, N, R, 0, 0)
        FNAME = FNAME(ISTART:ISTOP)
        if (FNAME .EQ. ' ') goto 15
    end IF
    inquire (file=FNAME, exist=EXISTS)
    if (.NOT. EXISTS) then
        NC = INDEX(FNAME, ' ')
        FNAME(NC:NC + 3) = '.nam'
        inquire (file=FNAME, exist=EXISTS)
        if (.NOT. EXISTS) then
            write (*, 480) FNAME(1:NC - 1), FNAME(1:NC + 3)
480         format(1X, 'Can''t find name file ', A, ' or ', A)
            call USTOP(' ')
        end if
    end if

    return
end subroutine


SUBROUTINE DEBUGFORSBERG()
    use FORSBERG
    doubleprecision :: hi, hf, GR, xc, yc, slab
    doubleprecision :: sy, posX, posY, posZ, dc, dr
    doubleprecision :: G = 6.673E-11
    doubleprecision :: density = 1000.d0
    doubleprecision :: pi = 3.1415926536
    doubleprecision :: uGal = 1e8

    sy = 1.d0
    hi = -10.0
    hf = -9.0
    posX = 0
    posY = 0
    posZ = 1
    dc = 1000000
    dr = 1000000
    xc = posX
    yc = posY

    slab = 2 * pi * G * density * sy * uGal

    GR = (GRAVITY(sy, hi, hf, posX, posY, posZ, xc, yc, dr, dc))
    write (*, 505) hf-hi, GR
505 format(' Prism gravity, single 100000 x 100000-m layer, ', F5.1, '-m thick prism, g = ', F7.3, ' uGal')
    write (*, 508) slab
508 format(' Horizontal infinite slab (2*pi*G), g = ', F0.3, ' uGal')

END SUBROUTINE

SUBROUTINE DEBUGFORSBERGGRID()

    use FORSBERG
    doubleprecision :: hi, hf, GR, xc, yc
    doubleprecision :: sy, posX, posY, posZ, dc, dr

    sy = 1.d0
    hi = -10.0
    hf = -9
    posX = 0
    posY = 0
    posZ = 1
    dc = 10
    dr = 10
    xc = posX
    yc = posY
    GR = 0.0

    do i = -50000, 49900, 10
        do j = -50000,49900, 10
            xc = posX + i + dc/2
            yc = posY + j +  dc/2
            GR = GR + GRAVITY(sy, hi, hf, posX, posY, posZ, xc, yc, dc, dr)
        end do
    end do

    write (*, 510) GR
510 format(' Forward gravity, single layer 100000 x 100000-m layer composed of 10x10x1-m prisms, g = ', F0.3, ' uGal')

END SUBROUTINE
!
!
SUBROUTINE DEBUGPOINTMASS()
    use POINTMASS
    doubleprecision :: hi, hf, GR
    doubleprecision :: sy, posX, posY, posZ, dc, dr
    doubleprecision :: xc, yc, m, ref_g, r, sr
    doubleprecision :: G = 6.673E-11
    doubleprecision :: density = 1000.d0
    doubleprecision :: pi = 3.1415
    doubleprecision :: uGal = 1e8

    ! A cube 20 x 20 x 20 , 100 m deep
    sy = 1.d0
    hi = -110.d0
    hf = -90.0d0
    posX = 0.0
    posY = 0.0
    posZ = 0.0
    dc = 20.0d0
    dr = 20.0d0
    xc = posX
    yc = posY

    ! Equivalent to a sphere of
    sr = 12.407009818
    m = (4.0/3.0) * pi * (sr ** 3) * density * sy
    r = (9.5 ** 2)

    ! g = Gm/r^2
    ref_g = uGal * G * m * 1.0/r**2
    GR = PMGRAVITY(sy, hi, hf, posX, posY, posZ, xc, yc, dc, dr)
    write (*, 515)
515      format( /1X, 'Gravitational attraction of a 20x20x20-m prism, 100 m deep, vs. an equivalent mass sphere:')
    write (*, 520) GR, ref_g
520 format(' Heavy point-mass gravity, g = ', F0.6, ' uGal; Newton''s Law (G*m/r^2) = ', F0.6, ' uGal')

END SUBROUTINE
