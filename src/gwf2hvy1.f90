module HVYGLOBAL
   integer, save, pointer ::NFHVY, IHVYUN
   real, save, pointer ::CFACT, HVYNOH
   doubleprecision, save, dimension(:, :, :), pointer ::HEAD
   real, save, dimension(:, :, :), pointer ::JBUFF
   integer, save, dimension(:, :, :), pointer ::IBUFF
   doubleprecision, save, dimension(:, :, :), pointer ::HVYBND
   character(len=20), save, dimension(:), pointer ::PRLBL
   character(len=20), save, dimension(:), pointer ::PMLBL
   integer, save, dimension(:), pointer ::KHVY
   integer, save, dimension(:), pointer ::KPM
   doubleprecision, save, dimension(:, :), pointer ::GRIDARRAY
   doubleprecision, save, dimension(:, :), pointer ::XYZHVY
   doubleprecision, save, dimension(:, :), pointer ::XYZPM
   doubleprecision, save, dimension(:, :, :), pointer ::SYHVY
   doubleprecision, save, dimension(:, :, :), pointer ::SSHVY
   doubleprecision, save, dimension(:), pointer ::HVYDELC
   doubleprecision, save, dimension(:), pointer ::HVYDELR
   doubleprecision, save, dimension(:), pointer ::XC
   doubleprecision, save, dimension(:), pointer ::YC
   doubleprecision, save, dimension(:), pointer ::XARR
   doubleprecision, save, dimension(:), pointer ::YARR
   doubleprecision, save, dimension(:, :, :), pointer ::HVYSTRT
   real, save, dimension(:), pointer :: GARR

        TYPE HVYGLOBALTYPE
        integer, pointer ::NFHVY,IHVYUN,IPMUN,NPMHVY
        real,    pointer ::CFACT,HVYNOH
                doubleprecision,    dimension(:),     pointer ::IPMGRAV
                doubleprecision,    dimension(:,:,:), pointer ::HEAD
                real,    dimension(:,:,:), pointer ::JBUFF
        integer, dimension(:,:,:), pointer ::IBUFF
        real,    dimension(:,:,:), pointer ::HVYBND
        character(LEN=20), dimension(:), pointer ::PRLBL
        character(LEN=20), dimension(:), pointer ::PMLBL
        integer, dimension(:),     pointer ::KHVY
        integer, dimension(:),     pointer ::KPM
        doubleprecision,    dimension(:,:),   pointer ::XYZHVY
        doubleprecision,    dimension(:,:),   pointer ::XYZPM
        doubleprecision,    dimension(:,:,:), pointer ::SYHVY
        doubleprecision,    dimension(:),     pointer ::HVYDELC
        doubleprecision,    dimension(:),     pointer ::HVYDELR
        doubleprecision,    dimension(:),     pointer ::XC
        doubleprecision,    dimension(:),     pointer ::YC
        doubleprecision,    dimension(:),     pointer ::XARR
        doubleprecision,    dimension(:),     pointer ::YARR
        doubleprecision,    dimension(:,:,:), pointer ::HVYSTRT
        end TYPE HVYGLOBALTYPE
        TYPE(HVYGLOBALTYPE), save ::HVYGLOBALDAT(10)
end module HVYGLOBAL
!
!
module HVYdata
   integer, save ::IOLIST, IUHEAD
   double precision, save, pointer, dimension(:) ::SIMTIME
   integer, save, pointer, dimension(:) ::NTOFF
   integer, save, pointer, dimension(:) ::NTSTP
end module
!
!
subroutine GWF2HVYAR(IU, GGRID, GRIDSTR)
!******************************************************
! Allocate and Read module for HEAVY microgravity
! module. This module reads input from the HVY file
! and stores the appropriate data for use in gravity
! calculations
!******************************************************
   use GLOBAL, ONLY: LENUNI, IOUT

   use HVYGLOBAL, ONLY: XC, YC, CFACT, NFHVY, IHVYUN, HVYNOH

   character*200 LINE
   character*200 MESSAGE
   character*24 ANAME(1)
   integer :: LLOC, ISTRT, ISTOP, GGRID, IU
   real :: RR
   logical HVYOK
   character*32 GRIDSTR
   HVYOK = .FALSE.

   data ANAME(1)/'          SPECIFIC YIELD'/
   allocate (CFACT, NFHVY, IHVYUN, HVYNOH)
!2----Identify the package
   write (IOUT, 1) IU
1  format(1X, /1X, 'HEAVY MICROGRAVITY FORWARD MODEL ', &
   & 'VERSION 1.0.0, 2/25/2023', /, 9X, 'INPUT read FROM UNIT ', I4)
!
!2----set conversion factor based on LENUNI
   if (LENUNI .eq. 0) then
      CFACT = 1.0
   elseif (LENUNI .eq. 1) then
      CFACT = 0.3048
   elseif (LENUNI .eq. 2) then
      CFACT = 1.0
   elseif (LENUNI .eq. 3) then
      CFACT = 100.0
   end if
   call URDCOM(IU, IOUT, LINE)
   LLOC = 1
   call URWORD(LINE, LLOC, ISTRT, ISTOP, 2, NFHVY, RR, IOUT, IU)
   call URWORD(LINE, LLOC, ISTRT, ISTOP, 2, IHVYUN, RR, IOUT, IU)

   write (IOUT, 2) IHVYUN
2  format(1X, 'MICROGRAVITY OUTPUT WILL BE SAVED ON UNIT ', I4)

!4a---Base Array Allocation
   call HVYALARR()

! If -g flag was specified, generate station locations:
   if (GGRID .gt. 0) then
      write(IOUT, 4) GRIDSTR
4     format(1X, "GENERATING STATIONS ON A GRID. COMMAND STRING: ", A)
      call HVYGRID(GGRID)
      HVYOK = .TRUE.
! Otherwise, read in station locations from .hvy file:
   else
!---Read comments and item 1

      write (IOUT, 5) NFHVY
5     format(1X, 'NUMBER OF OBSERVATION LOCATIONS=', I4)
      call HVYREADOL(IU)
      HVYOK = .TRUE.
   end if

   if (.NOT. HVYOK) then
      write (IOUT, 7)
7     format(1X, 'NO GRAVITY OBSERVATIONS FOUND')
      MESSAGE = 'NO GRAVITY OBSERVATIONS FOUND'
      call USTOP(MESSAGE)
   end if


!
!5--IMPORTANT! Keep XC, YC conversion at the end of this
!5--Subroutine, do not move!!!
   XC = XC*CFACT
   YC = YC*CFACT

end subroutine

subroutine HVYGRID(GGRID)
! ********************************************************************
! If command-line flag -g is used, create a grid of stations at which
!   to calculate gravity.
! ********************************************************************
!   use LINSPACE
   use GLOBAL, ONLY: IOUT, BOTM, IBOUND

   use HVYGLOBAL, ONLY: XC, YC, NFHVY, &
   &                   KHVY, XYZHVY, PRLBL, &
   &                   XARR, YARR, GARR, CFACT

   character*20 ctrcomb
   integer :: GN, GNN, KLAY, GGRID, CTR, R_idx, C_idx
   doubleprecision XRANGE, YRANGE
   doubleprecision, dimension(GGRID) :: xarray
   doubleprecision, dimension(GGRID) :: yarray
   NFHVY = GGRID*GGRID
   allocate (KHVY(NFHVY))
   allocate (XYZHVY(3, NFHVY))
   allocate (PRLBL(NFHVY))
   allocate (XARR(GGRID))
   allocate (YARR(GGRID))
   write (*, 520) GGRID, GGRID
520 format(' CREATING GRAVITY-STATION GRID, ', I0, 'x', I0)

   ZLOC = 50.0
   KLAY = 1
   CTR = 1
   ! if no subgrid coordinates are specified, all GARR elements will be zero
   if (sum(GARR) .gt. 0.2) then
   XRANGE = XC(size(XC)) - XC(1)
   YRANGE = YC(1) - YC(size(YC))  ! YC is decreasing - always the case?
   call LINSPACE(XC(1) + XRANGE*GARR(1), XC(1) + XRANGE*GARR(2), GGRID, xarray)
   call LINSPACE(YC(size(YC)) + YRANGE*GARR(3), YC(size(YC)) + YRANGE*GARR(4), GGRID, yarray)
   else
      call LINSPACE(XC(1), XC(size(XC)), GGRID, xarray)
      call LINSPACE(YC(1), YC(size(YC)), GGRID, yarray)
   end if
   do GN = 1, GGRID
      do GNN = 1, GGRID
         write (ctrcomb, 400) GN, GNN
400      format('G_', I0, '_', I0)
         PRLBL(CTR) = ctrcomb
         KHVY(CTR) = 1

         XYZHVY(1, CTR) = xarray(GN)*CFACT
         XYZHVY(2, CTR) = yarray(GNN)*CFACT

         ! Figure out elevations
         R_idx = minloc(abs(XC - xarray(GN)), DIM=1)
         C_idx = minloc(abs(YC - yarray(GNN)), DIM=1)
         if (IBOUND(R_idx, C_idx, 1) .ne. 0) then
            ZLOC = BOTM(R_idx, C_idx, 0) * CFACT
         else
            ZLOC = 0.0
         end if
         !print *, ZLOC
         XYZHVY(3, CTR) = ZLOC
         CTR = CTR + 1
      end do
   end do

22 format()
   write (IOUT, 23)
23 format(1X, /3X, 'HEAVY-GENERATED GRAVITY OBSERVATION LOCATIONS:', /1X, &
    & 'PRLBL', 19X, 'ZLOC', 10X, 'XLOC', 10X, 'YLOC', /1X, 56('-'))
   do N = 1, NFHVY
      write (IOUT, 24) PRLBL(N), XYZHVY(3, N), &
      & XYZHVY(1, N), XYZHVY(2, N)
24    format(1X, A14, 3F14.3)
   end do
end subroutine

subroutine HVYREADOL(IU)
   !*********************************************
   ! Read observation locations from .hvy file
   !*********************************************
   use GLOBAL, ONLY: IOUT

   use HVYGLOBAL, ONLY: NFHVY, IHVYUN, &
   &                   KHVY, XYZHVY, PRLBL

   character*200 LINE
   integer :: N, LLOC, NN, ISTRT, ISTOP
   real :: RR, ZLOC, XLOC, YLOC

   allocate (KHVY(NFHVY))
   allocate (XYZHVY(3, NFHVY))
   allocate (PRLBL(NFHVY))
!
!5----Read data set 2
   do N = 1, NFHVY
20    read (IU, '(A)') LINE
      if (LINE .EQ. ' ') GO TO 20
      LLOC = 1
      call URWORD(LINE, LLOC, ISTRT, ISTOP, 3, NN, ZLOC, IOUT, IU)
      call URWORD(LINE, LLOC, ISTRT, ISTOP, 3, NN, XLOC, IOUT, IU)
      call URWORD(LINE, LLOC, ISTRT, ISTOP, 3, NN, YLOC, IOUT, IU)
      call URWORD(LINE, LLOC, ISTRT, ISTOP, 1, NN, RR, IOUT, IU)
      PRLBL(N) = LINE(ISTRT:ISTOP)
      XYZHVY(1, N) = XLOC
      XYZHVY(2, N) = YLOC
      XYZHVY(3, N) = ZLOC

   end do
!
!7----print a table of values used in the gravity calculation
   write (IOUT, 23)
23 format(1X, /3X, 'GRAVITY OBSERVATION LOCATIONS:', /1X, &
    & 'PRLBL', 19X, 'ZLOC', 10X, 'XLOC', 10X, 'YLOC', /1X, 56('-'))
   do N = 1, NFHVY
      write (IOUT, 24) PRLBL(N), XYZHVY(3, N), &
      & XYZHVY(1, N), XYZHVY(2, N)
24    format(1X, A14, 3F14.3)
   end do

!2021-10-28: assuming gravity station coordinates are in meters!
!
!10--write header to output file
   write (IHVYUN, 25)
25 format(1X, 'PRLBL', 15X, 'KPER', 4X, 'KSTP', 9X, 'TOTIM', 8X, 'd-uGAL')

end subroutine
!
!
subroutine HVYALARR()
!*************************************************
! Create and allocate arrays that are used in both
! point mass and prism formulations of the Gravity
! forward model
!*************************************************
   use GLOBAL, ONLY: IBOUND, NCOL, NROW, NLAY, DELR, DELC, &
                     IUNIT
   use HVYGLOBAL, ONLY: XC, YC, SYHVY, SSHVY, HVYBND, HVYDELC, HVYDELR, &
   &                   HVYSTRT,  CFACT, NFHVY

   use GWFLPFmodule, ONLY: SC1, SC2, LAYTYP
   use GWFUPWmodule, ONLY: SC1UPW, SC2UPW, LAYTYPUPW

   character*200 MESSAGE
   integer :: N, J, I, K, S
   real, dimension(NCOL, NROW, NLAY) ::SCSY

   real, dimension(NCOL, NROW, NLAY) ::SCSS
   real, dimension(NLAY) ::LTYP
   real :: YLEN

   allocate (XC(NCOL))
   allocate (YC(NROW))
   allocate (HVYDELC(NROW))
   allocate (HVYDELR(NCOL))
   allocate (HVYSTRT(NCOL, NROW, NLAY))

   allocate (HVYBND(NCOL, NROW, NLAY))
   allocate (SYHVY(NCOL, NROW, NLAY))
   allocate (SSHVY(NCOL, NROW, NLAY))

!1----calculate cell centers for each direction
   do N = 1, NCOL
      if (N .eq. 1) then
         XC(N) = DELR(N)/2.0
      else
         XC(N) = XC(N - 1) + DELR(N - 1)/2.0 + DELR(N)/2.0
      end if
   end do

   YLEN = SUM(DELC)
   do N = 1, NROW
      if (N .eq. 1) then
         YC(N) = YLEN - DELC(N)/2.0
      else
         YC(N) = YC(N - 1) - DELC(N - 1)/2.0 - DELC(N)/2.0
      end if
   end do
!
!2--Get storage capacity from the appropriate flow pkg
! LPF
   if (IUNIT(3) .GT. 0) then
      if (LAYTYP(1) .eq. 0) then
         MESSAGE = "LAYER 1 IS CONFINED, CANNOT CALCULATE GRAVITY"
         call USTOP(MESSAGE)
      end if
      S = size(SC2,3)
      SCSY(:,:,1:S) = SC2
      SCSS = SC1
      LTYP = LAYTYP
! UPW
   else if (IUNIT(4) .GT. 0) then
      if (LAYTYPUPW(1) .eq. 0) then
         MESSAGE = "LAYER 1 IS CONFINED, CANNOT CALCULATE GRAVITY"
         call USTOP(MESSAGE)
      end if
      S = size(SC2,3)
      SCSS = SC1UPW
      SCSY(:,:,1:S) = SC2UPW
      LTYP = LAYTYPUPW
   else
      MESSAGE = 'NO COMPATIBLE FLOW PACKAGE FOUND'
      call USTOP(MESSAGE)
   end if
!
!---Prep arrays for later gravity calculations

   do K = 1, NLAY
      do J = 1, NCOL
         do I = 1, NROW
            if (LTYP(K) .eq. 0) then
               SYHVY(J, I, K) = 0
               SSHVY(J, I, K) = 0
            else
               ! SC = Storage capacity = (SY or SS) * delr * delc
               SYHVY(J, I, K) = SCSY(J, I, K)/(DELC(I)*DELR(J))
               SSHVY(J, I, K) = SCSS(J, I, K)/(DELC(I)*DELR(J))
            end if

            if (IBOUND(J, I, K) .gt. 0) then
               HVYBND(J, I, K) = 1
            else
               HVYBND(J, I, K) = 0
            end if
         end do
      end do
   end do

   HVYDELR = DELR*CFACT
   HVYDELC = DELC*CFACT

end subroutine
!
!
subroutine GWF2HVYHDLP(LAYER_TO_CALC)
! *******************************************
! Loop to read in the head file and call
! gravity subroutines for each observation.
! *******************************************
   use GLOBAL, ONLY: NPER, NSTP
   use HVYGLOBAL, ONLY: IHVYUN, HEAD, HVYSTRT

   integer ::KSTP, KPER, KS, KP, NN, PP, LAYER_TO_CALC
   double precision ::PERTIM, TOTIM
   LOGICAL ISOK
   NN = 1
   PP = 1
   if (LAYER_TO_CALC.eq.0) then
      print *, 'CALCULATING GRAVITY FOR ALL LAYERS'
   else
      print *, 'CALCULATING GRAVITY FOR LAYER ', LAYER_TO_CALC
   end if
   call readTIMESTEPHEAD(1, 1, KS, KP, PERTIM, TOTIM, ISOK)
   HVYSTRT = HEAD

   do KPER = 1, NPER
      do KSTP = 1, NSTP(KPER)
         ! Puts current timestep head into HEAD var
         call readTIMESTEPHEAD(KSTP, KPER, KS, KP, PERTIM, TOTIM, ISOK)
         if (ISOK) then
            if (IHVYUN .gt. 0) then
               call HVYPRISMLP(KP, KS, TOTIM, LAYER_TO_CALC)
               NN = NN + 1
            end if
         end if

      end do
   end do

end subroutine
!
!
subroutine HVYPRISMLP(KP, KS, TOTIM, LAYER_TO_CALC)
!***************************************************
! Loop which sets up and calls the forward gravity
! calculation
!***************************************************
   use GLOBAL, ONLY: NLAY, NROW, NCOL, BOTM
   use GWFBASmodule, ONLY: HDRY, HNOFLO
   use HVYGLOBAL, ONLY: NFHVY, HEAD, CFACT, &
       &                  PRLBL, IHVYUN, &
       &                  SYHVY, SSHVY, HVYSTRT, XYZHVY, XC, YC, HVYDELC, &
       &                  HVYDELR, HVYBND, HVYNOH

   use FORSBERG
   use POINTMASS
   character*10 :: formula
   integer ::KP, KS, N, J, I, L, LAYER_TO_CALC
   doubleprecision :: IBND, dh, TOTIM
   doubleprecision :: HF
   doubleprecision, dimension(NLAY) :: GRAV
   doubleprecision, dimension(NCOL, NROW) ::GVTEMP
   doubleprecision :: x, y, z, r2, dr2, prism_top, prism_bottom, d_rho
   doubleprecision :: g
!   DEBUG: write calculated g for every cell
!   OPEN(UNIT=66,FILE='prism.csv',ACTION='WRITE',position='append')
   do N = 1, NFHVY

!      if ( abs(XYZHVY(3, N)) < 0.01 ) cycle  ! Skip stations at z = 0

      do L = 1,NLAY
         if (LAYER_TO_CALC.ne.0) then
            if (L.ne.LAYER_TO_CALC) cycle
         end if
         do I = 1, NROW
            do J = 1, NCOL
                GVTEMP(J, I) = 0.d0

                HF = HEAD(J, I, L)
!                if (BOTM(J, I, LBOTM(L)-1).eq.BOTM(J, I, LBOTM(L))) cycle
                if (BOTM(J, I, L-1).eq.BOTM(J, I, L)) cycle
                if (HF.lt.BOTM(J, I, L)) cycle
                if (HVYSTRT(J, I, L).eq.HF) cycle
                if (HF .eq. HDRY) cycle
                if (HF .eq. HNOFLO) cycle
                if (HF .eq. HVYNOH) cycle

                HF = HF*CFACT
                IBND = HVYBND(J, I, L)

                if (IBND .ne. 0) then
                    dh = HF - HVYSTRT(J, I, L)*CFACT
                   ! Head higher than top of prism = confined

                   if ( HF .gt. BOTM(J, I, L-1) * CFACT) then
                      prism_bottom = BOTM(J, I, L) * CFACT
                      prism_top = BOTM(J, I, L-1) * CFACT

                      d_rho = SSHVY(J, I, L) * dh
                   else
                      HS = HVYSTRT(J, I, L)*CFACT
                      prism_bottom = min(HS, HF)
                      prism_top = max(HS, HF)

                      if (dh.gt.0) then

                        d_rho = SYHVY(J, I, L)

                      else
                        d_rho = (-1 * SYHVY(J, I, L) )
                      end if
                   end if
    !
                   x = XYZHVY(1, N) - XC(J)
                   y = XYZHVY(2, N) - YC(I)
                   z = XYZHVY(3, N) - (prism_top+prism_bottom)/2

                   r2 = x**2 + y**2 + z**2
                   formula="N"
                   dr2 = HVYDELR(J)**2 + HVYDELC(I)**2 + (prism_top-prism_bottom)**2
                   if ( r2/dr2 > 81) then
                      formula = "PM"
                      g = PMGRAVITY(d_rho, prism_bottom, prism_top, XYZHVY(1, N), &
                      &                 XYZHVY(2, N), XYZHVY(3, N), XC(J), YC(I), &
                      &                 HVYDELR(J), HVYDELC(I))
                      GVTEMP(J, I) = g
                   else
                      formula = "F"
                      g = GRAVITY(d_rho, prism_bottom, prism_top, XYZHVY(1, N), &
                      &                 XYZHVY(2, N), XYZHVY(3, N), XC(J), YC(I), &
                      &                 HVYDELR(J), HVYDELC(I))
                      GVTEMP(J, I) = g
                   end if

                   ! DEBUG
!                   if (abs(g) .gt. 0.001) then
!                     write(66,101) L, TOTIM, g, XC(J), YC(I), prism_bottom, &
!                     &       prism_top, HF, HVYSTRT(J, I, L)*CFACT, &
!                     &       BOTM(J, I, L) * CFACT, BOTM(J, I, L-1) * CFACT, &
!                     &       d_rho, XYZHVY(1, N), XYZHVY(2, N), XYZHVY(3, N), formula
!101                  format(I1,F15.5,E15.5,F15.5,F15.5,F15.5,F15.5,F15.5,F15.5,F15.5,F15.5,E15.5,F15.5,F15.5,F15.5,A5)
!                   end if
                end if
            end do
         end do
         GRAV(L) = SUM(GVTEMP)
      end do

!-----------Write output to heavy output file
    if (LAYER_TO_CALC.eq.0) then
      write (IHVYUN, 2) PRLBL(N), KP, KS, TOTIM, GRAV(:), SUM(GRAV)
2     format(1X, A16, 2I8, F20.3, *(F14.3))
    else
      write (IHVYUN, 2) PRLBL(N), KP, KS, TOTIM, GRAV(LAYER_TO_CALC)
    end if

   end do
end subroutine

subroutine LINSPACE(from, to, n, v_out)

      doubleprecision, intent(in) :: from, to
      doubleprecision, intent(inout), dimension(n) :: v_out
      doubleprecision :: steplen
      integer, intent(in) :: n
      integer :: i

      if (n > 1) then
        steplen = dble(to - from)/dble(n - 1)
        do i = 1, n
            v_out(i) = from + steplen*dble(i-1)
        end do
      elseif (n == 1) then
        v_out(1) = (from + to)/2
      end if

end subroutine

