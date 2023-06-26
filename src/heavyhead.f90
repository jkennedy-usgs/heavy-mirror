! ********************************************************************
! Module for reading MODFLOW head.
! ********************************************************************
MODULE DOUBLEBUDGET
   INTEGER, SAVE ::IPRBUD, IPRHEAD
   DOUBLE PRECISION, POINTER, SAVE, DIMENSION(:, :, :) ::DBLBUFF
   REAL(KIND=4), POINTER, SAVE, DIMENSION(:, :, :) ::SNGBUFF
END MODULE
!
!
SUBROUTINE HEADPRECISION(IU, IOUT, NCOL, NROW, NLAY)
! ********************************************************************
! Determine single or double precision file type for a MODFLOW binary
! head file:  0=unrecognized, 1=single, 2=double.
! ********************************************************************
   USE DOUBLEBUDGET, ONLY: IPRHEAD, DBLBUFF
!
   DOUBLE PRECISION PERTIMD, TOTIMD
   REAL(KIND=4) PERTIMS, TOTIMS
   CHARACTER*16 TEXT
   INTEGER IOUT
!
   IF (IPRHEAD .EQ. 1 .OR. IPRHEAD .EQ. 2) RETURN
!
!  SINGLE precision check
   READ (IU, ERR=50, END=50) KSTP, KPER, PERTIMS, TOTIMS, TEXT
   IF (TEXT .EQ. '            HEAD') THEN
      IPRHEAD = 1
      GO TO 100
   END IF
!
!  DOUBLE precision check
50 REWIND (IU)
   READ (IU, ERR=100, END=100) KSTP, KPER, PERTIMD, TOTIMD, TEXT
   IF (TEXT .EQ. '            HEAD') THEN
      IPRHEAD = 2
      ALLOCATE (DBLBUFF(NCOL, NROW, NLAY))
   END IF
!
100 IF (IPRHEAD .EQ. 1) THEN
      WRITE (IOUT, *) &
   &   'Single Precision Binary Head file'
   ELSE IF (IPRHEAD .EQ. 2) THEN
      WRITE (IOUT, *) &
   &   'Double Precision Binary Head file'
   ELSE
      CALL USTOP('Unable to determine the precision of the Head File')
   END IF
   REWIND (IU)
   RETURN
END subroutine
!
!
SUBROUTINE CHECKHEADFILE(IGRID, IOTYPE)
!***************************************************
! Method to check Headfile borrowed from MODPATH 6
!***************************************************
   USE GLOBAL, ONLY: NPER, NSTP
   USE HVYDATA, ONLY: SIMTIME, IOLIST

   INTEGER ::KSTP, KPER, KS, KP, IOTYPE, N
   DOUBLE PRECISION ::PERTIM, TOTIM
   LOGICAL ISOK

   ISOK = .FALSE.

   IF (IOTYPE .EQ. 1) THEN
      WRITE (IOLIST, *)
      WRITE (IOLIST, '(1X,A,I3)') 'CHECKING THE HEAD FILE FOR GRID ', IGRID
   END IF

   N = 0
   DO KPER = 1, NPER
      DO KSTP = 1, NSTP(KPER)
         CALL READTIMESTEPHEAD(KSTP, KPER, KS, KP, PERTIM, TOTIM, ISOK)
!            IF(KS.NE.KSTP .OR. KP.NE.KPER) THEN
!                WRITE(IOLIST, '(1X,A)') &
!     &          'FOUND UNEXPECTED TIME STEP IN HEAD FILE. STOP.'
!                ISOK=.FALSE.
!                RETURN
!            END IF
         IF (.NOT. ISOK) THEN
            CONTINUE
!                WRITE(MESSAGE,'(1X,A,I3,A)') 'THE HEAD FILE FOR GRID ', &
!     &                                 IGRID,' IS NOT COMPLETE. STOP'
!                CALL USTOP(MESSAGE)
         END IF
         N = N + 1
!   Store the simulation time for the end of this time step in the SIMTIME array
         IF (IGRID .EQ. 1) THEN
            SIMTIME(N) = SIMTIME(0) + TOTIM
         ELSE
            IF ((SIMTIME(0) + TOTIM) .NE. SIMTIME(N)) THEN
               WRITE (IOLIST, '(1X,A,I2,A)')  &
&               'INCONSISTENT TIME STEP STRUCTURE FOR GRID ', IGRID, &
&               ' STOP.'
               CALL USTOP('INCONSISTENT TIME STEP STRUCTURE. STOP.')
            END IF
         END IF
      END DO
   END DO

   ISOK = .TRUE.
   IF (IOTYPE .EQ. 1) THEN
      WRITE (IOLIST, '(1X,A,I3,A)') 'THE HEAD FILE FOR GRID ', &
   &                                   IGRID, ' IS COMPLETE.'
   END IF

   RETURN
END subroutine
!
!
SUBROUTINE READTIMESTEPHEAD(KSTP, KPER, KS, KP, PERTIM, TOTIM, ISOK)
   USE GLOBAL, ONLY: NIUNIT, NCOL, NROW, NLAY,  &
  &                 IUNIT

   USE DOUBLEBUDGET, ONLY: IPRHEAD
   USE HVYDATA, ONLY: IOLIST, IUHEAD
   DOUBLE PRECISION ::PERTIM, TOTIM
   INTEGER(8) ::IPOS
   INTEGER ::ISTEP, IU, KS, KP
   LOGICAL ISOK

   ISOK = .FALSE.

   IF (IUNIT(IUHEAD) .GT. 0) THEN
      IU = IUNIT(IUHEAD)
   ELSE
      WRITE (IOLIST, *) 'NO MODFLOW HEAD FILE IS OPENED. STOP.'
      ISOK = .FALSE.
      RETURN
   END IF

   IF (IPRHEAD .EQ. 0) THEN
      CALL HEADPRECISION(IU, IOLIST, NCOL, NROW, NLAY)
   END IF

   IF (IPRHEAD .EQ. 2) THEN
      IRECL = NLAY*(52 + 8*NCOL*NROW)
   ELSE IF (IPRHEAD .EQ. 1) THEN
      IRECL = NLAY*(44 + 4*NCOL*NROW)
   ELSE
      CALL USTOP( &
   &  'THE PRECISION OF THE HEAD FILE COULD NOT BE DETERMINED. STOP.')
   END IF

   CALL GETGLOBALTIMESTEP(ISTEP, KPER, KSTP)
   IPOS = (ISTEP - 1)*IRECL + 1

!     Set the file position to the beginning of the time step
   READ (IU, POS=IPOS)

   CALL READHEAD(KS, KP, PERTIM, TOTIM, ISOK)
!---- Remove the debug statement for now! IF(IOTYPE.EQ.1)
   IF (ISOK) then
      WRITE (IOLIST, '(1X,A,I5,A,I5,A,1PE14.6,A,1PE14.6)') &
     & 'READ HEAD FOR PERIOD ', KP, ' STEP ', KS, '  PERTIM = ', PERTIM, &
     & '  TOTIM = ', TOTIM
   END IF

   RETURN
END subroutine

SUBROUTINE READHEAD(KSTP, KPER, PERTIM, TOTIM, ISOK)
   USE GLOBAL, ONLY: NIUNIT, NCOL, NROW, NLAY, IUNIT
   USE HVYGLOBAL, ONLY: HEAD
   USE DOUBLEBUDGET, ONLY: IPRHEAD, DBLBUFF, SNGBUFF
   USE HVYDATA, ONLY: IOLIST, IUHEAD
   DOUBLE PRECISION    ::PERTIMD, TOTIMD
   CHARACTER(LEN=16)  ::STRING, TEXT
   INTEGER :: IU, IEND, NC, NR, K, L
   DOUBLE PRECISION ::PERTIM, TOTIM
   REAL(KIND=4) ::PERTIMS, TOTIMS
   LOGICAL ISOK

   ISOK = .FALSE.
   TEXT = '            HEAD'

   IF (IUNIT(IUHEAD) .GT. 0) THEN
      IU = IUNIT(IUHEAD)
   ELSE
      WRITE (IOLIST, *) 'NO MODFLOW HEAD FILE IS OPENED. STOP.'
      ISOK = .FALSE.
      RETURN
   END IF

   IF (IPRHEAD .EQ. 1) THEN
      IF (.NOT. ASSOCIATED(SNGBUFF)) ALLOCATE (SNGBUFF(NCOL, NROW, NLAY))
   ELSE IF (IPRHEAD .EQ. 2) THEN
      IF (.NOT. ASSOCIATED(DBLBUFF)) ALLOCATE (DBLBUFF(NCOL, NROW, NLAY))
   END IF

   DO L = 1, NLAY
      IF (IPRHEAD .EQ. 2) THEN
         READ (IU, END=10, ERR=30) KSTP, KPER, PERTIMD, TOTIMD, STRING, &
  &                             NC, NR, K
         PERTIM = PERTIMD
         TOTIM = TOTIMD
      ELSE
         READ (IU, END=10, ERR=30) KSTP, KPER, PERTIMS, TOTIMS, STRING, &
  &                             NC, NR, K
         PERTIM = PERTIMS
         TOTIM = TOTIMS
      END IF
      IF (K .NE. L) THEN
         WRITE (IOLIST, 1400) L, K
         ISOK = .FALSE.
         RETURN
      END IF

      IF (STRING .NE. TEXT) THEN
         WRITE (IOLIST, 1000) TEXT
         ISOK = .FALSE.
         RETURN
      END IF

      IF (IPRHEAD .EQ. 2) THEN
         READ (IU, ERR=60, END=70) ((DBLBUFF(J, I, 1), J=1, NCOL), I=1, NROW)
         DO I = 1, NROW
            DO J = 1, NCOL
               HEAD(J, I, K) = DBLBUFF(J, I, 1)
            END DO
         END DO
      ELSE
         READ (IU, ERR=60, END=70) ((SNGBUFF(J, I, 1), J=1, NCOL), I=1, NROW)
         DO I = 1, NROW
            DO J = 1, NCOL
               HEAD(J, I, K) = SNGBUFF(J, I, 1)
            END DO
         END DO
      END IF
   END DO

   IEND = 0
   ISOK = .TRUE.
   RETURN

!...  END OF FILE REACHED
10 IEND = 1
   ISOK = .FALSE.
   RETURN

!... PROBLEM READING HEADER RECORD
30 WRITE (IOLIST, 1100) TEXT
   ISOK = .FALSE.
   RETURN
40 WRITE (IOLIST, 1200) TEXT
   ISOK = .FALSE.
   RETURN
!... END OF FILE REACHED
50 WRITE (IOLIST, 1300) IU
   ISOK = .FALSE.
   RETURN
!... PROBLEM READING ARRAY
60 WRITE (IOLIST, 1500)
   ISOK = .FALSE.
   RETURN
70 WRITE (IOLIST, 1300) IU
   ISOK = .FALSE.
   RETURN

1000 FORMAT(1X, A, ' FILE DOES NOT HAVE VALID HEADER RECORD. STOP.')
1100 FORMAT(' ERROR READING HEADER RECORD IN UNFORMATTED ', A, &
   & ' FILE. STOP.')
1200 FORMAT(' ERROR READING HEADER RECORD IN FORMATTED ', A, &
   & ' FILE. STOP.')
1300 FORMAT(1X, 'END-OF-FILE ON UNIT ', I3, '. STOP.')
1400 FORMAT(1X, 'LAYER ', I3, 'READING HEAD FILE. ', &
   & 'EXPECTED HEAD FOR LAYER ', I3, '. FOUND HEAD FOR LAYER ', I3, &
   & '. STOP.')
1500 FORMAT(1X, 'ERROR READING HEAD ARRAY. STOP.')
END subroutine
!
!
SUBROUTINE GETGLOBALTIMESTEP(ISTEP, KPER, KSTP)
   USE GLOBAL, ONLY: NSTP, NPER
   USE HVYDATA, ONLY: NTOFF
   INTEGER ::KSTP, KPER, ISTEP
   ISTEP = 0
   IF (KPER .LT. 1 .OR. KPER .GT. NPER) THEN
      CALL USTOP('STRESS PERIOD NUMBER IS OUT OF RANGE. STOP.')
   END IF
   IF (KSTP .LT. 1 .OR. KSTP .GT. NSTP(KPER)) THEN
      CALL USTOP('TIME STEP NUMBER IS OUT OF RANGE. STOP.')
   END IF
   ISTEP = NTOFF(KPER) + KSTP
   RETURN
END

