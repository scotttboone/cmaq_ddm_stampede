C***********************************************************************
C   Portions of Models-3/CMAQ software were developed or based on      *
C   information from various groups: Federal Government employees,     *
C   contractors working on a United States Government contract, and    *
C   non-Federal sources (including research institutions).  These      *
C   research institutions have given the Government permission to      *
C   use, prepare derivative works, and distribute copies of their      *
C   work in Models-3/CMAQ to the public and to permit others to do     *
C   so.  EPA therefore grants similar permissions for use of the       *
C   Models-3/CMAQ software, but users are requested to provide copies  *
C   of derivative works to the Government without restrictions as to   *
C   use by others.  Users are responsible for acquiring their own      *
C   copies of commercial software associated with Models-3/CMAQ and    *
C   for complying with vendor requirements.  Software copyrights by    *
C   the MCNC Environmental Modeling Center are used with their         *
C   permissions subject to the above restrictions.                     *
C***********************************************************************

C RCS file, release, date & time of last delta, author, state, [and locker]
C $Header: /project/work/rep/CCTM/src/vdiff/acm2_inline/EMIS_DEFN.F,v 1.10 2010/02/24 18:13:07 yoj Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      MODULE S_EMIS_DEFN

#ifdef sens

C-----------------------------------------------------------------------
C Function: emissions interface for DDM-3D

C Revision History:
C   17 Aug 11 S.Napelenok: Initial implementation


C-----------------------------------------------------------------------

      USE EMIS_DEFN   
      USE AERO_EMIS           ! inherits GRID_CONF
      USE BEIS_DEFN           ! biogenic emissions
      USE BIOG_EMIS, ONLY: MSPCS
      USE PT3D_DEFN 
      USE PTMAP, ONLY: N_SPC_PTEM, PTEM_MAP, N_SPC_PTPM, PTPM_MAP
      USE PMEM_DEFN           ! defines emissions file(s) PM mapping

      IMPLICIT NONE

      REAL, ALLOCATABLE, SAVE :: SVDEMIS( :,:,:,:,: ) ! total emissions array

      INTEGER, PRIVATE, SAVE :: STRT_GC_EMIS, FINI_GC_EMIS, 
     &                          STRT_AE_EMIS, FINI_AE_EMIS,
     &                          STRT_NR_EMIS, FINI_NR_EMIS, 
     &                          STRT_TR_EMIS, FINI_TR_EMIS

      INTEGER, PRIVATE, SAVE :: STRT_GC_SPC, FINI_GC_SPC,
     &                          STRT_AE_SPC, FINI_AE_SPC,
     &                          STRT_NR_SPC, FINI_NR_SPC, 
     &                          STRT_TR_SPC, FINI_TR_SPC

      REAL,    PRIVATE, ALLOCATABLE, SAVE :: SBUFF( :,:,: )        ! read buffer
      REAL,    PRIVATE, ALLOCATABLE, SAVE :: VDEMIS_TEMP( :,:,:,: ) 
      REAL,    PRIVATE, ALLOCATABLE, SAVE :: SVDEMIS_AE( :,:,:,: )  ! aerosol emis
c     REAL,    PRIVATE, ALLOCATABLE, SAVE :: SVDEMIS_PT( :,:,:,: ) ! 3D pt src non-PM emis
c     REAL,    PRIVATE, ALLOCATABLE, SAVE :: SPMEMIS_PT( :,:,:,: ) ! 3D pt src PM emis

C     Emission rate of all aerosol species interpolated to current time
      REAL,    PRIVATE, ALLOCATABLE, SAVE :: S_EMBUFF( :,:,:,: )  ! in all grid cells
c     REAL,    PRIVATE, ALLOCATABLE, SAVE :: S_PM_EM( : )         ! in one grid cell


      LOGICAL, SAVE          :: ES
      CHARACTER( 16 ), SAVE  :: DDM3D_ES = 'DDM3D_ES'

      CONTAINS

C-----------------------------------------------------------------------

         FUNCTION S_EMIS_INIT( JDATE, JTIME, TSTEP, NPMAX ) RESULT ( SUCCESS )

         IMPLICIT NONE
C Includes:
         INCLUDE SUBST_GC_SPC
         INCLUDE SUBST_AE_SPC    ! aerosol species table
         INCLUDE SUBST_NR_SPC
         INCLUDE SUBST_TR_SPC
         INCLUDE SUBST_GC_EMIS   ! gas chem emis surrogate names & map table
         INCLUDE SUBST_AE_EMIS   ! aerosol emission surrogate names & map table
         INCLUDE SUBST_NR_EMIS   ! non-react emis surrogate names & map table
         INCLUDE SUBST_TR_EMIS   ! tracer emis surrogate names & map table

         INCLUDE SUBST_IOPARMS   ! I/O parameters definitions
#include         SUBST_IODECL    # I/O definitions and declarations

C Arguments:

         INTEGER JDATE, JTIME, TSTEP
         LOGICAL :: SUCCESS

         CHARACTER( 120 ) :: XMSG = ' '
         INTEGER  LOGDEV, STATUS

         INTEGER NPMAX

         CHARACTER( 16 ) :: PNAME = 'S_EMIS_INIT'

         CHARACTER( 80 ) :: VARDESC                ! environment variable description
         LOGICAL, EXTERNAL :: ENVYN


C-----------------------------------------------------------------------

         SUCCESS = .TRUE.

         LOGDEV = INIT3()


         ALLOCATE ( SVDEMIS ( N_SPC_EMIS+1,EMLAYS,MY_NCOLS,MY_NROWS,
     &              NPMAX ), STAT = STATUS )
         IF ( STATUS .NE. 0 ) THEN
            XMSG = 'SVDEMIS memory allocation failed'
            CALL M3WARN ( PNAME, JDATE, JTIME, XMSG )
            SUCCESS = .FALSE.; RETURN
         END IF

         STRT_GC_EMIS = 1
         FINI_GC_EMIS = N_GC_EMIS
         STRT_AE_EMIS = N_GC_EMIS + 1
         FINI_AE_EMIS = N_GC_EMIS + N_AE_EMIS
         STRT_NR_EMIS = N_GC_EMIS + N_AE_EMIS + 1
         FINI_NR_EMIS = N_GC_EMIS + N_AE_EMIS + N_NR_EMIS
         STRT_TR_EMIS = N_GC_EMIS + N_AE_EMIS + N_NR_EMIS + 1
         FINI_TR_EMIS = N_GC_EMIS + N_AE_EMIS + N_NR_EMIS + N_TR_EMIS

         STRT_GC_SPC = 1
         FINI_GC_SPC = N_GC_SPC
         STRT_AE_SPC = N_GC_SPC + 1
         FINI_AE_SPC = N_GC_SPC + N_AE_SPC
         STRT_NR_SPC = N_GC_SPC + N_AE_SPC + 1
         FINI_NR_SPC = N_GC_SPC + N_AE_SPC + N_NR_SPC
         STRT_TR_SPC = N_GC_SPC + N_AE_SPC + N_NR_SPC + 1
         FINI_TR_SPC = N_GC_SPC + N_AE_SPC + N_NR_SPC + N_TR_SPC

         ES = .FALSE.         ! default
         VARDESC = 'Flag for using split emissions files'
         ES = ENVYN( DDM3D_ES, VARDESC, ES, STATUS )
         IF ( STATUS .NE. 0 ) WRITE( LOGDEV, '(5X, A)' ) VARDESC
         IF ( STATUS .EQ. 1 ) THEN
            XMSG = 'Environment variable improperly formatted'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
         ELSE IF ( STATUS .EQ. -1 ) THEN
            XMSG =
     &       'Environment variable set, but empty ... Using default:'
            WRITE( LOGDEV, '(5X, A, I9)' ) XMSG, JTIME
         ELSE IF ( STATUS .EQ. -2 ) THEN
            XMSG = 'Environment variable not set ... Using default:'
            WRITE( LOGDEV, '(5X, A, I9)' ) XMSG, JTIME
         END IF

         ALLOCATE ( S_EMBUFF( MY_NCOLS,MY_NROWS,EMLAYS,N_EMIS_PM ),
     &                 STAT = STATUS )
         IF ( STATUS .NE. 0 ) THEN
            XMSG = '*** EMBUFF memory allocation failed'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF
         S_EMBUFF = 0.0

         ALLOCATE ( VDEMIS_TEMP( N_SPC_EMIS+1,EMLAYS,NCOLS,NROWS ),
     &                          STAT = STATUS )
         IF ( STATUS .NE. 0 ) THEN
            XMSG = 'VDEMIS_TEMP memory allocation failed'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG )
         END IF
         VDEMIS_TEMP = 0.0   ! array assignment

         IF ( N_AE_SPC .GT. 0 ) THEN   ! all EMLAYS
            ALLOCATE ( SVDEMIS_AE( N_AE_EMIS,EMLAYS,MY_NCOLS,MY_NROWS ),
     &                    STAT = STATUS )
            IF ( STATUS .NE. 0 ) THEN
               XMSG = 'SVDEMIS_AE memory allocation failed'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG )
            END IF
            SVDEMIS_AE = 0.0   ! array assignment
         ENDIF


         SUCCESS = .TRUE.; RETURN

         END FUNCTION S_EMIS_INIT

C-----------------------------------------------------------------------
         SUBROUTINE GET_SEMIS ( MDATE, MTIME, JDATE, JTIME, NPMAX,
     &                         TSTEP, CONVPA, RDX3F, RRHOJ, RJACM,
     &                         IPT, IPARM, ISPLIT,
     &                         GRID_FILE, GRID_NUM, PT3D_FILE, PT3D_NUM,
     &                         IAMOUNT,
     &                         ILAYER, STDATE, IDATE, ITIME, IREGION,
     &                         SEN_PAR )

         IMPLICIT NONE

C Includes:
         INCLUDE SUBST_AE_SPC    ! aerosol species table
         INCLUDE SUBST_GC_EMIS   ! gas chem emis surrogate names & map table
         INCLUDE SUBST_AE_EMIS   ! aerosol emission surrogate names & map table
         INCLUDE SUBST_NR_EMIS   ! non-react emis surrogate names & map table
         INCLUDE SUBST_TR_EMIS   ! tracer emis surrogate names & map table
         INCLUDE SUBST_IOPARMS   ! I/O parameters definitions
#include         SUBST_IODECL    # I/O definitions and declarations

C Arguments:

         INTEGER MDATE, MTIME, JDATE, JTIME, NPMAX
         INTEGER TSTEP( 2 )                 ! time step vector (HHMMSS)
         REAL    CONVPA             ! conversion factor to pressure in Pascals
         REAL    RDX3F( NLAYS )             ! reciprocal layer thickness
         REAL    RRHOJ( NCOLS,NROWS,NLAYS ) ! reciprocal density X Jacobian
         REAL    RJACM( NCOLS,NROWS,NLAYS ) ! reciprocal mid-layer Jacobian


         INTEGER, POINTER :: IPT( : )
         INTEGER, POINTER :: IPARM( :,: )
         INTEGER, POINTER :: ISPLIT( :,: )
         CHARACTER( 16 ), POINTER :: GRID_FILE( :,: )  ! ( NPMAX, 9 )  
         INTEGER, POINTER ::  GRID_NUM( : )
         CHARACTER( 16 ), POINTER :: PT3D_FILE( :,:,: )  ! ( NPMAX, 2, 9 )  
         INTEGER, POINTER ::  PT3D_NUM( : )
         REAL, POINTER :: IAMOUNT( :,:,: )
         INTEGER, POINTER :: ILAYER( :,: )
         INTEGER, POINTER :: IDATE( :,: )
         INTEGER, POINTER :: ITIME( :,: )
         REAL, POINTER ::  IREGION( :,:,:,: )
         INTEGER      STDATE       ! starting date,    format YYYYDDD
         CHARACTER( 8 ), POINTER :: SEN_PAR ( : )

         LOGICAL, SAVE :: FIRSTIME = .TRUE.
         INTEGER  LOGDEV

         INTEGER NP
         INTEGER COL, ROW, LAY, SPC, EMIS_SPC
         INTEGER N, V, S
         INTEGER FILENUM, LAYSNUM

         LOGICAL :: TIMEFLAG = .TRUE.
         INTEGER DATENUM, HOURNUM

         REAL             CNVTC         ! combined conversion factor
         REAL             CNVTI         ! intermediate combined conv. factor

         CHARACTER( 16 ) :: PNAME = 'GET_SEMIS'
         CHARACTER( 120 ) :: XMSG = ' '

         REAL    CONVEM          ! conversion for emissions rates to Kg/s
         INTEGER STATUS
         LOGICAL :: SUCCESS

         INTEGER, EXTERNAL :: TRIMLEN

c        LOGICAL :: PT3D_SENS = .FALSE.



         INTERFACE
           SUBROUTINE CKTIME(IDATE,JDATE,STDATE,ITIME,JTIME,TIMEFLAG,P )
             IMPLICIT NONE
             INTEGER,POINTER:: IDATE ( :,: )
             INTEGER,INTENT(IN):: JDATE
             INTEGER,INTENT(IN):: STDATE
             INTEGER,INTENT(IN):: JTIME
             INTEGER,POINTER:: ITIME (:,: )
             INTEGER,INTENT(IN)::  P
             LOGICAL,INTENT(OUT) :: TIMEFLAG
            END SUBROUTINE CKTIME
         END INTERFACE


C-----------------------------------------------------------------------

         IF ( FIRSTIME ) THEN

            FIRSTIME = .FALSE.
            LOGDEV = INIT3()

         END IF

C Set the date and hour counters used in sensitivity calls
         DATENUM = 1 + JDATE - STDATE
         HOURNUM = 1 + FLOOR ( REAL ( JTIME / 10000.0 ) )

         SVDEMIS = 0.0

c loop through all sensitivity parameters 

         DO NP = 1, NPMAX

            IF ( IPT( NP ) .NE. 3 ) CYCLE   ! do nothing if not an emissions sensitivity

c total emissions
            IF ( ISPLIT( NP, 1 ) .EQ. 1 ) THEN ! use all emissions

               SVDEMIS( :,:,:,:,NP ) = VDEMIS

            END IF

c gridded file emissions
            IF ( ISPLIT( NP, 2 ) .EQ. 1 ) THEN ! add emissions from gridded files

               IF ( ES ) THEN

                  DO FILENUM = 1, GRID_NUM( NP )  ! cycle through all files supplied

                     CALL S_OPEMIS( MDATE, MTIME, N_GC_EMIS, EM_TRAC, 
     &                               CONVEM, LAYSNUM, GRID_FILE( NP, FILENUM ) )

                     ALLOCATE ( SBUFF( NCOLS,NROWS,LAYSNUM ), STAT = STATUS )
                     IF ( STATUS .NE. 0 ) THEN
                        XMSG = 'BUFF memory allocation failed'
                        CALL M3WARN ( PNAME, MDATE, MTIME, XMSG )
                        SUCCESS = .FALSE.; RETURN
                     END IF
                     SBUFF = 0.0   ! array assignment

c                    ALLOCATE ( VDEMIS_TEMP( N_SPC_EMIS+1,LAYSNUM,NCOLS,NROWS ),
c    &                          STAT = STATUS )
c                    IF ( STATUS .NE. 0 ) THEN
c                       XMSG = 'VDEMIS_TEMP memory allocation failed'
c                       CALL M3WARN ( PNAME, MDATE, MTIME, XMSG )
c                       SUCCESS = .FALSE.; RETURN
c                    END IF
c                    VDEMIS_TEMP = 0.0   ! array assignment


c add the emissions from this file

C reactive gases (conversion to ppmv/s) 
                     DO V = STRT_GC_EMIS, FINI_GC_EMIS

                        S = V + 1 - STRT_GC_EMIS
                        SBUFF = 0.0   ! array assignment
                        IF ( .NOT. INTERPX( GRID_FILE( NP, FILENUM ), 
     &                          GC_EMIS( S ), PNAME,
     &                          STARTCOL,ENDCOL, STARTROW,ENDROW, 1,LAYSNUM,
     &                          MDATE, MTIME, SBUFF ) ) THEN
                           XMSG = 'Could not read  '
     &                            // GC_EMIS( S )( 1:TRIMLEN( GC_EMIS( S ) ) )
     &                            // ' from ' // GRID_FILE( NP, FILENUM )
                           CALL M3WARN( PNAME, MDATE, MTIME, XMSG )
                        END IF

                        CNVTI = CNVTP * GC_EMIS_FAC( S )
 
                        DO ROW = 1, MY_NROWS
                           DO COL = 1, MY_NCOLS
                              DO LAY = 1, LAYSNUM 
                                 CNVTC = CNVTI * RDX3F( LAY ) * RRHOJ( COL,ROW,LAY )

                                 VDEMIS_TEMP( V,LAY,COL,ROW ) = 
     &                                    CNVTC * SBUFF( COL,ROW,LAY )

                                 SVDEMIS( V,LAY,COL,ROW,NP ) = 
     &                                    SVDEMIS( V,LAY,COL,ROW,NP )
     &                                    + VDEMIS_TEMP( V,LAY,COL,ROW )
                              END DO
                           END DO
                        END DO

                     END DO   ! V

C aerosol emissions

                     IF ( N_AE_SPC .GT. 0 ) THEN   ! all EMLAYS

                        IF ( .NOT. S_AERO_EMIS_INIT ( MDATE, MTIME, 
     &                    TSTEP, LAYSNUM, GRID_FILE( NP, FILENUM ) ) ) THEN
                           XMSG = 'Failure in S_AERO_EMIS_INIT'
                           CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
                        END IF

                        SVDEMIS_AE = 0.0
                        CALL GET_S_AERO_EMIS ( MDATE, MTIME, TSTEP, LAYSNUM, RJACM )

                        DO ROW = 1, MY_NROWS
                           DO COL = 1, MY_NCOLS
                              DO LAY = 1, LAYSNUM
                                 DO V = STRT_AE_EMIS, FINI_AE_EMIS
                                    S = V + 1 - STRT_AE_EMIS
                                    SVDEMIS( V,LAY,COL,ROW,NP ) = SVDEMIS( V,LAY,COL,ROW,NP )
     &                                    + SVDEMIS_AE( S,LAY,COL,ROW )
                                 END DO
                              END DO
                           END DO
                        END DO

                     END IF


C non-reactive gases (conversion to ppmv/s) 

                     DO V = STRT_NR_EMIS, FINI_NR_EMIS
                        S = V + 1 - STRT_NR_EMIS
                        SBUFF = 0.0   ! array assignment
                        IF ( .NOT. INTERPX( GRID_FILE( NP, FILENUM ), 
     &                             NR_EMIS( S ), PNAME,
     &                             STARTCOL,ENDCOL, STARTROW,ENDROW, 1,LAYSNUM,
     &                             MDATE, MTIME, SBUFF ) ) THEN
                           XMSG = 'Could not read  '
     &                            // NR_EMIS( S )( 1:TRIMLEN( NR_EMIS( S ) ) )
     &                            // ' from ' // GRID_FILE( NP, FILENUM )
                           CALL M3WARN( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
                        END IF
                        CNVTI = CNVTP * NR_EMIS_FAC( S )
                        DO ROW = 1, MY_NROWS
                           DO COL = 1, MY_NCOLS
                              DO LAY = 1, LAYSNUM
                                 CNVTC = CNVTI * RDX3F( LAY ) * RRHOJ( COL,ROW,LAY )

                                 VDEMIS_TEMP( V,LAY,COL,ROW ) =
     &                                    CNVTC * SBUFF( COL,ROW,LAY )

                                 SVDEMIS( V,LAY,COL,ROW,NP ) =
     &                                    SVDEMIS( V,LAY,COL,ROW,NP )
     &                                    + VDEMIS_TEMP( V,LAY,COL,ROW )
                              END DO
                           END DO
                        END DO
                     END DO

                  DEALLOCATE ( SBUFF )

                  END DO

               ELSE
                  XMSG = 'Split emissions  sensitivity (ES) not set'
                  CALL M3EXIT( PNAME, 0, 0, XMSG, XSTAT1 )
               END IF

            END IF

c inline biogenics emissions
            IF ( ISPLIT( NP, 3 ) .EQ. 1 ) THEN ! add emissions from inline biogenics

               IF ( BIOGEMIS ) THEN
                  LAY = 1
                  CNVTI = CONVPA * CNVTE * RDX3F( LAY )  ! per sec
                  DO ROW = 1, MY_NROWS
                     DO COL = 1, MY_NCOLS
                        CNVTC = CNVTI * RRHOJ( COL,ROW,LAY )
                        N = 0
                        DO S = 1, MSPCS
                           V = BEIS_MAP( S )
                           IF ( V .GT. 0 ) THEN
                              N = N + 1
                              SVDEMIS( V,LAY,COL,ROW,NP ) = SVDEMIS( V,LAY,COL,ROW,NP )
     &                                    + CNVTC * VDEMIS_BI( N,COL,ROW )
                           END IF
                        END DO
                    END DO
                  END DO
               ELSE
                  XMSG = 'Can not specify sensitivity to inline beis' 
                  CALL M3EXIT( PNAME, 0, 0, XMSG, XSTAT1 )
               END IF

            END IF

c inline point sources
            IF ( ISPLIT( NP, 4 ) .EQ. 1 ) THEN  ! add emissions from inline sources

               IF ( PT3DEMIS ) THEN
c gas phase point sources
                  DO ROW = 1, MY_NROWS
                     DO COL = 1, MY_NCOLS
                        DO LAY = 1, EMLAYS
                           CNVTC = CNVTP * RDX3F( LAY ) * RRHOJ( COL,ROW,LAY )
                           DO S = 1, N_SPC_PTEM
                              V = PTEM_MAP( S )  ! <- pointer to non-PM (non-aerosol) species
                              SVDEMIS( V,LAY,COL,ROW,NP ) = SVDEMIS( V,LAY,COL,ROW,NP )
     &                                    + CNVTC * SVDEMIS_PT( COL,ROW,LAY,S,NP )
                           END DO
                        END DO
                     END DO
                  END DO

c aerosol point sources

                  IF ( N_AE_SPC .GT. 0 ) THEN   

c need only SULF out of the gas array, but it's easier to get it all
                     DO V = STRT_GC_EMIS, FINI_GC_EMIS
                        DO ROW = 1, MY_NROWS
                           DO COL = 1, MY_NCOLS
                              DO LAY = 1, EMLAYS
                                 VDEMIS_TEMP( V,LAY,COL,ROW ) = SVDEMIS( V,LAY,COL,ROW,NP )
                              END DO
                           END DO
                        END DO
                     END DO

                     DO N = 1, N_SPC_PTPM   ! merge inline point src PM species

                        V = PTPM_MAP( N )   ! set in PT3D_DEFN
                        DO LAY = 1, EMLAYS
                           DO ROW = 1, MY_NROWS
                              DO COL = 1, MY_NCOLS
                                 S_EMBUFF( COL,ROW,LAY,N ) = SPMEMIS_PT( COL,ROW,LAY,V,NP )
                              END DO
                           END DO
                        END DO
                     END DO

                     SVDEMIS_AE = 0.0
                     CALL GET_S_AERO_EMIS ( MDATE, MTIME, TSTEP, EMLYRS, RJACM )
  
                     DO ROW = 1, MY_NROWS
                        DO COL = 1, MY_NCOLS
                           DO LAY = 1, EMLAYS
                              DO V = STRT_AE_EMIS, FINI_AE_EMIS
                                 S = V + 1 - STRT_AE_EMIS
                                 SVDEMIS( V,LAY,COL,ROW,NP ) = SVDEMIS( V,LAY,COL,ROW,NP )
     &                              + SVDEMIS_AE( S,LAY,COL,ROW )
                              END DO
                           END DO
                        END DO
                     END DO

                  END IF

               ELSE
                  XMSG = 'PT3DEMIS not set to TRUE'
                  CALL M3EXIT( PNAME, 0, 0, XMSG, XSTAT1 )
               END IF

            END IF

            VDEMIS( N_SPC_EMIS+1,:,:,: ) = 0.0

c subset by IREGION, IDATE, and ITIME
            CALL CKTIME (IDATE, JDATE, STDATE, 
     &                   ITIME, JTIME, TIMEFLAG, NP)

            IF ( TIMEFLAG ) THEN
              DO ROW = 1, MY_NROWS
                 DO COL = 1, MY_NCOLS
                    DO LAY = 1, EMLAYS
                       DO SPC = 1, N_SPC_EMIS
                          SVDEMIS( SPC,LAY,COL,ROW,NP ) = 
     &                                 SVDEMIS( SPC,LAY,COL,ROW,NP ) * 
     &                                 IREGION( COL,ROW,LAY,NP ) *
     &                                 REAL( IDATE( NP,DATENUM ) ) 
                       END DO
                    END DO
                 END DO
              END DO

            ELSE
              SVDEMIS( :,:,:,:,NP ) = 0.0
            END IF

c subset by species (IPARM) and apply AMOUNT if specified

c gas species
            DO SPC = STRT_GC_EMIS, FINI_GC_EMIS
               IF ( IAMOUNT( NP, GC_EMIS_MAP( SPC ), 25 ) .EQ. 0.0 ) THEN
                  DO ROW = 1, MY_NROWS
                     DO COL = 1, MY_NCOLS
                        DO LAY = 1, EMLAYS
                           SVDEMIS( SPC,LAY,COL,ROW,NP ) = 
     &                     SVDEMIS( SPC,LAY,COL,ROW,NP ) *
     &                     REAL( IPARM( NP, GC_EMIS_MAP( SPC ) ) )
                        END DO
                     END DO
                  END DO
               ELSE
                  DO ROW = 1, MY_NROWS
                     DO COL = 1, MY_NCOLS
                        DO LAY = 1, EMLAYS
                           SVDEMIS( SPC,LAY,COL,ROW,NP ) =
     &                     IAMOUNT( NP, GC_EMIS_MAP( SPC ), HOURNUM ) *
     &                     CNVTE * RDX3F( LAY ) * RRHOJ( COL,ROW,LAY ) *
     &                     REAL( ILAYER ( NP, LAY ) ) *
     &                     IREGION( COL, ROW, LAY, NP ) *
     &                     REAL( IDATE ( NP, DATENUM ) ) *
     &                     REAL( IPARM ( NP, GC_EMIS_MAP( SPC ) ) )
                        END DO
                     END DO
                  END DO
               END IF
            END DO

c aerosol species
            EMIS_SPC = 0
            DO SPC = STRT_AE_EMIS, FINI_AE_EMIS
               EMIS_SPC = EMIS_SPC + 1

               IF ( IAMOUNT( NP, STRT_AE_SPC + AE_EMIS_MAP( EMIS_SPC ), 
     &              25 ) .EQ. 0.0 ) THEN
                  DO ROW = 1, MY_NROWS
                     DO COL = 1, MY_NCOLS
                        DO LAY = 1, EMLAYS
                           SVDEMIS( SPC,LAY,COL,ROW,NP ) =
     &                     SVDEMIS( SPC,LAY,COL,ROW,NP ) *
     &                     REAL( IPARM( NP, STRT_AE_SPC + AE_EMIS_MAP( EMIS_SPC ) ) )
                        END DO
                     END DO
                  END DO
               ELSE
                  DO ROW = 1, MY_NROWS
                     DO COL = 1, MY_NCOLS
                        DO LAY = 1, EMLAYS
                           SVDEMIS( SPC,LAY,COL,ROW,NP ) =
     &                     IAMOUNT( NP, STRT_AE_SPC + AE_EMIS_MAP( EMIS_SPC ), HOURNUM ) *
     &                     CNVTE * RDX3F( LAY ) * RRHOJ( COL,ROW,LAY ) *
     &                     REAL( ILAYER ( NP, LAY ) ) *
     &                     IREGION( COL, ROW, LAY, NP ) *
     &                     REAL( IDATE ( NP, DATENUM ) ) *
     &                     REAL( IPARM ( NP, STRT_AE_SPC + AE_EMIS_MAP( EMIS_SPC ) ) )
                        END DO
                     END DO
                  END DO
               END IF
            END DO

c nonreactive species

            IF ( N_TR_EMIS .GE. 1 ) THEN

              EMIS_SPC = 0
              DO SPC = STRT_NR_EMIS, FINI_NR_EMIS
                 EMIS_SPC = EMIS_SPC + 1
                 IF ( IAMOUNT( NP, STRT_NR_SPC + NR_EMIS_MAP( EMIS_SPC ), 
     &                25 ) .EQ. 0.0 ) THEN
                    DO ROW = 1, MY_NROWS
                       DO COL = 1, MY_NCOLS
                          DO LAY = 1, EMLAYS
                             SVDEMIS( SPC,LAY,COL,ROW,NP ) =
     &                       SVDEMIS( SPC,LAY,COL,ROW,NP ) *
     &                       REAL( IPARM( NP, STRT_NR_SPC + NR_EMIS_MAP( EMIS_SPC ) ) )
                          END DO
                       END DO
                    END DO
                 ELSE
                    DO ROW = 1, MY_NROWS
                       DO COL = 1, MY_NCOLS
                          DO LAY = 1, EMLAYS
                             SVDEMIS( SPC,LAY,COL,ROW,NP ) = 
     &                       IAMOUNT( NP, STRT_NR_SPC + NR_EMIS_MAP( EMIS_SPC ), HOURNUM ) *
     &                       CNVTE * RDX3F( LAY ) * RRHOJ( COL,ROW,LAY ) *
     &                       REAL( ILAYER ( NP, LAY ) ) * 
     &                       IREGION( COL, ROW, LAY, NP ) *
     &                       REAL( IDATE ( NP, DATENUM ) ) *
     &                       REAL( IPARM ( NP, STRT_NR_SPC + NR_EMIS_MAP( EMIS_SPC ) ) )
                          END DO 
                       END DO 
                    END DO
                 END IF
              END DO 

            END IF

c tracer species
            EMIS_SPC = 0
            DO SPC = STRT_TR_EMIS, FINI_TR_EMIS
               EMIS_SPC = EMIS_SPC + 1
               IF ( IAMOUNT( NP, STRT_TR_SPC + TR_EMIS_MAP( EMIS_SPC ), 
     &              25 ) .EQ. 0.0 ) THEN
                  DO ROW = 1, MY_NROWS
                     DO COL = 1, MY_NCOLS
                        DO LAY = 1, EMLAYS
                           SVDEMIS( SPC,LAY,COL,ROW,NP ) =
     &                     SVDEMIS( SPC,LAY,COL,ROW,NP ) *
     &                     REAL( IPARM( NP, STRT_TR_SPC + TR_EMIS_MAP( EMIS_SPC ) ) )
                        END DO 
                     END DO
                  END DO
               ELSE
                  DO ROW = 1, MY_NROWS
                     DO COL = 1, MY_NCOLS
                        DO LAY = 1, EMLAYS
                           SVDEMIS( SPC,LAY,COL,ROW,NP ) =
     &                     IAMOUNT( NP, STRT_TR_SPC + TR_EMIS_MAP( EMIS_SPC ), HOURNUM ) *
     &                     CNVTE * RDX3F( LAY ) * RRHOJ( COL,ROW,LAY ) *
     &                     REAL( ILAYER ( NP, LAY ) ) *
     &                     IREGION( COL, ROW, LAY, NP ) *
     &                     REAL( IDATE ( NP, DATENUM ) ) *
     &                     REAL( IPARM ( NP, STRT_TR_SPC + TR_EMIS_MAP( EMIS_SPC ) ) )
                        END DO
                     END DO
                  END DO
               END IF
            END DO

         END DO

         END SUBROUTINE GET_SEMIS

C-----------------------------------------------------------------------

         FUNCTION S_AERO_EMIS_INIT( MDATE, MTIME, TSTEP, LAYSNUM, EMISFILE )
     &                   RESULT ( SUCCESS )

         IMPLICIT NONE

         INCLUDE SUBST_IOPARMS   ! I/O parameters definitions
#include         SUBST_IODECL    # I/O definitions and declarations


         INTEGER      MDATE                      ! current model date, coded YYYYDDD
         INTEGER      MTIME                      ! current model time, coded HHMMSS
         INTEGER      TSTEP( 2 )                 ! time step vector (HHMMSS)
         INTEGER      LAYSNUM
         CHARACTER( 16 ) :: EMISFILE             ! emissions file to read

         LOGICAL :: SUCCESS 

         CHARACTER( 16 ) :: PNAME = 'S_AERO_EMIS_INIT'

         CHARACTER( 120 ) :: XMSG = ' '
         LOGICAL, SAVE     :: FIRSTIME = .TRUE.
         INTEGER, SAVE     :: LOGDEV
         INTEGER STATUS                   ! ENV..., ALLOCATE status
         INTEGER N
         INTEGER, EXTERNAL :: TRIMLEN

c        REAL,    PRIVATE, ALLOCATABLE, SAVE :: S_PM_EM( : )         ! in one grid cell


C-----------------------------------------------------------------------

         IF ( FIRSTIME ) THEN
            FIRSTIME = .FALSE.
            LOGDEV = INIT3()

c           ALLOCATE ( S_EMBUFF( MY_NCOLS,MY_NROWS,LAYSNUM,N_EMIS_PM ),
c    &                 STAT = STATUS )
c           IF ( STATUS .NE. 0 ) THEN
c              XMSG = '*** EMBUFF memory allocation failed'
c              CALL M3EXIT ( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
c           END IF
c           S_EMBUFF = 0.0

c           ALLOCATE ( S_PM_EM( N_EMIS_PM ), STAT = STATUS )
c           IF ( STATUS .NE. 0 ) THEN
c              XMSG = '*** S_PM_EM memory allocation failed'
c              CALL M3EXIT ( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
c           END IF
c           S_PM_EM = 0.0

         END IF    ! FIRSTIME

         DO N = 1, N_EMIS_PM   ! read 3d file
            IF ( .NOT. INTERPX( EMISFILE, PM_EMIS( N ), PNAME,
     &                             STARTCOL,ENDCOL, STARTROW,ENDROW, 1,LAYSNUM,
     &                             MDATE, MTIME, S_EMBUFF( :,:,:,N ) ) ) THEN
               XMSG = 'Could not read '
     &                 // PM_EMIS( N )( 1:TRIMLEN( PM_EMIS( N ) ) )
     &                 // ' from ' // EMISFILE
               CALL M3WARN ( PNAME, MDATE, MTIME, XMSG  )
            END IF
         END DO


         SUCCESS = .TRUE.; RETURN

         END FUNCTION S_AERO_EMIS_INIT


c-----------------------------------------------------------------------

         SUBROUTINE GET_S_AERO_EMIS ( MDATE, MTIME, TSTEP, LAYSNUM, RJACM )

         IMPLICIT NONE

         INCLUDE SUBST_CONST     ! physical and mathematical constants

         INCLUDE SUBST_IOPARMS   ! I/O parameters definitions
#include         SUBST_IODECL    # I/O definitions and declarations

         INTEGER      MDATE                      ! current model date, coded YYYYDDD
         INTEGER      MTIME                      ! current model time, coded HHMMSS
         INTEGER      TSTEP( 2 )                 ! time step vector (HHMMSS)
         INTEGER      LAYSNUM
         REAL         RJACM( NCOLS,NROWS,NLAYS ) ! reciprocal mid-layer Jacobian
c        CHARACTER( 16 ) :: EMISFILE             ! emissions file to read

         CHARACTER( 16 ) :: PNAME = 'GET_S_AERO_EMIS'

C     Emission rate of all aerosol species interpolated to current time
c        REAL, ALLOCATABLE, SAVE :: S_EMBUFF( :,:,:,: )  ! in all grid cells
         REAL, ALLOCATABLE, SAVE :: S_PM_EM( : )         ! in one grid cell

         INTEGER N
         INTEGER COL, ROW, LAY
         CHARACTER( 120 ) :: XMSG = ' '

c        INTEGER, EXTERNAL :: TRIMLEN

         LOGICAL, SAVE     :: FIRSTIME = .TRUE.

         REAL  GSFAC                                ! into [ug/m3/sec]
         REAL  CONVM     ! conversion factor for mass emissions [m3/mol]
         REAL  CONVN     ! conversion factor for number emissions [1e6*m3]
         REAL  CONVS     ! conversion factor for surface area emissions [m3/mol]

C     Aerosol mass-emission rates [ug/m3/s]
         REAL EPOC          ! fine primary organic carbon
         REAL EPEC          ! fine primary elemental carbon
         REAL EPNO3         ! fine primary nitrate
         REAL EPSO4         ! fine primary sulfate
         REAL EPMF          ! fine primary unspeciated pm
         REAL EPMCO         ! coarse anthropogenic aerosol
         REAL ESOILCO       ! coarse soil-derived aerosol
         REAL ESEASCO       ! coarse marine aerosol

C     Variables for handling vapor-phase sulfuric acid emissions
         REAL EMSULF                          ! emission rate [ppmv/s]
         REAL, PARAMETER :: MWH2SO4 = 98.0    ! molar mass [g/mol]

C     Speciation factors for coarse mode emissions, from paragraph #15 of
C     Binkowski & Roselle (2003)
         REAL, PARAMETER :: FAC_DUST  = 0.90  ! (fugitive dust)/PMC
         REAL, PARAMETER :: FAC_OTHER = 0.10  ! (non-fugitive dust)/PMC

C     Factors for splitting primary carbon emissions into Aitken and
C     accumulation modes, from paragraph #12 of Binkowski & Roselle (2003)
         REAL, PARAMETER :: FACEM25_ACC  = 0.999 ! accumulation mode
         REAL, PARAMETER :: FACEM25_ATKN = 0.001 ! Aitken mode

C     Mode-specific mass-emission rates [ug/m3/s]
         REAL EPM25AT       ! Aitken mode unspeciated aerosol
         REAL EPM25AC       ! accumulation mode unspeciated aerosol
         REAL EPORGAT       ! Aitken mode primary organic aerosol
         REAL EPORGAC       ! accumulation mode primary organic aerosol
         REAL EPECAT        ! Aitken mode elemental carbon
         REAL EPECAC        ! accumulation mode elemental carbon
         REAL EPSO4AT       ! Aitken mode primary sulfate
         REAL EPSO4AC       ! accumulation mode primary sulfate
         REAL EPSO4CO       ! coarse mode primary sulfate
         REAL EPNO3AT       ! Aitken mode primary nitrate
         REAL EPNO3AC       ! accumulation mode primary nitrate
c        REAL EPNAAT        ! Aitken mode sodium
c        REAL EPNAAC        ! accumulation mode sodium
c        REAL EPNACO        ! coarse mode sodium
c        REAL EPCLAT        ! Aitken mode chloride
c        REAL EPCLAC        ! accumulation mode chloride
c        REAL EPCLCO        ! coarse mode chloride
c        REAL EPH2OAT       ! Aitken mode primary water
c        REAL EPH2OAC       ! accumulation mode primary water
c        REAL EPH2OCO       ! coarse mode primary water

C     Geometric Constants
         REAL, PARAMETER :: F6DPI = 6.0 / PI
         REAL, PARAMETER :: F6DPIM9 = 1.0E-9 * F6DPI

C     Variables for converting emission rates into molar-mixing-ratio units
         REAL, PARAMETER :: GPKG = 1.0E+03     ! g/kg
         REAL, PARAMETER :: RAVO = 1.0 / AVO   ! reciprocal of Avogadro

C     Factors for converting aerosol mass concentration [ug/m3] to 3rd
C     moment concentration [m3/m3]
         REAL, PARAMETER :: SO4FAC  = F6DPIM9 / RHOSO4
         REAL, PARAMETER :: NO3FAC  = F6DPIM9 / RHONO3
         REAL, PARAMETER :: ORGFAC  = F6DPIM9 / RHOORG
         REAL, PARAMETER :: SOILFAC = F6DPIM9 / RHOSOIL
         REAL, PARAMETER :: SEASFAC = F6DPIM9 / RHOSEAS
         REAL, PARAMETER :: ANTHFAC = F6DPIM9 / RHOANTH

C     Third moment emissions rates [m3/m3/s]
         REAL EMISM3AT      ! Aitken mode
         REAL EMISM3AC      ! accumulation mode
         REAL EMISM3CO      ! coarse mode

C     Number emissions rates [1/m3/s]
         REAL EM_NUMATKN    ! Aitken mode
         REAL EM_NUMACC     ! accumulation mode
         REAL EM_NUMCOR     ! coarse mode

C     Surface area emission rates [m2/m3/s]
         REAL EM_SRFATKN    ! Aitken mode
         REAL EM_SRFACC     ! accumulation mode
         REAL EM_SRFCOR     ! coarse mode

         INTEGER STATUS                   ! ENV..., ALLOCATE status
         INTEGER, SAVE     :: LOGDEV


C-----------------------------------------------------------------------

         IF ( FIRSTIME ) THEN
            FIRSTIME = .FALSE.
            LOGDEV = INIT3()

            ALLOCATE ( S_PM_EM( N_EMIS_PM ), STAT = STATUS )
            IF ( STATUS .NE. 0 ) THEN
               XMSG = '*** S_PM_EM memory allocation failed'
               CALL M3EXIT ( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
            END IF
            S_PM_EM = 0.0

         END IF    ! FIRSTIME

         DO LAY = 1, LAYSNUM
            DO ROW = 1, MY_NROWS
               DO COL = 1, MY_NCOLS

C *** Store aerosol emission rates for this grid cell in a 1D array
                  DO N = 1, N_EMIS_PM
                     S_PM_EM( N ) = S_EMBUFF( COL,ROW,LAY,N )
                  END DO

C *** Calculate scaling factor for converting mass emissions into [ug/m3/s]
C     note: RJACM converts grid heights from sigma coordinates to meters
C     Also calculate scaling factors for converting to molar-mixing-ratio units

                  GSFAC = CONVEM_PM( LAY ) * RJACM( COL,ROW,LAY )
                  CONVM = MWAIR / GPKG / DENS( COL,ROW,LAY )   ! [m3/mol]
                  CONVN = CONVM * RAVO * 1.0E+06         ! 10^6 * [m3]
                  CONVS = CONVM                          ! [m3/mol]

C *** Calculate speciated mass emission rates for fine aerosol [ug/m3/s]

                  EPOC  = S_PM_EM( VPOC  ) * GSFAC
                  EPSO4 = S_PM_EM( VPSO4 ) * GSFAC
                  EPNO3 = S_PM_EM( VPNO3 ) * GSFAC
                  EPMF  = S_PM_EM( VPMF  ) * GSFAC
                  EPEC  = S_PM_EM( VPEC  ) * GSFAC

C *** If INDX10 is non-zero, the emission file did not contain PMC so 
C     PM_EM(VPMC) holds PM10 emissions.  Substract the PM2.5 species
C     from this so PM_EM(VPMC) holds only the coarse emission rate.

                  IF ( INDX10 .NE. 0 ) THEN
                     S_PM_EM( VPMC ) = MAX( 0.0,
     &                                    S_PM_EM( VPMC )
     &                                    - ( S_PM_EM( VPOC )
     &                                      + S_PM_EM( VPSO4 )
     &                                      + S_PM_EM( VPNO3 )
     &                                      + S_PM_EM( VPEC )
     &                                      + S_PM_EM( VPMF ) )
     &                                   )
                  END IF

C *** Assign coarse PM emissions (excluding sea salt) to the coarse mode.
C     Split coarse material between SOIL and UNSPECIATED material, as
C     described in paragraph #15 of Binkowski & Roselle (2003).

                  ESOILCO = FAC_DUST  * S_PM_EM( VPMC ) * GSFAC
                  EPMCO   = FAC_OTHER * S_PM_EM( VPMC ) * GSFAC
                  ESEASCO = 0.0

C *** Extract H2SO4 vapor emission rate from the VDEMIS array, add it to
C     the fine-PM sulfate emissions, and remove it from the gas emissions

                  EMSULF = VDEMIS_TEMP( VSULF,LAY,COL,ROW )
                  EPSO4 = EPSO4 + EMSULF * MWH2SO4 / CONVM
                  VDEMIS_TEMP( VSULF,LAY,COL,ROW ) = 0.0

C *** Assign fine-particle emission rates to the fine modes.
C     Assume all non-seasalt-emissions of sulfate, nitrate, and unspeciated
C     fine PM are in the accumulation mode.  Split the carbon emissions 
C     between the Aitken and accumulation modes, as described in paragraph 
C     #12 of Binkowski & Roselle (2003).

                  EPM25AT = 0.0
                  EPM25AC = EPMF
                  EPSO4AT = 0.0
                  EPSO4AC = EPSO4
                  EPNO3AT = 0.0
                  EPNO3AC = EPNO3
                  EPORGAT = FACEM25_ATKN * EPOC
                  EPORGAC = FACEM25_ACC  * EPOC
                  EPECAT  = FACEM25_ATKN * EPEC
                  EPECAC  = FACEM25_ACC  * EPEC

C *** Calculate emissions rate for third moments [m3/m3/s] of each mode 
C     (excluding sea salt), as in Equation 7a of Binkowski & Roselle (2003).

                  EMISM3AT =  ANTHFAC * ( EPM25AT + EPECAT )
     &                     +  ORGFAC  * EPORGAT
     &                     +  SO4FAC  * EPSO4AT
     &                     +  NO3FAC  * EPNO3AT
                  EMISM3AC =  ANTHFAC * ( EPM25AC + EPECAC )
     &                     +  ORGFAC  * EPORGAC
     &                     +  SO4FAC  * EPSO4AC
     &                     +  NO3FAC  * EPNO3AC
                  EMISM3CO =  SOILFAC * ESOILCO
     &                     +  SEASFAC * ESEASCO
     &                     +  ANTHFAC * EPMCO

C *** Calculate the number emissions rate for each mode [1/m3/s], using 
C     Equation 7b of Binkowski & Roselle (2003).

                  EM_NUMATKN = FACTNUMAT * EMISM3AT
                  EM_NUMACC  = FACTNUMAC * EMISM3AC
                  EM_NUMCOR  = FACTNUMC  * EMISM3CO

C *** Calculate the surface area emissions rate for the fine modes [m2/m3/s],
C     using Equation 7c of Binkowski & Roselle (2003).  Multiplying by PI 
C     converts 2nd moment to surface area.

                  EM_SRFATKN = PI * FACTM2AT * EMISM3AT
                  EM_SRFACC  = PI * FACTM2AC * EMISM3AC
                  EM_SRFCOR  = PI * FACTM2CO * EMISM3CO

                  SVDEMIS_AE( VSO4AI, LAY,COL,ROW ) = EPSO4AT * CONVM * R_AE_EM_MW( VSO4AI )
                  SVDEMIS_AE( VSO4AJ, LAY,COL,ROW ) = EPSO4AC * CONVM * R_AE_EM_MW( VSO4AJ )
                  SVDEMIS_AE( VNO3AI, LAY,COL,ROW ) = EPNO3AT * CONVM * R_AE_EM_MW( VNO3AI )
                  SVDEMIS_AE( VNO3AJ, LAY,COL,ROW ) = EPNO3AC * CONVM * R_AE_EM_MW( VNO3AJ )
                  SVDEMIS_AE( VORGPAI,LAY,COL,ROW ) = EPORGAT * CONVM * R_AE_EM_MW( VORGPAI )
                  SVDEMIS_AE( VORGPAJ,LAY,COL,ROW ) = EPORGAC * CONVM * R_AE_EM_MW( VORGPAJ )
                  SVDEMIS_AE( VECI,   LAY,COL,ROW ) = EPECAT  * CONVM * R_AE_EM_MW( VECI )
                  SVDEMIS_AE( VECJ,   LAY,COL,ROW ) = EPECAC  * CONVM * R_AE_EM_MW( VECJ )
                  SVDEMIS_AE( VP25AI, LAY,COL,ROW ) = EPM25AT * CONVM * R_AE_EM_MW( VP25AI )
                  SVDEMIS_AE( VP25AJ, LAY,COL,ROW ) = EPM25AC * CONVM * R_AE_EM_MW( VP25AJ )
                  SVDEMIS_AE( VANTHA, LAY,COL,ROW ) = EPMCO   * CONVM * R_AE_EM_MW( VANTHA )
                  SVDEMIS_AE( VSOILA, LAY,COL,ROW ) = ESOILCO * CONVM * R_AE_EM_MW( VSOILA )
                  SVDEMIS_AE( VAT0,   LAY,COL,ROW ) = EM_NUMATKN * CONVN
                  SVDEMIS_AE( VAC0,   LAY,COL,ROW ) = EM_NUMACC  * CONVN
                  SVDEMIS_AE( VCOR0,  LAY,COL,ROW ) = EM_NUMCOR  * CONVN
                  SVDEMIS_AE( VSURFAT,LAY,COL,ROW ) = EM_SRFATKN * CONVS
                  SVDEMIS_AE( VSURFAC,LAY,COL,ROW ) = EM_SRFACC  * CONVS
                  SVDEMIS_AE( VSURFCO,LAY,COL,ROW ) = EM_SRFCOR  * CONVS

               END DO
            END DO
         END DO

         END SUBROUTINE GET_S_AERO_EMIS

#endif


      END MODULE S_EMIS_DEFN

