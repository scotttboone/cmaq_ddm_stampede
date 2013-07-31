C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      MODULE S_PT3D_DEFN

C-----------------------------------------------------------------------
C Function: 3d point source emissions interface to DDM3D 

C Revision History:
C     25 Aug 2011 S.Napelenok: initial implementation from PT3D_DEFN

C-----------------------------------------------------------------------

c     USE UDTYPES           ! implementation for multiple groups of stack input files
c     USE PT3D_DEFN, ONLY: PT3DEMIS

      IMPLICIT NONE

c     LOGICAL,                       SAVE :: PT3DEMIS   ! flag in-lining plume rise
c     INTEGER,                       SAVE :: NPTGRPS    ! no. pt src input file groups
c     REAL,             ALLOCATABLE, SAVE :: VDEMIS_PT( :,:,:,: ) ! 3D pt src non-PM emis
c     REAL,             ALLOCATABLE, SAVE :: PMEMIS_PT( :,:,:,: ) ! 3D pt src PM emis

C Scenario start
      INTEGER, PRIVATE,              SAVE :: SDATE = 0   ! Julian start date (YYYYDDD)
      INTEGER, PRIVATE,              SAVE :: STIME = 0   ! start time (HHMMSS)

C Output layer fractions, dimensioned NSRC, emlays, and for diagnostic file
c     REAL,    PRIVATE, ALLOCATABLE, SAVE :: LFRAC( :,: ) ! per source local layer fractn
c     REAL,    PRIVATE, ALLOCATABLE, SAVE :: TFRAC( : )   ! Temporary LFRAC

C debugging zplume
c     REAL,    PRIVATE, ALLOCATABLE, SAVE :: ZPLUM( :,: ) ! has to be 2D for play file

C Stack emissions for sources within domain
c     REAL,    PRIVATE, ALLOCATABLE, SAVE :: STKEM( : )   ! per source local stk emis
      INTEGER, PRIVATE,              SAVE :: EMLYRS       ! no. of emis layers

C Wildfire emissions source parameters
c     REAL,    PRIVATE, ALLOCATABLE, SAVE :: ACRES( : )   ! acres burned
c     REAL,    PRIVATE, ALLOCATABLE, SAVE :: BFLUX( : )   ! Briggs buoyancy flux

C Vertical coord values
c     REAL,    PRIVATE, ALLOCATABLE, SAVE :: VGLVSXG( : )

C Play report log
c     INTEGER, PRIVATE, ALLOCATABLE, SAVE :: RDEV( : )
C Play files source id file
c#ifdef srcid
c     INTEGER, PRIVATE, ALLOCATABLE, SAVE :: LDEV( : )
c#endif
C Plume rise info report files
c     CHARACTER( 13 ), PRIVATE, ALLOCATABLE, SAVE :: REPTNAME( : )
c     CHARACTER( 11 ), PRIVATE                    :: REPSTR
C Minimum layer for reporting srcs w/ high plumes
c     INTEGER, PRIVATE,              SAVE :: REP_LAYR

C Layer fractions matrix output files
c     CHARACTER(  16 ), PRIVATE, ALLOCATABLE, SAVE :: PLAYNAME( : )
C Layer fractions source id output files
c     CHARACTER(  16 ), PRIVATE, ALLOCATABLE, SAVE :: PLAY_SRCID_NAME( : )

C Get heights from GRID_CRO file, if true
c     LOGICAL, PRIVATE,              SAVE :: ZSTATIC = .TRUE.
C Write point source 3d emis diagnostic file, if true
c     LOGICAL, PRIVATE,              SAVE :: PT3DDIAG = .FALSE.
C Point source 3d emis diagnostic file names
!     CHARACTER( 16 ),               SAVE :: PT3DNAME
C Point source 3d emis diagnostic integral average write buffers
c     REAL,    PRIVATE, ALLOCATABLE, SAVE :: VDEMIS_BUF( :,:,:,: )
c     REAL,    PRIVATE, ALLOCATABLE, SAVE :: PMEMIS_BUF( :,:,:,: )
C Diagnostic file avg factor
c     REAL,    PRIVATE                    :: DIVFAC
C Write layer fractions diagnostic file, if true
c     LOGICAL, PRIVATE,              SAVE :: PT3DFRAC = .FALSE.

C Numerical flag for plume vertical spread method
      INTEGER, PRIVATE,              SAVE :: IPVERT

C Stack parameters
c     TYPE( RARRY1 ), PRIVATE, ALLOCATABLE, SAVE :: STK_DM( : )
c     TYPE( RARRY1 ), PRIVATE, ALLOCATABLE, SAVE :: STK_HT( : )
c     TYPE( RARRY1 ), PRIVATE, ALLOCATABLE, SAVE :: STK_TK( : )
c     TYPE( RARRY1 ), PRIVATE, ALLOCATABLE, SAVE :: STK_VEL( : )

C Private Parameters:
      CHARACTER(  10 ), PRIVATE, PARAMETER :: BLANK10 = ' '

C Private External Functions:
      CHARACTER(   2 ), PRIVATE,  EXTERNAL :: CRLF

C cross-point layered grid file name
c     CHARACTER(  16 ), PRIVATE,     SAVE :: GC3NAME

C hourly stack emissions file names
c     CHARACTER(  16 ), PRIVATE, ALLOCATABLE, SAVE :: STKENAME( : )

      CHARACTER( 120 ), PRIVATE           :: XMSG = ' '

      CONTAINS

C=======================================================================

         FUNCTION S_PT3D_INIT ( N_SPC_EMIS, EMLAYS, JDATE, JTIME, TSTEP,
     &                          NPMAX, PT3D_FILE, PT3D_NUM,
     &                          ISPLIT )
     &                          RESULT ( SUCCESS )

C use N_SPC_EMIS+1 for non-accessed emissions

         USE PTMET, ONLY: MC2NAME, MC3NAME, MD3NAME
         USE PT3D_DEFN, ONLY: ZSTATIC, EMLYRS, IPVERT, VGLVSXG
         USE S_STK_EMIS

c        USE GRID_CONF ! horizontal & vertical domain specifications
c        USE STK_PRMS  ! stack groups file
c        USE STK_EMIS  ! hourly point source emissions
c        USE PMEM_DEFN ! defines emissions file(s) PM mapping
c        USE PTMET     ! defines met variables used in plume rise
         USE PTBILIN   ! defines bi-linear interpolation from grid to pt src locations
c        USE PTMAP     ! defines pt src species mapping to VDEMIS* arrays

         IMPLICIT NONE

C Includes:
c        INCLUDE SUBST_IOPARMS   ! I/O parameters definitions
c        INCLUDE SUBST_IOFDESC   ! file header data structure
!include         SUBST_IODECL    # I/O definitions and declarations
c        INCLUDE SUBST_IODECL    ! I/O definitions and declarations
c        INCLUDE SUBST_GC_EMIS   ! gas chem emis surrogate names & map table
c        INCLUDE SUBST_AE_EMIS   ! aerosol emission surrogate names & map table
c        INCLUDE SUBST_NR_EMIS   ! non-react emis surrogate names & map table
c        INCLUDE SUBST_TR_EMIS   ! tracer emis surrogate names & map table
c        INCLUDE SUBST_FILES_ID  ! file name parameters (for CTM_PT3D_DIAG)

C Arguments:
         INTEGER          N_SPC_EMIS ! total no. of model emissions species
         INTEGER          EMLAYS     ! number of emissions layers
         INTEGER          JDATE      ! Julian date (YYYYDDD)
         INTEGER          JTIME      ! time (HHMMSS)
         INTEGER          TSTEP      ! output time step
         INTEGER          NPMAX
         CHARACTER( 16 )  PT3D_FILE( NPMAX,2,9 )
         INTEGER          PT3D_NUM( NPMAX )
         INTEGER          ISPLIT( NPMAX,4 )
         LOGICAL          SUCCESS

C External Functions:
         INTEGER,         EXTERNAL :: SETUP_LOGDEV
         INTEGER,         EXTERNAL :: ENVINT
c        INTEGER,         EXTERNAL :: INDEX1, ENVINT, PROMPTFFILE
c        LOGICAL,         EXTERNAL :: ENVYN
c        CHARACTER( 16 ), EXTERNAL :: PROMPTMFILE

C Parameters:

C Local Variables:
c        CHARACTER( 16 ), SAVE :: CTM_PT3DEMIS = 'CTM_PT3DEMIS' ! env var for in-line
c                                                               ! 3d pt src emis
         CHARACTER( 16 )       :: PNAME = 'S_PT3D_INIT'  ! procedure name
c        CHARACTER( 16 )       :: VNAME    ! variable name buffer
c        CHARACTER( 16 ), SAVE, ALLOCATABLE :: STKGNAME( : ) ! stack groups file name

c        INTEGER          METCCOLS       ! cross point grid number of columns
c        INTEGER          METCROWS       ! cross point grid number of rows
c        INTEGER          METDCOLS       ! dot point grid number of columns
c        INTEGER          METDROWS       ! dot point grid number of rows
         INTEGER          IOS            ! i/o and allocate memory status
         INTEGER, SAVE :: NSTEPS = 1
c        INTEGER          NDATE          ! next timestep file Julian date (YYYYDDD)
c        INTEGER          NTIME          ! next timestep file time (HHMMSS)

         INTEGER          LOGDEV

c        INTEGER I, J, K, L, M, N, S, SRC, V  ! counters and indices

c        INTERFACE
c           SUBROUTINE DELTA_ZS( EMLYRS, MY_NSRC, SRC_MAP, STKHT, ZF, ZSTK, DDZF )
c              IMPLICIT NONE
c              INTEGER, INTENT( IN )  :: EMLYRS, MY_NSRC
c              INTEGER, INTENT( IN )  :: SRC_MAP( : )
c              REAL,    INTENT( IN )  :: STKHT( : )
c              REAL,    INTENT( IN )  :: ZF  ( EMLYRS,MY_NSRC )
c              REAL,    INTENT( OUT ) :: ZSTK( EMLYRS,MY_NSRC )
c              REAL,    INTENT( OUT ) :: DDZF( EMLYRS,MY_NSRC )
c           END SUBROUTINE DELTA_ZS
c        END INTERFACE

C-----------------------------------------------------------------------

         SUCCESS = .TRUE.

         LOGDEV = SETUP_LOGDEV()


C get scenario timing

         SDATE = ENVINT( 'LAYP_STDATE',
     &                   'Scenario Starting Date (YYYYDDD)',
     &                    2001188, IOS )

         STIME = ENVINT( 'LAYP_STTIME',
     &                   'Scenario Starting Time (HHMMSS)',
     &                    0, IOS )

         NSTEPS = ENVINT( 'LAYP_NSTEPS',
     &                    'Scenario Run Duration (integer steps)',
     &                     24, IOS )

         JDATE = SDATE
         JTIME = STIME

C Initialize stack emissions data
         print*, "BEFORE S_STK_EMIS_INIT"

         IF ( .NOT. S_STK_EMIS_INIT( JDATE, JTIME, NPMAX, 
     &                               PT3D_FILE, PT3D_NUM,
     &                               ISPLIT ) ) THEN
            XMSG = 'Could not initialize stack parameters'
            CALL M3WARN ( PNAME, JDATE, JTIME, XMSG )
            SUCCESS = .FALSE.; RETURN
         END IF



         print*,"AFTER  S_STK_EMIS_INIT =-=-=-=-=-=-=-="









         END FUNCTION S_PT3D_INIT













C-----------------------------------------------------------------------

      END MODULE S_PT3D_DEFN

