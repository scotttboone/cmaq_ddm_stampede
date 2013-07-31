
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
C $Header: /opt/home/wt3/CMAQ4.7/models/CCTM/src/vdiff/acm2/vdiffacm2.F,v 1.1.1.1 2008/11/18 20:07:38 sjr Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C @(#)vdiffim.F 1.8 /project/mod3/CMAQ/src/vdiff/eddy/SCCS/s.vdiffim.F 25 Jul 1997 12:57:45

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#ifdef sens
      SUBROUTINE VDIFF ( CGRID, SENGRID, JDATE, JTIME, TSTEP,
     &                   NPMAX, IPT, IPARM, ISPLIT, IAMOUNT,
     &                   ILAYER, STDATE, IDATE, ITIME, IREGION,
     &                   SEN_PAR )
#else
      SUBROUTINE VDIFF ( CGRID,JDATE,JTIME,TSTEP )
#endif


C-----------------------------------------------------------------------
C Asymmetric Convective Model v2 (ACM2) -- Pleim(2006)
C Function:
C   calculates and writes dry deposition.
C   calculates vertical diffusion


C Subroutines and Functions Called:
C   INIT3, SEC2TIME, TIME2SEC, WRITE3, NEXTIME,
C   M3EXIT, EDDYX, TRI, MATRIX, PA_UPDATE_EMIS, PA_UPDATE_DDEP

C Revision History:
C   Analogous to VDIFFIM (Eddy diffusion PBL scheme)

C   31 Jan 05 J.Young: dyn alloc - establish both horizontal & vertical
C                      domain specifications in one module (GRID_CONF)
C    7 Jun 05 P.Bhave: added call to OPSSEMIS if MECHNAME='AE4';
C                      added TSTEP to RDEMIS_AE call vector
C    Aug 05 J. Pleim Update to v4.5
C
C    Jan 06 J. Pleim ACM2 implimentation

C    Jul 09 W. Tang & D. Cohan: DDM implimentation 
C-----------------------------------------------------------------------

      USE AERO_EMIS           ! inherits GRID_CONF

      USE SUBST_MODULES              ! stenex
!     USE SUBST_GLOBAL_SUM_MODULE    ! stenex

      IMPLICIT NONE

!     INCLUDE SUBST_HGRD_ID   ! horizontal dimensioning parameters
!     INCLUDE SUBST_VGRD_ID   ! vertical dimensioning parameters
      INCLUDE SUBST_RXCMMN    ! model mechanism name

      INCLUDE SUBST_GC_SPC    ! gas chemistry species table
      INCLUDE SUBST_GC_EMIS   ! gas chem emis surrogate names and map table
      INCLUDE SUBST_GC_DEPV   ! gas chem dep vel surrogate names and map table
      INCLUDE SUBST_GC_DDEP   ! gas chem dry dep species and map table
      INCLUDE SUBST_GC_DIFF   ! gas chem diffusion species and map table

      INCLUDE SUBST_AE_SPC    ! aerosol species table
!     INCLUDE SUBST_AE_EMIS   ! aerosol emis surrogate names and map table
      INCLUDE SUBST_AE_DEPV   ! aerosol dep vel surrogate names and map table
      INCLUDE SUBST_AE_DDEP   ! aerosol dry dep species and map table
      INCLUDE SUBST_AE_DIFF   ! aerosol diffusion species and map table

      INCLUDE SUBST_NR_SPC    ! non-reactive species table
      INCLUDE SUBST_NR_EMIS   ! non-react emis surrogate names and map table
      INCLUDE SUBST_NR_DEPV   ! non-react dep vel surrogate names and map table
      INCLUDE SUBST_NR_DDEP   ! non-react dry dep species and map table
      INCLUDE SUBST_NR_DIFF   ! non-react diffusion species and map table

      INCLUDE SUBST_TR_SPC    ! tracer species table
      INCLUDE SUBST_TR_EMIS   ! tracer emis surrogate names and map table
      INCLUDE SUBST_TR_DEPV   ! tracer dep vel surrogate names and map table
      INCLUDE SUBST_TR_DDEP   ! tracer dry dep species and map table
      INCLUDE SUBST_TR_DIFF   ! tracer diffusion species and map table

!     INCLUDE SUBST_EMLYRS_ID ! emissions layers parameter
#ifdef emis_chem
      INCLUDE SUBST_EMPR_CH   ! emissions processing in chem
#else
      INCLUDE SUBST_EMPR_VD   ! emissions processing in vdif 
#endif
      INCLUDE SUBST_PACTL_ID  ! PA control parameters
      INCLUDE SUBST_CONST     ! constants
      INCLUDE SUBST_FILES_ID  ! file name parameters
      INCLUDE SUBST_IOPARMS   ! I/O parameters definitions
#include      SUBST_IODECL    # I/O definitions and declarations
!     INCLUDE SUBST_COORD_ID  ! coordinate and domain definitions (req IOPARMS)

      CHARACTER( 120 ) :: XMSG = ' '

C Arguments:

!     REAL         CGRID( NCOLS,NROWS,NLAYS,* )  !  concentrations
!     REAL      :: CGRID( :,:,:,: )              !  concentrations
      REAL, POINTER :: CGRID( :,:,:,: )              !  concentrations
      INTEGER      JDATE        ! current model date, coded YYYYDDD
      INTEGER      JTIME        ! current model time, coded HHMMSS
      INTEGER      TSTEP( 2 )   ! time step vector (HHMMSS)
                                ! TSTEP(1) = local output step
                                ! TSTEP(2) = sciproc sync. step (chem)

C Parameters:

C explicit, THETA = 0, implicit, THETA = 1
      REAL, PARAMETER :: THETA = 0.5,  ! For dry deposition term
     &                   THBAR = 1.0 - THETA
      REAL THRAT  ! THBAR/THETA

      INTEGER, PARAMETER :: N_SPC_DDEP = N_GC_DDEP
     &                                 + N_AE_DDEP
     &                                 + N_NR_DDEP
     &                                 + N_TR_DDEP

C global dep vel species
      INTEGER, PARAMETER :: N_SPC_DEPV = N_GC_DEPV
     &                                 + N_AE_DEPV
     &                                 + N_NR_DEPV
     &                                 + N_TR_DEPV

C global diffusion species
      INTEGER, PARAMETER :: N_SPC_DIFF = N_GC_DIFF
     &                                 + N_AE_DIFF
     &                                 + N_NR_DIFF
     &                                 + N_TR_DIFF

C global emissions species

      INTEGER, SAVE      :: N_SPC_EMIS
!     INTEGER, PARAMETER :: N_SPC_EMIS = NEMIS
!!   &                                 + N_AE_EMIS
!    &                                 + N_NR_EMIS
!    &                                 + N_TR_EMIS

C ACM parameters

      REAL, PARAMETER :: M2PHA = 1.0E+04       ! 1 hectare = 1.0e4 m**2
      REAL, PARAMETER :: CMLMR = 1.0E+06       ! ppmV/Molar Mixing Ratio
      REAL, PARAMETER :: CNVTD = M2PHA / CMLMR / MWAIR ! combined ddep
                                                       ! conversion factor
      REAL, PARAMETER :: GPKG = 1.0E+03        ! g/Kg
      REAL, PARAMETER :: MGPG = 1.0E+06        ! micro-g/g

      REAL, PARAMETER :: CRANKP = 0.5
      REAL, PARAMETER :: CRANKQ = 1.0 - CRANKP
      REAL, PARAMETER :: KARMAN = 0.4
      REAL, PARAMETER :: EPS = 1.0E-06

!     INTEGER, PARAMETER :: IFACM2 = 1    ! 1 = acm2, 0 = acm1
      INTEGER, PARAMETER :: IFACM = 1     ! 1 = acm, 0 = no acm

C External Functions not previously declared in IODECL3.EXT:

      INTEGER, EXTERNAL :: SECSDIFF, SEC2TIME, TIME2SEC
      LOGICAL, EXTERNAL :: ENVYN

C File variables:

      REAL         RDEPVHT( NCOLS,NROWS )        ! air dens / dep vel height
      REAL         RJACM  ( NCOLS,NROWS,NLAYS )  ! reciprocal mid-layer Jacobian
      REAL         RVJACMF( NCOLS,NROWS,NLAYS )  ! 1/ mid-full layer vert Jac
      REAL         RRHOJ  ( NCOLS,NROWS,NLAYS )  ! reciprocal density X Jacobian
      REAL         DENS1  ( NCOLS,NROWS )        ! layer 1 air density
!     REAL         DEPV   ( NCOLS,NROWS,N_SPC_DEPV+1 ) ! deposition velocities
!     REAL         MDEPV  ( NCOLS,NROWS,N_SPC_DEPV+1 ) ! deposition velocities
                                                       ! X air density for all
                                                       ! but aerosol species
      REAL         DEPV  ( N_SPC_DEPV+1,NCOLS,NROWS ) ! dep vel X dens/msfx2

C Local Variables:

      CHARACTER( 16 ), SAVE :: PNAME = 'VDIFFIM'
      CHARACTER( 16 ), SAVE :: DDEP_SPC( N_SPC_DDEP + 1 )
      CHARACTER( 16 ), SAVE :: CTM_SSEMDIAG = 'CTM_SSEMDIAG'  ! env var for SSEMDIAG file
      CHARACTER( 80 ) :: VARDESC                ! environment variable description

      LOGICAL, SAVE :: FIRSTIME = .TRUE.
      INTEGER, SAVE :: WSTEP  = 0               ! local write counter
      INTEGER  STATUS                           ! ENV... status

      REAL          DX1, DX2                    ! CX x1- and x2-cell widths
      REAL, ALLOCATABLE, SAVE :: DX3F ( : )
      REAL, ALLOCATABLE, SAVE :: RDX3F( : )     ! reciprocal layer thickness
      REAL          X3M  ( NLAYS )              ! middle layer heigth
      REAL          DX3M ( NLAYS )              ! layer thickness at middle
      REAL, ALLOCATABLE, SAVE :: RDX3M( : )     ! reciprocal layer thickness
      REAL          CONVPA           ! conversion factor to pressure in Pascals
      REAL          CONVEM           ! conversion for emissions rates to Kg/s
      REAL, SAVE :: CNVTE                       ! combined conversion factor
      REAL       :: CNVTR                       ! combined conversion factor

      REAL, ALLOCATABLE, SAVE :: CNGRD( :,:,:,: )  ! cgrid replacement

      REAL, ALLOCATABLE, SAVE :: DDEP( :,:,: )  ! ddep accumulator
      REAL, ALLOCATABLE, SAVE :: DDEP_PA( :,:,: )! ddep for process analysis
      REAL          WRDD( NCOLS,NROWS )         ! ddep write buffer

      INTEGER, SAVE :: DEPV_MAP( N_SPC_DEPV+1 ) ! global depv map to CGRID
      INTEGER, SAVE :: DIFF_MAP( N_SPC_DIFF+1 ) ! global diff map to CGRID
      INTEGER, SAVE :: DF2DV   ( N_SPC_DIFF+1 ) ! map from diff spc to depv spc
      INTEGER, SAVE :: DF2EM   ( N_SPC_DIFF+1 ) ! map from diff spc to emis spc
      INTEGER, SAVE :: DD2DV   ( N_SPC_DDEP+1 ) ! map from ddep spc to depv spc
      INTEGER, SAVE :: DV2DF   ( N_SPC_DEPV )   ! map from depv spc to diff spc

      INTEGER, SAVE :: ELAYS                    ! no. of emis integration layers
                                                ! ELAYS must be .LT. NLAYS
      INTEGER          EMISLYRS                 ! no. of file emissions layers

      REAL, ALLOCATABLE, SAVE :: VDEMIS( :,:,:,: ) ! total emissions array
      REAL, ALLOCATABLE, SAVE :: EMIS_PA( :,:,:,: ) ! emis for process analysis
      REAL, ALLOCATABLE, SAVE :: VDEMIS_AE( :,:,:,: ) ! aerosol emissions
      REAL, ALLOCATABLE, SAVE :: VDEMIS_NR( :,:,:,: ) ! nonreactive gas emis
      REAL, ALLOCATABLE, SAVE :: VDEMIS_TR( :,:,:,: ) ! tracer emissions

      LOGICAL, SAVE :: EM_TRAC = .FALSE.        ! do tracer emissions?
      LOGICAL, SAVE :: SSEMDIAG                 ! flag for creating SSEMIS
                                                ! output file
      INTEGER, SAVE :: NEMIS_AE                 ! no. of aero emis species
      INTEGER, SAVE :: N_SPC_CGRID              ! no. of CGRID species
      INTEGER, SAVE :: NAESPCEMIS               ! no. of species on the PM
                                                ! emissions input file. Set in
                                                ! OPEMIS the value changes with
                                                ! the type of emissions file.

      REAL, SAVE :: DD_CONV( N_SPC_DEPV+1 )     ! ddep spc conversion factors

      REAL         DD_FAC( N_SPC_DEPV)          ! combined subexpression
      REAL           DDBF( N_SPC_DEPV)          ! secondary DDEP
      REAL           CONC( N_SPC_DIFF,NLAYS )   ! secondary CGRID expression
      REAL           EMIS( N_SPC_DIFF,NLAYS )   ! emissions subexpression
      REAL         EDDYV ( NCOLS,NROWS,NLAYS )  ! from EDYINTB
      REAL         SEDDY ( NLAYS,NCOLS,NROWS )  ! flipped EDDYV
      INTEGER      NSTEPS( NCOLS,NROWS )        ! diffusion time steps
!     REAL         DT    ( NCOLS,NROWS )        ! eddy diff. delta T
      REAL         DELT                         ! DT
      REAL         DTDENS1                      ! DT * layer 1 air density
      REAL         DTSEC                        ! model time step in seconds

C ACM Local Variables
      REAL        MBAR                          ! ACM2 mixing rate (S-1)
      REAL        HOL   ( NCOLS,NROWS )         ! PBL over Monin-Obukhov Len
      REAL        XPBL  ( NCOLS,NROWS )         ! PBL HT in gen coords
      INTEGER     LPBL  ( NCOLS,NROWS )         ! layer containing PBL HT
      LOGICAL     CONVCT( NCOLS,NROWS )         ! flag for ACM
      REAL        MEDDY
      REAL        EDDY  ( NLAYS )
      REAL        MBARKS( NLAYS )               ! by layer
      REAL        MDWN  ( NLAYS )               ! ACM down mix rate
      REAL        MFAC                          ! intermediate loop factor
      REAL        AA    ( NLAYS )               ! matrix column one
      REAL        BB    ( NLAYS )               ! diagonal
      REAL        CC    ( NLAYS )               ! subdiagonal
      REAL        EE    ( NLAYS )               ! superdiagonal
      REAL        DD    ( N_SPC_DIFF,NLAYS )    ! R.H.S
      REAL        UU    ( N_SPC_DIFF,NLAYS )    ! returned solution
      REAL        XPLUS
      REAL        XMINUS
      REAL        EFAC1 ( N_SPC_DEPV )
      REAL        EFAC2 ( N_SPC_DEPV )
      REAL        DVS1  ( N_SPC_DEPV )
      REAL        FNL
      INTEGER     NLP, NL, LCBL
      REAL        DTLIM, DTS, DTACM, RZ, DELC, LFAC1, LFAC2

      INTEGER, SAVE :: LOGDEV
 
      INTEGER      ALLOCSTAT
      INTEGER      C, R, L, S, V, N             ! loop induction variables
      INTEGER      STRT, FINI                   ! loop induction variables
      INTEGER      MDATE, MTIME, MSTEP          ! internal simulation date&time
                                                ! array in vert. mixing
!     LOGICAL, SAVE :: EDDY_STATS = .FALSE.
!     REAL DT_AVG                               ! avg eddy delta T 
!     REAL NSTP_AVG                             ! avg no. of integration steps

!     CHARACTER( 16 ) :: VNAME

#ifdef sens
c  VARIABLES SPECIFIC TO SENSITIVITY
!     REAL         SENGRID( NCOLS,NROWS,NLAYS,NPMAX,* ) ! sensitivities
      REAL, POINTER :: SENGRID( :,:,:,:,: )             ! sensitivities
      CHARACTER( 8 ), POINTER :: SEN_PAR( : )           ! names for sens parameter (eg,'_EMIS_NO2')
      REAL SEN_TEMP( NPMAX, N_SPC_DIFF, NLAYS )
      REAL, ALLOCATABLE, SAVE :: SVDEMIS   ( :,:,:,:,: ) ! total sens emis array
      REAL SVDEMIS_TEMP( NPMAX, N_SPC_DIFF, NLAYS )  ! stores SVDEMIS*DTS
c SLN - dry deposition sensitivity
      REAL, ALLOCATABLE, SAVE :: SENDDEP( :,:,:,: )     ! ddep sens accumulator 
      CHARACTER( 7 ) :: SPCNAME
      INTEGER, EXTERNAL :: TRIMLEN            ! externally defined function to trim string  
      CHARACTER( 16 ) :: PARNAME              ! full name of sensitivity paramber (eg 'NO2_EMIS_NOX')
      REAL SENDDBF( N_SPC_DEPV, NPMAX )       ! incremental SENDDEP

      REAL, ALLOCATABLE, SAVE :: SNGRD( :,:,:,:,: )  ! sengrid replacement
c SLN - dry deposition velocity sensitivity
      REAL SDEPVCON( N_SPC_DEPV, NLAYS, NPMAX )       ! dry deposition velocity sens
c SLN - eddy diffusion coefficient sensitivity
      REAL SDIFFCON( N_SPC_DIFF, NLAYS, NPMAX )       ! eddy diffusion coefficient sens
c Variables used for split emissions files
      CHARACTER( 16 ), SAVE  :: DDM3D_ES = 'DDM3D_ES'
      LOGICAL, SAVE          :: ES
      REAL, ALLOCATABLE, SAVE :: VDEMIS_A ( :,:,:,: )
      REAL, ALLOCATABLE, SAVE :: VDEMIS_B ( :,:,:,: )
      REAL, ALLOCATABLE, SAVE :: VDEMIS_M ( :,:,:,: )
      REAL, ALLOCATABLE, SAVE :: VDEMIS_P ( :,:,:,: )
      REAL, ALLOCATABLE, SAVE :: VDEMIS_N ( :,:,:,: )
! ythu 11202006 added
      REAL, ALLOCATABLE, SAVE :: VDEMIS_AE_A ( :,:,:,: )
      REAL, ALLOCATABLE, SAVE :: VDEMIS_AE_B ( :,:,:,: )
      REAL, ALLOCATABLE, SAVE :: VDEMIS_AE_M ( :,:,:,: )
      REAL, ALLOCATABLE, SAVE :: VDEMIS_AE_P ( :,:,:,: )
      REAL, ALLOCATABLE, SAVE :: VDEMIS_AE_N ( :,:,:,: )
      REAL, ALLOCATABLE, SAVE :: VDEMIS_NR_A ( :,:,:,: )
      REAL, ALLOCATABLE, SAVE :: VDEMIS_NR_B ( :,:,:,: )
      REAL, ALLOCATABLE, SAVE :: VDEMIS_NR_M ( :,:,:,: )
      REAL, ALLOCATABLE, SAVE :: VDEMIS_NR_P ( :,:,:,: )
      REAL, ALLOCATABLE, SAVE :: VDEMIS_NR_N ( :,:,:,: )
!     INTEGER, PARAMETER :: NEMIS_B = N_GC_EMIS_B
!!   &                                N_AE_EMIS_B + 
!!   &                                N_NR_EMIS_B
!     INTEGER, PARAMETER :: NEMIS_M = N_GC_EMIS_M +
!    &                                N_AE_EMIS_M
!!   &                                N_NR_EMIS_M


      INTEGER  NPMAX                    ! number of sens parameters
      INTEGER, POINTER :: IPT( : )
!      INTEGER IPT( NPMAX )              ! Type of sensitivity
                                        ! 1 = INIT; 2 = BOUN; 3 = EMIS
                                        ! 4 = HIGH; 5 = RATE; 6 = DDPV
                                        ! 7 = DIFC
      INTEGER, POINTER :: IPARM( :,: )
!      INTEGER IPARM( NPMAX, MXSPCS )    ! set to 1 if sens to that species
                                        ! set to 0 otherwise
      INTEGER, POINTER :: ISPLIT( :,: )
!      INTEGER ISPLIT( NPMAX, 6 )        ! 1 = total; 2 = area; 3 = bio
                                        ! 4 = mobile; 5 = point;
                                        ! 6 = non-road
      REAL, POINTER :: IAMOUNT( :,:,: )
!      REAL IAMOUNT ( NPMAX, MXSPCS, 25 )  ! emission rate during the hour

      INTEGER, POINTER :: ILAYER( :,: )
!      INTEGER ILAYER ( NPMAX, NLAYS )   ! 1 if that layer's emissions included, 
                                        ! 0 else
                                        ! only used if emis amount specified;
                                        ! otherwise set to 1 for all layers
      INTEGER      STDATE       ! starting date,    format YYYYDDD

      INTEGER, POINTER :: IDATE( :,: )
!      INTEGER IDATE ( NPMAX, NDAYS ) ! set to 1 if that date is included

      INTEGER, POINTER :: ITIME( :,: )
!      INTEGER ITIME ( NPMAX, 2 )        ! begin time & end time for each sens
                                        ! default is 0, 240000
                                        ! Note: to span midnight, begin time
                                        ! may be larger than endtime

      REAL, POINTER :: IREGION( :,:,:,: )
!      REAL IREGION( SENCOLS, SENROWS, NLAYS, NPMAX )


      REAL         SDD ( NPMAX, N_SPC_DIFF, NLAYS ) ! DD for sens
      REAL         SUU ( NPMAX, N_SPC_DIFF, NLAYS ) ! UU for sens

      LOGICAL TIMEFLAG
      INTEGER DATENUM, HOURNUM
      INTEGER  P                                ! loop counter for NPMAX
      INTEGER STRT_B, FINI_B, STRT_M, FINI_M
      INTEGER S2
      REAL SDELC

C END SENSITIVITY VARIABLES
#endif

      INTERFACE
         SUBROUTINE RDMET( MDATE, MTIME, RDEPVHT, RJACM, RVJACMF, RRHOJ,
     &                     DENS1 )
            IMPLICIT NONE
            INTEGER, INTENT( IN )       :: MDATE, MTIME
            REAL, INTENT( OUT )         :: RDEPVHT( :,: )
            REAL, INTENT( OUT )         :: RJACM  ( :,:,: )
            REAL, INTENT( OUT )         :: RVJACMF( :,:,: )
            REAL, INTENT( OUT )         :: RRHOJ  ( :,:,: )
            REAL, INTENT( OUT )         :: DENS1  ( :,: )
         END SUBROUTINE RDMET
         SUBROUTINE RDDEPV ( MDATE, MTIME, MSTEP, CGRID, DEPV )
            IMPLICIT NONE
            INTEGER, INTENT( IN )       :: MDATE, MTIME, MSTEP
            REAL, POINTER               :: CGRID( :,:,:,: )
            REAL, INTENT( OUT )         :: DEPV( :,:,: )
         END SUBROUTINE RDDEPV
         SUBROUTINE RDEMIS_GC ( MDATE, MTIME, EMISLYRS, NSPC_EMIS, VDEMIS )
            IMPLICIT NONE
            INTEGER, INTENT( IN )       :: MDATE, MTIME, EMISLYRS, NSPC_EMIS
            REAL, INTENT( OUT )         :: VDEMIS( :,:,:,: )
         END SUBROUTINE RDEMIS_GC
         SUBROUTINE RDEMIS_NR ( MDATE, MTIME, EMISLYRS, NSPC_EMIS, VDEMIS )
            IMPLICIT NONE
            INTEGER, INTENT( IN )       :: MDATE, MTIME, EMISLYRS, NSPC_EMIS
            REAL, INTENT( OUT )         :: VDEMIS( :,:,:,: )
         END SUBROUTINE RDEMIS_NR
         SUBROUTINE RDEMIS_TR ( MDATE, MTIME, EMISLYRS, NSPC_EMIS, VDEMIS )
            IMPLICIT NONE
            INTEGER, INTENT( IN )       :: MDATE, MTIME, EMISLYRS, NSPC_EMIS
            REAL, INTENT( OUT )         :: VDEMIS( :,:,:,: )
         END SUBROUTINE RDEMIS_TR
#ifdef sens
         SUBROUTINE S_RDEMIS_GC ( MDATE, MTIME, EMISLYRS, NEMIS,
     &                            VDEMIS, EMISFILE )
            IMPLICIT NONE
            INTEGER, INTENT( IN )        :: MDATE, MTIME, EMISLYRS, NEMIS
            REAL, INTENT( OUT )          :: VDEMIS( :,:,:,: )
            CHARACTER( * ), INTENT( IN ) :: EMISFILE
         END SUBROUTINE S_RDEMIS_GC
         SUBROUTINE S_RDEMIS_NR ( MDATE, MTIME, EMISLYRS, NEMIS,
     &                            VDEMIS, EMISFILE )
            IMPLICIT NONE
            INTEGER, INTENT( IN )        :: MDATE, MTIME, EMISLYRS, NEMIS
            REAL, INTENT( OUT )          :: VDEMIS( :,:,:,: )
            CHARACTER( * ), INTENT( IN ) :: EMISFILE
         END SUBROUTINE S_RDEMIS_NR
#endif

         SUBROUTINE PA_UPDATE_EMIS ( PNAME, VDEMIS, JDATE, JTIME, TSTEP )
            IMPLICIT NONE
            CHARACTER( * ), INTENT( IN ) :: PNAME
            REAL, INTENT( IN )           :: VDEMIS( :,:,:,: )
            INTEGER, INTENT( IN )        :: JDATE, JTIME
            INTEGER, INTENT( IN )        :: TSTEP( 2 )
         END SUBROUTINE PA_UPDATE_EMIS
         SUBROUTINE PA_UPDATE_DDEP ( PNAME, DDEP, JDATE, JTIME, TSTEP )
            IMPLICIT NONE
            CHARACTER( * ), INTENT( IN ) :: PNAME
            REAL, INTENT( IN )           :: DDEP( :,:,: )
            INTEGER, INTENT( IN )        :: JDATE, JTIME
            INTEGER, INTENT( IN )        :: TSTEP( 2 )
         END SUBROUTINE PA_UPDATE_DDEP
         SUBROUTINE CONV_CGRID ( CGRID, JDATE, JTIME, CNGRD )
            IMPLICIT NONE
            REAL, POINTER :: CGRID( :,:,:,: )
            INTEGER, INTENT( IN )        :: JDATE, JTIME
            REAL, INTENT( OUT ) :: CNGRD( :,:,:,: )
         END SUBROUTINE CONV_CGRID
         SUBROUTINE REV_CGRID ( CNGRD, JDATE, JTIME, CGRID )
            IMPLICIT NONE
            REAL, INTENT( IN ) :: CNGRD( :,:,:,: )
            INTEGER, INTENT( IN )        :: JDATE, JTIME
            REAL, POINTER :: CGRID( :,:,:,: )
         END SUBROUTINE REV_CGRID
#ifdef sens
         SUBROUTINE CONV_SENGRID ( SENGRID, JDATE, JTIME, SNGRD, NPMAX )
            IMPLICIT NONE
            REAL, POINTER :: SENGRID( :,:,:,:,: )
            INTEGER, INTENT( IN )        :: JDATE, JTIME
            REAL, INTENT( OUT ) :: SNGRD( :,:,:,:,: )
            INTEGER, INTENT( IN )        :: NPMAX
         END SUBROUTINE CONV_SENGRID
         SUBROUTINE REV_SENGRID ( SNGRD, JDATE, JTIME, SENGRID, NPMAX )
            IMPLICIT NONE
            REAL, INTENT( IN ) :: SNGRD( :,:,:,:,: )
            INTEGER, INTENT( IN )        :: JDATE, JTIME
            REAL, POINTER :: SENGRID( :,:,:,:,: )
            INTEGER, INTENT( IN )        :: NPMAX
         END SUBROUTINE REV_SENGRID
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
         SUBROUTINE OPDDEP(JDATE, JTIME, TSTEP, N_SPC_DDEP, NPMAX, SEN_PAR )
           IMPLICIT NONE
           INTEGER, INTENT( IN ) :: JDATE, JTIME, TSTEP, N_SPC_DDEP
           INTEGER, INTENT( IN ) :: NPMAX
           CHARACTER( * ), POINTER :: SEN_PAR( : )
         END SUBROUTINE OPDDEP
#endif

         SUBROUTINE EDDYX ( JDATE, JTIME, TSTEP,
     &                      EDDYV, HOL, XPBL, LPBL, CONVCT )
            IMPLICIT NONE
            INTEGER, INTENT( IN )       :: JDATE, JTIME, TSTEP
            REAL, INTENT( OUT )         :: EDDYV ( :,:,: )
            REAL, INTENT( OUT )         :: HOL   ( :,: )
            REAL, INTENT( OUT )         :: XPBL  ( :,: )
            INTEGER, INTENT( OUT )      :: LPBL  ( :,: )
            LOGICAL, INTENT( OUT )      :: CONVCT( :,: )
         END SUBROUTINE EDDYX
      END INTERFACE

C-----------------------------------------------------------------------

      IF ( FIRSTIME ) THEN

         FIRSTIME = .FALSE.
#ifdef sens
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
#endif

         LOGDEV = INIT3()

C for emissions (form COORD.EXT) .......................................

         IF ( GDTYP_GD .EQ. LATGRD3 ) THEN
            DX1 = DG2M * XCELL_GD ! in m.
            DX2 = DG2M * YCELL_GD
     &          * COS( PI180*( YORIG_GD + YCELL_GD * FLOAT( GL_NROWS/2 ))) ! in m.
         ELSE
            DX1 = XCELL_GD        ! in m.
            DX2 = YCELL_GD        ! in m.
         END IF

C create global maps

         CALL VDIFF_MAP ( DF2EM, DF2DV, DD2DV, DEPV_MAP, DIFF_MAP, DDEP_SPC,
     &                    DV2DF )

C set vertical layer definitions from COORD.EXT

         ALLOCATE ( DX3F( NLAYS ),
     &              RDX3F( NLAYS ),
     &              RDX3M( NLAYS ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating DX3F, RDX3F or RDX3M'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

         DO L = 1, NLAYS
            DX3F( L )  = X3FACE_GD( L ) - X3FACE_GD( L-1 )
            RDX3F( L ) = 1.0 / DX3F( L )
            X3M( L ) = 0.5 * ( X3FACE_GD( L ) + X3FACE_GD( L-1 ) )
         END DO

         DO L = 1, NLAYS - 1
            RDX3M( L ) = 1.0 / ( X3M( L+1 ) - X3M( L ) )
         END DO
         RDX3M( NLAYS ) = 0.0

C set molecular weights

         S = 0
         DO V = 1, N_GC_DEPV
            S = S + 1
            DD_CONV( S ) = CNVTD * GC_MOLWT( GC_DEPV_MAP( V ) )
         END DO 

         DO V = 1, N_AE_DEPV
            S = S + 1
            IF ( AE_SPC( AE_DEPV_MAP( V ) )( 1:3 ) .EQ. 'NUM' ) THEN
!              DD_CONV( S ) = M2PHA ! irrelevant, since not deposited
               DD_CONV( S ) = CNVTD * AVO * 1.0E+03    ! --> #/Ha
            ELSE IF ( AE_SPC( AE_DEPV_MAP( V ) )( 1:3 ) .EQ. 'SRF' ) THEN
!              DD_CONV( S ) = M2PHA ! irrelevant, since not deposited
               DD_CONV( S ) = M2PHA * 1.0E+03 / MWAIR  ! --> M**2/Ha
            ELSE
!              DD_CONV( S ) = M2PHA / GPKG / MGPG
               DD_CONV( S ) = CNVTD * AE_MOLWT( AE_DEPV_MAP( V ) )
            END IF
         END DO

         DO V = 1, N_NR_DEPV
            S = S + 1
            DD_CONV( S ) = CNVTD * NR_MOLWT( NR_DEPV_MAP( V ) )
         END DO

         DO V = 1, N_TR_DEPV
            S = S + 1
            DD_CONV( S ) = CNVTD * TR_MOLWT( TR_DEPV_MAP( V ) )
         END DO

C Open the met files

         CALL OPMET ( JDATE, JTIME, CONVPA )
 
C Open Emissions files

         CALL OPEMIS ( JDATE, JTIME, NEMIS, EM_TRAC, CONVEM, EMISLYRS )


         ELAYS = MIN ( EMISLYRS, NLAYS - 1 )


#ifdef sens
           IF ( ES ) THEN
              CALL S_OPEMIS( JDATE, JTIME, NEMIS, EM_TRAC, CONVEM,
     &                       EMISLYRS, EMIS_A )
              CALL S_OPEMIS( JDATE, JTIME, NEMIS, EM_TRAC, CONVEM,
     &                       EMISLYRS, EMIS_B )
              CALL S_OPEMIS( JDATE, JTIME, NEMIS, EM_TRAC, CONVEM,
     &                       EMISLYRS, EMIS_M )
              CALL S_OPEMIS( JDATE, JTIME, NEMIS, EM_TRAC, CONVEM,
     &                       EMISLYRS, EMIS_P )
              CALL S_OPEMIS( JDATE, JTIME, NEMIS, EM_TRAC, CONVEM,
     &                       EMISLYRS, EMIS_N )
           END IF

         IF ( MYPE .EQ. 0 ) CALL OPDDEP ( JDATE, JTIME, TSTEP( 1 ), N_SPC_DDEP,
     &                                    NPMAX, SEN_PAR )

#else
C Set output file characteristics based on COORD.EXT and open the dry dep file

         IF ( MYPE .EQ. 0 ) CALL OPDDEP ( JDATE, JTIME, TSTEP( 1 ), N_SPC_DDEP )
#endif


C Get sea-salt-emission diagnostic file flag 

         SSEMDIAG = .FALSE.         ! default
      VARDESC = 'Flag for writing the sea-salt-emission diagnostic file'
         SSEMDIAG = ENVYN( CTM_SSEMDIAG, VARDESC, SSEMDIAG, STATUS )
         IF ( STATUS .NE. 0 ) WRITE( LOGDEV, '(5X, A)' ) VARDESC
         IF ( STATUS .EQ. 1 ) THEN
            XMSG = 'Environment variable improperly formatted'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
         ELSE IF ( STATUS .EQ. -1 ) THEN
            XMSG =
     &          'Environment variable set, but empty ... Using default:'
            WRITE( LOGDEV, '(5X, A, I9)' ) XMSG, JTIME
         ELSE IF ( STATUS .EQ. -2 ) THEN
            XMSG = 'Environment variable not set ... Using default:'
            WRITE( LOGDEV, '(5X, A, I9)' ) XMSG, JTIME
         END IF

C Open sea-salt emission file if running the AE4 or AE5 aerosol mechanism

         IF ( ( INDEX ( MECHNAME, 'AE3' ) .EQ. 0 ) .AND. SSEMDIAG ) THEN
            IF ( MYPE .EQ. 0 ) CALL OPSSEMIS ( JDATE, JTIME, TSTEP( 1 ) )
         END IF

C Allocate and initialize dry deposition array

!        ALLOCATE ( DDEP( MY_NCOLS,MY_NROWS,N_SPC_DEPV ), STAT = ALLOCSTAT )
         ALLOCATE ( DDEP( N_SPC_DEPV,MY_NCOLS,MY_NROWS ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating DDEP'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF
            
         DDEP = 0.0

#ifdef sens
C SLN - allocate sensitivity dry deposition array
         ALLOCATE ( SENDDEP( N_SPC_DEPV,MY_NCOLS,MY_NROWS,NPMAX ), STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating SENDDEP'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF

         SENDDEP = 0.0
#endif


         IF ( N_AE_SPC .GT. 0 ) THEN
         WRITE( LOGDEV,'( /5X, A )' ) 'Aerosol Emissions Processing in '
     &                                // 'Vertical diffusion ...'
            NEMIS_AE = N_AE_EMIS  ! from AE_EMIS.EXT
         ELSE
            NEMIS_AE = 0
         END IF

         N_SPC_EMIS = NEMIS
     &              + NEMIS_AE
     &              + N_NR_EMIS
     &              + N_TR_EMIS

         ALLOCATE ( VDEMIS ( N_SPC_EMIS+1,ELAYS,MY_NCOLS,MY_NROWS ),
     &              STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'VDEMIS memory allocation failed'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

#ifdef sens
C Create SVDEMIS, the sensitivity equivalent of VDEMIS, and init to 0
      ALLOCATE ( SVDEMIS ( N_SPC_EMIS+1,ELAYS,MY_NCOLS,MY_NROWS,NPMAX ),
     &           STAT = ALLOCSTAT )
      IF ( ALLOCSTAT .NE. 0 ) THEN
         XMSG = 'SVDEMIS memory allocation failed'
         CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
      END IF

      IF ( ES ) THEN

          ALLOCATE ( VDEMIS_A( N_SPC_EMIS+1,ELAYS,MY_NCOLS,MY_NROWS ),
     &               STAT = ALLOCSTAT )
          IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'VDEMIS_A memory allocation failed'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
          END IF

          ALLOCATE ( VDEMIS_B( N_SPC_EMIS+1,ELAYS,MY_NCOLS,MY_NROWS ),
     &               STAT = ALLOCSTAT )
          IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'VDEMIS_B memory allocation failed'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
          END IF

          ALLOCATE ( VDEMIS_M( N_SPC_EMIS+1,ELAYS,MY_NCOLS,MY_NROWS ),
     &               STAT = ALLOCSTAT )
          IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'VDEMIS_M memory allocation failed'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
          END IF

          ALLOCATE ( VDEMIS_P( N_SPC_EMIS+1,ELAYS,MY_NCOLS,MY_NROWS ),
     &               STAT = ALLOCSTAT )
          IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'VDEMIS_P memory allocation failed'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
          END IF

          ALLOCATE ( VDEMIS_N( N_SPC_EMIS+1,ELAYS,MY_NCOLS,MY_NROWS ),
     &               STAT = ALLOCSTAT )
          IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'VDEMIS_N memory allocation failed'
            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
          END IF

C ythu 11202006 added 
          IF ( N_NR_EMIS .GT. 0 ) THEN
            ALLOCATE ( VDEMIS_NR_A( N_NR_EMIS,ELAYS,MY_NCOLS,MY_NROWS ),
     &              STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
              XMSG = 'VDEMIS_NR_A memory allocation failed'
              CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF

            ALLOCATE ( VDEMIS_NR_B( N_NR_EMIS,ELAYS,MY_NCOLS,MY_NROWS ),
     &              STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
              XMSG = 'VDEMIS_NR_B memory allocation failed'
              CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF

            ALLOCATE ( VDEMIS_NR_M( N_NR_EMIS,ELAYS,MY_NCOLS,MY_NROWS ),
     &              STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
              XMSG = 'VDEMIS_NR_M memory allocation failed'
              CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF

            ALLOCATE ( VDEMIS_NR_P( N_NR_EMIS,ELAYS,MY_NCOLS,MY_NROWS ),
     &              STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
              XMSG = 'VDEMIS_NR_P memory allocation failed'
              CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF

            ALLOCATE ( VDEMIS_NR_N( N_NR_EMIS,ELAYS,MY_NCOLS,MY_NROWS ),
     &              STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
              XMSG = 'VDEMIS_NR_N memory allocation failed'
              CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF
          END IF  ! N_NR_EMIS .GT. 0

C ythu added for AE Emissions Splits 

          IF (  N_AE_SPC .GT. 0 ) THEN
            ALLOCATE ( VDEMIS_AE_A( NEMIS_AE,ELAYS,MY_NCOLS,MY_NROWS ),
     &              STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
              XMSG = 'VDEMIS_AE_A memory allocation failed'
              CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF

            ALLOCATE ( VDEMIS_AE_B( NEMIS_AE,ELAYS,MY_NCOLS,MY_NROWS ),
     &              STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
              XMSG = 'VDEMIS_AE_A memory allocation failed'
              CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF

            ALLOCATE ( VDEMIS_AE_M( NEMIS_AE,ELAYS,MY_NCOLS,MY_NROWS ),
     &              STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
              XMSG = 'VDEMIS_AE_A memory allocation failed'
              CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF

            ALLOCATE ( VDEMIS_AE_P( NEMIS_AE,ELAYS,MY_NCOLS,MY_NROWS ),
     &              STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
              XMSG = 'VDEMIS_AE_P memory allocation failed'
              CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF

            ALLOCATE ( VDEMIS_AE_N( NEMIS_AE,ELAYS,MY_NCOLS,MY_NROWS ),
     &              STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
              XMSG = 'VDEMIS_AE_N memory allocation failed'
              CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF
          END IF  !  N_AE_SPC .GT. 0
C End of calls specific to split
         END IF
#endif

         IF ( N_AE_SPC .GT. 0 ) THEN
            ALLOCATE ( VDEMIS_AE( NEMIS_AE,ELAYS,MY_NCOLS,MY_NROWS ),
     &                 STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = 'VDEMIS_AE memory allocation failed'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF
         END IF

         IF ( N_NR_EMIS .GT. 0 ) THEN
            ALLOCATE ( VDEMIS_NR( N_NR_EMIS,ELAYS,MY_NCOLS,MY_NROWS ),
     &                 STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = 'VDEMIS_NR memory allocation failed'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF
         END IF

         IF ( EM_TRAC .AND. N_TR_EMIS .GT. 0 ) THEN
            ALLOCATE ( VDEMIS_TR( N_TR_EMIS,ELAYS,MY_NCOLS,MY_NROWS ),
     &                 STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = 'VDEMIS_TR memory allocation failed'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF
         END IF

         IF ( LIPR ) THEN
            ALLOCATE ( EMIS_PA( MY_NCOLS,MY_NROWS,ELAYS,N_SPC_EMIS+1 ),
     &                 STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = 'EMIS_PA memory allocation failed'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF
!           ALLOCATE ( DDEP_PA( MY_NCOLS,MY_NROWS,N_SPC_DDEP ),
            ALLOCATE ( DDEP_PA( MY_NCOLS,MY_NROWS,N_SPC_DEPV ),
     &                 STAT = ALLOCSTAT )
            IF ( ALLOCSTAT .NE. 0 ) THEN
               XMSG = 'DDEP_PA memory allocation failed'
               CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF
         END IF

C combined gas emssions conversion factor

         CNVTE = CMLMR * CONVPA * CONVEM * MWAIR / ( DX1 * DX2 )

         N_SPC_CGRID = SIZE ( CGRID,4 )

         ALLOCATE ( CNGRD( N_SPC_CGRID,NLAYS,MY_NCOLS,MY_NROWS ),
     &              STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating CNGRD'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF

         CNGRD = 0.0

#ifdef sens
         ALLOCATE ( SNGRD( NPMAX,N_SPC_CGRID,NLAYS,MY_NCOLS,MY_NROWS ),
     &              STAT = ALLOCSTAT )
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating SNGRD'
            CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
            END IF

         SNGRD = 0.0
#endif

      END IF          !  if Firstime

      MDATE = JDATE
      MTIME = JTIME
      MSTEP = TIME2SEC( TSTEP( 2 ) )
      DTSEC = FLOAT( MSTEP )
      CALL NEXTIME ( MDATE, MTIME, SEC2TIME( MSTEP / 2 ) )

#ifdef sens
C Set the date and hour counters used in sensitivity calls
      DATENUM = 1 + JDATE - STDATE
      HOURNUM = 1 + FLOOR ( REAL ( JTIME / 10000.0 ) )
#endif

C read & interpolate met data

      CALL RDMET ( MDATE, MTIME, RDEPVHT, RJACM, RVJACMF, RRHOJ, DENS1 )

C read & interpolate deposition velocities

      CALL RDDEPV ( MDATE, MTIME, TSTEP( 2 ), CGRID, DEPV )

C Initialize deposition velocities for nondeposited species to zero

      DO R = 1, MY_NROWS
         DO C = 1, MY_NCOLS
            DEPV( N_SPC_DEPV+1,C,R ) = 0.0 ! accounts for dry dep. species
         END DO                            ! names as a subset of the
      END DO                               ! vert. diffused species list

C read & interpolate emissions (create VDEMIS in the species class order)

      VDEMIS = 0.0

      CALL RDEMIS_GC ( MDATE, MTIME, ELAYS, NEMIS, VDEMIS )

C reactive gases (conversion to ppmv/s) VDEMIS in this order from RDEMIS

      STRT = 1
      FINI = NEMIS
      DO R = 1, MY_NROWS
         DO C = 1, MY_NCOLS
            DO L = 1, ELAYS
               CNVTR = CNVTE * RDX3F( L ) * RRHOJ( C,R,L )
               DO V = STRT, FINI
                  VDEMIS( V,L,C,R ) = VDEMIS( V,L,C,R ) * CNVTR
               END DO
            END DO
         END DO
      END DO

#ifdef sens
      ! Initialize SVDEMIS to 0.0
      SVDEMIS = 0.0
      VDEMIS_A = 0.0
      VDEMIS_B = 0.0
      VDEMIS_M = 0.0
      VDEMIS_P = 0.0
      VDEMIS_N = 0.0

      IF ( ES ) THEN

         CALL S_RDEMIS_GC( MDATE, MTIME, ELAYS, NEMIS,
     &                     VDEMIS_A, EMIS_A )
         CALL S_RDEMIS_GC( MDATE, MTIME, ELAYS, NEMIS,
     &                     VDEMIS_B, EMIS_B )
         CALL S_RDEMIS_GC( MDATE, MTIME, ELAYS, NEMIS,
     &                     VDEMIS_M, EMIS_M )
         CALL S_RDEMIS_GC( MDATE, MTIME, ELAYS, NEMIS,
     &                     VDEMIS_P, EMIS_P )
         CALL S_RDEMIS_GC( MDATE, MTIME, ELAYS, NEMIS,
     &                     VDEMIS_N, EMIS_N )

         STRT = 1
         FINI = NEMIS
         DO R = 1, MY_NROWS
           DO C = 1, MY_NCOLS
             DO L = 1, ELAYS
               CNVTR = CNVTE * RDX3F( L ) * RRHOJ( C,R,L )
               DO V = STRT, FINI
                 VDEMIS_A( V,L,C,R ) = VDEMIS_A( V,L,C,R ) * CNVTR
                 VDEMIS_B( V,L,C,R ) = VDEMIS_B( V,L,C,R ) * CNVTR
                 VDEMIS_M( V,L,C,R ) = VDEMIS_M( V,L,C,R ) * CNVTR
                 VDEMIS_N( V,L,C,R ) = VDEMIS_N( V,L,C,R ) * CNVTR
                 VDEMIS_P( V,L,C,R ) = VDEMIS_P( V,L,C,R ) * CNVTR
               END DO ! V
             END DO ! L
           END DO  ! C
         END DO  ! R
      END IF
#endif

C aerosol emissions - all units conversions done in RDEMIS_AE for aerosols

      IF ( N_AE_SPC .GT. 0 ) THEN
#ifdef sens
C        RDEMIS_AE in f90 module AERO_EMIS
         VDEMIS_AE_A = 0.0
         VDEMIS_AE_B = 0.0
         VDEMIS_AE_M = 0.0
         VDEMIS_AE_P = 0.0
         VDEMIS_AE_N = 0.0

         CALL RDEMIS_AE ( MDATE, MTIME, TSTEP, ELAYS, RJACM, VDEMIS, VDEMIS_AE, 
     &                    VDEMIS_A, VDEMIS_AE_A, VDEMIS_B, VDEMIS_AE_B,
     &                    VDEMIS_M, VDEMIS_AE_M, VDEMIS_P, VDEMIS_AE_P,
     &                    VDEMIS_N,VDEMIS_AE_N )

#else
         CALL RDEMIS_AE ( MDATE, MTIME, TSTEP, ELAYS, RJACM, VDEMIS, VDEMIS_AE )
#endif
      END IF

      STRT = NEMIS + 1
      FINI = NEMIS + NEMIS_AE
      DO R = 1, MY_NROWS
         DO C = 1, MY_NCOLS
            DO L = 1, ELAYS
               DO V = STRT, FINI
                  S = V + 1 - STRT
                  VDEMIS( V,L,C,R ) = VDEMIS_AE( S,L,C,R )
               END DO
            END DO
         END DO
      END DO

#ifdef sens
      IF ( ES ) THEN
C aerosol emissions - all units conversions done in S_RDEMIS_AE for aerosols
      STRT = NEMIS + 1
      FINI = NEMIS + NEMIS_AE
      DO R = 1, MY_NROWS
         DO C = 1, MY_NCOLS
            DO L = 1, ELAYS
               DO V = STRT, FINI
                  S = V + 1 - STRT
                  VDEMIS_A( V,L,C,R ) = VDEMIS_AE_A( S,L,C,R )
                  VDEMIS_B( V,L,C,R ) = VDEMIS_AE_B( S,L,C,R )
                  VDEMIS_M( V,L,C,R ) = VDEMIS_AE_M( S,L,C,R )
                  VDEMIS_P( V,L,C,R ) = VDEMIS_AE_P( S,L,C,R )
                  VDEMIS_N( V,L,C,R ) = VDEMIS_AE_N( S,L,C,R )
               END DO
            END DO
         END DO
      END DO
      END IF
#endif

C non-reactive gases (conversion to ppmv/s) VDEMIS in this order from RDEMIS

      IF ( N_NR_EMIS .GT. 0 ) THEN
         CALL RDEMIS_NR ( MDATE, MTIME, ELAYS, N_NR_EMIS, VDEMIS_NR )
      END IF

      STRT = NEMIS + NEMIS_AE + 1
      FINI = NEMIS + NEMIS_AE + N_NR_EMIS
      DO R = 1, MY_NROWS
         DO C = 1, MY_NCOLS
            DO L = 1, ELAYS
               CNVTR = CNVTE * RDX3F( L ) * RRHOJ( C,R,L )
               DO V = STRT, FINI
                  S = V + 1 - STRT
                  VDEMIS( V,L,C,R ) = VDEMIS_NR( S,L,C,R ) * CNVTR
               END DO
            END DO
         END DO
      END DO

#ifdef sens
C non-reactive gases 

      IF ( ES ) THEN

      IF ( N_NR_EMIS .GT. 0 ) THEN
        ! Initialize to 0.0
        VDEMIS_NR_A = 0.0
        VDEMIS_NR_B = 0.0
        VDEMIS_NR_M = 0.0
        VDEMIS_NR_P = 0.0
        VDEMIS_NR_N = 0.0
        CALL S_RDEMIS_NR ( MDATE,MTIME,ELAYS,N_NR_EMIS,VDEMIS_NR_A, EMIS_A )
        CALL S_RDEMIS_NR ( MDATE,MTIME,ELAYS,N_NR_EMIS,VDEMIS_NR_B, EMIS_B )
        CALL S_RDEMIS_NR ( MDATE,MTIME,ELAYS,N_NR_EMIS,VDEMIS_NR_M, EMIS_M )
        CALL S_RDEMIS_NR ( MDATE,MTIME,ELAYS,N_NR_EMIS,VDEMIS_NR_P, EMIS_P )
        CALL S_RDEMIS_NR ( MDATE,MTIME,ELAYS,N_NR_EMIS,VDEMIS_NR_N, EMIS_N )
      END IF

      STRT = NEMIS + NEMIS_AE + 1
      FINI = NEMIS + NEMIS_AE + N_NR_EMIS
      DO R = 1, MY_NROWS
         DO C = 1, MY_NCOLS
            DO L = 1, ELAYS
               CNVTR = CNVTE * RDX3F( L ) * RRHOJ( C,R,L )
               DO V = STRT, FINI
                  S = V + 1 - STRT
                  VDEMIS_A( V,L,C,R ) = VDEMIS_NR_A( S,L,C,R ) * CNVTR
                  VDEMIS_B( V,L,C,R ) = VDEMIS_NR_B( S,L,C,R ) * CNVTR
                  VDEMIS_M( V,L,C,R ) = VDEMIS_NR_M( S,L,C,R ) * CNVTR
                  VDEMIS_P( V,L,C,R ) = VDEMIS_NR_P( S,L,C,R ) * CNVTR
                  VDEMIS_N( V,L,C,R ) = VDEMIS_NR_N( S,L,C,R ) * CNVTR
               END DO
            END DO
         END DO
      END DO
C End block of calls specific to split emissions

      END IF
#endif

C tracer gases (conversion to ppmv/s)

      IF ( EM_TRAC .AND. N_TR_EMIS .GT. 0 ) THEN
         CALL RDEMIS_TR ( MDATE, MTIME, ELAYS, N_TR_EMIS, VDEMIS_TR )
      END IF

      STRT = NEMIS + NEMIS_AE + N_NR_EMIS + 1
      FINI = NEMIS + NEMIS_AE + N_NR_EMIS + N_TR_EMIS
      DO R = 1, MY_NROWS
         DO C = 1, MY_NCOLS
            DO L = 1, ELAYS
               CNVTR = CNVTE * RDX3F( L ) * RRHOJ( C,R,L )
               DO V = STRT, FINI
                  S = V + 1 - STRT
                  VDEMIS( V,L,C,R ) = VDEMIS_TR( S,L,C,R ) * CNVTR
               END DO
            END DO
         END DO
      END DO

C zero out emissions values for species not included in diffused list

      DO R = 1, MY_NROWS
         DO C = 1, MY_NCOLS
            DO L = 1, ELAYS
               VDEMIS( N_SPC_EMIS+1,L,C,R ) = 0.0  ! accounts for emissions
            END DO                                 ! species names as a subset
         END DO                                    ! of the vert. diffused
      END DO                                       ! species list

#ifdef sens
        DO P = 1, NPMAX
          DO R = 1, MY_NROWS
            DO C = 1, MY_NCOLS
              DO L = 1, ELAYS
                SVDEMIS( N_SPC_EMIS+1,L,C,R,P ) = 0.0  ! accounts for emissions
              END DO                         ! species names as a subset
            END DO                            ! of the vert. diffused
          END DO                               ! species list
        END DO

C If IPT = 3, then put the emissions into SVDEMIS for the relevant species
C Otherwise, SVDEMIS = 0
        DO P = 1, NPMAX
          IF ( IPT( P ) .EQ. 3 ) THEN
            CALL CKTIME (IDATE,JDATE,STDATE, ITIME,
     &      JTIME, TIMEFLAG , P)
            IF ( TIMEFLAG ) THEN
              DO R = 1, MY_NROWS
                DO C = 1, MY_NCOLS
                  ! Total emissions (not split as area, bio, ...)
                  IF ( ( ISPLIT( P, 1 ) .EQ. 1 ) .OR.
     &                 ( ES .EQ. .FALSE. ) ) THEN
                    DO L = 1, ELAYS
                      DO V = 1, NEMIS
                        IF ( IAMOUNT( P, GC_EMIS_MAP( V ), 25 )
     &                       .EQ. 0.0 ) THEN
                          SVDEMIS( V, L, C, R, P ) =
     &                      VDEMIS( V, L, C, R ) *
     &                      IREGION( C, R, L, P ) *
     &                      REAL( IDATE ( P, DATENUM ) ) *
     &                      REAL( IPARM( P, GC_EMIS_MAP( V ) ) )
                        ELSE
                          SVDEMIS( V, L, C, R, P ) =
     &                      IAMOUNT( P, GC_EMIS_MAP( V ), HOURNUM ) *
     &                      CNVTE * RDX3F( L ) * RRHOJ( C,R,L ) *
     &                      REAL( ILAYER ( P, L ) ) *
     &                      IREGION( C, R, L, P ) *
     &                      REAL( IDATE ( P, DATENUM ) ) *
     &                      REAL( IPARM( P, GC_EMIS_MAP( V ) ) )
                        END IF
                      END DO
                    END DO
c aerosol emissions
                    S = 0
                    STRT = NEMIS + 1
                    FINI = NEMIS + NEMIS_AE
                    DO V = STRT, FINI
                    S = S + 1
                      DO L = 1, ELAYS
                        IF ( IAMOUNT( P,
     &                         N_GC_SPC + 1 + AE_EMIS_MAP( S ),
     &                         25 ) .EQ. 0.0 ) THEN
                          SVDEMIS( V, L, C, R, P ) =
     &                      VDEMIS( V, L, C, R ) *
     &                      IREGION( C, R, L, P ) *
     &                      REAL( IDATE ( P, DATENUM ) ) *
     &                      REAL( IPARM( P,
     &                        N_GC_SPC + 1 + AE_EMIS_MAP( S ) ))
                        ELSE
                          SVDEMIS( V, L, C, R, P ) =
     &                      IAMOUNT( P,
     &                        N_GC_SPC +1 + AE_EMIS_MAP( S ),
     &                        HOURNUM ) *
     &                      CNVTE * RDX3F( L ) * RRHOJ( C,R,L ) *
     &                      REAL( ILAYER ( P, L ) ) *
     &                      IREGION( C, R, L, P ) *
     &                      REAL( IDATE ( P, DATENUM ) ) *
     &                      REAL( IPARM( P,
     &                        N_GC_SPC+1 + AE_EMIS_MAP( S ) ) )
                        END IF
                      END DO
                    END DO
c non-reactive emissions (NH3)
                    S = 0
                    STRT = NEMIS + NEMIS_AE + 1
                    FINI = NEMIS + NEMIS_AE + N_NR_EMIS
                    DO V = STRT, FINI
                    S = S + 1
                      DO L = 1, ELAYS
                        IF ( IAMOUNT( P,
     &                         N_GC_SPC + 1+N_AE_SPC + NR_EMIS_MAP( S ),
     &                         25 ) .EQ. 0.0 ) THEN
                          SVDEMIS( V, L, C, R, P ) =
     &                      VDEMIS( V, L, C, R ) *
     &                      IREGION( C, R, L, P ) *
     &                      REAL( IDATE ( P, DATENUM ) ) *
     &                      REAL( IPARM( P,
     &                        N_GC_SPC +1+N_AE_SPC + NR_EMIS_MAP( S ) ))
                        ELSE
                          SVDEMIS( V, L, C, R, P ) =
     &                      IAMOUNT( P,
     &                        N_GC_SPC +1+ N_AE_SPC + NR_EMIS_MAP( S ),
     &                        HOURNUM ) *
     &                      CNVTE * RDX3F( L ) * RRHOJ( C,R,L ) *
     &                      REAL( ILAYER ( P, L ) ) *
     &                      IREGION( C, R, L, P ) *
     &                      REAL( IDATE ( P, DATENUM ) ) *
     &                      REAL( IPARM( P,
     &                        N_GC_SPC+1+N_AE_SPC + NR_EMIS_MAP( S ) ) )
                        END IF
                      END DO
                    END DO
                  ELSE  ! Emissions are split
C ythu 11202006 added
C gas emissions
                    DO L = 1, ELAYS
                      DO V = 1, NEMIS
                          SVDEMIS( V, L, C, R, P ) =
     &                      VDEMIS_A( V, L, C, R ) *
     &                      REAL( ISPLIT( P, 2 ) ) *
     &                      IREGION( C, R, L, P ) *
     &                      REAL( IDATE ( P, DATENUM ) ) *
     &                      REAL( IPARM( P, GC_EMIS_MAP( V ) ) )
     &                     +VDEMIS_B( V, L, C, R ) *
     &                      REAL( ISPLIT( P, 3 ) ) *
     &                      IREGION( C, R, L, P ) *
     &                      REAL( IDATE ( P, DATENUM ) ) *
     &                      REAL( IPARM( P, GC_EMIS_MAP( V ) ) )
     &                     +VDEMIS_M( V, L, C, R ) *
     &                      REAL( ISPLIT( P, 4 ) ) *
     &                      IREGION( C, R, L, P ) *
     &                      REAL( IDATE ( P, DATENUM ) ) *
     &                      REAL( IPARM( P, GC_EMIS_MAP( V ) ) )
     &                     +VDEMIS_P( V, L, C, R ) *
     &                      REAL( ISPLIT( P, 5 ) ) *
     &                      IREGION( C, R, L, P ) *
     &                      REAL( IDATE ( P, DATENUM ) ) *
     &                      REAL( IPARM( P, GC_EMIS_MAP( V ) ) )
     &                     +VDEMIS_N( V, L, C, R ) *
     &                      REAL( ISPLIT( P, 6 ) ) *
     &                      IREGION( C, R, L, P ) *
     &                      REAL( IDATE ( P, DATENUM ) ) *
     &                      REAL( IPARM( P, GC_EMIS_MAP( V ) ) )
                      END DO
                    END DO
c aerosol emissions
                    S = 0
                    STRT = NEMIS + 1
                    FINI = NEMIS + NEMIS_AE
                    DO V = STRT, FINI
                    S = S + 1
                      DO L = 1, ELAYS
                          SVDEMIS( V, L, C, R, P ) =
     &                      VDEMIS_A( V, L, C, R ) *
     &                      REAL( ISPLIT( P, 2 ) ) *
     &                      IREGION( C, R, L, P ) *
     &                      REAL( IDATE ( P, DATENUM ) ) *
     &                      REAL( IPARM( P,
     &                        N_GC_SPC + 1 + AE_EMIS_MAP( S ) ))
     &                     +VDEMIS_B( V, L, C, R ) *
     &                      REAL( ISPLIT( P, 3 ) ) *
     &                      IREGION( C, R, L, P ) *
     &                      REAL( IDATE ( P, DATENUM ) ) *
     &                      REAL( IPARM( P,
     &                        N_GC_SPC + 1 + AE_EMIS_MAP( S ) ))
     &                     +VDEMIS_M( V, L, C, R ) *
     &                      REAL( ISPLIT( P, 4 ) ) *
     &                      IREGION( C, R, L, P ) *
     &                      REAL( IDATE ( P, DATENUM ) ) *
     &                      REAL( IPARM( P,
     &                        N_GC_SPC + 1 + AE_EMIS_MAP( S ) ))
     &                     +VDEMIS_P( V, L, C, R ) *
     &                      REAL( ISPLIT( P, 5 ) ) *
     &                      IREGION( C, R, L, P ) *
     &                      REAL( IDATE ( P, DATENUM ) ) *
     &                      REAL( IPARM( P,
     &                        N_GC_SPC + 1 + AE_EMIS_MAP( S ) ))
     &                     +VDEMIS_N( V, L, C, R ) *
     &                      REAL( ISPLIT( P, 6 ) ) *
     &                      IREGION( C, R, L, P ) *
     &                      REAL( IDATE ( P, DATENUM ) ) *
     &                      REAL( IPARM( P,
     &                        N_GC_SPC + 1 + AE_EMIS_MAP( S ) ))
                      END DO
                    END DO
c non-reactive emissions (NH3)
                    S = 0
                    STRT = NEMIS + NEMIS_AE + 1
                    FINI = NEMIS + NEMIS_AE + N_NR_EMIS
                    DO V = STRT, FINI
                    S = S + 1
                      DO L = 1, ELAYS
                          SVDEMIS( V, L, C, R, P ) =
     &                      VDEMIS_A( V, L, C, R ) *
     &                      REAL( ISPLIT( P, 2 ) ) *
     &                      IREGION( C, R, L, P ) *
     &                      REAL( IDATE ( P, DATENUM ) ) *
     &                      REAL( IPARM( P,
     &                        N_GC_SPC +1+N_AE_SPC + NR_EMIS_MAP( S ) ))
     &                     +VDEMIS_B( V, L, C, R ) *
     &                      REAL( ISPLIT( P, 3 ) ) *
     &                      IREGION( C, R, L, P ) *
     &                      REAL( IDATE ( P, DATENUM ) ) *
     &                      REAL( IPARM( P,
     &                        N_GC_SPC +1+N_AE_SPC + NR_EMIS_MAP( S ) ))
     &                     +VDEMIS_M( V, L, C, R ) *
     &                      REAL( ISPLIT( P, 4 ) ) *
     &                      IREGION( C, R, L, P ) *
     &                      REAL( IDATE ( P, DATENUM ) ) *
     &                      REAL( IPARM( P,
     &                        N_GC_SPC +1+N_AE_SPC + NR_EMIS_MAP( S ) ))
     &                     +VDEMIS_P( V, L, C, R ) *
     &                      REAL( ISPLIT( P, 5 ) ) *
     &                      IREGION( C, R, L, P ) *
     &                      REAL( IDATE ( P, DATENUM ) ) *
     &                      REAL( IPARM( P,
     &                        N_GC_SPC +1+N_AE_SPC + NR_EMIS_MAP( S ) ))
     &                     +VDEMIS_N( V, L, C, R ) *
     &                      REAL( ISPLIT( P, 6 ) ) *
     &                      IREGION( C, R, L, P ) *
     &                      REAL( IDATE ( P, DATENUM ) ) *
     &                      REAL( IPARM( P,
     &                        N_GC_SPC +1+N_AE_SPC + NR_EMIS_MAP( S ) ))
                      END DO
                    END DO
                  END IF   ! if split emissions or total
                END DO     ! end loop over C
              END DO       ! end loop over R
            END IF         ! TIMEFLAG
          END IF           ! if IPT = 3
        END DO             ! end loop over P
C END SENSITIVITY CALLS
#endif

      IF ( LIPR ) THEN
         DO S = 1, N_SPC_EMIS+1
            DO L = 1, ELAYS
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     EMIS_PA( C,R,L,S ) = VDEMIS( S,L,C,R )
                  END DO
               END DO
            END DO
         END DO
         CALL PA_UPDATE_EMIS ( 'VDIF', EMIS_PA, JDATE, JTIME, TSTEP )
      END IF

      CALL EDDYX ( JDATE, JTIME, TSTEP( 2 ),
     &             EDDYV, HOL, XPBL, LPBL, CONVCT )

C EDDYV returned = Kz, where Kz is in m**2/sec

      DO R = 1, MY_NROWS
         DO C = 1, MY_NCOLS
            DO L = 1, NLAYS
               SEDDY( L,C,R ) = EDDYV( C,R,L )
     &                        * RVJACMF( C,R,L ) * RDX3M( L ) ! * DELT
            END DO
         END DO
      END DO

      IF ( IFACM .EQ. 0 ) CONVCT = .FALSE.   ! no ACM

C Convert non-molar mixing ratio species and re-order CGRID

      CALL CONV_CGRID ( CGRID, MDATE, MTIME, CNGRD )
#ifdef sens

C Convert non-molar mixing ratio species and re-order SENGRID

      CALL CONV_SENGRID ( SENGRID, MDATE, MTIME, SNGRD, NPMAX )

#endif


      IF ( WSTEP .EQ. 0 ) DDEP = 0.0

C ------------------------------------------- Row, Col LOOPS -----------

      DO 345 R = 1, MY_NROWS
      DO 344 C = 1, MY_NCOLS

C ACM insert

         DTLIM = DTSEC

C Note: DT has been moved from EDDY to here, dt = .75 dzf dzh / Kz

         DO L = 1, NLAYS - 1
            DTLIM = MIN( DTLIM, 0.75 / ( SEDDY( L,C,R ) * RDX3F( L ) ) )
         END DO
         MBARKS = 0.0
         MDWN = 0.0

C New couple ACM & EDDY ------------------------------------------------

         MBAR = 0.0

         IF ( CONVCT( C,R ) ) THEN   ! Do ACM for this column
            LCBL = LPBL( C,R )
            MEDDY = SEDDY( 1,C,R ) / ( XPBL( C,R ) - X3FACE_GD( 1 ) )
            FNL = 1.0 / ( 1.0 + ( ( KARMAN / ( -HOL( C,R ) ) ) ** 0.3333 )
     &                / ( 0.72 * KARMAN ) )

!           IF ( FNL .GT. 1.0 ) WRITE( LOGDEV,* ) ' FNL= ', FNL

            MBAR = MEDDY * FNL
            DO L = 1, LCBL - 1
               SEDDY( L,C,R ) = SEDDY( L,C,R  ) * ( 1.0 - FNL )
            END DO

            IF ( MBAR .LT. EPS ) THEN
               WRITE( LOGDEV,* ) ' EDDYV, MBAR, FNL, HOL = ',
     &                             EDDYV( C,R,1 ), MBAR, FNL, HOL( C,R )
               CONVCT( C,R ) = .FALSE.
               LCBL = 1
               XMSG = '*** ACM fails ***'
               CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT2 )
            END IF

            IF ( ( FNL .LE. 0.0 ) .OR.   ! never gonna happen for CONVCT
     &           ( LCBL .GE. NLAYS-1 ) .OR.    ! .GT. never gonna happen
     &           ( HOL( C,R ) .GT. -0.00001 ) )   ! never gonna happen
     &         WRITE( LOGDEV,1015 ) LCBL, MBAR, FNL, EDDYV( C,R,1 ),
     &                              SEDDY( 1,C,R ), HOL( C,R )
1015           FORMAT( ' LCBL, MBAR, FNL, SEDDY1, HOL:', I3, 1X, 5(1PE13.5) )

            DO L = 1, LCBL - 1
               MBARKS( L ) = MBAR
               MDWN( L ) = MBAR * ( XPBL( C,R ) - X3FACE_GD( L-1 ) )
     &                   * RDX3F( L )
            END DO

            MBARKS( LCBL ) = MBAR * ( XPBL( C,R ) - X3FACE_GD( LCBL-1 ) )
     &                     * RDX3F( LCBL )
            MDWN( LCBL ) = MBARKS( LCBL )

C Modify Timestep for ACM

            RZ     = ( X3FACE_GD( LCBL ) - X3FACE_GD( 1 ) ) * RDX3F( 1 )
            DTACM  = 1.0 / ( MBAR * RZ )
            DTLIM  = MIN( 0.75 * DTACM, DTLIM )
         ELSE
            LCBL = 1
         END IF

C-----------------------------------------------------------------------

         NLP = INT( DTSEC / DTLIM + 0.99 )
         DTS = DTSEC / NLP
         DTDENS1 = DTS * DENS1( C,R )

!         IF ( R .EQ. MY_NROWS / 2 .AND. C .EQ. MY_NCOLS / 2 )
!     &      WRITE( LOGDEV,1021 ) CONVCT( C,R ), DTS, EDDYV( C,R,1 ),
!     &                           MBAR, FNL
!1021        FORMAT( ' CONVCT, DTS, EDDYV, MBAR, FNL: ',
!     &              L3, 1X, 4(1PE13.5) )

C End ACM insert

         DO L = 1, NLAYS
            DO V = 1, N_SPC_DIFF
               CONC( V,L ) = CNGRD( DIFF_MAP( V ),L,C,R )
!              CONC( V,L ) = CGRID( C,R,L,DIFF_MAP( V ) )
            END DO
         END DO

         DO V = 1, N_SPC_DEPV
!           DDBF( V ) = DDEP( C,R,V )
            DDBF( V ) = DDEP( V,C,R )
            DD_FAC( V ) = DTDENS1 * DD_CONV( V ) * DEPV( V,C,R )
         END DO
#ifdef sens
            DO V = 1, N_SPC_DEPV
               DO P = 1, NPMAX
                  SENDDBF( V, P ) = SENDDEP( V, C, R, P )
               ENDDO
            END DO
#endif

         EMIS = 0.0
         DO L = 1, ELAYS
            DO V = 1, N_SPC_DIFF
               EMIS( V,L ) = VDEMIS( DF2EM( V ),L,C,R ) * DTS
            END DO
         END DO
#ifdef sens
           ! Store sengrid and svdemis of a single col and row
           ! into smaller arrays for computational efficiency
           SEN_TEMP = 0.0
           SVDEMIS_TEMP = 0.0
           DO L = 1, NLAYS
             DO V = 1, N_SPC_DIFF
               S = DIFF_MAP( V )
               DO P = 1, NPMAX
!                SEN_TEMP( P, V, L ) = SENGRID( C, R, L, P, S )
                 SEN_TEMP( P, V, L ) = SNGRD( P, S, L, C, R )
               END DO
             END DO
           END DO

           DO L = 1, ELAYS
             DO V = 1, N_SPC_DIFF
               S2 = DF2EM( V )
               DO P = 1, NPMAX
                 SVDEMIS_TEMP( P, V, L ) =
     &               SVDEMIS( S2, L, C, R, P ) * DTS
               END DO
             END DO
           END DO
#endif


C-----------------------------------------------------------------------

         DO L = 1, NLAYS
            EDDY( L ) = SEDDY( L,C,R )
         END DO

         DO V = 1, N_SPC_DEPV
         EFAC1( V ) = EXP( -DEPV( V,C,R ) * RDEPVHT( C,R ) * THBAR * DTS )
         EFAC2( V ) = EXP( -DEPV( V,C,R ) * RDEPVHT( C,R ) * THETA * DTS )
         DVS1( V ) = DEPV( V,C,R ) * DTS * RDEPVHT( C,R )
         END DO

         DO 301 NL = 1, NLP      ! loop over sub time

            DO V = 1, N_SPC_DEPV
               DDBF( V ) = DDBF( V )
     &                   + THBAR * DD_FAC( V ) * CONC( DV2DF( V ),1 )
               CONC( DV2DF( V ),1 ) = EFAC1( V ) * CONC( DV2DF( V ),1 )
            END DO
#ifdef sens
            DO V = 1, N_SPC_DEPV
              DO P = 1, NPMAX
               SENDDBF( V, P ) = SENDDBF( V, P ) + THBAR * DD_FAC( V )
     &                      * SEN_TEMP( P, DV2DF( V ), 1 )
               SEN_TEMP( P,DV2DF( V ),1 ) = EFAC1( V ) *  
     &                                      SEN_TEMP( P,DV2DF( V ),1 )
              ENDDO
            END DO
#endif

#ifdef sens
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Begin code of sensitivity specific to dry deposition velocity & eddy diffusion
C coefficient 
C Added by W.Tang & D.Cohan 07102009
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C If IPT = 6, then put the dry deposition into SDEPVCON for the relevant species
C Otherwise, SDEPVCON = 0
        SDEPVCON = 0.0                  !Initialize SDEPVCON
        DO P = 1, NPMAX
          IF ( IPT( P ) .EQ. 6 ) THEN
            CALL CKTIME (IDATE,JDATE,STDATE, ITIME,
     &      JTIME, TIMEFLAG , P)
            IF ( TIMEFLAG ) THEN
                    L = 1
c gas species deposition
                    DO V = 1, N_GC_DEPV
                        SDEPVCON( V, L, P ) =
     &                      DVS1( V ) * CONC( DV2DF( V ),L ) *
     &                      IREGION( C, R, L, P ) *
     &                      REAL( IPARM( P, GC_DEPV_MAP( V ) ) )
                    END DO
c aerosol deposition
                    S = 0
                    STRT = N_GC_DEPV + 1
                    FINI = N_GC_DEPV + N_AE_DEPV
                    DO V = STRT, FINI
                    S = S + 1
                          SDEPVCON( V, L, P ) =
     &                      DVS1( V ) * CONC( DV2DF( V ),L ) *
     &                      IREGION( C, R, L, P ) *
     &                      REAL( IPARM( P, N_GC_SPC + 1 +
     &                                      AE_DEPV_MAP( S ) ) )
                    END DO
c non-reactive deposition
                    S = 0
                    STRT = N_GC_DEPV + N_AE_DEPV + 1
                    FINI = N_GC_DEPV + N_AE_DEPV + N_NR_DEPV
                    DO V = STRT, FINI
                    S = S + 1
                          SDEPVCON( V, L, P ) =
     &                      DVS1( V ) * CONC( DV2DF( V ), L ) *
     &                      IREGION( C, R, L, P ) *
     &                      REAL( IPARM( P, N_GC_SPC + 1 + N_AE_SPC +
     &                                      NR_DEPV_MAP( S ) ) )
                    END DO
            END IF !TIMEFLAG
          END IF !IPT=6
        END DO !P

C If IPT = 7, then put the vertical diffusion into SDIFFCON
C Otherwise, SDIFFCON = 0
        SDIFFCON = 0.0                  !Initialize SDIFFCON
        DO P = 1, NPMAX
          IF ( IPT( P ) .EQ. 7 ) THEN
            CALL CKTIME (IDATE,JDATE,STDATE, ITIME,
     &      JTIME, TIMEFLAG , P)
            IF ( TIMEFLAG ) THEN
             DO L = 1, NLAYS
              DO V = 1 , N_SPC_DIFF
               IF ( L.GE.2.AND.L.LE.( NLAYS-1 ) ) THEN
                 SDIFFCON( V,  L,  P ) =
     &                          ( ( EDDY( L ) * RDX3F( L ) * DTS )
     &                         * ( CONC( V,L+1 ) - CONC( V,L ) )
     &                         -  ( EDDY( L-1 ) * RDX3F( L ) * DTS )
     &                         * ( CONC( V,L ) - CONC( V,L-1 ) ) )
     &                         * IREGION( C, R, L, P )
     &                         * REAL( ILAYER( P,L ) )
                 IF ( CONVCT(C,R).AND.L.LE.LCBL ) THEN
                    SDIFFCON(V,L,P) = SDIFFCON(V,L,P) +  
     &                              ( DTS * ( MBARKS( L ) * CONC( V,1 )
     &                              - MDWN( L ) * CONC( V,L )
     &                              + (DX3F( L+1 ) * RDX3F( L ) 
     &                              * MDWN( L+1 )) * CONC( V,L+1 ) ) )
     &                              * IREGION( C, R, L, P )
     &                              * REAL( ILAYER( P,L ) )
                 ELSE
                    SDIFFCON(V,L,P) = SDIFFCON(V,L,P)
                 END IF

               ELSE IF ( L.EQ.1 ) THEN
                 SDIFFCON( V,  L,  P ) =
     &                          ( EDDY( L ) * RDX3F( L ) * DTS ) *
     &                          ( CONC( V,L+1 ) - CONC( V,L ) ) *
     &                          IREGION( C, R, L, P ) *
     &                          REAL( ILAYER( P,L ) )
                 IF ( CONVCT(C,R) ) THEN 
                    SDIFFCON(V,L,P) = SDIFFCON(V,L,P) +
     &                                (-( MBARKS( 1 ) * ( XPBL( C,R ) 
     &                                - X3FACE_GD( 1 ) ) * RDX3F( 1 ) 
     &                                * DTS)* CONC( V,1 )
     &                                + (MDWN( 2 ) * DX3F( 2 ) 
     &                                * RDX3F( 1 ) * DTS)* CONC( V,2 ))
     &                                * IREGION( C, R, L, P )
     &                                * REAL( ILAYER( P,L ) )
                 ELSE
                    SDIFFCON(V,L,P) = SDIFFCON(V,L,P)
                 END IF

               ELSE
                 SDIFFCON( V,  L,  P ) =
     &                          ( EDDY( L-1 ) * RDX3F( L ) * DTS ) *
     &                          ( CONC( V,L ) - CONC( V,L-1 ) ) *
     &                          IREGION( C, R, L, P ) *
     &                          REAL( ILAYER( P,L ) )
 
               END IF !LAYER
              END DO !SPECIES 
             END DO !LAYER
            END IF !TIMEFLAG
          END IF !IPT
         END DO !P                          
#endif

C Init variables for use below

           DO L = 1, NLAYS
               AA( L ) = 0.0
               BB( L ) = 0.0
               CC( L ) = 0.0
               EE( L ) = 0.0
             DO V = 1, N_SPC_DIFF
                  DD( V,L ) = 0.0
                  UU( V,L ) = 0.0
#ifdef sens
               DO P = 1, NPMAX
                  SDD( P,V,L) = 0.0
                  SUU( P,V,L) = 0.0
               END DO
#endif
             END DO
           END DO

C Compute tendency of CBL concentrations - semi-implicit solution
C Define arrays A,B,E which make up MATRIX and D which is RHS

            IF ( CONVCT( C,R ) ) THEN
               DO L = 2, LCBL
                  AA( L )   = -CRANKP * MBARKS( L ) * DTS
                  BB( L )   = 1.0 + CRANKP * MDWN( L ) * DTS
                  EE( L-1 ) = -CRANKP * MDWN( L ) * DTS * DX3F( L )
     &                      * RDX3F( L-1 )
                  MFAC = DX3F( L+1 ) * RDX3F( L ) * MDWN( L+1 )
                  DO V = 1, N_SPC_DIFF
                     DELC = DTS * ( MBARKS( L ) * CONC( V,1 )
     &                              - MDWN( L ) * CONC( V,L )
     &                              + MFAC      * CONC( V,L+1 ) )
                     DD( V,L ) = CONC( V,L ) + CRANKQ * DELC
#ifdef sens
                  DO P = 1, NPMAX
                     SDELC = DTS * ( MBARKS( L ) * SEN_TEMP( P,V,1 )
     &                              - MDWN( L ) * SEN_TEMP( P,V,L )
     &                              + MFAC      * SEN_TEMP( P,V,L+1 ) )
                     SDD( P,V,L ) = SEN_TEMP( P,V,L ) + CRANKQ * SDELC
                  END DO
#endif
                  END DO
               END DO
            END IF

            AA( 2 ) = AA( 2 ) - EDDY( 1 ) * CRANKP * RDX3F( 2 ) * DTS
            EE( 1 ) = EE( 1 ) - EDDY( 1 ) * CRANKP * RDX3F( 1 ) * DTS

            DO L = 2, NLAYS
               IF ( L .GT. LCBL ) THEN
                  BB( L ) = 1.0
                  DO V = 1, N_SPC_DIFF
                     DD( V,L ) = CONC( V,L )
#ifdef sens
                  DO P = 1, NPMAX
                     SDD( P,V,L ) = SEN_TEMP( P,V,L )
                  END DO
#endif
                  END DO
               END IF
               XPLUS  = EDDY( L )   * RDX3F( L ) * DTS
               XMINUS = EDDY( L-1 ) * RDX3F( L ) * DTS
               BB( L ) = BB( L ) + ( XPLUS + XMINUS ) * CRANKP
               CC( L ) = - XMINUS * CRANKP
               EE( L ) = EE( L ) - XPLUS * CRANKP
               IF ( L .EQ. NLAYS ) THEN
                  DO V = 1, N_SPC_DIFF
                     DD( V,L ) = DD( V,L )
     &                         - CRANKQ * XMINUS
     &                         * ( CONC( V,L ) - CONC( V,L-1 ) )
#ifdef sens                 
                  DO P = 1, NPMAX
                     SDD( P,V,L ) = SDD( P,V,L )
     &                     - CRANKQ * XMINUS
     &                     * ( SEN_TEMP( P,V,L ) - SEN_TEMP( P,V,L-1 ) )
     &                     - SDIFFCON ( V, L ,P )
                  END DO
#endif
                  END DO
               ELSE
                  LFAC1 = CRANKQ * XPLUS
                  LFAC2 = CRANKQ * XMINUS
                  DO V = 1, N_SPC_DIFF
                     DD( V,L ) = DD( V,L )
     &                         + LFAC1
     &                         * ( CONC( V,L+1 ) - CONC( V,L ) )
     &                         - LFAC2
     &                         * ( CONC( V,L ) - CONC( V,L-1 ) )
                 IF ( L .LE. ELAYS ) DD( V,L ) = DD( V,L ) + EMIS( V,L )
#ifdef sens
                  DO P = 1, NPMAX
                     SDD( P,V,L ) = SDD( P,V,L )
     &                     + LFAC1
     &                     * ( SEN_TEMP( P,V,L+1 ) - SEN_TEMP( P,V,L ) )
     &                     - LFAC2
     &                     * ( SEN_TEMP( P,V,L ) - SEN_TEMP( P,V,L-1 ) )
     &                     + SDIFFCON ( V, L, P )
                 IF ( L .LE. ELAYS )  SDD( P,V,L ) = SDD( P,V,L ) 
     &                                           + SVDEMIS_TEMP( P,V,L )
                  END DO
#endif
                  END DO
               END IF
            END DO

            BB( 1 ) = 1.0
            DO V = 1, N_SPC_DIFF
               DD( V,1 ) = CONC( V,1 )
#ifdef sens
            DO P = 1, NPMAX
               SDD( P,V,1 ) = SEN_TEMP( P,V,1 )
            END DO
#endif
            END DO

            IF ( CONVCT( C,R ) ) THEN
               LFAC1 = ( XPBL( C,R ) - X3FACE_GD( 1 ) ) * RDX3F( 1 ) * DTS
               LFAC2 = CRANKQ * MDWN( 2 ) * DX3F( 2 )
     &               * RDX3F( 1 ) * DTS
               BB( 1 ) = BB( 1 ) + CRANKP * MBARKS( 1 ) * LFAC1
               LFAC1 = CRANKQ * MBARKS( 1 ) * LFAC1
               DO V = 1, N_SPC_DIFF
                  DD( V,1 ) = DD( V,1 )
     &                      - LFAC1 * CONC( V,1 )
     &                      + LFAC2 * CONC( V,2 ) ! net mixing above
#ifdef sens
               DO P = 1, NPMAX
                  SDD( P,V,1 ) = SDD( P,V,1 )
     &                      - LFAC1 * SEN_TEMP( P,V,1 )
     &                      + LFAC2 * SEN_TEMP( P,V,2 )
               END DO
#endif
               END DO
            END IF

            BB( 1 ) = BB( 1 ) + CRANKP * EDDY( 1 ) * RDX3F( 1 ) * DTS
            LFAC1 = CRANKQ * EDDY( 1 ) * RDX3F( 1 ) * DTS
            DO V = 1, N_SPC_DIFF
               DD( V,1 ) = DD( V,1 )
     &                   + LFAC1 * ( CONC( V,2 ) - CONC( V,1 ) )
     &                   + EMIS( V,1 )
#ifdef sens
            DO P = 1, NPMAX
               SDD( P,V,1 ) = SDD( P,V,1 )
     &               + LFAC1 * ( SEN_TEMP( P,V,2 ) - SEN_TEMP( P,V,1 ) )
     &               + SVDEMIS_TEMP( P,V,1 )
     &               - SDEPVCON ( DF2DV(V), 1, P )
     &               + SDIFFCON ( V, 1, P )
            END DO
#endif 
            END DO

C Subroutine MATRIX then solves for U if ACM2, else TRI solves for U

            IF ( CONVCT( C,R ) ) THEN
               CALL MATRIX ( AA, BB, CC, DD, EE, UU )
            ELSE
               CALL TRI ( CC, BB, EE, DD, UU )
            END IF
#ifdef sens
            IF ( CONVCT( C,R ) ) THEN
               CALL SMATRIX ( AA, BB, CC, SDD, EE, SUU, NPMAX )
            ELSE
               CALL STRI ( CC, BB, EE, SDD, SUU, NPMAX )
            END IF
#endif

C Load into CGRID
            DO L = 1, NLAYS
               DO V = 1, N_SPC_DIFF
                  CONC( V,L ) = UU( V,L )
#ifdef sens
               DO P = 1, NPMAX
                  SEN_TEMP( P,V,L ) = SUU ( P,V,L )
               END DO
#endif
               END DO
           END DO
            
            DO V = 1, N_SPC_DEPV
               DDBF( V ) = DDBF( V )
     &                   + THETA * DD_FAC( V ) * CONC( DV2DF( V ),1 )
               CONC( DV2DF( V ),1 ) = EFAC2( V ) * CONC( DV2DF( V ),1 )
            END DO
#ifdef sens
            DO V = 1, N_SPC_DEPV
              DO P = 1, NPMAX
               SENDDBF( V, P ) = SENDDBF( V, P ) + THETA * DD_FAC( V )
     &                      * SEN_TEMP( P, DV2DF( V ), 1 )
               SEN_TEMP( P,DV2DF( V ),1 ) = EFAC2( V ) *
     &                                      SEN_TEMP( P,DV2DF( V ),1 )
              ENDDO
            END DO
#endif

301      CONTINUE                 ! end sub time loop

         DO L = 1, NLAYS
            DO V = 1, N_SPC_DIFF
               CNGRD( DIFF_MAP( V ),L,C,R ) = CONC( V,L )
!              CGRID( C,R,L,DIFF_MAP( V ) ) = CONC( V,L )
            END DO
         END DO

            DO V = 1, N_SPC_DEPV
!              DDEP( C,R,V ) = DDBF( V )
               DDEP( V,C,R ) = DDBF( V )
               END DO 

#ifdef sens
           ! Upload SEN_TEMP back in to SNGRD
           DO L = 1, NLAYS
             DO V = 1, N_SPC_DIFF
               S = DIFF_MAP( V )
               DO P = 1, NPMAX
                 SNGRD( P, S, L, C, R ) = SEN_TEMP( P, V, L )
               END DO
             END DO
           END DO
C SLN - update dry deposition sensitivities
           DO V = 1, N_SPC_DEPV
              DO P = 1, NPMAX
                 SENDDEP( V, C, R, P ) = SENDDBF( V, P )
              END DO
           END DO
#endif


344   CONTINUE         !  end loop on col C
345   CONTINUE         !  end loop on row R

C Revert non-molar mixing ratio species and re-order CGRID

      CALL REV_CGRID ( CNGRD, MDATE, MTIME, CGRID )
#ifdef sens

C Revert non-molar mixing ratio species and re-order SENGRID

      CALL REV_SENGRID ( SNGRD, MDATE, MTIME, SENGRID, NPMAX )

#endif

C If last call this hour:  write accumulated depositions:

      WSTEP = WSTEP + TIME2SEC( TSTEP( 2 ) )
      IF ( WSTEP .GE. TIME2SEC( TSTEP( 1 ) ) ) THEN
         MDATE = JDATE
         MTIME = JTIME
         CALL NEXTIME( MDATE, MTIME, TSTEP( 2 ) )
         WSTEP = 0

         DO V = 1, N_SPC_DDEP
            S = DD2DV( V )
            DO R = 1, MY_NROWS
               DO C = 1, MY_NCOLS
                  WRDD( C,R ) = DDEP( S,C,R )
               END DO
            END DO

            IF ( .NOT. WRITE3( CTM_DRY_DEP_1, DDEP_SPC( V ),
     &                 MDATE, MTIME, WRDD ) ) THEN
               XMSG = 'Could not write ' // CTM_DRY_DEP_1 // ' file'
               CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
            END IF

         END DO

         WRITE( LOGDEV, '( /5X, 3( A, :, 1X ), I8, ":", I6.6 )' )
     &         'Timestep written to', CTM_DRY_DEP_1,
     &         'for date and time', MDATE, MTIME

#ifdef sens
C SLN - write out dry deposition sensitivities
         DO V = 1, N_SPC_DDEP
            S = DD2DV( V )
            SPCNAME = DDEP_SPC( V )                              ! <- these lines have to
            DO P = 1, NPMAX
               PARNAME = SPCNAME( 1:TRIMLEN( SPCNAME )) //       ! <- match what is in
     &                              '_' // SEN_PAR( P )          ! <- opddep.F
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     WRDD( C,R ) = SENDDEP( S, C, R, P )
                  ENDDO
               ENDDO

               IF ( .NOT. WRITE3( CTM_SDRYDEP_1, PARNAME,
     &                 MDATE, MTIME, WRDD ) ) THEN
                  XMSG = 'Could not write ' // CTM_SDRYDEP_1 // ' file'
                  CALL M3EXIT( PNAME, MDATE, MTIME, XMSG, XSTAT1 )
               END IF

            ENDDO
         ENDDO
         WRITE( LOGDEV, '( /5X, 3( A, :, 1X ), I8, ":", I6.6 )' )
     &         'Timestep written to', CTM_SDRYDEP_1,
     &         'for date and time', MDATE, MTIME

#endif

         IF ( LIPR ) THEN
!           DO V = 1, N_SPC_DDEP
            DO V = 1, N_SPC_DEPV
               DO R = 1, MY_NROWS
                  DO C = 1, MY_NCOLS
                     DDEP_PA( C,R,V ) = DDEP( V,C,R )
                  END DO
               END DO
            END DO
            CALL PA_UPDATE_DDEP ( 'VDIF', DDEP_PA, JDATE, JTIME, TSTEP )
         END IF

C re-set dry deposition array to zero

         DDEP = 0.0
#ifdef sens
         SENDDEP = 0.0
#endif

      END IF

      RETURN
      END
