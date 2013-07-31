
C RCS file, release, date & time of last delta, author, state, [and locker]
C $Header: /project/work/rep/CCTM/src/vdiff/acm2_inline/STK_EMIS.F,v 1.1 2008/05/28 13:38:17 yoj Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      MODULE S_STK_EMIS

C-----------------------------------------------------------------------
C Function: stack emissions by source from stack emissions file and
C           3d point source emissions

C Revision History:
C     16 Jan 2007 J.Young: initial implementation

C-----------------------------------------------------------------------

      USE UDTYPES, ONLY: RARRY2, CARRY1
      USE STK_EMIS, ONLY: NVARS_FIRE, FGRPS

      IMPLICIT NONE

      INTEGER,         ALLOCATABLE, SAVE :: S_NSRC_EMIS( :,: )   ! no. of pt sources
      INTEGER,         ALLOCATABLE, SAVE :: S_NVARS_EMIS( :,: )  ! no. of species

      LOGICAL,         ALLOCATABLE, SAVE :: S_FIREFLAG( :,: )    ! which is the fire file?

      TYPE( RARRY2 ), ALLOCATABLE, SAVE :: S_STKEMIS( :,: )   ! stack emissions [mol/s]
      TYPE( CARRY1 ), ALLOCATABLE, SAVE :: S_STKSPC( :,: )    ! emitted species names
      TYPE( RARRY2 ), ALLOCATABLE, SAVE :: S_FIREVAR( :,: )   ! fire variables, 2 assumed
      TYPE( CARRY1 ), ALLOCATABLE, SAVE :: S_FIRENAM( :,: )   ! fire variable names

      INTEGER,        ALLOCATABLE, SAVE :: S_STKDATE( :,: )   ! src grp timestamp date
      INTEGER,        ALLOCATABLE, SAVE :: S_SYNCOFF( :,: )   ! src grp timestamp date offset

      CONTAINS

         FUNCTION S_STK_EMIS_INIT ( JDATE, JTIME, NPMAX,
     &                              PT3D_FILE, PT3D_NUM,
     &                              ISPLIT ) RESULT ( SUCCESS )

         IMPLICIT NONE

         INCLUDE SUBST_IOPARMS   ! I/O parameters definitions
         INCLUDE SUBST_IODECL    ! I/O definitions and declarations
         INCLUDE SUBST_IOFDESC   ! file header data structure

         INTEGER          JDATE, JTIME
         INTEGER          NPMAX
         CHARACTER( 16 )  PT3D_FILE( NPMAX,2,9 )
         INTEGER          PT3D_NUM( NPMAX )
         INTEGER          ISPLIT( NPMAX,4 )
         LOGICAL SUCCESS

         INTEGER, EXTERNAL :: SECSDIFF

         INTEGER NP
         INTEGER I, J, N, NGRPS, LEN1, LEN2, V

         INTEGER ASTAT
         CHARACTER( 96 ) :: XMSG = ' '
         CHARACTER( 16 ) :: PNAME = 'S_STK_EMIS_INIT'

C-----------------------------------------------------------------------

         SUCCESS = .TRUE.

c        NGRPS = SIZE( STKENAME )
         NGRPS = 9 ! maximum of 9 groups per sensitivity

         ALLOCATE ( S_NSRC_EMIS( NGRPS, NPMAX ), STAT = ASTAT )
         CALL CHECKMEM( ASTAT, 'S_NSRC_EMIS', PNAME )

         ALLOCATE ( S_NVARS_EMIS( NGRPS, NPMAX ), STAT = ASTAT )
         CALL CHECKMEM( ASTAT, 'NVARS_EMIS', PNAME )

         ALLOCATE ( S_STKDATE( NGRPS, NPMAX ), STAT = ASTAT )
         CALL CHECKMEM( ASTAT, 'S_STKDATE', PNAME )

         ALLOCATE ( S_SYNCOFF( NGRPS, NPMAX ), STAT = ASTAT )
         CALL CHECKMEM( ASTAT, 'S_SYNCOFF', PNAME )

         ALLOCATE ( S_STKEMIS( NGRPS, NPMAX ), STAT = ASTAT )
         CALL CHECKMEM( ASTAT, 'S_STKEMIS', PNAME )

         ALLOCATE ( S_STKSPC( NGRPS, NPMAX ), STAT = ASTAT )
         CALL CHECKMEM( ASTAT, 'S_STKSPC', PNAME )

         ALLOCATE ( S_FIREFLAG( NGRPS, NPMAX ), STAT = ASTAT )
         CALL CHECKMEM( ASTAT, 'S_FIREFLAG', PNAME )
         S_FIREFLAG = .FALSE.   ! array

         ALLOCATE ( S_FIREVAR( FGRPS, NPMAX ), STAT = ASTAT )
         CALL CHECKMEM( ASTAT, 'S_FIREVAR', PNAME )

         ALLOCATE ( S_FIRENAM( FGRPS, NPMAX ), STAT = ASTAT )
         CALL CHECKMEM( ASTAT, 'S_FIRENAM', PNAME )

         DO NP = 1, NPMAX

         IF ( ISPLIT( NP, 4 ) .NE. 1 ) CYCLE

            DO N = 1, PT3D_NUM( NP )

               IF ( .NOT. OPEN3( PT3D_FILE( NP, 2, N ), FSREAD3, PNAME ) ) THEN
                  XMSG = 'Could not open '// 
     &            TRIM( PT3D_FILE( NP, 2, N ) ) // ' file'
                  CALL M3MESG( XMSG )
                  SUCCESS = .FALSE.; RETURN
               END IF

               IF ( .NOT. DESC3( PT3D_FILE( NP, 2, N ) ) ) THEN
                 XMSG = 'Could not get ' 
     &           // TRIM( PT3D_FILE( NP, 2, N ) ) // ' file description'
                 CALL M3MESG( XMSG )
                 SUCCESS = .FALSE.; RETURN
               END IF

               S_STKDATE( N,NP ) = SDATE3D
               S_SYNCOFF( N,NP ) = SECSDIFF( JDATE,JTIME, S_STKDATE( N,NP ),JTIME )
               S_NSRC_EMIS( N,NP ) = NROWS3D
               S_NVARS_EMIS( N,NP ) = NVARS3D
               LEN1 = S_NSRC_EMIS( N,NP )
               LEN2 = S_NVARS_EMIS( N,NP )
               DO I = 1, LEN2
                  IF ( VNAME3D( I ) .EQ. 'HFLUX' ) S_FIREFLAG( N,NPMAX ) = .TRUE.
               END DO
               IF ( S_FIREFLAG( N,NP ) ) THEN   ! subtract off 'HFLUX'
                  LEN2 = LEN2 - NVARS_FIRE
                  S_NVARS_EMIS( N,NP ) = LEN2
               END IF

               S_STKEMIS( N,NP )%LEN1 = LEN1
               S_STKEMIS( N,NP )%LEN2 = LEN2
               ALLOCATE ( S_STKEMIS( N,NP )%ARRY( LEN1,LEN2 ), STAT = ASTAT )
               CALL CHECKMEM( ASTAT, 'S_STKEMIS', PNAME )
               S_STKEMIS( N,NP )%ARRY = 0.0   ! array

               S_STKSPC( N,NP )%LEN = LEN2
               ALLOCATE ( S_STKSPC( N,NP )%ARRY( LEN2 ), STAT = ASTAT )
               CALL CHECKMEM( ASTAT, 'S_STKSPC', PNAME )

               IF ( S_FIREFLAG( N,NP ) ) THEN   ! subtract and save off 'HFLUX'
                  ALLOCATE ( S_FIREVAR( 1,NP )%ARRY( LEN1,NVARS_FIRE ), STAT = ASTAT )
                  CALL CHECKMEM( ASTAT, 'FIREVAR', PNAME )
                  ALLOCATE ( S_FIRENAM( 1,NP )%ARRY( NVARS_FIRE ), STAT = ASTAT )
                  CALL CHECKMEM( ASTAT, 'FIRENAM', PNAME )
                  I = 0
                  DO J = 1, NVARS3D
C The FIRENAM_FIREVAR to match use in PT3D_DEFN
                     IF ( TRIM( VNAME3D( J ) ) .EQ. 'HFLUX' ) THEN
                        S_FIRENAM( 1,NP )%ARRY( 1 ) = VNAME3D( J )
                        CYCLE
                     END IF
                     I = I + 1
                     S_STKSPC( N,NP )%ARRY( I ) = VNAME3D( J )
                  END DO
                  S_FIRENAM( 1,NP )%LEN = NVARS_FIRE
                  S_FIREVAR( 1,NP )%LEN1 = LEN1
                  S_FIREVAR( 1,NP )%LEN2 = NVARS_FIRE
               ELSE
                  S_STKSPC( N,NP )%ARRY = VNAME3D   ! array
               END IF

            END DO ! N

         END DO ! NP

         RETURN

         END FUNCTION S_STK_EMIS_INIT

      END MODULE S_STK_EMIS
