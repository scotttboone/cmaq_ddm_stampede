
C RCS file, release, date & time of last delta, author, state, [and locker]
C $Header: /project/work/rep/CCTM/src/vdiff/acm2_inline/preplm.f,v 1.2 2008/08/30 13:32:48 yoj Exp $

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      SUBROUTINE PREPLM( FIREFLG, EMLAYS, HMIX, HTS, PSFC, TS, DDZF, QV,
     &                   TA, UW, VW, ZH, ZF, PRES, LSTK, LPBL, TSTK, 
     &                   WSTK, DTHDZ, WSPD )

C-----------------------------------------------------------------------
C Description:
C   Computes the values needed for the PLMRIS subroutine from the 
C   meteorology data.
 
C Preconditions:
C   Interpolated (to the location of a source) meteorology data as input,
C   vertical grid structure.
 
C Subroutines and Functions Called:
C    I/O API 
 
C Revision History:
C    Copied from preplm.f v 1.2 in DAQM-V2 Emissions Preprocessor by
C        M. Houyoux 3/99
 
C-----------------------------------------------------------------------
C Modified from:
 
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling System
C File: @(#)$Id: preplm.f,v 1.2 2008/08/30 13:32:48 yoj Exp $
C COPYRIGHT (C) 2002, MCNC Environmental Modeling Center
C All Rights Reserved
C See file COPYRIGHT for conditions of use.
C Environmental Modeling Center
C MCNC
C P.O. Box 12889
C Research Triangle Park, NC  27709-2889
C smoke@emc.mcnc.org
C Pathname: $Source: /project/work/rep/CCTM/src/vdiff/acm2_inline/preplm.f,v $
C Last updated: $Date: 2008/08/30 13:32:48 $ 
C-----------------------------------------------------------------------
 
      IMPLICIT NONE
 
C Includes:

C External Functions:
      REAL,    EXTERNAL :: POLY

C Parameters:
      INTEGER, PARAMETER :: DEG = 3       ! degree of interpolationg polynomial
      REAL,    PARAMETER :: CTOK = 273.15 ! conversion from deg. C to deg. K

C Arguments:
C (Note: All met parms are per-source)
      LOGICAL, INTENT( IN )  :: FIREFLG         ! .true. => processing fire source
      INTEGER, INTENT( IN )  :: EMLAYS          ! no. emissions layers
      REAL,    INTENT( IN )  :: HMIX            ! mixing height
      REAL,    INTENT( IN )  :: HTS             ! stack height
      REAL,    INTENT( IN )  :: PSFC            ! surface pressure
      REAL,    INTENT( IN )  :: TS              ! surface temperature
      REAL,    INTENT( IN )  :: DDZF( EMLAYS )  ! 1/( zf(l) - zf(l-1) )
      REAL,    INTENT( IN )  :: QV  ( EMLAYS )  ! mixing ratio
      REAL,    INTENT( IN )  :: TA  ( EMLAYS )  ! absolute temperature
      REAL,    INTENT( IN )  :: UW  ( EMLAYS )  ! x-direction winds
      REAL,    INTENT( IN )  :: VW  ( EMLAYS )  ! y-direction winds
      REAL,    INTENT( IN )  :: ZH  ( EMLAYS )  ! layer center  height [m]
      REAL,    INTENT( IN )  :: ZF  ( EMLAYS )  ! layer surface height [m]
      REAL,    INTENT( IN )  :: PRES( EMLAYS+1 ) ! pres at full layer hts (mod by YOJ)
      INTEGER, INTENT( OUT ) :: LSTK            ! first L: ZF(L) > STKHT
      INTEGER, INTENT( OUT ) :: LPBL            ! first L: ZF(L) > mixing layer
      REAL,    INTENT( OUT ) :: TSTK            ! temperature @ top of stack [K]
      REAL,    INTENT( OUT ) :: WSTK            ! wind speed @ top of stack [m/s]
      REAL,    INTENT( OUT ) :: DTHDZ( EMLAYS ) ! potential temp. grad. 
      REAL,    INTENT( OUT ) :: WSPD ( EMLAYS ) ! wind speed [m/s]

C Local Variables:
      INTEGER      L, M
      REAL         ES
      REAL         QSFC
      REAL         TVSFC
      REAL         THETG
      REAL         THV1
      REAL         THVK
      REAL         TV( EMLAYS )   ! Virtual temperature
      REAL         TF( EMLAYS )   ! Full-layer height temperatures
      REAL         P, Q
      REAL         DZZ
      REAL         DELZ

C-----------------------------------------------------------------------

C Convert pressure to millibars from pascals, compute wind speed,
C and virtual temperature

      DO L = 1, EMLAYS
         P = UW( L )
         Q = VW( L )
         WSPD( L ) = SQRT( P * P + Q * Q )
         TV( L ) = TA( L ) * ( 1.0 + 0.622 * ( QV( L ) / ( 1.0 + QV( L ) ) ) )
      END DO

      ES    = 6.1078 * EXP( 5384.21 / CTOK - 5384.21 / TS )
      QSFC  = 0.622 * ES / ( PSFC - ES )
      TVSFC = TS * ( 1.0 + 0.6077 * QSFC )
      THETG = TVSFC * ( 1000.0 / PSFC ) ** 0.286
      IF ( HMIX .LE. ZF( 1 ) ) LPBL = 1
      IF ( HTS .LE. ZF( 1 ) ) LSTK = 1

C Interpolate the virtual temperatures at the full-layer face heights (at ZFs)
      DO L = 1, EMLAYS - 1
         DELZ = ZH( L+1 ) - ZH( L )
         TF( L ) = TV( L ) + ( TV( L+1 ) - TV( L ) ) * ( ZF( L ) - ZH( L ) ) / DELZ
      END DO
      L = EMLAYS
      DELZ = ZH( L ) - ZH( L-1 )
      TF( L ) = TV( L ) + ( TV( L ) - TV( L-1 ) ) * ( ZF( L ) - ZH( L ) ) / DELZ

      THV1  = TF( 1 ) * ( 1000.0 / PRES( 2 ) ) ** 0.286
!     DTHDZ( 1 ) = ( THV1 - THETG ) / ZF( 1 )

      DO L = 2, EMLAYS
 
         IF ( HMIX .GT. ZF( L-1 ) ) LPBL = L
         IF ( HTS .GT. ZF( L-1 ) ) LSTK = L
 
         THVK = TF( L ) * ( 1000.0 / PRES( L+1 ) ) ** 0.286
         DTHDZ( L ) = DDZF( L ) * ( THVK - THV1 )
         THV1 = THVK
 
      END DO

C Set the 1st level vertical THETV gradient to the 2nd layer value -
C overrides the layer 1 gradient determined above
      DTHDZ( 1 ) = DTHDZ( 2 )

      IF ( .NOT. FIREFLG ) THEN
C Interpolate ambient temp. and windspeed to top of stack using DEG deg. polynomial
         M    = MAX( 1, LSTK - DEG - 1 )
         TSTK =      POLY( HTS, ZH( M ), TA( M ), DEG )
         WSTK = MAX( POLY( HTS, ZH( M ), WSPD( M ), DEG ), 0.1 )
      ELSE
         TSTK = TS
         WSTK = WSPD( 1 )
      END IF

      RETURN

      END SUBROUTINE PREPLM
