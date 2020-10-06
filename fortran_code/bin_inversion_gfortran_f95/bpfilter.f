CC FILE: bpfilter.f
      subroutine BandPassFilter
     2 (data1,NP,DT,n_but,i_pass,lowcorner,highcorner)
c Apply a simple band pass filter to data
      real*4 data1(*)
      real*4 lowcorner,highcorner
Cf2py intent(out) data1
Cf2py intent(in) data1,NP,DT,n_but,i_pass,lowcorner,highcorner
      if (lowcorner.gt.1.e-10) then
         if (highcorner.lt.2.0) then
            call IIRFILT( DATA1, NP, 'BUTTER  ',n_but,
     +                   'BP      ', lowcorner, highcorner, DT, i_pass)
         elseif (highcorner.ge.2.0) then
            call IIRFILT( DATA1, NP, 'BUTTER  ',n_but,
     +                   'HP      ', lowcorner, highcorner, DT, i_pass)
         endif
      elseif (lowcorner.le.1.e-10) then
         call IIRFILT( DATA1, NP, 'BUTTER  ',n_but,
     +                 'LP      ', lowcorner, highcorner, DT, i_pass)
      endif
      return
      end

C  Copyright 1990  Regents of the University of California
C
C
C  Author:  Dave Harris
C
C           Lawrence Livermore National Laboratory
C           L-205
C           P.O. Box 808
C           Livermore, CA  94550
C           USA
C
C           (415) 423-0617
C                                                               APPLY
C  Subroutine to apply an iir filter to a DATA1 sequence.
C    The filter is assumed to be stored as second order sections.
C    Filtering is in-place.
C    Zero-phase (forward and reverse) is an option.
C
C  Input Arguments:
C  ----------------
C
C    DATA1                           Array containing DATA1
C
C    NSAMPS                         Number of DATA1 samples
C
C    ZP                             Logical variable, true for
C                                     zero phase filtering, false
C                                     for single pass filtering
C
C    SN                             Numerator polynomials for second
C                                     order sections.
C
C    SD                             Denominator polynomials for second
C                                     order sections.
C
C    NSECTS                         Number of second-order sections
C
C  Output Arguments:
C  -----------------
C
C    DATA1                          DATA1 array (same as input)
C
C
      SUBROUTINE APPLY( DATA1, NSAMPS, ZP, SN, SD, NSECTS )
C
        REAL*4 SN(1), SD(1)
        REAL*4 OUTPUT,DATA1(*)
        LOGICAL ZP
C
        JPTR = 1
        DO    1 J = 1, NSECTS
C
          X1 = 0.0
          X2 = 0.0
          Y1 = 0.0
          Y2 = 0.0
          B0 = SN(JPTR)
          B1 = SN(JPTR+1)
          B2 = SN(JPTR+2)
          A1 = SD(JPTR+1)
          A2 = SD(JPTR+2)
C
          DO    2 I = 1, NSAMPS
C
            OUTPUT = B0*DATA1(I) + B1*X1 + B2*X2
            OUTPUT = OUTPUT - ( A1*Y1 + A2*Y2 )
            Y2 = Y1
            Y1 = OUTPUT
            X2 = X1
            X1 = DATA1(I)
            DATA1(I) = OUTPUT
C
    2     CONTINUE
C
          JPTR = JPTR + 3
C
    1   CONTINUE
C
        IF (   ZP ) THEN
C
          JPTR = 1
          DO    3 J = 1, NSECTS
C
            X1 = 0.0
            X2 = 0.0
            Y1 = 0.0
            Y2 = 0.0
            B0 = SN(JPTR)
            B1 = SN(JPTR+1)
            B2 = SN(JPTR+2)
            A1 = SD(JPTR+1)
            A2 = SD(JPTR+2)
C
            DO    4 I = NSAMPS, 1, -1
C
              OUTPUT = B0*DATA1(I) + B1*X1 + B2*X2
              OUTPUT = OUTPUT - ( A1*Y1 + A2*Y2 )
              Y2 = Y1
              Y1 = OUTPUT
              X2 = X1
              X1 = DATA1(I)
              DATA1(I) = OUTPUT
C
    4       CONTINUE
C
            JPTR = JPTR + 3
C
    3     CONTINUE
C
        END IF
C
      RETURN
      END
C
C
C BEROOTS -- SUBROUTINE TO RETURN BESSEL POLES FOR
C   NORMALIZED LOWPASS FILTER
C
C LAST MODIFIED:  April 15, 1992. Changed P and RTYPE to adjustable
C                 array by using an "*" rather than a "1".
C
C  OUTPUT ARGUMENTS:
C  -----------------
C      P              COMPLEX ARRAY CONTAINING POLES
C                       CONTAINS ONLY ONE FROM EACH
C                       COMPLEX CONJUGATE PAIR, AND
C                       ALL REAL POLES
C
C      RTYPE          CHARACTER ARRAY INDICATING 2ND ORDER SECTION
C                       TYPE:
C                         (SP)  SINGLE REAL POLE
C                         (CP)  COMPLEX CONJUGATE POLE PAIR
C                         (CPZ) COMPLEX CONJUGATE POLE-ZERO PAIRS
C
C      DCVALUE        MAGNITUDE OF FILTER AT ZERO FREQUENCY
C
C      NSECTS         NUMBER OF SECOND ORDER SECTIONS
C
C  INPUT ARGUMENTS:
C  ----------------
C
C      IORD           DESIRED FILTER ORDER
C
C
      SUBROUTINE BEROOTS( P, RTYPE, DCVALUE, NSECTS, IORD )
C
        COMPLEX P(*)
        INTEGER NSECTS, IORD
        CHARACTER*3 RTYPE(*)
C       !-- LUDAN TEST IORD
C
        IF (   IORD .EQ. 1 ) THEN
C
          P(1) = CMPLX( -1.0, 0.0 )
          RTYPE(1) = 'SP'
C
        ELSE IF (  IORD .EQ. 2 ) THEN
C
          P(1) = CMPLX( -1.1016013,  0.6360098 )
          RTYPE(1) = 'CP'
C
        ELSE IF (  IORD .EQ. 3 ) THEN
C
          P(1) = CMPLX( -1.0474091, 0.9992645 )
          RTYPE(1) = 'CP'
          P(2) = CMPLX( -1.3226758, 0.0 )
          RTYPE(2) = 'SP'
C
        ELSE IF (  IORD .EQ. 4 ) THEN
C
          P(1) = CMPLX( -0.9952088,  1.2571058 )
          RTYPE(1) = 'CP'
          P(2) = CMPLX( -1.3700679, 0.4102497 )
          RTYPE(2) = 'CP'
C
        ELSE IF (  IORD .EQ. 5 ) THEN
C
          P(1) = CMPLX( -0.9576766,  1.4711244 )
          RTYPE(1) = 'CP'
          P(2) = CMPLX( -1.3808774,  0.7179096 )
          RTYPE(2) = 'CP'
          P(3) = CMPLX( -1.5023160, 0.0 )
          RTYPE(3) = 'SP'
C
        ELSE IF (  IORD .EQ. 6 ) THEN
C
          P(1) = CMPLX( -0.9306565,  1.6618633 )
          RTYPE(1) = 'CP'
          P(2) = CMPLX( -1.3818581,  0.9714719 )
          RTYPE(2) = 'CP'
          P(3) = CMPLX( -1.5714904,  0.3208964 )
          RTYPE(3) = 'CP'
C
        ELSE IF (  IORD .EQ. 7 ) THEN
C
          P(1) = CMPLX( -0.9098678,  1.8364514 )
          RTYPE(1) = 'CP'
          P(2) = CMPLX( -1.3789032,  1.1915667 )
          RTYPE(2) = 'CP'
          P(3) = CMPLX( -1.6120388,  0.5892445 )
          RTYPE(3) = 'CP'
          P(4) = CMPLX( -1.6843682, 0.0 )
          RTYPE(4) = 'SP'
C
        ELSE IF (  IORD .EQ. 8 ) THEN
C
          P(1) = CMPLX( -0.8928710,  1.9983286 )
          RTYPE(1) = 'CP'
          P(2) = CMPLX( -1.3738431,  1.3883585 )
          RTYPE(2) = 'CP'
          P(3) = CMPLX( -1.6369417,  0.8227968 )
          RTYPE(3) = 'CP'
          P(4) = CMPLX( -1.7574108,  0.2728679 )
          RTYPE(4) = 'CP'
C
        END IF
C
         NSECTS = IORD - IORD/2
C        LUDAN TEST NSECTS

C
        DCVALUE = 1.0
C
C  DONE
C
      RETURN
      END
C
C
C  Transforms an analog filter to a digital filter via the bilinear transformati
C    Assumes both are stored as second order sections.  The transformation is
C    done in-place.
C
C  Input Arguments:
C  ----------------
C
C    SN                   Array containing numerator polynomial coefficients for
C                           second order sections.  Packed head-to-tail.
C
C    SD                   Array containing denominator polynomial coefficients f
C                           second order sections.  Packed head-to-tail.
C
C    NSECTS               Number of second order sections.
C
C
      SUBROUTINE BILIN2( SN, SD, NSECTS )
C
        REAL*4 SN(1), SD(1)
C
        IPTR = 1
        DO    1 I = 1, NSECTS
C
          A0 = SD(IPTR)
          A1 = SD(IPTR+1)
          A2 = SD(IPTR+2)
C
          SCALE = A2 + A1 + A0
          SD(IPTR)   = 1.
          SD(IPTR+1) = (2.*(A0 - A2)) / SCALE
          SD(IPTR+2) = (A2 - A1 + A0) / SCALE
C
          A0 = SN(IPTR)
          A1 = SN(IPTR+1)
          A2 = SN(IPTR+2)
C
          SN(IPTR)   = (A2 + A1 + A0) / SCALE
          SN(IPTR+1) = (2.*(A0 - A2)) / SCALE
          SN(IPTR+2) = (A2 - A1 + A0) / SCALE
C
          IPTR = IPTR + 3
C
    1   CONTINUE
C
      RETURN
      END
C
C
C BUROOTS -- SUBROUTINE TO COMPUTE BUTTERWORTH POLES FOR
C   NORMALIZED LOWPASS FILTER
C
C LAST MODIFIED:  SEPTEMBER 7, 1990
C
C  OUTPUT ARGUMENTS:
C  -----------------
C      P              COMPLEX ARRAY CONTAINING POLES
C                       CONTAINS ONLY ONE FROM EACH
C                       COMPLEX CONJUGATE PAIR, AND
C                       ALL REAL POLES
C
C      RTYPE          CHARACTER ARRAY INDICATING 2ND ORDER SECTION
C                       TYPE:
C                         (SP)  SINGLE REAL POLE
C                         (CP)  COMPLEX CONJUGATE POLE PAIR
C                         (CPZ) COMPLEX CONJUGATE POLE-ZERO PAIRS
C
C      DCVALUE        MAGNITUDE OF FILTER AT ZERO FREQUENCY
C
C      NSECTS         NUMBER OF SECOND ORDER SECTIONS
C
C  INPUT ARGUMENTS:
C  ----------------
C
C      IORD           DESIRED FILTER ORDER
C
C
      SUBROUTINE BUROOTS( P, RTYPE, DCVALUE, NSECTS, IORD )
C
        COMPLEX P(1)
        INTEGER HALF
        CHARACTER*3 RTYPE(1)
C
        PI=3.14159265
C
        HALF = IORD/2
C
C TEST FOR ODD ORDER, AND ADD POLE AT -1
C
        NSECTS = 0
        IF (    2*HALF .LT. IORD ) THEN
          P(1) = CMPLX( -1., 0. )
          RTYPE(1) = 'SP'
          NSECTS = 1
        END IF
C
        DO    1  K = 1, HALF
          ANGLE = PI * ( .5 + FLOAT(2*K-1) / FLOAT(2*IORD) )
          NSECTS = NSECTS + 1
          P(NSECTS) = CMPLX( COS(ANGLE), SIN(ANGLE) )
          RTYPE(NSECTS) = 'CP'
    1   CONTINUE
C
        DCVALUE = 1.0
C
      RETURN
      END



C
C  Subroutine to alter the cutoff of a filter.  Assumes that the
C    filter is structured as second order sections.  Changes
C    the cutoffs of normalized lowpass and highpass filters through
C    a simple polynomial transformation.
C
C  Input Arguments:
C  ----------------
C
C    F                       New cutoff frequency
C
C  Input/Output Arguments:
C  -----------------------
C
C    SN                      Numerator polynomials for second order
C                              sections.
C
C    SD                      Denominator polynomials for second order
C                              sections.
C
C    NSECTS                  Number of second order sectionsects
C
C
      SUBROUTINE CUTOFFS( SN, SD, NSECTS, F )
C
        REAL*4 SN(1), SD(1)
C
        SCALE = 2.*3.14159265*F
C
        IPTR = 1
        DO    1 I = 1, NSECTS
C
          SN( IPTR + 1 ) = SN( IPTR + 1 ) / SCALE
          SN( IPTR + 2 ) = SN( IPTR + 2 ) / (SCALE*SCALE)
          SD( IPTR + 1 ) = SD( IPTR + 1 ) / SCALE
          SD( IPTR + 2 ) = SD( IPTR + 2 ) / (SCALE*SCALE)
          IPTR = IPTR + 3
C
    1   CONTINUE
C
      RETURN
      END
C
C
C  Subroutine to design IIR digital filters from analog prototypes.
C
C  Input Arguments:
C  ----------------
C
C    IORD                Filter order (10 MAXIMUM)
C
C    TYPE                Character*2 variable containing filter type
C                          LOWPASS (LP)
C                          HIGHPASS (HP)
C                          BANDPASS (BP)
C                          BANDREJECT (BR)
C
C    APROTO              Character*2 variable designating analog prototype
C                          Butterworth (BU)
C                          Bessel (BE)

C    FL                  Low-frequency cutoff
C
C    FH                  High-frequency cutoff
C
C    TS                  Sampling interval (in seconds)
C
C  Output Arguments:
C  -----------------
C
C    SN                  Array containing numerator coefficients of
C                        second-order sections packed head-to-tail.
C
C    SD                  Array containing denominator coefficients
C                        of second-order sections packed head-to-tail.
C
C    NSECTS              Number of second-order sections.
C
C
      SUBROUTINE DESIGN( IORD, TYPE, APROTO,
     &                   FL, FH, TS, SN, SD, NSECTS )
C
        COMPLEX P(10), Z(10)
        CHARACTER*2 TYPE, APROTO
        CHARACTER*3 STYPE(10)
        REAL*4 SN(1), SD(1)
C
C  Analog prototype selection
C
        IF (     APROTO .EQ. 'BU' ) THEN
C
          CALL BUROOTS( P, STYPE, DCVALUE, NSECTS, IORD )
C
        ELSE IF (    APROTO .EQ. 'BE' ) THEN
C
          CALL BEROOTS( P, STYPE, DCVALUE, NSECTS, IORD )
C
        END IF
C
C  Analog mapping selection
C
        IF (     TYPE .EQ. 'BP' ) THEN
C
          FLW = WARP( FL*TS/2., 2. )
          FHW = WARP( FH*TS/2., 2. )
          CALL LPTBP( P, Z, STYPE, DCVALUE, NSECTS, FLW, FHW, SN, SD )
C
        ELSE IF (   TYPE .EQ. 'BR' ) THEN
C
          FLW = WARP( FL*TS/2., 2. )
          FHW = WARP( FH*TS/2., 2. )
          CALL LPTBR( P, Z, STYPE, DCVALUE, NSECTS, FLW, FHW, SN, SD )
C
        ELSE IF (    TYPE .EQ. 'LP' ) THEN
C
          FHW = WARP( FH*TS/2., 2. )
          CALL LP( P, Z, STYPE, DCVALUE, NSECTS, SN, SD )
          CALL CUTOFFS( SN, SD, NSECTS, FHW )
C
        ELSE IF (    TYPE .EQ. 'HP' ) THEN
C
          FLW = WARP( FL*TS/2., 2. )
          CALL LPTHP( P, Z, STYPE, DCVALUE, NSECTS, SN, SD )
          CALL CUTOFFS( SN, SD, NSECTS, FLW )
C
        END IF
C
C  Bilinear analog to digital transformation
C
        CALL BILIN2( SN, SD, NSECTS )
C
      RETURN
      END








C  ARGUMENTS:
C  ----------
C
C    DATA1           REAL ARRAY CONTAINING SEQUENCE TO BE FILTERED
C                     ORIGINAL DATA1 DESTROYED, REPLACED BY FILTERED DATA1
C
C    NSAMPS         NUMBER OF SAMPLES IN DATA1
C
C
C    APROTO         CHARACTER*8 VARIABLE, CONTAINS TYPE OF ANALOG
C                     PROTOTYPE FILTER
C                     '(BU)TTER  ' -- BUTTERWORTH FILTER
C                     '(BE)SSEL  ' -- BESSEL FILTER
C
C    IORD           ORDER (#POLES) OF ANALOG PROTOTYPE
C                   NOT TO EXCEED 10 IN THIS CONFIGURATION.  4 - 5
C                   SHOULD BE AMPLE.
C
C    TYPE           CHARACTER*8 VARIABLE CONTAINING FILTER TYPE
C                     'LP' -- LOW PASS
C                     'HP' -- HIGH PASS
C                     'BP' -- BAND PASS
C                     'BR' -- BAND REJECT
C
C    FLO            LOW FREQUENCY CUTOFF OF FILTER (HERTZ)
C                   IGNORED IF TYPE = 'LP'
C
C    FHI            HIGH FREQUENCY CUTOFF OF FILTER (HERTZ)
C                   IGNORED IF TYPE = 'HP'
C
C    TS             SAMPLING INTERVAL (SECONDS)
C
C    PASSES           INTEGER VARIABLE CONTAINING THE NUMBER OF PASSES
C                   1 -- FORWARD FILTERING ONLY
C                   2 -- FORWARD AND REVERSE (I.E. ZERO PHASE) FILTERING
C
C
C  SUBPROGRAMS REFERENCED:  BILIN2, BUROOTS, WARP, CUTOFFS, LPTHP, LPTBP,
C    LP, LPTBR, BEROOTS, DESIGN, APPLY
C
      SUBROUTINE IIRFILT( DATA1, NSAMPS, APROTO, IORD,
     +                   TYPE, FLO, FHI, TS, PASSES )
C

        REAL*4 DATA1(*)
        CHARACTER*8 TYPE, APROTO
        INTEGER NSAMPS, PASSES, IORD
        REAL*4  FLO, FHI, TS, SN(30), SD(30)
        LOGICAL ZP
C
C  Filter designed
        CALL DESIGN( IORD, TYPE(1:2), APROTO(1:2),
     &               FLO, FHI, TS, SN, SD, NSECTS )
C
C  Filter DATA1
C
        IF (   PASSES .EQ. 1 ) THEN
          ZP = .FALSE.
        ELSE
          ZP = .TRUE.
        END IF
        CALL APPLY( DATA1, NSAMPS, ZP, SN, SD, NSECTS )
C
      RETURN
      END
C
C  Copyright 1990  Regents of the University of California
C
C
C
C  Subroutine to generate second order section parameterization
C    from an pole-zero description for lowpass filters.
C
C  Input Arguments:
C  ----------------
C
C    P                       Array containing poles
C
C    Z                       Array containing zeros
C
C    RTYPE                   Character array containing root type information
C                              (SP)  Single real pole or
C                              (CP)  Complex conjugate pole pair
C                              (CPZ) Complex conjugate pole and zero pairs
C
C    DCVALUE                 Zero-frequency value of prototype filter
C
C    NSECTS                  Number of second-order sections
C
C  Output Arguments:
C  -----------------
C
C    SN                      Numerator polynomials for second order
C                              sections.
C
C    SD                      Denominator polynomials for second order
C                              sections.
C
C
      SUBROUTINE LP( P, Z, RTYPE, DCVALUE, NSECTS, SN, SD )
C
        COMPLEX P(*), Z(*)
        CHARACTER*3 RTYPE(*)
        REAL*4 SN(*), SD(*), DCVALUE
C
        IPTR = 1
        DO    1 I = 1, NSECTS
C
          IF (   RTYPE(I) .EQ. 'CPZ' ) THEN
C
            SCALE = REAL( P(I) * CONJG( P(I) ) )
     &            / REAL( Z(I) * CONJG( Z(I) ) )
            SN( IPTR )     = REAL( Z(I) * CONJG( Z(I) ) ) * SCALE
            SN( IPTR + 1 ) = -2. * REAL( Z(I) ) * SCALE
            SN( IPTR + 2 ) = 1. * SCALE
            SD( IPTR )     = REAL( P(I) * CONJG( P(I) ) )
            SD( IPTR + 1 ) = -2. * REAL( P(I) )
            SD( IPTR + 2 ) = 1.
            IPTR = IPTR + 3
C
          ELSE IF (   RTYPE(I) .EQ. 'CP' ) THEN
C
            SCALE = REAL( P(I) * CONJG( P(I) ) )
            SN( IPTR )     = SCALE
            SN( IPTR + 1 ) = 0.
            SN( IPTR + 2 ) = 0.
            SD( IPTR )     = REAL( P(I) * CONJG( P(I) ) )
            SD( IPTR + 1 ) = -2. * REAL( P(I) )
            SD( IPTR + 2 ) = 1.
            IPTR = IPTR + 3
C
          ELSE IF (  RTYPE(I) .EQ. 'SP' ) THEN
C
            SCALE = -REAL( P(I) )
            SN( IPTR )     = SCALE
            SN( IPTR + 1 ) = 0.
            SN( IPTR + 2 ) = 0.
            SD( IPTR )     = -REAL( P(I) )
            SD( IPTR + 1 ) = 1.
            SD( IPTR + 2 ) = 0.
            IPTR = IPTR + 3
C
          END IF
C
    1   CONTINUE
C
        SN(1) = DCVALUE * SN(1)
        SN(2) = DCVALUE * SN(2)
        SN(3) = DCVALUE * SN(3)

C
      RETURN
      END


C
C  Subroutine to convert an prototype lowpass filter to a bandpass filter via
C    the analog polynomial transformation.  The lowpass filter is
C    described in terms of its poles and zeros (as input to this routine).
C    The output consists of the parameters for second order sections.
C
C  Input Arguments:
C  ----------------
C
C    P                       Array containing poles
C
C    Z                       Array containing zeros
C
C    RTYPE                   Character array containing type information
C                              (SP) single real pole  or
C                              (CP) complex conjugate pole pair  or
C                              (CPZ) complex conjugate pole/zero pairs
C
C    DCVALUE                 Zero frequency value of filter
C
C    NSECTS                  Number of second-order sections upon input
C
C    FL                      Low-frequency cutoff
C
C    FH                      High-frequency cutoff
C
C  Output Arguments:
C  -----------------
C
C    SN                      Numerator polynomials for second order
C                              sections.
C
C    SD                      Denominator polynomials for second order
C                              sections.
C
C    NSECTS                  Number of second order sections upon output
C                              This subroutine doubles the number of
C                              sections.
C
C
      SUBROUTINE LPTBP( P, Z, RTYPE, DCVALUE, NSECTS, FL, FH, SN, SD )
C
        COMPLEX P(*), Z(*), CTEMP, P1, P2, Z1, Z2, S, H
        CHARACTER*3 RTYPE(*)
        REAL*4 SN(*), SD(*), DCVALUE
C
        PI = 3.14159265
        TWOPI = 2.*PI
        A = TWOPI*TWOPI*FL*FH
        B = TWOPI*( FH - FL )
C
        N = NSECTS
        NSECTS = 0
        IPTR = 1
        DO    1 I = 1, N
C
          IF (    RTYPE(I) .EQ. 'CPZ' ) THEN
C
            CTEMP = ( B*Z(I) )**2 - 4.*A
            CTEMP = CSQRT( CTEMP )
            Z1 = 0.5*( B*Z(I) + CTEMP )
            Z2 = 0.5*( B*Z(I) - CTEMP )
            CTEMP = ( B*P(I) )**2 - 4.*A
            CTEMP = CSQRT( CTEMP )
            P1 = 0.5*( B*P(I) + CTEMP )
            P2 = 0.5*( B*P(I) - CTEMP )
            SN( IPTR )     = REAL( Z1 * CONJG( Z1 ) )
            SN( IPTR + 1 ) = -2. * REAL( Z1 )
            SN( IPTR + 2 ) = 1.
            SD( IPTR )     = REAL( P1 * CONJG( P1 ) )
            SD( IPTR + 1 ) = -2. * REAL( P1 )
            SD( IPTR + 2 ) = 1.
            IPTR = IPTR + 3
            SN( IPTR )     = REAL( Z2 * CONJG( Z2 ) )
            SN( IPTR + 1 ) = -2. * REAL( Z2 )
            SN( IPTR + 2 ) = 1.
            SD( IPTR )     = REAL( P2 * CONJG( P2 ) )
            SD( IPTR + 1 ) = -2. * REAL( P2 )
            SD( IPTR + 2 ) = 1.
            IPTR = IPTR + 3
C
            NSECTS = NSECTS + 2
C
          ELSE IF (   RTYPE(I) .EQ. 'CP' ) THEN
C
            CTEMP = ( B*P(I) )**2 - 4.*A
            CTEMP = CSQRT( CTEMP )
            P1 = 0.5*( B*P(I) + CTEMP )
            P2 = 0.5*( B*P(I) - CTEMP )
            SN( IPTR )     = 0.
            SN( IPTR + 1 ) = B
            SN( IPTR + 2 ) = 0.
            SD( IPTR )     = REAL( P1 * CONJG( P1 ) )
            SD( IPTR + 1 ) = -2. * REAL( P1 )
            SD( IPTR + 2 ) = 1.
            IPTR = IPTR + 3
            SN( IPTR )     = 0.
            SN( IPTR + 1 ) = B
            SN( IPTR + 2 ) = 0.
            SD( IPTR )     = REAL( P2 * CONJG( P2 ) )
            SD( IPTR + 1 ) = -2. * REAL( P2 )
            SD( IPTR + 2 ) = 1.
            IPTR = IPTR + 3
C
            NSECTS = NSECTS + 2
C
          ELSE IF (  RTYPE(I) .EQ. 'SP' ) THEN
C
            SN( IPTR )     = 0.
            SN( IPTR + 1 ) = B
            SN( IPTR + 2 ) = 0.
            SD( IPTR )     = A
            SD( IPTR + 1 ) = -B*REAL( P(I) )
            SD( IPTR + 2 ) = 1.
            IPTR = IPTR + 3
C
            NSECTS = NSECTS + 1
C
          END IF
C
    1   CONTINUE
C
C  Scaling - use the fact that the bandpass filter amplitude at sqrt( omega_l *
C            equals the amplitude of the lowpass prototype at d.c.
C
        S = CMPLX( 0., SQRT(A) )
        H = CMPLX( 1., 0. )
C
        IPTR = 1
        DO    2 I = 1, NSECTS
          H = H * ( ( SN(IPTR+2)*S + SN(IPTR+1) )*S + SN(IPTR) )
     &      / ( ( SD(IPTR+2)*S + SD(IPTR+1) )*S + SD(IPTR) )
          IPTR = IPTR + 3
    2   CONTINUE
        SCALE = DCVALUE / SQRT( REAL( H ) * CONJG( H ) )

        SN(1) = SN(1) * SCALE
        SN(2) = SN(2) * SCALE
        SN(3) = SN(3) * SCALE
C
      RETURN
      END
C
C
C  Subroutine to convert a lowpass filter to a band reject filter
C    via an analog polynomial transformation.  The lowpass filter is
C    described in terms of its poles and zeros (as input to this routine).
C    The output consists of the parameters for second order sections.
C
C  Input Arguments:
C  ----------------
C
C    P                       Array containing poles
C
C    Z                       Array containing zeros
C
C    RTYPE                   Character array containing type information
C                              (SP)  single real pole or
C                              (CP)  complex conjugate pole pair
C                              (CPZ) complex conjugate pole/zero pairs
C
C    DCVALUE                 Zero-frequency value of prototype filter
C
C    NSECTS                  Number of second-order sections
C                              prior to transformation
C
C    FL                      Low-frequency cutoff
C
C    FH                      High-frequency cutoff
C
C  Output Arguments:
C  -----------------
C
C    SN                      Numerator polynomials for second order
C                              sections.
C
C    SD                      Denominator polynomials for second order
C                              sections.
C
C    NSECTS                  Number of second order sections following
C                              transformation.  The number is doubled.
C
C
      SUBROUTINE LPTBR( P, Z, RTYPE, DCVALUE, NSECTS, FL, FH, SN, SD )
C
        COMPLEX P(*), Z(*), CINV, CTEMP, P1, P2, Z1, Z2
        CHARACTER*3 RTYPE(*)
        REAL*4 SN(*), SD(*)
C
        PI = 3.14159265
        TWOPI = 2.*PI
        A = TWOPI*TWOPI*FL*FH
        B = TWOPI*( FH - FL )
C
        N = NSECTS
        NSECTS = 0
        IPTR = 1
        DO    1 I = 1, N
C
          IF (    RTYPE(I) .EQ. 'CPZ' ) THEN
C
            CINV = 1./Z(I)
            CTEMP = ( B*CINV )**2 - 4.*A
            CTEMP = CSQRT( CTEMP )
            Z1 = 0.5*( B*CINV + CTEMP )
            Z2 = 0.5*( B*CINV - CTEMP )
            CINV = 1./P(I)
            CTEMP = ( B*CINV )**2 - 4.*A
            CTEMP = CSQRT( CTEMP )
            P1 = 0.5*( B*CINV + CTEMP )
            P2 = 0.5*( B*CINV - CTEMP )
            SN( IPTR )     = REAL( Z1 * CONJG( Z1 ) )
            SN( IPTR + 1 ) = -2. * REAL( Z1 )
            SN( IPTR + 2 ) = 1.
            SD( IPTR )     = REAL( P1 * CONJG( P1 ) )
            SD( IPTR + 1 ) = -2. * REAL( P1 )
            SD( IPTR + 2 ) = 1.
            IPTR = IPTR + 3
            SN( IPTR )     = REAL( Z2 * CONJG( Z2 ) )
            SN( IPTR + 1 ) = -2. * REAL( Z2 )
            SN( IPTR + 2 ) = 1.
            SD( IPTR )     = REAL( P2 * CONJG( P2 ) )
            SD( IPTR + 1 ) = -2. * REAL( P2 )
            SD( IPTR + 2 ) = 1.
            IPTR = IPTR + 3
C
            NSECTS = NSECTS + 2
C
          ELSE IF (   RTYPE(I) .EQ. 'CP' ) THEN
C
            CINV = 1./P(I)
            CTEMP = ( B*CINV )**2 - 4.*A
            CTEMP = CSQRT( CTEMP )
            P1 = 0.5*( B*CINV + CTEMP )
            P2 = 0.5*( B*CINV - CTEMP )
            SN( IPTR )     = A
            SN( IPTR + 1 ) = 0.
            SN( IPTR + 2 ) = 1.
            SD( IPTR )     = REAL( P1 * CONJG( P1 ) )
            SD( IPTR + 1 ) = -2. * REAL( P1 )
            SD( IPTR + 2 ) = 1.
            IPTR = IPTR + 3
            SN( IPTR )     = A
            SN( IPTR + 1 ) = 0.
            SN( IPTR + 2 ) = 1.
            SD( IPTR )     = REAL( P2 * CONJG( P2 ) )
            SD( IPTR + 1 ) = -2. * REAL( P2 )
            SD( IPTR + 2 ) = 1.
            IPTR = IPTR + 3
C
            NSECTS = NSECTS + 2
C
          ELSE IF (  RTYPE(I) .EQ. 'SP' ) THEN
C
            SN( IPTR )     = A
            SN( IPTR + 1 ) = 0.
            SN( IPTR + 2 ) = 1.
            SD( IPTR )     = -A*REAL( P(I) )
            SD( IPTR + 1 ) = B
            SD( IPTR + 2 ) = -REAL( P(I) )
            IPTR = IPTR + 3
C
            NSECTS = NSECTS + 1
C
          END IF
C
    1   CONTINUE
C
C  Scaling - use the fact that the bandreject filter amplitude  at d.c.
C            equals the lowpass prototype amplitude at d.c.
C
        H = 1.0
C
        IPTR = 1
        DO    2 I = 1, NSECTS
          H = H * SN(IPTR) / SD(IPTR)
          IPTR = IPTR + 3
    2   CONTINUE
        SCALE = DCVALUE / ABS(H)
        SN(1) = SN(1) * SCALE
        SN(2) = SN(2) * SCALE
        SN(3) = SN(3) * SCALE
C
      RETURN
      END
C
C
C  Subroutine to convert a lowpass filter to a highpass filter via
C    an analog polynomial transformation.  The lowpass filter is
C    described in terms of its poles and zeroes (as input to this routine).
C    The output consists of the parameters for second order sections.
C
C  Input Arguments:
C  ----------------
C
C    P                       Array containing poles
C
C    Z                       Array containing zeroes
C
C    RTYPE                   Character array containing root type information
C                              (SP) single real pole or
C                              (CP)  complex conjugate pair
C                              (CPZ) complex pole/zero pairs
C
C    DCVALUE                 Zero-frequency value of prototype filter
C
C    NSECTS                  Number of second-order sections
C
C  Output Arguments:
C  -----------------
C
C    SN                      Numerator polynomials for second order
C                              sections.
C
C    SD                      Denominator polynomials for second order
C                              sections.
C
C
      SUBROUTINE LPTHP( P, Z, RTYPE, DCVALUE, NSECTS, SN, SD )
C
        COMPLEX P(*), Z(*)
        CHARACTER*3 RTYPE(*)
        REAL*4 SN(*), SD(*), DCVALUE
C
        IPTR = 1
        DO    1 I = 1, NSECTS
C
          IF (     RTYPE(I) .EQ. 'CPZ' ) THEN
C
            SCALE = REAL( P(I) * CONJG( P(I) ) )
     &            / REAL( Z(I) * CONJG( Z(I) ) )
            SN( IPTR )     = 1.  *  SCALE
            SN( IPTR + 1 ) = -2. * REAL( Z(I) )  *  SCALE
            SN( IPTR + 2 ) = REAL( Z(I) * CONJG( Z(I) ) )  *  SCALE
            SD( IPTR )     = 1.
            SD( IPTR + 1 ) = -2. * REAL( P(I) )
            SD( IPTR + 2 ) = REAL( P(I) * CONJG( P(I) ) )
            IPTR = IPTR + 3
C
          ELSE IF (   RTYPE(I) .EQ. 'CP' ) THEN
C
            SCALE = REAL( P(I) * CONJG( P(I) ) )
            SN( IPTR )     = 0.
            SN( IPTR + 1 ) = 0.
            SN( IPTR + 2 ) = SCALE
            SD( IPTR )     = 1.
            SD( IPTR + 1 ) = -2. * REAL( P(I) )
            SD( IPTR + 2 ) = REAL( P(I) * CONJG( P(I) ) )
            IPTR = IPTR + 3
C
          ELSE IF (  RTYPE(I) .EQ. 'SP' ) THEN
C
            SCALE = -REAL( P(I) )
            SN( IPTR )     = 0.
            SN( IPTR + 1 ) = SCALE
            SN( IPTR + 2 ) = 0.
            SD( IPTR )     = 1.
            SD( IPTR + 1 ) = -REAL( P(I) )
            SD( IPTR + 2 ) = 0.
            IPTR = IPTR + 3
C
          END IF
C
    1   CONTINUE
C
        SN(1) = SN(1) * DCVALUE
        SN(2) = SN(2) * DCVALUE
        SN(3) = SN(3) * DCVALUE
C
      RETURN
      END
C
C
C WARP -- FUNCTION, APPLIES TANGENT FREQUENCY WARPING TO COMPENSATE
C         FOR BILINEAR ANALOG -> DIGITAL TRANSFORMATION
C
C ARGUMENTS:
C ----------
C
C      F       ORIGINAL DESIGN FREQUENCY SPECIFICATION (HERTZ)
C      TS      SAMPLING INTERVAL (SECONDS)
C
C  LAST MODIFIED:  SEPTEMBER 20, 1990
C
      REAL FUNCTION WARP( F , TS )
C
        TWOPI = 6.2831853
        ANGLE = TWOPI*F*TS/2.
        WARP = 2.*TAN(ANGLE)/TS
        WARP = WARP/TWOPI
C
      RETURN
      END

Cf2py END FILE bpfilter.f      
