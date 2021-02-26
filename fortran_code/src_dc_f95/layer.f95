!***************************************************************
!  F-K: @(#) layer.h			1.0 4/29/2000
!
!  Copyright (c) 2000 by L. Zhu
!  See README file for copying and redistribution conditions.
!
!  include header for transfering parameters of any given layer
!  from one subroutines to another
!***************************************************************
module layer


   use omp_lib
   use constants, only : nlay
   implicit none
   real*8 kd,exa,exb,mu2,y,y1
   complex*16 ra,rb,Ca,Cb,Ya,Xa,Yb,Xb,r,r1,ka(nlay), kb(nlay)
!$omp threadprivate(kd, mu2, exa, exb, ra, rb, Ca, Cb, Ya, Xa, Yb, Xb) 
!$omp threadprivate(r, r1, y, y1, ka, kb)


end module layer
