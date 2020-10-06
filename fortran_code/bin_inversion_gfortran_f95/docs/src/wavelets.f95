module wavelets


   use constants, only : pi, twopi, inptd, nnsta
   use wavelet_param, only : jmin, jmax, lnpt, nlen
   implicit none
   complex :: rwt1(inptd, 12), rwt2(inptd, 12), c1, c2
   integer :: kkk(4200, 15)
   real :: cos_fft(4200, 15), sin_fft(4200, 15)


contains


   subroutine wavelet_obs(cr, cz, u)
   !!
   !!   Compute wavelet transform. Input data is in time domain
   !!
   implicit none
!       ******* Explaination ********
!       In order to normalized the amplitude, we move the coeffficient
!       sqrt(T/kmax)
!                                       Jichen 1997, 1, 20
!
! we split between this routine and wavelet_syn routine to hopefully reduce
! computation time
!
   real :: cr(inptd), cz(inptd), u(inptd)
   integer lcc, kmax, is, j, k, i
   
   complex fre(inptd)
   complex cc
!
   do i = 1, nlen
      cr(i) = u(i)
      cz(i) = 0.0
   end do
   call cfft(cr, cz, lnpt)

   call realtr(cr, cz, lnpt)
   do j = 1, nlen
      fre(j) = cmplx(cr(j), cz(j))
      cr(j) = 0.0
      cz(j) = 0.0
      u(j) = 0.0
   end do
!
! c1 = wave(pi2, 2), c2 = wave(pi, 2).
!
   u(1) = 2.0*real(fre(2)*c1)
   u(2) = 2.0*real(fre(2)*c2+fre(3)*c1)
   u(3) = 2.0*real(-fre(2)*c2+fre(3)*c1)
   j = max(3, jmin) - 1
   kmax = 2 ** (j-1)
   do j = max(3, jmin), jmax
      kmax = 2*kmax
      do is = 1, kmax
         cc = fre(is)*rwt1(is, j)+fre(is+kmax)*rwt2(is, j)
         cr(is) = real(cc)
         cz(is) = aimag(cc)
      end do
      lcc = j-1
      call cifft(cr, cz, Lcc)
      do k = 1, kmax
         u(kmax+k-1) = 2.0*cr(k)
      end do
   end do
   end subroutine wavelet_obs


   pure subroutine wavelet_syn(cr, cz, u)
   !!
   !!   Compute wavelet transform. Input data is in frequency domain
   !!
   implicit none
!       ******* Explaination ********
!       In order to normalized the amplitude, we move the coeffficient
!       sqrt(T/kmax)
!                         Jichen 1997, 1, 20
   real, intent(out) :: u(inptd)
   real, intent(inout) :: cr(inptd), cz(inptd)
   integer lcc, kmax, is, lb, i, i1, i2, j, k
        
   complex fre(inptd)
   complex cc

   lcc = 2**(lnpt-1)
   lb = lcc-1
   lcc = lcc+1
   do i = 1, lb
      i1 = lcc+i
      i2 = lcc-i
      cr(i1) = cr(i2)
      cz(i1) = -cz(i2)
   end do
   cz(lcc) = 0.0
   do j = 1, nlen
      fre(j) = cmplx(cr(j)/nlen, cz(j)/nlen)
   end do
!
! c1 = wave(pi2, 2), c2 = wave(pi, 2).
!
   u(1) = 2.0*real(fre(2)*c1)
   u(2) = 2.0*real(fre(2)*c2+fre(3)*c1)
   u(3) = 2.0*real(-fre(2)*c2+fre(3)*c1)
        
   j = max(3, jmin) - 1
   kmax = 2 ** (j-1)
   do j = max(3, jmin), jmax
      kmax = 2*kmax  ! ** (j - 1)   ! 2*kmax
      do is = 1, kmax
         cc = fre(is)*rwt1(is, j)+fre(is+kmax)*rwt2(is, j)
         cr(is) = real(cc)
         cz(is) = aimag(cc)
      end do
      lcc = j-1
      call cifft(cr, cz, Lcc)
      do k = 1, kmax
         u(kmax+k-1) = 2.0*cr(k)
      end do
   end do

   end subroutine wavelet_syn


   pure subroutine cfft(xr, xi, n)
   !!
   !!   Fast Fourier transform
   !!
   implicit none
   integer, intent(in) :: n
   real, intent(inout) :: xr(*), xi(*)
   integer k, i, ib, nb, lx, l, lb, lbh, ist, jh, j1, j
   real wkr, wki, qr, qi, holdr, holdi, flx
!
! we store in memory the values of 'k', so as to not need to initialize them time and again.
!
   LX = 2**N
   FLX = real(LX)
   DO L = 1, N
      NB = 2**(L-1)
      LB = LX/NB
      LBH = LB/2
      DO IB = 1, NB           ! 2 ** (l-1) operaciones
         WKR = cos_fft(ib, n)
         WKI = -sin_fft(ib, n)
         IST = LB*(IB-1)
         DO J = IST+1, IST+LBH
            JH = J+LBH
            QR = XR(JH)*WKR-XI(JH)*WKI
            QI = XR(JH)*WKI+XI(JH)*WKR
            XR(JH) = (XR(J)-QR)
            XI(JH) = (XI(J)-QI)
            XR(J) = (XR(J)+QR)
            XI(J) = (XI(J)+QI)
         end do
      end do
   end do
   K = 0
   DO J = 1, LX
      K = KKK(J, N)
      IF(K.LT.J) cycle
      HOLDR = XR(J)
      HOLDI = XI(J)
      J1 = K+1
      XR(J) = XR(J1)
      XI(J) = XI(J1)
      XR(J1) = HOLDR
      XI(J1) = HOLDI
   end do
   DO I = 1, LX
      XR(I) = XR(I)/FLX
      XI(I) = XI(I)/FLX
   ENDDO
   END subroutine cfft


   pure subroutine cifft(xr, xi, n)
   !!
   !!   Inverse Fast Fourier transform
   !!
   implicit none
   integer, intent(in) :: n
   real, intent(inout) :: xr(*), xi(*)
   integer k, ib, nb, lx, l, lb, lbh, ist, jh, j1, j
   real wkr, wki, qr, qi, holdr
   LX = 2**N
   DO L = 1, N
      NB = 2**(L-1)
      LB = LX/NB
      LBH = LB/2
      DO IB = 1, NB           ! 2 ** (l-1) operaciones
         WKR = cos_fft(ib, n)
         WKI = sin_fft(ib, n)
         IST = LB*(IB-1)
         DO J = IST+1, IST+LBH
            JH = J+LBH
            QR = XR(JH)*WKR-XI(JH)*WKI
            QI = XR(JH)*WKI+XI(JH)*WKR
            XR(JH) = XR(J)-QR
            XI(JH) = XI(J)-QI
            XR(J) = XR(J)+QR
            XI(J) = XI(J)+QI
         ENDDO
      ENDDO
   ENDDO
   K = 0
   DO J = 1, LX
      K = KKK(j, N)
      IF(K.LT.J) cycle
      HOLDR = XR(J)
      J1 = K+1
      XR(J) = XR(J1)
      XR(J1) = HOLDR
   end do
   end subroutine cifft


   subroutine fourier_coefs()
   implicit none
   !!
   !!  We load into memory certain values of sine and cosine, frequently used in computing the
   !!  Fast Fourier Transform.
   !!
   integer i, n, nb, nnb, ib, k, kk(4200, 15)
   real*8 :: omega

   do n = 2, 12
      k = 0
      nb = 2 ** (n-1)
      nnb = 2 ** n
      do ib = 1, nb
         kk(ib, n) = k
         omega = twopi*k/nnb
         cos_fft(ib, n) = cos(omega)
         sin_fft(ib, n) = sin(omega)
         do i = 2, n
            if (k .lt. 2**(n-i)) exit
            k = k-2**(n-i)
         end do
         k = k+2**(n-i)         
      end do
      k = 0
      do ib = 1, nnb
         kkk(ib, n) = k
         do i = 1, n
            if (k .lt. 2**(n-i)) exit
            k = k-2**(n-i)
         end do
         k = k+2**(n-i)         
      end do
   end do

   end subroutine fourier_coefs


   subroutine meyer_yamada()
   !!
   !!   We load into memory the values of the Meyer wavelet function, 
   !!   which are used when computing the wavelet transform of an input data.
   !!
   implicit none
   real*8 :: omega1, omega2
   integer j, is, kmax
   c1 = wave(twopi, 2)
   c2 = wave(pi, 2)
   do j = 1, 12
      kmax = 2**(j-1)
      do is = 1, kmax
         omega1 = twopi*(is-1)/(kmax)
         omega2 = omega1+twopi
         rwt1(is, j) = wave(omega1, 2)
         rwt2(is, j) = wave(omega2, 2)
      end do
   end do
   end subroutine meyer_yamada


   function wave(w, ic) result(wave1)
   !!
   !!  Meyer wavelet function in frequency domain.
   !!
   implicit none
   real*8 :: w
   integer :: ic
   integer :: i, j, k
   real*8, parameter :: p43 = 4.d0*pi/3.d0, p23 = 2.d0*pi/3.d0
   real*8 :: g(2), p(2), fw, wf, wd, wp, wg, ww, hw, cct
   complex*16 :: pp, a
   complex :: wave1
   hw = w*0.5d0
   a = cmplx(cos(hw),-sin(hw), kind(1.d0))
   ww = abs(w)
   if (ww .gt. (2.d0*p43) .or. ww .lt. p23) then
      wave1 = cmplx(0.0, 0.0)
      return
   end if
   do i = 1, 2
      wp = hw*i
      fw = 1.d0
      do j = 1, 2
         if (j .eq. 1) wf = -wp
         if (j .eq. 2) wf = wp
         do k = 1, 2
            if (k .eq. 1) wd = p43-wf
            if (k .eq. 2) wd = wf-p23
            if (wd .le. 0) then
               g(k) = 0.d0
            else
               g(k) = exp(-1.d0/wd/wd)
            end if
         end do
         wg = g(1)/(g(1)+g(2))
         fw = fw*wg
      end do
      p(i) = fw
   end do
   cct = p(1)-p(2)
   pp = a*dsqrt(cct)
   if (ic .eq. 1) then
      wave1 = pp 
   else
      wave1 = cmplx(real(pp),-aimag(pp))    
   end if
   end function wave


   subroutine realtr(xr, xi, n)
   implicit none
   real xr(*), xi(*)
   integer :: n
   integer :: i, i1, i2, lh, lb
   lh = 2**(n-1)
   lb = lh-1
   lh = lh+1
   do i = 1, lb
      i1 = lh+I
      i2 = lh-I
      xr(i1) = xr(i2)
      xi(i1) = -xi(i2)
   end do
   xi(lh) = 0.0
   end subroutine realtr


end module wavelets
