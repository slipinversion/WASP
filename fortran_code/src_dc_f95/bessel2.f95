module bessel2


   use constants, only : pi, twopi, ndis
   use retrieve_gf, only : dist_max, lnpt, dt
   implicit none
   integer, parameter :: dmax=3000, fmax=300000
   real :: aj0s(dmax, fmax), aj1s(dmax, fmax), aj2s(dmax, fmax)


contains


   subroutine load_bessel(x,nd_max)
   implicit none
   real :: aj0, aj1, aj2, x(ndis)
   real :: k, dk, z, pmax, omega, dw
   integer :: i, n, ix, kc, nfft2, nfft, nd_max
   dk = 0.2
   dk = dk*pi/dist_max
   if(dk.gt.0.03)dk=0.03
   k=0.5*dk 
   kc = 30 / 0.5
   nfft = 2 ** lnpt
   nfft2 = nfft
   dw = twopi/(nfft*dt)
   omega=(nfft2-1)*dw
   pmax = 1.11
   n=(sqrt(kc*kc+(pmax*omega)**2)-k)/dk
   write(*,*)"maximum n is ",n
   if(n.gt.fmax)then
      write(*,*)n,fmax
      write(*,*)"You might want to increase the size of nkk_max"
      n=fmax
   endif
   do i=1,n
      do ix=1,nd_max
         z=x(ix)*k
         call besselFn(z,aj0,aj1,aj2)
         aj0s(ix,i)=aj0
         aj1s(ix,i)=aj1
         aj2s(ix,i)=aj2
      enddo
      k=k+dk
   enddo

   end subroutine load_bessel


end module bessel2
