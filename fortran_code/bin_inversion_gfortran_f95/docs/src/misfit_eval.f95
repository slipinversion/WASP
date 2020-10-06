module misfit_eval


   use constants, only : inptd
   use wavelet_param, only : lnpt, jmin, jmax, jfmax, nlen
   use get_stations_data, only : weight, wave_obs, wmax, misfit_type, t_max, wavelet_weight
   implicit none


contains


   pure subroutine misfit(nf, wave_syn, error)
   !!  
   !!   Misfit between observed and synthetic waveforms, in wavelet domain.
   !!
   real, intent(inout) :: wave_syn(inptd)
   integer, intent(in) :: nf
   real, intent(out) :: error
   real*8 :: aa, ab, bb, er, er1, er2, erp, error2
   real :: ramp
   integer :: i, j, k, n1, n2, n_be, n_begin, n_delt, n_ed
!     
   error2 = 0.d0
   ramp = wmax(nf)
   n1 = 2**(jmin-1)
   n2 = 2**(jmax)-1
   do i = n1, n2
      wave_syn(i) = wave_syn(i)/ramp
   end do
!   ramp = 1
   do j = jmin, jmax
      er = 0.d0
      if (misfit_type(j, nf) .eq. 0) cycle
      if (wavelet_weight(j, nf) .lt. 1.0e-5) cycle
      n_begin = 2**(j-1)
      n_delt = nlen/n_begin
      n_be = n_begin
      n_ed = n_begin+int(t_max(nf)/n_delt+0.5)-1
      if (n_ed .lt. n_be) n_ed = n_be
!  j = 1 L1 Norm
      if (misfit_type(j, nf) .eq. 1) then
         do k = n_be, n_ed
            er = er+abs(wave_syn(k)-wave_obs(k, nf))
         end do
         er = wavelet_weight(j, nf)*(er/(n_ed-n_be+1))
         error2 = error2 + er
      end if
!       j = 2 L2 norm
      if (misfit_type(j, nf) .eq. 2) then
         do k = n_be, n_ed
            erp = wave_syn(k)-wave_obs(k, nf)
            er = er+erp*erp
         end do
         er = wavelet_weight(j, nf)*sqrt(er/(n_ed-n_be+1))
         error2 = error2 + er
      end if
!       j = 3 L1+L2 Norm
      if (misfit_type(j, nf) .eq. 3) then
         er1 = 0.d0
         er2 = 0.d0
         do k = n_be, n_ed
            erp = abs(wave_syn(k)-wave_obs(k, nf))
            er1 = er1+erp 
            er2 = er2+erp*erp
         end do
         er = wavelet_weight(j, nf) &
      &  *(0.5d0*sqrt(er2 / (n_ed-n_be+1))+0.5d0*er1 / (n_ed-n_be+1))
         error2 = error2 + er
      end if
!     j = 4 correlation
      if (misfit_type(j, nf) .eq. 4) then
         ab = 0.d0
         aa = 0.d0
         bb = 0.d0
         do k = n_be, n_ed
            ab = wave_syn(k)*wave_obs(k, nf)+ab
            aa = wave_syn(k)*wave_syn(k)+aa
            bb = wave_obs(k, nf)*wave_obs(k, nf)+bb
         end do
         er = wavelet_weight(j, nf)*(1.d0-2.d0*ab/(aa+bb))
         error2 = error2 + er
      end if
   end do
   error2 = error2/(jmax-jmin+1)
   error2 = error2*weight(nf)
   error = real(error2)
   end subroutine misfit


end module misfit_eval
