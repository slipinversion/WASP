module wavelet_param


   implicit none
   integer :: lnpt, nlen, jmin, jmax, max_freq


contains

   
   subroutine get_data_param(lnpt0, jmin0, jmax0, nlen0, max_freq0)
   implicit none
   integer :: lnpt0, nlen0, jmin0, jmax0, max_freq0
   lnpt0 = lnpt
   jmin0 = jmin
   jmax0 = jmax
   nlen0 = nlen
   max_freq0 = max_freq
   end subroutine get_data_param


   subroutine set_params(lnpt0, jmin0, jmax0, nlen0, max_freq0)
   implicit none
   integer :: lnpt0, nlen0, jmin0, jmax0, max_freq0
   lnpt = lnpt0
   jmin = jmin0
   jmax = jmax0
   nlen = nlen0
   max_freq = max_freq0
   end subroutine set_params


end module wavelet_param
