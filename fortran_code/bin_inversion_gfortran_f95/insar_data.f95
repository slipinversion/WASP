!     This program load the static subfault GFs and generate the matrix
!     for baseline correction.         Ji, Chen, 2003
!     Reviewed. The previous version have an error in loading GFs.
!                      
!                                   Ji, Chen, 2004
!  Fixed
!
!
!
module insar_data


   use constants, only : max_subf, max_seg, dpi
   use model_parameters, only : nxs_sub, nys_sub, nx_p, ny_p, segments
   implicit none
   integer, parameter, private :: n_stations = 2500
   integer, private :: n_chan, lines_asc, lines_desc
   integer :: ramp_length1, ramp_length2, ramp_length
   real*8 :: synm_whole(n_stations), weight_sum
   real*8 :: ramp_gf(18, n_stations)
   real :: lat(n_stations), lon(n_stations), max_los, syn_disp(n_stations)
   real, private, allocatable :: green(:, :, :, :)
   real :: obse(n_stations), look(3, n_stations), weight(n_stations)
   character(len=6) :: sta_name(n_stations)
   character(len=10) :: ramp_type1, ramp_type2
  

contains


   subroutine get_insar_gf()
   implicit none
   allocate(green(6, n_stations, max_subf, max_seg))
   end subroutine get_insar_gf
   

   subroutine get_insar_data()
   implicit none
   integer i, no, k
   real :: weight_asc, weight_desc
!
   open(9, file='insar_data.txt', status='old')
   read(9,*) n_chan
   read(9,*) lines_asc, lines_desc, weight_asc, weight_desc
   do i = 1, n_chan
      read(9,*) no, sta_name(i), lat(i), lon(i), obse(i), (look(k, i), k = 1, 3)
      if (i .le. lines_asc) weight(i) = weight_asc 
      if (i .gt. lines_asc) weight(i) = weight_desc 
   end do
   weight_sum = lines_asc*weight_asc + lines_desc*weight_desc
   max_los = maxval(abs(obse(:)))
   close(9) 
   end subroutine get_insar_data


   subroutine is_ramp(ramp_gf_file)
   implicit none
   integer :: channel
   logical :: ramp_gf_file
   inquire( file = 'ramp_gf.txt', exist = ramp_gf_file )
   if (ramp_gf_file) then
      open(23, file='ramp_gf.txt', status='old')
      read(23,*) ramp_type1, ramp_type2
      ramp_length1 = 0
      ramp_length2 = 0
      select case (ramp_type1)
      case ('linear')
         if (lines_asc .gt. 0) ramp_length1 = ramp_length1 + 3
      case ('bilinear')
         if (lines_asc .gt. 0) ramp_length1 = ramp_length1 + 6
      case ('quadratic')
         if (lines_asc .gt. 0) ramp_length1 = ramp_length1 + 5
      end select
      select case (ramp_type2)
      case ('linear')
         if (lines_desc .gt. 0) ramp_length2 = ramp_length2 + 3
      case ('bilinear')
         if (lines_desc .gt. 0) ramp_length2 = ramp_length2 + 6
      case ('quadratic')
         if (lines_desc .gt. 0) ramp_length2 = ramp_length2 + 5
      end select
      ramp_length = ramp_length1 + ramp_length2
      do channel = 1, n_chan
         ramp_gf(:, channel) = 0.d0
         read(23,*)ramp_gf(:ramp_length, channel)
      end do
      close(23)
   end if
   end subroutine is_ramp


   subroutine initial_ramp(ramp)
   implicit none
   real*8 :: ramp(18)
   integer :: i
   ramp(:) = 0.d0
   do i=1, ramp_length, 2
      ramp(i) = 20
   end do
   do i=2, ramp_length, 2
      ramp(i) = -20
   end do
   end subroutine initial_ramp


   subroutine initial_insar(slip, rake, ramp)
   implicit none
   real*8, optional :: ramp(18)
   integer k, j, segment, channel, nxy
   integer iys, ixs, n_tt!, lines_asc, lines_desc
   real slip(max_subf, max_seg), rake(max_subf, max_seg)
   real :: cosal, sinal, angle
   real*8 :: disp, ramp2
   logical is_file
!

   open(33, file='Insar_static_subfault.txt', status='old')
   read(33,*) n_tt
   do channel = 1, n_chan
      read(33,*)
      do segment = 1, segments
         read(33,*)
         nxy = 0
         do iys = 1, nys_sub(segment)
            do ixs = 1, nxs_sub(segment)
               nxy = nxy + 1
               read(33,*)(green(k, channel, nxy, segment), k = 1, 6)
            end do
         end do
      end do
   end do
   close(33)

   do channel = 1, n_chan
      syn_disp(channel) = 0.d0
      do k = 1, 3
         j = 2 * k - 1
         disp = 0.d0
         do segment = 1, segments
            nxy = 0
            do iys = 1, nys_sub(segment)
               do ixs = 1, nxs_sub(segment)
                  nxy = nxy + 1
                  angle = rake(nxy, segment)*dpi 
                  sinal = sin(angle)
                  cosal = cos(angle)
                  disp = disp + slip(nxy, segment) &
     &  *(sinal*green(j, channel, nxy, segment)+cosal*green(j+1, channel, nxy, segment))
               end do
            end do
         end do
         disp = disp * look(k, channel)
         syn_disp(channel) = syn_disp(channel) + disp
      end do
      if (present(ramp)) then
         ramp2 = 0.d0 
         do k = 1, ramp_length
            ramp2 = ramp2 + ramp_gf(k, channel)*ramp(k)
         end do 
         syn_disp(channel) = syn_disp(channel) + ramp2
      endif
   end do
   open(10,file='insar_synthetics.txt')
   write(10,*) n_chan
   do channel = 1, n_chan
      write(10,*) channel, sta_name(channel), lat(channel), lon(channel), syn_disp(channel)
   end do
   close(10)
   
   if (present(ramp)) then
      open(11,file='ramp_coefficients.txt')
      write(11,*) ramp(:)
      close(11)
   
      open(12,file='insar_ramp.txt')
      write(12,*) n_chan, lines_asc, lines_desc
      do channel = 1, n_chan
         ramp2 = 0.d0 
         do k = 1, ramp_length
            ramp2 = ramp2 + ramp_gf(k, channel)*ramp(k)
         end do 
         write(12,*) channel, sta_name(channel), lat(channel), lon(channel), ramp2
      end do
      close(12)
   end if
   end subroutine initial_insar
  

!
! routine for loading static synthetic seismograms, given a rupture model
!
   subroutine insar_synthetic(slip, rake, subfaults_segment, err, ramp)
   implicit none
   real*8, optional :: ramp(18)
   real slip(max_subf, max_seg), rake(max_subf, max_seg)!, err
   integer k, j, segment, channel, subfaults_segment(max_seg), nxy
   real err, dif, angle, sinal, cosal
   real*8 :: disp, err2, ramp2
      
   err2 = 0.d0
   synm_whole(:) = 0.d0
   do channel = 1, n_chan
      do k = 1, 3
         j = 2 * k - 1
         disp = 0.d0
         do segment = 1, segments
            do nxy = 1, subfaults_segment(segment)
               angle = rake(nxy, segment)*dpi 
               sinal = sin(angle)
               cosal = cos(angle)
               disp = disp + slip(nxy, segment) &
       &  *(sinal*green(j, channel, nxy, segment)+cosal*green(j+1, channel, nxy, segment))
            end do
         end do
         disp = disp * look(k, channel)
         synm_whole(channel) = synm_whole(channel) + disp
      end do
      if (present(ramp)) then
         ramp2 = 0.d0 
         do k = 1, ramp_length
            ramp2 = ramp2 + ramp_gf(k, channel)*ramp(k)
         end do 
         synm_whole(channel) = synm_whole(channel) + ramp2
      endif
      dif = synm_whole(channel) - obse(channel)
      err2 = err2 + weight(channel) * dif * dif / max_los / max_los!100.0
   end do
   err2 = sqrt(err2/weight_sum)
   err = real(err2)
   end subroutine insar_synthetic
            

!  
! subroutine for removing the static response of current subfault for all static stations
!  
   subroutine insar_remove_subfault(slip, rake, n_s, n_sub)
   implicit none
   real, intent(in) :: slip, rake
   integer, intent(in) :: n_s, n_sub
   integer k, j, channel
   real disp, angle, sinal, cosal
!
   angle = rake*dpi
   sinal = sin(angle)
   cosal = cos(angle)

   do channel = 1, n_chan
      disp = 0.0
      do k = 1, 3
         j = 2*k-1
         disp = disp + slip*look(k, channel) &
     &   *(sinal*green(j, channel, n_sub, n_s)+cosal*green(j+1, channel, n_sub, n_s))
      end do
      synm_whole(channel) = synm_whole(channel)-disp
   end do

   end subroutine insar_remove_subfault


!
! subroutine for testing the new response of the current subfault, for all stations
! we also give the misfit error of static data
!
   pure subroutine insar_modify_subfault(slip, rake, n_s, n_sub, err)
   implicit none
   real, intent(in) :: slip, rake
   real, intent(out) :: err
   integer, intent(in) :: n_s, n_sub
   integer k, j, channel
   real disp, angle, dif, sinal, cosal
   real*8 :: err2
!
   angle = rake*dpi
   sinal = sin(angle)
   cosal = cos(angle)

   err2 = 0.d0
   do channel = 1, n_chan
      disp = 0.0
      do k = 1, 3
         j = 2*k-1
         disp = disp + slip*look(k, channel) &
       & *(sinal*green(j, channel, n_sub, n_s)+cosal*green(j+1, channel, n_sub, n_s))
      end do
      dif = (synm_whole(channel) + disp) - obse(channel)
      err2 = err2 + weight(channel) * dif * dif / max_los / max_los!100.0
   end do
   err2 = sqrt(err2/weight_sum)
   err = real(err2)

   end subroutine insar_modify_subfault


!
! subroutine for asliping response of current subfault, for all stations.
! we also give the misfit error of static data
!
   subroutine insar_add_subfault(slip, rake, n_s, n_sub)
   implicit none
   real, intent(in) :: slip, rake
   integer, intent(in) :: n_s, n_sub
   integer k, j, channel
   real disp, angle, dif, sinal, cosal

   angle = rake*dpi
   sinal = sin(angle)
   cosal = cos(angle)

   do channel = 1, n_chan
      disp = 0.0
      do k = 1, 3
         j = 2*k-1
         disp = disp + look(k, channel)*slip &
       & *(sinal*green(j, channel, n_sub, n_s)+cosal*green(j+1, channel, n_sub, n_s))
      end do
      synm_whole(channel) = synm_whole(channel)+disp
   end do

   end subroutine insar_add_subfault


!  
! subroutine for removing the static response of current subfault for all static stations
!  
   subroutine insar_remove_ramp(ramp)
   implicit none
   real*8, intent(in) :: ramp(18)
   integer k, j, channel
   real ramp2
!
   do channel = 1, n_chan
      ramp2 = 0.0
      do k = 1, ramp_length
         ramp2 = ramp2 + ramp_gf(k, channel)*ramp(k)
      end do
      synm_whole(channel) = synm_whole(channel)-ramp2
   end do

   end subroutine insar_remove_ramp


!
! subroutine for testing the new response of the current subfault, for all stations
! we also give the misfit error of static data
!
   pure subroutine insar_modify_ramp(ramp, err)
   implicit none
   real*8, intent(in) :: ramp(18)
   real, intent(out) :: err
   integer k, j, channel
   real*8 :: err2, ramp2, dif
!
   err2 = 0.d0
   do channel = 1, n_chan
      ramp2 = 0.0
      do k = 1, ramp_length
         ramp2 = ramp2 + ramp_gf(k, channel)*ramp(k)
      end do
      dif = (synm_whole(channel) + ramp2) - obse(channel)
      err2 = err2 + weight(channel) * dif * dif / max_los / max_los!100.0
   end do
   err2 = sqrt(err2/weight_sum)
   err = real(err2)

   end subroutine insar_modify_ramp


!
! subroutine for asliping response of current subfault, for all stations.
! we also give the misfit error of static data
!
   subroutine insar_add_ramp(ramp)
   implicit none
   real*8, intent(in) :: ramp(18)
   integer k, j, channel
   real*8 :: err2, ramp2, dif
   
   do channel = 1, n_chan
      ramp2 = 0.d0
      do k = 1, ramp_length
         ramp2 = ramp2 + ramp_gf(k, channel)*ramp(k)
      end do
      synm_whole(channel) = synm_whole(channel)+ramp2
   end do

   end subroutine insar_add_ramp


   subroutine deallocate_insar_gf()
   implicit none
   deallocate(green)
   end subroutine deallocate_insar_gf


end module insar_data
