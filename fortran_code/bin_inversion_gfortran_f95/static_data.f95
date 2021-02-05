!     This program load the static subfault GFs and generate the matrix
!     for baseline correction.         Ji, Chen, 2003
!     Reviewed. The previous version have an error in loading GFs.
!                      
!                                   Ji, Chen, 2004
!  Fixed
!
!
!
module static_data


   use constants, only : max_subf, max_seg, dpi
   use model_parameters, only : nxs_sub, nys_sub, nx_p, ny_p, segments
   implicit none
   integer, parameter, private :: n_stations = 500
   integer :: n_chan
   real*8 :: synm_whole(n_stations, 3), weight_sum
   real :: lat(n_stations), lon(n_stations)
   real :: green(n_stations, 6, max_subf, max_seg), syn_disp(n_stations, 3)
   real :: obse(n_stations, 3), weight(n_stations, 3)
   character(len=6) :: sta_name(n_stations)
  

contains


   subroutine initial_gps(slip, rake)
   integer k, j, segment, channel, no, nxy, i, segmenteg
   integer iys, ixs, n_tt
   real slip(max_subf, max_seg), rake(max_subf, max_seg)
   real :: cosal, sinal, angle
   real*8 :: disp
   logical is_file
!
   inquire( file = 'Readlp.static', exist = is_file )
   if (is_file) then 
      weight_sum = 0.d0
      open(9, file='Readlp.static', status='old')
      read(9,*) n_chan
      read(9,*)
      do i = 1, n_chan
         read(9,*) no, sta_name(i), lat(i), lon(i), (obse(i, k), k = 1, 3), (weight(i, k), k = 1, 3)
         do k = 1, 3
            weight_sum = weight_sum + weight(i, k)
         end do
      end do
      close(9) 

      open(33, file='Green_static_subfault', status='old')
      read(33,*) n_tt
      do channel = 1, n_chan
         read(33,*)
         do segmenteg = 1, segments
            read(33,*)
            nxy = 0
            do iys = 1, nys_sub(segmenteg)
               do ixs = 1, nxs_sub(segmenteg)
                  nxy = nxy + 1
                  read(33,*)(green(channel, k, nxy, segmenteg), k = 1, 6)
               end do
            end do
         end do
      end do
      close(33)
      
      do k = 1, 3
         j = 2 * k - 1
         do channel = 1, n_chan
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
     &  *(sinal*green(channel, j, nxy, segment)+cosal*green(channel, j+1, nxy, segment))
                  end do
               end do
            end do
            syn_disp(channel, k) = disp
         end do
      end do 
      open(10,file='synm.static')
      write(10,*) n_chan
      do channel = 1, n_chan
         write(10,*) channel, sta_name(channel), lat(channel), lon(channel),(syn_disp(channel, k), k=1, 3)
      end do
      close(10)
   end if
   end subroutine initial_gps
   
   
!
! routine for loading static synthetic seismograms, given a rupture model
!
   subroutine static_synthetic(slip, rake, subfaults_segment, err)
   real slip(max_subf, max_seg), rake(max_subf, max_seg)!, err
   integer k, j, segment, channel, subfaults_segment(max_seg), nxy
   real err, dif, angle, sinal, cosal
   real*8 :: disp, err2
      
   err2 = 0.d0
   synm_whole(:, :) = 0.d0
   do k = 1, 3
      j = 2 * k - 1
      do channel = 1, n_chan
         disp = 0.d0
         do segment = 1, segments
            do nxy = 1, subfaults_segment(segment)
               angle = rake(nxy, segment)*dpi 
               sinal = sin(angle)
               cosal = cos(angle)
               disp = disp + slip(nxy, segment) &
       &  *(sinal*green(channel, j, nxy, segment)+cosal*green(channel, j+1, nxy, segment))
            end do
         end do
         synm_whole(channel, k) = disp
         dif = synm_whole(channel, k) - obse(channel, k)
         err2 = err2 + weight(channel, k) * dif * dif / 100.0
      end do
   end do
   err2 = sqrt((err2/weight_sum))
   err = real(err2)
   end subroutine static_synthetic
            

!  
! subroutine for removing the static response of current subfault for all static stations
!  
   subroutine static_remove_subfault(slip, rake, n_s, n_sub)
   real, intent(in) :: slip(max_subf, max_seg), rake(max_subf, max_seg)
   integer, intent(in) :: n_s, n_sub
   integer k, j, channel
   real disp, angle, sinal, cosal
!
   angle = rake(n_sub, n_s)*dpi
   sinal = sin(angle)
   cosal = cos(angle)

   do k = 1, 3
      j = 2*k-1
      do channel = 1, n_chan
         disp = slip(n_sub, n_s) &
     &   *(sinal*green(channel, j, n_sub, n_s)+cosal*green(channel, j+1, n_sub, n_s))
         synm_whole(channel, k) = synm_whole(channel, k)-disp
      end do
   end do

   end subroutine static_remove_subfault


!
! subroutine for testing the new response of the current subfault, for all stations
! we also give the misfit error of static data
!
   pure subroutine static_modify_subfault(slip, rake, n_s, n_sub, err)
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
   do k = 1, 3
      j = 2*k-1
      do channel = 1, n_chan
         disp = slip &
       & *(sinal*green(channel, j, n_sub, n_s)+cosal*green(channel, j+1, n_sub, n_s))
         dif = (synm_whole(channel, k) + disp) - obse(channel, k)
         err2 = err2 + weight(channel, k) * dif * dif / 100.0
      end do
   end do
   err2 = sqrt(err2/weight_sum)
   err = real(err2)

   end subroutine static_modify_subfault


!
! subroutine for asliping response of current subfault, for all stations.
! we also give the misfit error of static data
!
   subroutine static_add_subfault(slip, rake, n_s, n_sub, err)
   real, intent(in) :: slip(max_subf, max_seg), rake(max_subf, max_seg)
   real, intent(out) :: err
   integer, intent(in) :: n_s, n_sub
   integer k, j, channel
   real disp, angle, dif, sinal, cosal
   real*8 :: err2
!
   angle = rake(n_sub, n_s)*dpi
   sinal = sin(angle)
   cosal = cos(angle)

   err2 = 0.d0
   do k = 1, 3
      j = 2*k-1
      do channel = 1, n_chan
         disp = slip(n_sub, n_s) &
       & *(sinal*green(channel, j, n_sub, n_s)+cosal*green(channel, j+1, n_sub, n_s))
         synm_whole(channel, k) = synm_whole(channel, k)+disp
         dif = synm_whole(channel, k) - obse(channel, k)
         err2 = err2 + weight(channel, k) * dif * dif / 100.0
      end do
   end do
   err2 = sqrt(err2/weight_sum)
   err = real(err2)

   end subroutine static_add_subfault


end module static_data
