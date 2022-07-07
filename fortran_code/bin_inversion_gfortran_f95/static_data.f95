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


   use constants, only : max_subf, max_seg, max_subfaults, dpi
   use model_parameters, only : nxs_sub, nys_sub, nx_p, ny_p, segments, subfaults
   implicit none
   integer, parameter, private :: n_stations = 500
   integer :: n_chan
   real*8 :: synm_whole(n_stations, 3), weight_sum
   real :: lat(n_stations), lon(n_stations), max_gps
   real :: green(n_stations, 6, max_subfaults), syn_disp(n_stations, 3)
   real :: obse(n_stations, 3), weight(n_stations, 3)
   character(len=6) :: sta_name(n_stations)
  

contains


   subroutine initial_gps(slip, rake)
   implicit none
   integer k, j, segment, channel, no, i
   integer iys, ixs, n_tt, subfault
   real slip(:), rake(:)
   real :: cosal, sinal, angle
   real*8 :: disp
   logical is_file
!
   inquire( file = 'static_data.txt', exist = is_file )
   if (is_file) then 
      weight_sum = 0.d0
      open(9, file='static_data.txt', status='old')
      read(9,*) n_chan
      read(9,*)
      do i = 1, n_chan
         read(9,*) no, sta_name(i), lat(i), lon(i), (obse(i, k), k = 1, 3), (weight(i, k), k = 1, 3)
         do k = 1, 3
            weight_sum = weight_sum + weight(i, k)
         end do
      end do
      max_gps = maxval(abs(obse(:, :)))
      close(9) 

      open(33, file='Green_static_subfault.txt', status='old')
      read(33,*) n_tt
      do channel = 1, n_chan
         read(33,*)
         subfault = 0
         do segment = 1, segments
            read(33,*)
            do j = 1, nys_sub(segment)*nxs_sub(segment)
               subfault = subfault + 1
               read(33,*)(green(channel, k, subfault), k = 1, 6)
            enddo
         end do
      end do
      close(33)
      
      do k = 1, 3
         j = 2 * k - 1
         do channel = 1, n_chan
            disp = 0.d0
            do subfault = 1, subfaults
               angle = rake(subfault)*dpi 
               sinal = sin(angle)
               cosal = cos(angle)
               disp = disp + slip(subfault) &
     &  *(sinal*green(channel, j, subfault)+cosal*green(channel, j+1, subfault))
            end do
            syn_disp(channel, k) = disp
         end do
      end do 
      open(10,file='static_synthetics.txt')
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
   subroutine static_synthetic(slip, rake, err)
   implicit none
   real slip(:), rake(:)!, err
   integer k, j, segment, channel, subfault
   real err, dif, angle, sinal, cosal
   real*8 :: disp, err2
   
   err2 = 0.d0
   synm_whole(:, :) = 0.d0
   do k = 1, 3
      j = 2 * k - 1
      do channel = 1, n_chan
         disp = 0.d0
         do subfault = 1, subfaults
            angle = rake(subfault)*dpi 
            sinal = sin(angle)
            cosal = cos(angle)
            disp = disp + slip(subfault) &
       &  *(sinal*green(channel, j, subfault)+cosal*green(channel, j+1, subfault))
         end do
         synm_whole(channel, k) = disp
         dif = synm_whole(channel, k) - obse(channel, k)
         err2 = err2 + weight(channel, k) * dif * dif / max_gps / max_gps!100.0
      end do
   end do
   err2 = sqrt((err2/weight_sum))
   err = real(err2)
   end subroutine static_synthetic
            

!  
! subroutine for removing the static response of current subfault for all static stations
!  
   subroutine static_remove_subfault(slip, rake, subfault)
   implicit none
   real, intent(in) :: slip, rake
   integer, intent(in) :: subfault
   integer k, j, channel
   real disp, angle, sinal, cosal
!
   angle = rake*dpi
   sinal = sin(angle)
   cosal = cos(angle)

   do k = 1, 3
      j = 2*k-1
      do channel = 1, n_chan
         disp = slip &
     &   *(sinal*green(channel, j, subfault)+cosal*green(channel, j+1, subfault))
         synm_whole(channel, k) = synm_whole(channel, k)-disp
      end do
   end do

   end subroutine static_remove_subfault


!
! subroutine for testing the new response of the current subfault, for all stations
! we also give the misfit error of static data
!
   pure subroutine static_modify_subfault(slip, rake, subfault, err)
   implicit none
   real, intent(in) :: slip, rake
   real, intent(out) :: err
   integer, intent(in) :: subfault
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
       & *(sinal*green(channel, j, subfault)+cosal*green(channel, j+1, subfault))
         dif = (synm_whole(channel, k) + disp) - obse(channel, k)
         err2 = err2 + weight(channel, k) * dif * dif / max_gps / max_gps!100.0
      end do
   end do
   err2 = sqrt(err2/weight_sum)
   err = real(err2)

   end subroutine static_modify_subfault


!
! subroutine for asliping response of current subfault, for all stations.
! we also give the misfit error of static data
!
   subroutine static_add_subfault(slip, rake, subfault)
   implicit none
   real, intent(in) :: slip, rake
   integer, intent(in) :: subfault
   integer k, j, channel
   real disp, angle, dif, sinal, cosal
   
   angle = rake*dpi
   sinal = sin(angle)
   cosal = cos(angle)

   do k = 1, 3
      j = 2*k-1
      do channel = 1, n_chan
         disp = slip &
       & *(sinal*green(channel, j, subfault)+cosal*green(channel, j+1, subfault))
         synm_whole(channel, k) = synm_whole(channel, k)+disp
      end do
   end do

   end subroutine static_add_subfault


end module static_data
