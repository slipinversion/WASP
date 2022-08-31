module get_stations_data


   use constants, only : max_stations, max_seg, npuse, wave_pts2, wave_pts, n_data
   use wavelet_param, only : lnpt, nlen, jmin, jmax, max_freq
   use wavelets, only : fourier_coefs, meyer_yamada, wavelet_obs
   implicit none
   integer, parameter :: nnsta_tele = 80
   real :: weight(max_stations), wave_obs(wave_pts2, max_stations), wmax(max_stations), dt_channel(max_stations)
   integer :: misfit_type(12, max_stations), t_max_val(max_stations)
   integer :: t_min(max_stations), t_max(max_stations), channels
   real :: wavelet_weight(12, max_stations)
   character(len=15) :: sta_name(max_stations)
   character(len=3) :: component(max_stations)
   integer :: mmm(max_stations), llove(max_stations), io_up(max_stations), idata(max_stations), &
         &     disp_or_vel(max_stations), used_data


contains


   subroutine get_data(strong, cgps, body, surf, dart)
   implicit none
!
!  Here, we load into memory, wavelet transform of observed data, and 
!       other properties of stations
!
   integer :: ll_in, ll_out
   logical :: strong, cgps, body, surf, dart

   write(*,*)'Get stations metadata and waveforms and store them in memory...'
   call fourier_coefs()
   call meyer_yamada()
   ll_in = 0
   ll_out = 0
   dt_channel(:) = 0.0
   wave_obs(:, :) = 0.0
   used_data = 0
   if (strong) then
      call get_strong_motion_stations(ll_in, ll_out)
      ll_in = ll_out
   end if
   if (cgps) then
      call get_cgps_stations(ll_in, ll_out)
      ll_in = ll_out
   end if
   if (body) then
      call get_body_waves_stations(ll_in, ll_out)
      ll_in = ll_out
   end if
   if (surf) then
      call get_surface_waves_stations(ll_in, ll_out)
      ll_in = ll_out
   end if
   if (dart) then
      call get_dart_stations(ll_in, ll_out)
      ll_in = ll_out
   end if
   channels = ll_out
   end subroutine get_data

   
   subroutine get_strong_motion_stations(ll_in, ll_out)
   implicit none
   integer ll_in, ll_out, ll_g, io_chan, j_con(11), k, ir, no, &
   &  n_wave_weight, ir_max, n_chan, &
   &  nos(max_stations), io_mod(max_stations), n_chan3
   real dt, weig(max_stations), j_wig(11), dt_sample, lat_s, lon_s
   logical, parameter :: cgps=.False.
   character(len=20) filename
!   character(len=6) sta_name(max_stations)
!   character(len=3) component(max_stations)

   write(*,*)'Get strong motion stations metadata and waveforms...'
   n_chan3 = 0
!
!  suppose the ni = u3e+11, then then moment of 1cm*1km^2 
!       subfault is 3e+21. The gfs is for Mo = 1e+20
!
   open(9,file='channels_strong.txt',status='old')
   open(15,file='wavelets_strong.txt', status='old')
   read(15,*) jmin, jmax, max_freq
   read(15,*)
   read(15,*) n_wave_weight

   read(9,*)
   read(9,*) 
   read(9,*) lnpt, dt_sample
   read(9,*)
   nlen = 2**lnpt
   dt = dt_sample

   read(9,*) ir_max, n_chan
   read(9,*)
   call error1(ll_in, n_chan)
   call error2(n_wave_weight, n_chan)
   ir = 0
   io_chan = 0
 
   do ir = 1, ir_max
      ll_g = ir+ll_in
      read(9,*) no, sta_name(ll_g), lat_s, lon_s, io_mod(ir), component(ll_g), weig(ir), nos(ir)
      if (weig(ir) .gt. 0) n_chan3 = n_chan3 + 1
   end do
   close(9)
!   write(*,*)'n_chan: ', n_chan3

   filename = 'waveforms_strong.txt'
   filename = trim(filename)
   call get_waveforms(filename, ir_max, ll_in, ll_out, cgps)

   do ir = 1, ir_max

      io_chan = io_chan+1
      ll_g = io_chan+ll_in
      weight(ll_g) = weig(ir) / n_chan3
      dt_channel(ll_g) = dt_sample
      
      do k = 1, 11
         misfit_type(k, ll_g) = 0
         wavelet_weight(k, ll_g) = 0
      end do
      read(15,*)
      read(15,*)(j_con(k), k = 1, jmax)
      read(15,*)(j_wig(k), k = 1, jmax)
      do k = jmin, jmax
         misfit_type(k, ll_g) = j_con(k)
         wavelet_weight(k, ll_g) = real(j_wig(k))
      end do
      do k = jmax+1, 11
         misfit_type(k, ll_g) = 0
         wavelet_weight(k, ll_g) = 0
      end do
   end do
   close(15)  
   ll_out = ll_in + n_chan
   call count_wavelets(ll_in, n_chan)
   end subroutine get_strong_motion_stations

   
   subroutine get_cgps_stations(ll_in, ll_out)
   implicit none
   integer ll_in, ll_out, ll_g, io_chan, j_con(11), k, ir, no, &
   &  n_wave_weight, ir_max, n_chan, &
   &  nos(max_stations), io_mod(max_stations), n_chan3
   real :: lat_s, lon_s, dt, weig(max_stations), j_wig(11), dt_sample
   logical, parameter :: cgps=.True.
   character(len=20) filename
!   character(len=6) sta_name(max_stations)
!   character(len=3) component(max_stations)

   write(*,*)'Get cGPS stations metadata and waveforms...'
   n_chan3 = 0
!
!  suppose the ni = u3e+11, then then moment of 1cm*1km^2 
!       subfault is 3e+21. The gfs is for Mo = 1e+20
!
   open(9, file='channels_cgps.txt', status='old')
   open(15, file='wavelets_cgps.txt', status='old')
   read(15,*) jmin, jmax, max_freq
   read(15,*)
   read(15,*) n_wave_weight

   read(9,*)
   read(9,*)
   read(9,*) lnpt, dt_sample
   read(9,*)
   nlen = 2**lnpt
   dt = dt_sample

   read(9,*) ir_max, n_chan
   read(9,*)
   call error1(ll_in, n_chan)
   call error2(n_wave_weight, n_chan)
   ir = 0
   io_chan = 0
 
   do ir = 1, ir_max
      ll_g = ir+ll_in
      read(9,*) no, sta_name(ll_g), lat_s, lon_s, io_mod(ir), component(ll_g), weig(ir), nos(ir)
      if (weig(ir) .gt. 0) n_chan3 = n_chan3 + 1
   end do
   close(9)

   filename = 'waveforms_cgps.txt'
   filename = trim(filename)
   call get_waveforms(filename, ir_max, ll_in, ll_out, cgps) 

   do ir = 1, ir_max

      io_chan = io_chan+1
      ll_g = io_chan+ll_in
      weight(ll_g) = weig(ir) / n_chan3
      dt_channel(ll_g) = dt_sample
      
      do k = 1, 11
         misfit_type(k, ll_g) = 0
         wavelet_weight(k, ll_g) = 0
      end do
      read(15,*)
      read(15,*)(j_con(k), k = 1, jmax)
      read(15,*)(j_wig(k), k = 1, jmax)
      do k = jmin, jmax
         misfit_type(k, ll_g) = j_con(k)
         wavelet_weight(k, ll_g) = real(j_wig(k))
      end do
      do k = jmax+1, 11
         misfit_type(k, ll_g) = 0
         wavelet_weight(k, ll_g) = 0
      end do
   end do
   close(15)
   ll_out = ll_in + n_chan
   call count_wavelets(ll_in, n_chan)
   end subroutine get_cgps_stations

   
   subroutine get_body_waves_stations(ll_in, ll_out)
   implicit none
   integer ll_in, ll_out, ll_g, j_con(11), k, ir, n_wave_weight, n_chan, &
   &  nos, n_chan3, idts, nstaon, love, int1
   real lat_sta, lon_sta, j_wig(11), dt, & 
   &  rang, az, earth_angle, float1, float2
   logical, parameter :: cgps=.False.
   character(len=20) filename
   character(len=6)  earth, sttyp 
   character(len=14) fname

   write(*,*)'Get body waves stations metadata and waveforms...'
   n_chan3 = 0
!
!  suppose the ni = u3e+11, then then moment of 1cm*1km^2 
!       subfault is 3e+21. The gfs is for Mo = 1e+20
!
   open(9, file='channels_body.txt', status='old')
   open(14, file='body_wave_weight.txt', status='old')
   if (max_freq .gt. npuse) then
      write(*,*)'You should stop and check dimension sww,sws'
      stop
   end if
   read(9,*)
   read(9,*) idts, lnpt, dt
   read(9,*)
   read(9,*) nstaon
   call error1(ll_in, nstaon)
   nlen = 2 ** lnpt

   open(15, file='wavelets_body.txt', status='old')
   read(15,*) jmin, jmax, max_freq
   read(15,*)
   read(15,*) n_wave_weight
   do ir = 1, nstaon
      ll_g = ll_in+ir
      read(14,*) weight(ll_g)
      if (weight(ll_g) .gt. 0) n_chan3 = n_chan3 + 1
   end do
   rewind 14
!   write(*,*)'n_chan: ', n_chan3
   n_chan = nstaon 
   call error2(n_wave_weight, nstaon)
   
   filename = 'waveforms_body.txt'
   filename = trim(filename)
   call get_waveforms(filename, nstaon, ll_in, ll_out, cgps)

   do ir = 1, nstaon
      ll_g = ir+ll_in
      read(9,*) nos, earth, sttyp, sta_name(ll_g), fname, &
      & rang, az, lat_sta, lon_sta, earth_angle, float1, &
      & mmm(ir), disp_or_vel(ir), float2, llove(ir), int1, idata(ir)
      if (idata(ir) .gt. 0 .or. mmm(ir) .eq. 3) cycle
      dt_channel(ll_g) = dt
      read(14,*) weight(ll_g)
      weight(ll_g) = weight(ll_g) / n_chan3
      love = llove(ir)
      component(ll_g) = 'P'
      if (love .gt. 0) component(ll_g) = 'SH'
      
      do k = 1, 11
         misfit_type(k, ll_g) = 0
         wavelet_weight(k, ll_g) = 0
      end do
      read(15,*)
      read(15,*)(j_con(k), k = 1, jmax)
      read(15,*)(j_wig(k), k = 1, jmax)
      do k = jmin, jmax
         misfit_type(k, ll_g) = j_con(k)
         wavelet_weight(k, ll_g) = real(j_wig(k))
      end do
      do k = jmax+1, 11
         misfit_type(k, ll_g) = 0
         wavelet_weight(k, ll_g) = 0
      end do
   end do
   close(9)
   close(14)
   close(15)
   ll_out = ll_in + n_chan
   call count_wavelets(ll_in, n_chan)
   end subroutine get_body_waves_stations


   subroutine get_surface_waves_stations(ll_in, ll_out)
   implicit none
   integer ll_in, ll_out, ll_g, io_chan, j_con(11), k, ir, no, &
   &  n_wave_weight, ir_max, n_chan, &
   &  nos(max_stations), io_mod(max_stations), n_chan3, &
   &  io_ns(max_stations), idts, io_ew(max_stations), io_str(max_stations)
   real lat_s(max_stations), lon_s(max_stations), dt, &
   &  ang_ns(max_stations), ang_ew(max_stations), weig(max_stations, 3), j_wig(11), dt_sample, &
   &  dip, rake, theta
   logical, parameter :: cgps=.False.
   character(len=20) filename
!   character(len=6) sta_name(max_stations)
   character(len=250) modes

   write(*,*)'Get stations metadata and waveforms for long period surface waves...'
   n_chan3 = 0
!
!  suppose the ni = u3e+11, then then moment of 1cm*1km^2 
!       subfault is 3e+21. The gfs is for Mo = 1e+20
!
 
   open(9, file='channels_surf.txt', status='old')

   open(15, file='wavelets_surf.txt', status='old')
   read(15,*) jmin, jmax, max_freq
   read(15,'(a)')modes
   read(15,*) n_wave_weight
   idts = 3
 
   read(9,*)
   read(9,*)
   read(9,*) dip, rake, theta, lnpt, dt_sample
   read(9,*)
   if (lnpt.ne.10) then
      write(*,*)"please check input LNPT"
   end if
   nlen = 2**lnpt
   dt = dt_sample
  
   read(9,*) ir_max, n_chan
   read(9,*)
   call error1(ll_in, n_chan)
   call error2(n_wave_weight, n_chan)
!       
!       Here we read the green functions of long period surface waves
!  
   do ir = 1, ir_max    
      ll_g = ir+ll_in
      read(9,*) no, sta_name(ll_g), lat_s(ir), lon_s(ir), io_mod(ir), &
      & io_up(ir), io_ns(ir), io_ew(ir), ang_ns(ir), ang_ew(ir), &
      & io_str(ir),(weig(ir, k), k = 1, 3), nos(ir)
      if ((weig(ir, 1) .gt. 0) .and. (io_up(ir) .eq. 1)) n_chan3 = n_chan3 + 1
      if ((weig(ir, 2) .gt. 0) .and. (io_ns(ir) .eq. 1)) n_chan3 = n_chan3 + 1
   end do
   close(9)

   io_chan = 0
   do ir = 1, n_chan
      if (io_up(ir) .eq. 1) then
         io_chan = io_chan + 1
         ll_g = io_chan + ll_in
         weight(ll_g) = weig(ir, 1) / n_chan3
         component(ll_g) = 'R'
      end if
      if (io_ns(ir) .eq. 1) then
         io_chan = io_chan + 1
         ll_g = io_chan + ll_in
         weight(ll_g) = weig(ir, 2) /n_chan3
         component(ll_g) = 'L'
      end if
   end do
   
   filename = 'waveforms_surf.txt'
   filename = trim(filename)
   call get_waveforms(filename, n_chan, ll_in, ll_out, cgps) 

   io_chan = 0
   do ir = 1, n_chan
      io_chan = io_chan+1
      ll_g = io_chan+ll_in
      dt_channel(ll_g) = dt_sample
      
      do k = 1, 11
         misfit_type(k, ll_g) = 0
         wavelet_weight(k, ll_g) = 0
      end do
      read(15,*)
      read(15,*)(j_con(k), k = 1, jmax)
      read(15,*)(j_wig(k), k = 1, jmax)
      do k = jmin, jmax
         misfit_type(k, ll_g) = j_con(k)
         wavelet_weight(k, ll_g) = real(j_wig(k))
      end do
      do k = jmax+1, 11
         misfit_type(k, ll_g) = 0
         wavelet_weight(k, ll_g) = 0
      end do
   end do
   close(15)   
   ll_out = ll_in+n_chan
   call count_wavelets(ll_in, n_chan)
   end subroutine get_surface_waves_stations


   subroutine get_dart_stations(ll_in, ll_out)
   implicit none
   integer ll_in, ll_out, ll_g, j_con(11), k, ir, no, &
   &  n_wave_weight, ir_max, n_chan, nos, io_mod, n_chan3
   real :: lat_s, lon_s, dt, weig(max_stations), j_wig(11), dt_sample
   logical, parameter :: cgps=.False.
   character(len=20) filename
!   character(len=6) sta_name(max_stations)
!   character(len=3) component(max_stations)

   n_chan3 = 0
!
!  suppose the ni = u3e+11, then then moment of 1cm*1km^2 
!       subfault is 3e+21. The gfs is for Mo = 1e+20
!
   open(9, file='channels_dart.txt', status='old')
   open(15, file='wavelets_dart.txt', status='old')
   read(15,*) jmin, jmax, max_freq
   read(15,*)
   read(15,*) n_wave_weight

   read(9,*)
   read(9,*)
   read(9,*) lnpt, dt_sample
   read(9,*)
   nlen = 2**lnpt
   dt = dt_sample

   read(9,*) ir_max, n_chan
   read(9,*)
   call error1(ll_in, n_chan)
   call error2(n_wave_weight, n_chan)
 
   do ir = 1, ir_max
      ll_g = ir+ll_in
      read(9,*) no, sta_name(ll_g), lat_s, lon_s, io_mod, component(ll_g), weig(ir), nos
      if (weig(ir) .gt. 0) n_chan3 = n_chan3 + 1
   end do
   close(9)
   
   filename = 'waveforms_dart.txt'
   filename = trim(filename)
   call get_waveforms(filename, ir_max, ll_in, ll_out, cgps) 

   do ir = 1, ir_max

      ll_g = ir+ll_in
      weight(ll_g) = weig(ir) / n_chan3
      dt_channel(ll_g) = dt_sample
      
      do k = 1, 11
         misfit_type(k, ll_g) = 0
         wavelet_weight(k, ll_g) = 0
      end do
      read(15,*)
      read(15,*)(j_con(k), k = 1, jmax)
      read(15,*)(j_wig(k), k = 1, jmax)
      do k = jmin, jmax
         misfit_type(k, ll_g) = j_con(k)
         wavelet_weight(k, ll_g) = real(j_wig(k))
      end do
      do k = jmax+1, 11
         misfit_type(k, ll_g) = 0
         wavelet_weight(k, ll_g) = 0
      end do
   end do
   close(15)
   ll_out = ll_in + n_chan
   call count_wavelets(ll_in, n_chan)
   end subroutine get_dart_stations


   subroutine get_waveforms(filename, ir_max, ll_in, ll_out, cgps)
   implicit none
   integer ll_in, ll_out, ll_g, nm, i, ir, no, &
   &  ir_max
   real cr(wave_pts2), cz(wave_pts2), obser(n_data), &
   &  dto, amp_max, obse(n_data, max_stations), mean
   logical :: cgps
   character(len=40) string
   character(len=20) filename
  
   filename = trim(filename)
   open(13, file=filename, status='old')
   do ir = 1, ir_max
      ll_g = ir+ll_in
      t_min(ll_g) = 0
      read(13,*)
      read(13,*)
      read(13,*) string, dto
      read(13,*) string, no
      read(13,*) string, t_max(ll_g)
      read(13,*)
      read(13,*)(obse(i, ll_g), i = 1, no)
      mean = sum(obse(no - 20:no, ll_g)) / 20.0
!      no = min(2**lnpt - 20, no)
      do i = 1, wave_pts2
         cz(i) = 0.0
         cr(i) = 0.0
         if (i .lt. no) then
            obser(i) = obse(i, ll_g)
         else
            obser(i) = 0.0
            if (cgps) obser(i) = mean 
         endif
         
!         if (cgps) then
!            if (i .lt. no) then
!               obser(i) = obse(i, ll_g)
!            elseif (i .lt. (2**lnpt - 20)) then
!               obser(i) = mean
!            else
!               obser(i) = obse(1, ll_g)
!            endif
!         else
!            if (i .lt. no) then
!               obser(i) = obse(i, ll_g)
!            else
!               obser(i) = obse(1, ll_g) 
!            endif
!         end if
      end do
      do i = 1,nlen
         cr(i) = obser(i)
         cz(i) = 0.0
      enddo

      call wavelet_obs(cr, cz, obser)
      amp_max = 0.0
      nm = 1
      if (cgps .eqv. .False.) then
         do i = 1, nlen
            if (amp_max .lt. abs(obser(i))) then
               amp_max = abs(obser(i))
               nm = i
            end if
         end do
      else
         do i = 1, nlen
            if (amp_max .lt. abs(obser(i))) then
               amp_max = abs(obser(i))
               nm = i
            end if
         end do
      endif

      wmax(ll_g) = amp_max
      t_max_val(ll_g) = nm
      do i = 1, nlen
         wave_obs(i, ll_g) = obser(i)/(amp_max)
      end do
   enddo
   close(13)
   ll_out = ll_in + ir_max
   end subroutine get_waveforms


   subroutine count_wavelets(ll_in, ir_max)
   implicit none
   integer :: k, j, ll_in, ir, ir_max, ll_g, n_begin, n_delt, length
   do ir = 1, ir_max
      ll_g = ir+ll_in
      if (weight(ll_g) .gt. (1e-3)) then
         do j = jmin, jmax
            n_begin = 2**(j-1)
            n_delt = nlen/n_begin
            length = int(t_max(ll_g)/n_delt+0.5)-1
            if (wavelet_weight(j, ll_g) .gt. (1e-3)) then
               used_data = used_data + length
            endif
         end do
      endif
   end do
   end subroutine count_wavelets


   subroutine error1(ll_in, n_chan)
   implicit none
   integer :: ll_in, n_chan
   if (n_chan + ll_in .gt. max_stations) then
      write(*,*)'Error: Maximum allowed number of channels: ', max_stations
      write(*,*)'But amount of channels used for modelling is at least', n_chan + ll_in
      write(*,*)'Please check maximum allowed number of channels'
      stop
   end if
   end subroutine error1
   

   subroutine error2(n_wave_weight, n_chan)
   implicit none
   integer :: n_wave_weight, n_chan
   if (n_wave_weight.ne.n_chan) then
      write(0,*)'n_wave_chan n_chan',n_wave_weight,n_chan
      write(0,*)'Mismatch in amount of stations betweeen wavelet file and stations file'
      stop
   end if
   end subroutine error2


end module get_stations_data      
