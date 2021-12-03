!
! TODO: change to velocity if io_v_d = 1
!
module retrieve_gf


   use constants, only : npuse, max_stations, max_subfaults, max_subf, nt, nny, ndis, &
   &  wave_pts2, max_seg, block_far, ltde, n_data, max_psources, max_dip_psources, &
   &  max_dip_subfaults, block_stg, pi, twopi, wave_pts
   use model_parameters, only : point_sources, segments, v_ref, nxs_sub, nys_sub, delay_seg, &
   &  dxs, dys, nx_p, ny_p, dip, strike, shear, t_latest
   use wavelets, only : fourier_coefs, meyer_yamada
   use wavelet_param, only : lnpt, nlen, jmin, jmax, max_freq
   use retrieve_surf_gf, only : get_surf_gf_data, interp_gf, get_surf_gf, &
   &  npt_bank, dt_bank, check_bounds
   use rad_pattern, only : rad_coef
   use geodesics, only : distaz
   use get_stations_data, only : sta_name1, sta_name2, sta_name3, disp_or_vel, &
          & component1, component2, idata, mmm, llove, dt_channel
   implicit none
   integer, parameter :: nnsta_tele = 80
   complex, allocatable :: green_dip(:, :, :), green_stk(:, :, :)


contains


   subroutine get_gf(strong, cgps, body, surf, dart)
   implicit none
!
   integer ll_in, ll_out
   logical :: strong, cgps, body, surf, dart
   allocate(green_dip(npuse, max_stations, max_subfaults))
   allocate(green_stk(npuse, max_stations, max_subfaults))

   write(*,*)'Store GF in memory...'
   call fourier_coefs()
   call meyer_yamada()
!
!  Here, we load into memory the green functions for each subfault, for every used station
!  
   ll_in = 0
   ll_out = 0
   if (strong) then
      call get_strong_motion_gf(ll_in, ll_out)
      ll_in = ll_out
   end if
   if (cgps) then
      call get_cgps_gf(ll_in, ll_out)
      ll_in = ll_out
   end if
   if (body) then
      call get_body_waves_gf(ll_in, ll_out)
      ll_in = ll_out
   end if
   if (surf) then
      call get_surface_waves_gf(ll_in, ll_out)
      ll_in = ll_out
   end if
!   if (dart) then
!      call get_dart_gf(ll_in, ll_out)
!      ll_in = ll_out
!   end if
   end subroutine get_gf

   
   subroutine get_strong_motion_gf(ll_in, ll_out)
   implicit none
   integer :: ll_in, ll_out, io_v_d, ll_g, ll, &
   &  io_chan, i, segment, channel, channel_max, n_chan, &
   &  ixs, iys
   real :: omega, block, dt, df, dt_sample, w, tlen
   complex :: z0, z
   character(len=80) filename
   character(len=3) comp!component(max_stations), comp
   character(len=1) channel2!component(max_stations), comp
   
   write(*,*)'Store strong motion GF in memory...'
   z0 = cmplx(0.0, 0.0)
!
!  suppose the ni = u3e+11, then then moment of 1cm*1km^2 
!       subfault is 3e+21. The gfs is for Mo = 1e+20
!
   block = dxs * dys * (1.e-10) 

   open(9, file='channels_strong.txt', status='old')

   read(9,*)
   read(9,*)
   read(9,*) lnpt, dt_sample
   read(9,*) io_v_d
   nlen = 2**lnpt
   dt = dt_sample
   tlen = dt*nlen
   df = 1.0/tlen

   read(9,*) channel_max, n_chan
   close(9)
   channel = 0
   io_chan = 0
!       
!       Here we read the green functions of strong motion waves
!
   df = 1. / ((2.0 ** lnpt) * dt)
   nlen = 2 ** lnpt
   tlen = nlen * dt

   do channel = 1, channel_max

      comp = component1(channel)
      channel2 = comp(3:3)

      io_chan = io_chan+1
      ll_g = io_chan+ll_in
      if (channel2 .eq. 'Z') filename = trim(sta_name1(channel))//'1'
      if (channel2 .eq. 'N') filename = trim(sta_name1(channel))//'2'
      if (channel2 .eq. 'E') filename = trim(sta_name1(channel))//'3'
      filename = trim(filename)

      open(12, file=filename, status='old', access='direct',recl=block_stg)
      ll = 0
!
!       Here, we read the green functions and derivate them
!
      ll_g = ll_in+io_chan
      do segment = 1, segments
         do iys = 1, nys_sub(segment)
            do ixs = 1, nxs_sub(segment)
               ll = ll+1
               read(12, rec = ll) &
               & (green_dip(i, ll_g, ll), i = 1, max_freq),(green_stk(i, ll_g, ll), i = 1, max_freq)
               do i = 1, max_freq
!
! we eventually shift synthetics in time, in case the fault plane used has a delay
!
                  omega = twopi*(i-1)*df
                  w = -omega*delay_seg(segment)
                  z = cmplx(0.0, w)
                  z = cexp(z)
                  green_dip(i, ll_g, ll) = green_dip(i, ll_g, ll)*z*block
                  green_stk(i, ll_g, ll) = green_stk(i, ll_g, ll)*z*block
               end do
            end do
         end do
      end do
      close(12)
   end do      
   ll_out = ll_in+n_chan
   end subroutine get_strong_motion_gf

   
   subroutine get_cgps_gf(ll_in, ll_out)
   implicit none
   integer ll_in, ll_out, ll_g, ll, &
   &  io_chan, i, segment, channel, &
   &  channel_max, n_chan, ixs, iys!, jmin, jmax, max_freq
   real omega, block, dt, df, dt_sample, w, low_freq, tlen
   complex z0, z
   character(len=80) filename
   character(len=3) comp
   character(len=1) channel2

   write(*,*)'Store cGPS GF in memory...'
   z0 = cmplx(0.0, 0.0)
!
!  suppose the ni = u3e+11, then then moment of 1cm*1km^2 
!       subfault is 3e+21. The gfs is for Mo = 1e+20
!
   block = dxs * dys * (1.e-10) 

   open(9, file='channels_cgps.txt', status='old')

   read(9,*)
   read(9,*) 
   read(9,*) lnpt, dt_sample
   read(9,*)
   nlen = 2**lnpt
   dt = dt_sample
   tlen = dt*nlen
   df = 1.0/tlen

   read(9,*) channel_max, n_chan
   close(9)
   channel = 0
   io_chan = 0
!       
!       Here we read the green functions of strong motion waves
!
   df = 1. / ((2.0 ** lnpt) * dt)
   nlen = 2 ** lnpt
   tlen = nlen * dt
   low_freq = 1.0 / 200.0

   do channel = 1, channel_max

      comp = component2(channel)
      channel2 = comp(3:3)

      io_chan = io_chan+1
      ll_g = io_chan+ll_in
      if (channel2 .eq. 'Z') filename = trim(sta_name2(channel))//'.cgps.1'
      if (channel2 .eq. 'N') filename = trim(sta_name2(channel))//'.cgps.2'
      if (channel2 .eq. 'E') filename = trim(sta_name2(channel))//'.cgps.3'
      filename = trim(filename)

      open(12, file=filename, status='old', access='direct', recl=block_stg)
      ll = 0
!
!       Here, we read the green functions and derivate them
!
      ll_g = ll_in+io_chan
      do segment = 1, segments
         do iys = 1, nys_sub(segment)
            do ixs = 1, nxs_sub(segment)
               ll = ll+1
               read(12, rec = ll)(green_dip(i, ll_g, ll), i = 1, max_freq),(green_stk(i, ll_g, ll), i = 1, max_freq)
               do i = 1, max_freq
                  omega = twopi*(i-1)*df
                  w = -omega*delay_seg(segment) ! we shift synthetics in time by delay_seg(segment) seconds
                  z = cmplx(0.0, w)
                  z = cexp(z)
                  green_dip(i, ll_g, ll) = green_dip(i, ll_g, ll)*z*block
                  green_stk(i, ll_g, ll) = green_stk(i, ll_g, ll)*z*block
               end do
            end do
         end do
      end do
      close(12)
   end do
   ll_out = ll_in+n_chan
   end subroutine get_cgps_gf


   subroutine get_body_waves_gf(ll_in, ll_out)
   implicit none
   integer nstaon, channel, ll_g, k, nsta, n_chan, &
   &  love, ll, segment, iys, iyp, io_seg, iys_c, iy_c, jf, i, ipy, npxy, &
   &  nkxy, nxs_c, nys_c, nxp_c, nyp_c, ll_s, kxy, ixs, kpxy, ixp, & 
!, iud(max_stations_tele), &
!   &  idata(max_stations_tele), nos(max_stations_tele), llove(max_stations_tele), mmm(max_stations_tele), &
   &  ll_in, ll_out
   real dt, df, ddelt, time, block, w, tlen!, &
!   &  earth_angle(max_stations_tele), hcru(max_stations_tele), &
!   &  disp_or_vel(max_stations_tele), rang(max_stations_tele), az(max_stations_tele), &
!   &  ttvl(max_stations_tele), &
!   &  lat_sta(max_stations_tele), lon_sta(max_stations_tele)
   real, allocatable :: tdel(:,:,:,:)
   complex :: z, z0, wsyn(wave_pts)
   complex :: kahan_y1(wave_pts2), kahan_t1(wave_pts2), kahan_c1(wave_pts2)
   complex :: kahan_y2(wave_pts2), kahan_t2(wave_pts2), kahan_c2(wave_pts2)
   complex, allocatable :: green_dip0(:,:,:,:,:)
   complex, allocatable :: green_stk0(:,:,:,:,:)
!
!   character(len=5) earth(max_stations_tele)
!   character(len=14) fname(max_stations_tele)
!   character(len=6) stname(max_stations_tele), sttyp(max_stations_tele)
   character(len=14) fname4, fname6
   allocate(green_dip0(npuse, max_dip_psources, max_dip_subfaults, max_seg, nnsta_tele))
   allocate(green_stk0(npuse, max_dip_psources, max_dip_subfaults, max_seg, nnsta_tele))
   allocate(tdel(max_psources, max_subf, max_seg, nnsta_tele))
   
   write(*,*)'Store body waves GF in memory...'
   z0 = cmplx(0.0, 0.0)

   open(9, file='channels_body.txt', status='old')
!   open(13,file='Obser.tele',status='unknown')
   read(9,*)
   read(9,*)
   read(9,*)
   read(9,*) nstaON
   block = dxs*dys*(1.0e-10)
   block = block*1.e+4
!  
!  Because the data is micrometer, so The amplitude of block should
!  time e+4 Only in this test
!
!   open(15,file='Wave.tele')
!   read(15,*)! jmin, jmax, max_freq
!   read(15,*)
!   read(15,*) n_wave_weight

   dt = dt_channel(ll_in + 1)
   tlen = nlen*dt
   npxy = nx_p*ny_p
   do channel = 1, nstaon
      ll_g = ll_in+channel

      if (idata(channel) .gt. 0 .or. mmm(channel) .eq. 3) cycle
      nsta = nstaON
      n_chan = nstaon
      love = llove(channel)
!
!  MM = 0 FOR FAR; MM = 1 FOR UPPER; MM = 3 FOR PNL
!  Here we read the green functions of the teleseismic body waves
!  
      if (love .eq. 0) then
         fname4 = trim(sta_name3(channel))//'.GRE'
         fname6 = trim(sta_name3(channel))//'.TDE'
      else
         fname4 = trim(sta_name3(channel))//'SH.GRE'
         fname6 = trim(sta_name3(channel))//'SH.TDE'
      end if
      open(12, file=fname4, status='old', access='direct', recl=block_far)
      open(32, file=fname6, status='old', access='direct', recl=ltde)
      LL = 0
      do segment = 1, segments
         if (segment .le. 10) then
            do iys = 1, nys_sub(segment)
               do IPY = 1, ny_p
                  LL = LL+1
                  read(12, rec = LL) io_seg, iys_c, iy_c, JF, DT, DF, &
     &  (green_dip0(i, ipy, iys, segment, channel), i = 1, 2*max_freq)
                  LL = LL+1
                  read(12, rec = LL) io_seg, iys_c, iy_c, JF, LNPT, NLEN, &
     &  (green_stk0(i, ipy, iys, segment, channel), I = 1, 2*max_freq)
!
!       Sanity checks
!
                  if ((io_seg.ne.segment) .or. (iys_c.ne.iys) .or. (iy_c.ne.ipy)) then
                     write(*,*)'io_seg vs segment: ', io_seg,segment
                     write(*,*)'iys_c, iys: ', iys_c,iys
                     write(*,*)'iy_c, ipy: ', iy_c,ipy
                     write(*,*)"Green function is not matched with fault model"
                     stop
                  end if
                  do i = 1, 2*max_freq
                     green_dip0(i, ipy, iys, segment, channel) = &
                     & block*green_dip0(i, ipy, iys, segment, channel)
                     green_stk0(i, ipy, iys, segment, channel) = &
                     & block*green_stk0(i, ipy, iys, segment, channel)
                  end do
               end do
            end do
            nkxy = nxs_sub(segment)*nys_sub(segment)
            read(32, rec = segment) nxs_c, nys_c, nxp_c, nyp_c, &
     &  ((tdel(k, ll_s, segment, channel), k = 1, npxy), ll_s = 1, nkxy)
            if ((nxs_c.ne.nxs_sub(segment)) .or. (nys_c.ne.nys_sub(segment)) &
     &      .or.(nxp_c.ne.nx_p) .or. (nyp_c.ne.ny_p)) then
                write(*,*)'nxs',nxs_c,nxs_sub(segment)
                write(*,*)'nys',nys_c,nys_sub(segment)
                write(*,*)'nxp',nxp_c,nx_p
                write(*,*)'nyp',nyp_c,ny_p
                write(*,'(a)')'Mismatch in amount of point sources or subfaults &
                &between the specified in Fault.time, and those used in the &
                &green functions.'
                stop
             end if
          else
             npxy = nx_p*ny_p
             nkxy = nxs_sub(segment)*nys_sub(segment)
             do iys = 1, nys_sub(segment)
                do ipy = 1, ny_p
                   do i = 1, 2*max_freq
                      green_dip0(i, ipy, iys, segment, channel) = z0
                      green_stk0(i, ipy, iys, segment, channel) = z0
                  end do
               end do
            end do
            do ll_s = 1, nkxy
               do k = 1, npxy
                  tdel(k, ll_s, segment, channel) = 0.0
               end do
            end do
         end if
      end do
      close(12)
      close(32)
   end do
   close(9)
   
   do channel = 1, nstaon
      ll_g = ll_in+channel

      do I = 1, wave_pts
         wsyn(i) = z0
      end do
      LL = 0
      ddelt = 0.0
      z = z0
      do segment = 1, segments
         kxy = 0
         do iys = 1, nys_sub(segment)
            do ixs = 1, nxs_sub(segment)
               kxy = kxy+1
               ll = ll+1
               do i = 1, 2*max_freq
                  green_dip(i, ll_g, ll) = z0
                  green_stk(i, ll_g, ll) = z0
               end do
               kahan_y1(:) = z0
               kahan_t1(:) = z0
               kahan_c1(:) = z0
               kahan_y2(:) = z0
               kahan_t2(:) = z0
               kahan_c2(:) = z0
               kpxy = 0
               do iyp = 1, ny_p
                  do ixp = 1, nx_p
                     kpxy = kpxy+1
                     time = point_sources(4, ixp, iyp, ixs, iys, segment)/v_ref+delay_seg(segment)
                     time = time + tdel(kpxy, kxy, segment, channel)
                     do i = 1, 2*max_freq
                        W = -twopi*(i-1)*df*time
                        z = cmplx(0.0, w)
                        z = cexp(z)
                        kahan_y1(i) = green_dip0(i, iyp, iys, segment, channel)*z - kahan_c1(i)
                        kahan_t1(i) = green_dip(i, ll_g, ll) + kahan_y1(i)
                        kahan_c1(i) = (kahan_t1(i) - green_dip(i, ll_g, ll)) - kahan_y1(i)
                        green_dip(i, ll_g, ll) = kahan_t1(i)
                        kahan_y2(i) = green_stk0(i, iyp, iys, segment, channel)*z - kahan_c2(i)
                        kahan_t2(i) = green_stk(i, ll_g, ll) + kahan_y2(i)
                        kahan_c2(i) = (kahan_t2(i) - green_stk(i, ll_g, ll)) - kahan_y2(i)
                        green_stk(i, ll_g, ll) = kahan_t2(i)
                     end do
                  end do
               end do
               do i = 1, 2*max_freq
                  if (disp_or_vel(channel) .eq. 0) then
                     z = cmplx(1.0, 0.0)
                  elseif (disp_or_vel(channel) .eq. 1) then
                     w = twopi*(i-1)*df
                     z = cmplx(0.0, w)
                  end if
                  green_dip(i, ll_g, ll) = z*green_dip(i, ll_g, ll)/npxy
                  green_stk(i, ll_g, ll) = z*green_stk(i, ll_g, ll)/npxy
               end do
            end do
         end do
      end do
   end do
   ll_out = ll_in+nstaon
   deallocate(green_stk0)
   deallocate(green_dip0)
   deallocate(tdel)
   end subroutine get_body_waves_gf


   subroutine get_surface_waves_gf(ll_in, ll_out)
   implicit none
   integer ll_in, ll_out, nf1, nf2, nf3, nf4, no, &
   &  npp, ix, iy, segment_subfault, j, ll_g, ll, io_chan, i, k, io_mod(max_stations), &
   &  segment, channel, i_ch, channel_max, n_chan, &
   &  iys, ixs, io_up(max_stations), io_ew(max_stations), io_ns(max_stations)

   real*8 :: dip_segment, theta, dis, az, baz, rad_c, coef_v(2, 3), coef_r(2, 5)
   real filter(wave_pts), f1, f2, f3, f4, const_c, tsub(max_psources), h, omega, &
   &  dist_min, dist_max, depth_sub, dep_min, dep_max, lat_p, lon_p, shear_subfault, &
   &  lat_sta, lon_sta, df_bank, tlen_bank, &
   &  time, a, block, lat_s(max_stations), lon_s(max_stations), dt, rake, &
   &  ang_ns(max_stations), ang_ew(max_stations), df, area, dt_sample, tlen
   
   complex :: kahan_y, kahan_t, kahan_c
   complex sour_sub(wave_pts), green_s(wave_pts2, 10), www, wss, z0

   character(len=250) modes
   character(len=100) surf_gf_bank
   character(len=6) sta_name(max_stations)

   write(*,*)'Store long period surface waves GF in memory...'
   z0 = cmplx(0.0, 0.0)
!
!  suppose the ni = u3e+11, then then moment of 1cm*1km^2 
!       subfault is 3e+21. The gf is for Mo = 1e+20 dyne.cm
!  the unit of surface wave is mm.
!
   area = dxs*dys
   block = 1000.0*area*(1.e-10)
  
   open(9, file='channels_surf.txt', status='old')
   open(15, file='wavelets_surf.txt', status='old')
   read(15,*)
   read(15,'(a)')modes
 
   read(9,*)
   read(9,*) 
   read(9,*) dip_segment, rake, theta, lnpt, dt_sample
   read(9,*)
   if (lnpt.ne.10) then
      write(*,*)"please check input LNPT"
   end if
   nlen = 2**lnpt
   dt = dt_channel(ll_in + 1)!dt_sample
   tlen = dt*nlen
   df = 1.0/tlen
  
   read(9,*) channel_max, n_chan
   read(9,*)
   do channel = 1, channel_max    
      read(9,*) no, sta_name(channel), lat_s(channel), lon_s(channel), io_mod(channel), &
      &  io_up(channel), io_ns(channel), io_ew(channel), ang_ns(channel), ang_ew(channel)
   end do
!  int
!  by default frequency window: 0.003 0.004 0.007 0.008
!
   f1 = 0.003
   f2 = 0.004
   f3 = 0.006
   f4 = 0.007
   nf1 = int(f1/df)
   nf2 = int(f2/df)
   nf3 = int(f3/df)
   nf4 = int(f4/df)+1
   do k = 1, nlen/2
      filter(k) = 1.0
      if (k .ge. nf4 .or. k .lt. nf1) filter(k) = 0.0
      if (k .ge. nf1 .and. k .lt. nf2 .and. ((nf2-nf1) .gt. 0)) then
         filter(k) = 0.5*(1.0-cos(pi*(k-nf1)/(nf2-nf1)))
      end if
      if (k .gt. nf3 .and. k .lt. nf4 .and. ((nf4-nf3) .gt. 0)) then
         filter(k) = 0.5*(1.0+cos(pi*(k-nf3)/(nf4-nf3)))
      end if
      filter(k) = filter(k)*block!/nlen !filter(k) = filter(k)*block
   end do

   io_chan = 0
   npp = nx_p*ny_p
!       
!       Here we read the green functions of long period surface waves
!
   call get_surf_gf_data(modes, surf_gf_bank)
   tlen_bank = npt_bank*dt_bank
   df_bank = 1.0/tlen_bank
   const_c = dt_bank/dt
   ll = 0
   dist_min = 20.0
   dist_max = 120.0
   do segment = 1, segments
      dip_segment = dip(segment)
      theta = strike(segment)
      dep_min = point_sources(3, 1, 1, 1, 1, segment)
      dep_max = point_sources(3, 1, ny_p, 1, nys_sub(segment), segment)
      call check_bounds(30.0, 90.0, dep_min, dep_max)
      do iys = 1, nys_sub(segment)
         segment_subfault = (iys-1)*nxs_sub(segment)+1
         depth_sub = point_sources(3, 1, int(ny_p/2) + 1, 1, iys, segment)
         dep_min = point_sources(3, 1, 1, 1, iys, segment) - 1.0
         dep_max = point_sources(3, 1, ny_p, 1, iys, segment) + 1.0
         if (dep_min .lt. 4) dep_min = 4
!  
!  Sanity check:
!
         if (dep_max .lt. dep_min) dep_max = 4
!
         call get_surf_gf(surf_gf_bank, dist_min, dist_max, dep_min, dep_max)
!
!       Note that if the tlen is same in both gf bank and used in here
!       values at each frequency are same. In this code, we will let dt_bank = 2 sec
!       and npt_bank = 2048, in contrast, dt_sample = 4 sec, npt_bank = 1024 
!
         do ixs = 1, nxs_sub(segment)
            ll = ll+1
            segment_subfault = (iys-1)*nxs_sub(segment)+ixs
            lat_p = point_sources(1, int(nx_p/2) + 1, int(ny_p/2) + 1, ixs, iys, segment)
            lon_p = point_sources(2, int(nx_p/2) + 1, int(ny_p/2) + 1, ixs, iys, segment)
            shear_subfault = shear(segment_subfault, segment)
            do i = 1, max_freq
               omega = (i-1)*df_bank*twopi
               sour_sub(i) = cmplx(0.0, 0.0)
               kahan_y = z0
               kahan_t = z0
               kahan_c = z0
               do iy = 1, ny_p
                  do ix = 1, nx_p
                     h = point_sources(4, ix, iy, ixs, iys, segment)
                     time = h/v_ref
                     tsub(1) = time+delay_seg(segment)
                     a = -omega*tsub(1)
                     kahan_y = cmplx(cos(a), sin(a)) - kahan_c
                     kahan_t = sour_sub(i) + kahan_y
                     kahan_c = (kahan_t - sour_sub(i)) - kahan_y
                     sour_sub(i) = kahan_t
!                     sour_sub(i) = sour_sub(i)+cmplx(cos(a), sin(a))
                  end do
               end do
               sour_sub(i) = sour_sub(i)*const_c/npp
            end do
            io_chan = 0
            do channel = 1, channel_max
               lat_sta = lat_s(channel)
               lon_sta = lon_s(channel)
               call distaz(lat_sta, lon_sta, lat_p, lon_p, dis, az, baz)
               dis = dis/111.32
               green_s = interp_gf(dis, depth_sub, dist_min, dist_max, dep_min, dep_max)
               do i_ch = 1, 3
                  if (i_ch .eq. 1.and.io_up(channel).ne.1) then
                     cycle
                  elseif (i_ch .eq. 1 .and. io_up(channel) .eq. 1) then
                     rad_c = 0.0
                  end if
                  if (i_ch .eq. 2.and.io_ns(channel).ne.1) then
                     cycle
                  elseif (i_ch .eq. 2 .and. io_ns(channel) .eq. 1) then 
                     rad_c = ang_ns(channel)
                  end if
                  if (i_ch .eq. 3.and.io_ew(channel).ne.1) then
                     cycle
                  elseif (i_ch .eq. 3 .and. io_ew(channel) .eq. 1) then
                     rad_c = ang_ew(channel)
                  end if

                  io_chan = io_chan+1
                  ll_g = io_chan+ll_in
                  call rad_coef(dip_segment, theta, az, rad_c, coef_v, coef_r)
                  do i = 1, max_freq
                     www = cmplx(0.0, 0.0)
                     wss = cmplx(0.0, 0.0)
                     if (i_ch .eq. 1) then
                        do j = 1, 3
                           www = www+coef_v(1, j)*green_s(i, j+5)
                           wss = wss+coef_v(2, j)*green_s(i, j+5)
                        end do
                     else
                        do j = 1, 5
                           www = www+coef_r(1, j)*green_s(i, j)
                           wss = wss+coef_r(2, j)*green_s(i, j)
                        end do
                     end if
                     if (i .le. max_freq .and. segment .le. 10) then
                        green_stk(i, ll_g, ll) = wss*filter(i)*sour_sub(i)*shear_subfault
                        green_dip(i, ll_g, ll) = www*filter(i)*sour_sub(i)*shear_subfault
                     else
                        green_dip(i, ll_g, ll) = cmplx(0.0, 0.0)
                        green_stk(i, ll_g, ll) = cmplx(0.0, 0.0)
                     end if
                  end do
               end do
            end do
         end do
      end do
   end do

   close(9)
   ll_out = ll_in+n_chan

   return
   end subroutine get_surface_waves_gf


   subroutine deallocate_gf()
   deallocate(green_stk)
   deallocate(green_dip)
   end subroutine deallocate_gf


end module retrieve_gf
