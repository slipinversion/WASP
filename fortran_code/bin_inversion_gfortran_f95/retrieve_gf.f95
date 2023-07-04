!
! TODO: change to velocity if io_v_d = 1
!
module retrieve_gf


   use constants, only : npuse, max_stations, max_subfaults, max_subf, nt, nny, ndis, &
   &  wave_pts2, max_seg, block_far, ltde, n_data, max_psources, max_dip_psources, &
   &  max_dip_subfaults, block_stg, pi, twopi, wave_pts
   use wavelet_param, only : get_data_param
   use model_parameters, only : point_sources
   use retrieve_surf_gf, only : get_surf_gf_data, interp_gf, get_surf_gf, npt_bank, &
   &  dt_bank, check_bounds
   use rad_pattern, only : rad_coef
   use geodesics, only : distaz
   use get_stations_data, only : get_properties, disp_or_vel, &
          & idata, mmm, llove, get_event_sta
   implicit none
   integer, parameter :: nnsta_tele = 150
   complex, allocatable :: green_dip(:, :, :), green_stk(:, :, :)
!   complex, allocatable :: green_dip2(:, :, :), green_stk2(:, :, :)
   integer :: segments, nxs_sub(max_seg), nys_sub(max_seg), nx_p, ny_p
   integer :: cum_subfaults(max_seg)
   real :: dip(max_seg), strike(max_seg), delay_seg(max_seg)
   real :: shear(max_subfaults)
   real :: dxs, dys, v_ref, dt_channel(max_stations)
   integer :: lnpt, nlen, jmin, jmax, max_freq, channels, event_sta(max_stations)
   logical :: segment_in_event(max_seg, 10)
   logical :: subfault_in_event(max_subfaults, 10)
   character(len=15) :: sta_name(max_stations)
   character(len=3) :: component(max_stations)
   

contains


   subroutine retrievegf_set_data_properties()
   implicit none
   call get_properties(sta_name, component, dt_channel, channels)
   call get_data_param(lnpt, jmin, jmax, nlen, max_freq)
   call get_event_sta(event_sta)
   end subroutine retrievegf_set_data_properties


   subroutine retrievegf_set_fault_parameters()
   use model_parameters, only : get_shear, get_segments, &
   &  get_subfaults, get_events_segments
   implicit none
   integer :: subfaults
   real :: v_min, v_max
   call get_shear(shear)
   call get_segments(nxs_sub, nys_sub, dip, strike, delay_seg, segments, subfaults, cum_subfaults)
   call get_subfaults(dxs, dys, nx_p, ny_p, v_min, v_max, v_ref)
   call get_events_segments(segment_in_event, subfault_in_event)
   end subroutine retrievegf_set_fault_parameters


   subroutine get_gf(strong, cgps, body, surf, dart, many_events)
   implicit none
!
   integer ll_in, ll_out
   logical :: strong, cgps, body, surf, dart, many_events
   allocate(green_dip(npuse, max_stations, max_subfaults))
   allocate(green_stk(npuse, max_stations, max_subfaults))

   write(*,*)'Store GF in memory...'
!
!  Here, we load into memory the green functions for each subfault, for every used station
!  
   ll_in = 0
   ll_out = 0
   if (strong) then
      call get_near_field_gf(ll_in, ll_out, many_events, strong, .False.)
      ll_in = ll_out
   end if
   if (cgps) then
      call get_near_field_gf(ll_in, ll_out, many_events, .False., cgps)
      ll_in = ll_out
   end if
   if (body) then
      call get_body_waves_gf(ll_in, ll_out, many_events)
      ll_in = ll_out
   end if
   if (surf) then
      call get_surface_waves_gf(ll_in, ll_out, many_events)
      ll_in = ll_out
   end if
   if (dart) then
      call get_dart_gf(ll_in, ll_out)
      ll_in = ll_out
   end if
   end subroutine get_gf


   subroutine get_near_field_gf(ll_in, ll_out, many_events, strong, cgps)
   implicit none
   integer :: ll_in, ll_out, io_v_d, ll_g, ll, &
   &  io_chan, i, segment, channel, channel_max, n_chan, &
   &  ixs, iys, event
   real :: omega, block, dt, df, dt_sample, w, tlen
   complex :: z0, z
   character(len=80) filename, filename2
   character(len=3) comp!component(max_stations), comp
   character(len=1) channel2!component(max_stations), comp
   character(len=2) event2!component(max_stations), comp
   logical :: many_events, strong, cgps
   
   if (strong) write(*,*)'Store strong motion GF in memory...'
   if (cgps) write(*,*)'Store cGPS GF in memory...'
   z0 = cmplx(0.0, 0.0)
!
!  suppose the ni = u3e+11, then then moment of 1cm*1km^2 
!       subfault is 3e+21. The gfs is for Mo = 1e+20
!
   block = dxs * dys * (1.e-10) 

   filename = 'channels_strong.txt'
   if (cgps) filename = 'channels_cgps.txt'
   filename = trim(filename)
   open(9, file=filename, status='old')

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

      io_chan = io_chan+1
      ll_g = io_chan+ll_in
      comp = component(ll_g)
      channel2 = comp(3:3)
      filename2 = trim(sta_name(ll_g))//'.'//comp
      filename2 = trim(filename2)

      open(12, file=filename2, status='old', access='direct',recl=block_stg)
      ll = 0
!
!       Here, we read the green functions and derivate them
!
      if (many_events) event = event_sta(ll_g)
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
               if ((many_events) .and. (segment_in_event(segment, event) .eqv. .False.)) then
                  green_dip(:max_freq, ll_g, ll) = 0.0
                  green_stk(:max_freq, ll_g, ll) = 0.0
               endif
            end do
         end do
      enddo
      close(12)
   end do      
   ll_out = ll_in+n_chan
   end subroutine get_near_field_gf


   subroutine get_body_waves_gf(ll_in, ll_out, many_events)
   implicit none
   integer nstaon, channel, ll_g, k, nsta, n_chan, subfault, psource, &
   &  love, ll, segment, iys, iyp, io_seg, iys_c, iy_c, jf, i, ipy, npxy, &
   &  nkxy, nxs_c, nys_c, nxp_c, nyp_c, ll_s, kxy, ixs, kpxy, ixp, & 
   &  ll_in, ll_out, event
   real dt, df, ddelt, time, block, w, tlen!, &
   real, allocatable :: tdel(:,:,:,:)
   complex :: z, z0, wsyn(wave_pts)
   complex :: kahan_y1(wave_pts2), kahan_t1(wave_pts2), kahan_c1(wave_pts2)
   complex :: kahan_y2(wave_pts2), kahan_t2(wave_pts2), kahan_c2(wave_pts2)
   complex, allocatable :: green_dip0(:,:,:,:,:)
   complex, allocatable :: green_stk0(:,:,:,:,:)
!
   character(len=14) fname4, fname6
   logical :: many_events
   allocate(green_dip0(npuse, max_dip_psources, max_dip_subfaults, max_seg, nnsta_tele))
   allocate(green_stk0(npuse, max_dip_psources, max_dip_subfaults, max_seg, nnsta_tele))
   allocate(tdel(max_psources, max_subf, max_seg, nnsta_tele))
   
   write(*,*)'Store body waves GF in memory...'
   z0 = cmplx(0.0, 0.0)

   open(9, file='channels_body.txt', status='old')
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
         fname4 = trim(sta_name(ll_g))//'.GRE'
         fname6 = trim(sta_name(ll_g))//'.TDE'
      else
         fname4 = trim(sta_name(ll_g))//'SH.GRE'
         fname6 = trim(sta_name(ll_g))//'SH.TDE'
      end if
      open(12, file=fname4, status='old', access='direct', recl=block_far)
      open(32, file=fname6, status='old', access='direct', recl=ltde)
      LL = 0
      if (many_events) event = event_sta(ll_g)
      do segment = 1, segments
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
     &   .or.(nxp_c.ne.nx_p) .or. (nyp_c.ne.ny_p)) then
            write(*,*)'nxs',nxs_c,nxs_sub(segment)
            write(*,*)'nys',nys_c,nys_sub(segment)
            write(*,*)'nxp',nxp_c,nx_p
            write(*,*)'nyp',nyp_c,ny_p
            write(*,'(a)')'Mismatch in amount of point sources or subfaults &
            &between the specified in Fault.time, and those used in the &
            &green functions.'
            stop
         end if
         if (many_events .and. (segment_in_event(segment, event) .eqv. .False.)) then
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
               do psource = 1, npxy
                  tdel(psource, ll_s, segment, channel) = 0.0
               end do
            end do
         endif
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
      subfault = 0
      do segment = 1, segments
         kxy = 0
         do iys = 1, nys_sub(segment)
            do ixs = 1, nxs_sub(segment)
               kxy = kxy+1
               subfault = subfault+1
               do i = 1, 2*max_freq
                  green_dip(i, ll_g, subfault) = z0
                  green_stk(i, ll_g, subfault) = z0
               end do
               kahan_y1(:) = z0
               kahan_t1(:) = z0
               kahan_c1(:) = z0
               kahan_y2(:) = z0
               kahan_t2(:) = z0
               kahan_c2(:) = z0
               psource = 0
               do iyp = 1, ny_p
                  do ixp = 1, nx_p
                     psource = psource+1
                     time = point_sources(4, psource, subfault)/v_ref+delay_seg(segment)
                     time = time + tdel(psource, kxy, segment, channel)
                     do i = 1, 2*max_freq
                        W = -twopi*(i-1)*df*time
                        z = cmplx(0.0, w)
                        z = cexp(z)
                        kahan_y1(i) = green_dip0(i, iyp, iys, segment, channel)*z - kahan_c1(i)
                        kahan_t1(i) = green_dip(i, ll_g, subfault) + kahan_y1(i)
                        kahan_c1(i) = (kahan_t1(i) - green_dip(i, ll_g, subfault)) - kahan_y1(i)
                        green_dip(i, ll_g, subfault) = kahan_t1(i)
                        kahan_y2(i) = green_stk0(i, iyp, iys, segment, channel)*z - kahan_c2(i)
                        kahan_t2(i) = green_stk(i, ll_g, subfault) + kahan_y2(i)
                        kahan_c2(i) = (kahan_t2(i) - green_stk(i, ll_g, subfault)) - kahan_y2(i)
                        green_stk(i, ll_g, subfault) = kahan_t2(i)
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
                  green_dip(i, ll_g, subfault) = z*green_dip(i, ll_g, subfault)/npxy
                  green_stk(i, ll_g, subfault) = z*green_stk(i, ll_g, subfault)/npxy
!                  green_dip2(i, ll, ll_g) = green_dip(i, ll_g, ll)
!                  green_stk2(i, ll, ll_g) = green_stk(i, ll_g, ll)
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


   subroutine get_surface_waves_gf(ll_in, ll_out, many_events)
   implicit none
   integer ll_in, ll_out, nf1, nf2, nf3, nf4, no, subfault, psource, &
   &  npp, ix, iy, segment_subfault, j, ll_g, ll, io_chan, i, k, io_mod(max_stations), &
   &  segment, channel, i_ch, channel_max, n_chan, event, &
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
   character(len=6) sta_name1
   logical :: many_events

   write(*,*)'Store long period surface waves GF in memory...'
   z0 = cmplx(0.0, 0.0)
!
!  suppose the ni = u3e+11, then then moment of 1cm*1km^2 
!       subfault is 3e+21. The gf is for Mo = 1e+20 dyne.cm
!  the unit of surface wave is mm.
!
   area = dxs*dys
   block = 1000.0*area*(1.e-10)

   open(9, file='surf_filter.txt', status='old')
   read(9,*)f1, f2, f3, f4
   close(9)
  
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
      read(9,*) no, sta_name1, lat_s(channel), lon_s(channel), io_mod(channel), &
      &  io_up(channel), io_ns(channel), io_ew(channel), ang_ns(channel), ang_ew(channel)
   end do
!  int
!  by default frequency window: 0.003 0.004 0.007 0.008
!
!   f1 = 0.003
!   f2 = 0.004
!   f3 = 0.006
!   f4 = 0.007
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
   subfault = 0
   dist_min = 20.0
   dist_max = 120.0
   do segment = 1, segments
      dip_segment = dip(segment)
      theta = strike(segment)
      subfault = cum_subfaults(segment) + 1
      dep_min = point_sources(3, 1, subfault)
      subfault = cum_subfaults(segment)
      subfault = subfault + (nys_sub(segment)-1)*nxs_sub(segment) + 1
      psource = (ny_p-1)*nx_p + 1
      dep_max = point_sources(3, psource, subfault)
      call check_bounds(30.0, 90.0, dep_min, dep_max)
      do iys = 1, nys_sub(segment)
         segment_subfault = (iys-1)*nxs_sub(segment)+1
         psource = int(ny_p/2)*nx_p + 1
         subfault = cum_subfaults(segment)
         subfault = subfault + (iys-1)*nxs_sub(segment) + 1
         depth_sub = point_sources(3, psource, subfault)
         dep_min = point_sources(3, 1, subfault) - 1.0
         psource = (ny_p-1)*nx_p + 1
         dep_max = point_sources(3, psource, subfault) + 1.0
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
         subfault = subfault - 1
         do ixs = 1, nxs_sub(segment)
            subfault = subfault + 1
            psource = int(ny_p/2)*nx_p + int(nx_p/2) + 1
            lat_p = point_sources(1, psource, subfault)
            lon_p = point_sources(2, psource, subfault)
            shear_subfault = shear(subfault)
            do i = 1, max_freq
               omega = (i-1)*df_bank*twopi
               sour_sub(i) = cmplx(0.0, 0.0)
               kahan_y = z0
               kahan_t = z0
               kahan_c = z0
               do psource = 1, ny_p*nx_p
                  h = point_sources(4, psource, subfault)
                  time = h/v_ref
                  tsub(1) = time+delay_seg(segment)
                  a = -omega*tsub(1)
                  kahan_y = cmplx(cos(a), sin(a)) - kahan_c
                  kahan_t = sour_sub(i) + kahan_y
                  kahan_c = (kahan_t - sour_sub(i)) - kahan_y
                  sour_sub(i) = kahan_t
!                  sour_sub(i) = sour_sub(i)+cmplx(cos(a), sin(a))
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
                  if (many_events) event = event_sta(ll_g)
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
                     if (many_events) then
                        if (i .le. max_freq .and. segment_in_event(segment, event)) then
                           green_stk(i, ll_g, subfault) = wss*filter(i)*sour_sub(i)*shear_subfault
                           green_dip(i, ll_g, subfault) = www*filter(i)*sour_sub(i)*shear_subfault
                        else
                           green_dip(i, ll_g, subfault) = cmplx(0.0, 0.0)
                           green_stk(i, ll_g, subfault) = cmplx(0.0, 0.0)
                        end if
                     else
                        if (i .le. max_freq .and. segment .le. 10) then
                           green_stk(i, ll_g, subfault) = wss*filter(i)*sour_sub(i)*shear_subfault
                           green_dip(i, ll_g, subfault) = www*filter(i)*sour_sub(i)*shear_subfault
                        else
                           green_dip(i, ll_g, subfault) = cmplx(0.0, 0.0)
                           green_stk(i, ll_g, subfault) = cmplx(0.0, 0.0)
                        end if
                     endif
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


   subroutine get_dart_gf(ll_in, ll_out)
   implicit none
   integer :: ll_in, ll_out, io_v_d, ll_g, subfault, psource, &
   &  io_chan, i, segment, channel, channel_max, n_chan, &
   &  ixs, iys, length, etc
   real :: omega, dt, df, dt_sample, w, tlen, real, imag, time
   complex :: z0, z
   character(len=80) filename
   
   write(*,*)'Store DART GF in memory...'
   z0 = cmplx(0.0, 0.0)
!
!  suppose the ni = u3e+11, then then moment of 1cm*1km^2 
!       subfault is 3e+21. The gfs is for Mo = 1e+20
!

   open(9, file='channels_dart.txt', status='old')

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
   io_chan = 0
!       
!       Here we read the green functions of strong motion waves
!
   df = 1. / ((2.0 ** lnpt) * dt)
   nlen = 2 ** lnpt
   tlen = nlen * dt
   psource = int(ny_p/2 + 1)*nx_p + int(nx_p/2 + 1)

   do channel = 1, channel_max

      ll_g = ll_in+channel
      filename = trim(sta_name(ll_g))//'_gf.txt'
      open(12, file=filename, status='old')
      io_chan = io_chan+1
!
!       Here, we read the green functions and derivate them
!
      subfault = 0
      do segment = 1, segments
         do iys = 1, nys_sub(segment)
            do ixs = 1, nxs_sub(segment)
               subfault = subfault+1
               time = point_sources(4, psource, subfault)/v_ref
               green_dip(:,ll_g,subfault) = z0
               green_stk(:,ll_g,subfault) = z0
               read(12, *)etc, max_freq
               do i = 1, max_freq
                  read(12, *)real, imag
                  green_dip(i, ll_g, subfault) = cmplx(real, imag)
               enddo
               do i = 1, max_freq
!
! we eventually shift synthetics in time, in case the fault plane used has a delay
!
                  omega = twopi*(i-1)*df
                  w = -omega*(time+delay_seg(segment))
                  z = cmplx(0.0, w)
                  z = cexp(z)
                  green_dip(i, ll_g, subfault) = green_dip(i, ll_g, subfault)*z
               end do
            end do
         end do
      end do
      close(12)
   end do      
   ll_out = ll_in+n_chan
   end subroutine get_dart_gf

   
   subroutine deallocate_gf()
   deallocate(green_stk)
   deallocate(green_dip)
   end subroutine deallocate_gf


end module retrieve_gf
