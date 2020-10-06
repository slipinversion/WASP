!c
!c	For strong motion , we need use geometry coordanation directly.
!c
program get_strong_motion


   use constants, only : dpi, twopi, nnpy, nnpx, nnys, nnxs, mpx, mpy, inptd, max_seg, block_stg, nt
   use vel_model_data, only : new_vel_s, new_dens, s_layer, read_vel_model, update_model
   use wave_travel, only : trav
   use retrieve_gf, only : dt_gfs, lnpt_gfs, t_cor, get_gf_data, get_gf, interp_gf, deallocate_gf
   use geodesics, only : distaz 
   use fast_fourier_trans, only : fft, fourier_coefs
   use bpfilter, only : bandpassfilter
   use rad_pattern, only : rad_coef
   implicit none
   character(len=100) :: vel_model, gf_file, gf_bank, wave_file, stat_file
   character(len=70) :: string1, input
   character(len=6) :: sta_name(200)
   character(len=12) :: filename
   character(len=3) :: comp
   real*8 :: dip, theta, az, baz, ang_d, coef_v(2,3), coef_r(2,5)
   real :: dep_p, dep_min, dep_max
   real :: lat_e, lon_e, dep_e
   real :: fau_mod(7, nnpx, nnpy, nnxs, nnys, max_seg)
   real :: lat_sta, lon_sta, lat_p, lon_p
   real*8 :: azi(mpx, mpy), bazi(mpx, mpy), distances(mpx, mpy), dis
   real :: time(mpx,mpy), dist(mpx), tmin(mpx)
   real :: dt, kahan_y1(inptd), kahan_y2(inptd), kahan_t1(inptd), kahan_t2(inptd), kahan_c1(inptd), kahan_c2(inptd)
   real*8 :: df, w
   real :: t1, real_dip(inptd), imag_dip(inptd), niu
   real :: real_stk(inptd), imag_stk(inptd)
   real :: niu_fault(nnpy, nnys, max_seg)
   real :: depth_sources(nnpy, nnys, max_seg)
   real :: low_freq, high_freq
   integer :: i, j, nsta
   integer :: ir_max, n_chan, ir, nx, ky, kx, LL, no
   integer :: iys, ixs, iyp, ixp, nxp, nyp
   integer :: npt, lnpt, k, io_mod(200), nn_comp
   integer :: jfmax
   integer :: n_pos, n_delt, nleft
   real :: green_s(nt, 8), v_ref, green
   real :: ww, ws, green_dip(inptd), green_stk(inptd)
   complex :: ssww(inptd), ssws(inptd)
   complex*16 :: z, z1
   integer :: n_seg, i_seg
   integer :: nxs_sub(max_seg), nys_sub(max_seg)
   real :: float1, float2
   integer :: integer1, integer2

   real :: dip_sub(max_seg), stk_sub(max_seg)
   real :: c_depth
   logical :: disp

   call getarg(1, input)
   disp = (input.eq.'cgps')
!
!	Input fault data
!
   call fourier_coefs()
   open(22, file='Fault.time', status='old')
   read(22, *)integer1, integer2, c_depth
   read(22, *)n_seg, float1, float2, nxp, nyp
   read(22, *)float1, float2, integer1, v_ref
   do i_seg = 1, n_seg
      read(22, *)integer1, dip_sub(i_seg),stk_sub(i_seg)
      read(22, *)nxs_sub(i_seg),nys_sub(i_seg)
      do iys = 1, nys_sub(i_seg)
         do ixs = 1, nxs_sub(i_seg)
            read(22,*)
         enddo
      enddo
   enddo
   close(22)
   open(22, file='Fault.pos', status='old')
   do i_seg=1,n_seg
      read(22,*)
      do iys = 1, nys_sub(i_seg)
         do ixs = 1, nxs_sub(i_seg)
            do iyp = 1, nyp
               do ixp = 1, nxp
                  read(22,*)(fau_mod(k, ixp, iyp, ixs, iys, i_seg),k=1,7)
               enddo
            enddo
         enddo
      enddo
   enddo
   close(22)

   open(1, file='filtro_strong', status='old')
   read(1,*)string1, low_freq, high_freq
   close(1)
   if (disp) low_freq = 0.0
!
!	Load vel. model into memory
! 
   gf_file = 'Green.in'
   if (disp) gf_file = 'Green_cgps.in'
   call get_gf_data(gf_file, vel_model, gf_bank)
   call read_vel_model(vel_model)
   dep_min = 1.e+10
   dep_max = 0.0
   do i_seg=1,n_seg
      do iys = 1, nys_sub(i_seg)
         do iyp = 1, nyp
            dep_p = fau_mod(3, 1, iyp, 1, iys, i_seg)
            if(dep_p.lt.dep_min)dep_min = dep_p
            if(dep_p.gt.dep_max)dep_max = dep_p
            call update_model(dep_p)
            niu = new_vel_s(s_layer)*new_vel_s(s_layer)*new_dens(s_layer) * (1.e+10)
            niu_fault(iyp, iys, i_seg) = niu
            depth_sources(iyp, iys, i_seg) = dep_p
         enddo
      enddo
   enddo
   call get_gf(gf_bank, dep_min, dep_max)
!
   wave_file = 'Wave.str'
   if (disp) wave_file = 'Wave.cgps'
   open(13, file=wave_file, status='old')
   read(13,*)integer1, integer2, jfmax
   read(13,*)
   read(13,*)nsta
   if(nsta.lt.1)then
      write(*,*)"In wavelet file, you have not specified any stations"
      stop
   endif
   close(13)
!
!	Load hypocenter location
! 
   stat_file = 'Readlp.inf'
   if (disp) stat_file = 'Readlp.cgps'
   open(12, file=stat_file, status='old')
   read(12,*)
   read(12,*)lat_e, lon_e, dep_e
   read(12,*)lnpt,dt
   read(12,*)
   if(abs(c_depth-dep_e).gt.1e-4)then
      write(*,*)"c_depth,dep_e"
      write(*,*)c_depth,dep_e
      write(*,*)'Mismatch between hypocenter depth in Readlp.inf and in &
     & Fault.time'
      stop
   endif
   npt = 2 ** lnpt
   nleft = npt - 80
   if(npt.gt.nt)then
      write(*,*)"The length of green's function must be shorten"
      stop
   endif
   df = 1.d0 / npt / dt
   if(abs(lnpt_gfs-lnpt).gt.1e-6.or.abs(dt-dt_gfs).gt.1.0e-6)then
      write(*,*)"Please keep the consistency between Readlp.inf and Gfs bank"
      write(*,*)"In Readlp.inf dt= ",dt, " Lnpt= ",lnpt
      write(*,*)"In Gfs        dt= ",dt_gfs," Lnpt= ",lnpt_gfs
      stop
   endif

   read(12,*)Ir_max,n_chan
   read(12,*)
   write(*,*)'Begin to calculate the response',ir_max,n_chan,ir
!
!	Load GF for each station and channel
!
   do ir = 1, ir_max
      read(12,*)no, sta_name(ir), lat_sta, lon_sta, io_mod(ir), comp
      call distaz(lat_sta,lon_sta,lat_e,lon_e,dis,az,baz)
      write(*,*)'az=',az,dis,baz,lat_sta,lon_sta,lat_e,lon_e
      if(comp.eq.'HNZ')then
         nn_comp = 14
         filename=trim(sta_name(ir))//'1'
      endif
      if(comp.eq.'LXZ')then
         nn_comp = 14
         filename=trim(sta_name(ir))//'.cgps.1'
      endif
      if(comp.eq.'HNN')then
         nn_comp = 15
         filename=trim(sta_name(ir))//'2'
      endif
      if(comp.eq.'LXN')then
         nn_comp = 15
         filename=trim(sta_name(ir))//'.cgps.2'
      endif
      if(comp.eq.'HNE')then
         nn_comp = 16
         filename=trim(sta_name(ir))//'3'
      endif
      if(comp.eq.'LXE')then
         nn_comp = 16
         filename=trim(sta_name(ir))//'.cgps.3'
      endif
      open(20,file=filename,status='unknown',access='direct',recl=block_stg)
      ll = 0
      ky = 0
      do i_seg = 1,n_seg
         dip = dip_sub(i_seg)
         theta = stk_sub(i_seg)
         do iys = 1, nys_sub(i_seg)
!
!	 Store azimuth, back-azimuth, distance, and arrivals, of point sources
!
            do iyp = 1, nyp
               ky = (iys - 1) * nyp + iyp
               kx = 0
               do ixs = 1, nxs_sub(i_seg)
                  do ixp = 1, nxp
                     kx = kx + 1
                     lat_p = fau_mod(1, ixp, iyp, ixs, iys, i_seg)
                     lon_p = fau_mod(2, ixp, iyp, ixs, iys, i_seg)
                     call distaz(lat_sta, lon_sta, lat_p, lon_p, dis, az, baz)
                     if (dis.lt.0.25) dis = 0.25d0
                     dist(kx) = dis
                     distances(kx, ky) = dis
                     azi(kx,ky) = az
                     bazi(kx,ky) = baz
                  enddo
               enddo
               nx = nxs_sub(i_seg) * nxp
               dep_p = depth_sources(iyp, iys, i_seg)
               call update_model(dep_p)
               call trav(dist, nx, tmin)
               time(:nx,ky) = tmin(:nx)
            enddo
!
!	Calculate subfault response
!	First save the file: vertical:  STA.1 - STA.2 - STA.3 
!
            do ixs = 1, nxs_sub(i_seg)
               LL = LL + 1
               green_dip(:) = 0.0
               green_stk(:) = 0.0
               kahan_y1(:) = 0.0
               kahan_y2(:) = 0.0
               kahan_t1(:) = 0.0
               kahan_t2(:) = 0.0
               kahan_c1(:) = 0.0
               kahan_c2(:) = 0.0
               do iyp = 1, nyp
                  niu = niu_fault(iyp, iys, i_seg)
                  dep_p = depth_sources(iyp, iys, i_seg)
                  ky = (iys - 1) * nyp + iyp
                  do ixp = 1, nxp
                     kx = (ixs - 1) * nxp + ixp
                     dis = distances(kx, ky)
                     az = azi(kx,ky)
                     baz = bazi(kx,ky)   
                     t1 = fau_mod(4, ixp, iyp, ixs, iys, i_seg) / v_ref
                     t1 = t1 + time(kx, ky)
                     n_delt = int(t1 / dt +0.5)
                     green_s = interp_gf(dis, dep_p, dep_min, dep_max) 
                     ang_d = (baz - 180.0d0)*dpi
                     call rad_coef(dip,theta,az,ang_d,nn_comp,coef_v,coef_r)
!
! n_delt: position of array for which the wave arrives to the station from the given subfault.
! n_pos = i - n_delt, means we initialize the green function when p-wave from the subfault arrives to the station
!
                     if(nn_comp.eq.14)then
                        do i = n_delt + 1, npt
                           ww = 0.0
                           ws = 0.0
                           n_pos = i - n_delt
                           do j = 1, 3
                              green = green_s(n_pos, j + 5) * niu
                              ww = ww + coef_v(1, j) * green
                              ws = ws + coef_v(2, j) * green
                           enddo
                           kahan_y1(i) = ww - kahan_c1(i)
                           kahan_t1(i) = green_dip(i) + kahan_y1(i)
                           kahan_c1(i) = (kahan_t1(i) - green_dip(i)) - kahan_y1(i)
                           green_dip(i) = kahan_t1(i)
                           kahan_y2(i) = ws - kahan_c2(i)
                           kahan_t2(i) = green_stk(i) + kahan_y2(i)
                           kahan_c2(i) = (kahan_t2(i) - green_stk(i)) - kahan_y2(i)
                           green_stk(i) = kahan_t2(i)
!                           green_dip(i) = green_dip(i) + ww
!                           green_stk(i) = green_stk(i) + ws
                        enddo
                     else
                        do i = n_delt + 1, npt
                           ww = 0.0
                           ws = 0.0
                           n_pos = i - n_delt
                           do j = 1,5
                              green = green_s(n_pos, j) * niu
                              ww = ww + coef_r(1, j) * green
                              ws = ws + coef_r(2, j) * green
                           enddo
                           kahan_y1(i) = ww - kahan_c1(i)
                           kahan_t1(i) = green_dip(i) + kahan_y1(i)
                           kahan_c1(i) = (kahan_t1(i) - green_dip(i)) - kahan_y1(i)
                           green_dip(i) = kahan_t1(i)
                           kahan_y2(i) = ws - kahan_c2(i)
                           kahan_t2(i) = green_stk(i) + kahan_y2(i)
                           kahan_c2(i) = (kahan_t2(i) - green_stk(i)) - kahan_y2(i)
                           green_stk(i) = kahan_t2(i)
!                           green_dip(i) = green_dip(i) + ww
!                           green_stk(i) = green_stk(i) + ws
                        enddo
                     endif
                  enddo
               enddo
!
! synthetic in dip and strike component
!
               imag_dip(:) = 0.0
               real_dip(:) = 0.0

               real_dip(:nleft) = green_dip(:nleft) / (nxp * nyp)
               call bandpassfilter(real_dip, inptd, dt, 4, 2, low_freq, high_freq)
               call fft(real_dip, imag_dip, lnpt, -1)
               imag_stk(:) = 0.0
               real_stk(:) = 0.0
               real_stk(:nleft) = green_stk(:nleft) / (nxp * nyp)
               call bandpassfilter(real_stk, inptd, dt, 4, 2, low_freq, high_freq)
               call fft(real_stk, imag_stk, lnpt, -1)
               w = t_cor*twopi*df
               z1 = cmplx(cos(w), sin(w), kind(1d0))
               z = cmplx(1.d0, 0.d0, kind(1d0))
               do i = 1, jfmax
!                  w = t_cor * twopi * (i - 1) * df  ! t_cor seconds backward shift
!                  z = cmplx(0.0, w)
!                  z = cexp(z)
                  ssww(i) = cmplx(real_dip(i), imag_dip(i)) * z
                  ssws(i) = cmplx(real_stk(i), imag_stk(i)) * z
                  z = z*z1
               enddo
               write(20, rec=ll)(ssww(i),i=1,jfmax),(ssws(i),i=1,jfmax)
            enddo
         enddo
      enddo
      close(20)
   enddo
   close(12)
   call deallocate_gf()

   open(1,file='green_functions_details')
   write(1,*)'dt: ', dt
   write(1,*)'Frequency values: ', jfmax
   write(1,*)'Frequency step: '
   close(1)


end program get_strong_motion
