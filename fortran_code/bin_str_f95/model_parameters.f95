module model_parameters


   use constants, only : max_seg, nnpx, nnpy, nnxs, nnys, nnxy, nt1
   use modelling_inputs, only : t_latest
   implicit none
   real :: slip0(nnxy, max_seg), rake0(nnxy, max_seg), rupt_time0(nnxy, max_seg)
   real :: tl0(nnxy, max_seg), tr0(nnxy, max_seg)
   integer :: n_seg, nxs_sub(max_seg), nys_sub(max_seg), nx_p, ny_p, nxs0, nys0
   real :: dip_seg(max_seg), stk_seg(max_seg), delay_seg(max_seg)
   real :: point_sources(7, nnpx, nnpy, nnxs, nnys, max_seg)
   real :: cniu(nnxy, max_seg)
   real :: c_depth, dxs, dys, v_ref, v_min, v_max, tbl, tbr
   real :: ta0, dta
   real :: time_min(nnxy, max_seg), time_max(nnxy, max_seg), rake_min
   integer :: msou, jfmax
   real :: beg(nt1), dp(nt1)
   integer :: np(nt1)
!
! for regularization
!
   integer :: nleft(3, nnys, nnxs, max_seg), nright(3, nnys, nnxs, max_seg), & 
   & nup(3, nnys, nnxs, max_seg), ndown(3, nnys, nnxs, max_seg)


contains


   subroutine get_faults_data()
   implicit none
   integer n_s, i_s, io_x, io_y, ixs, iys, io_v_d, nxy, k, nxys(max_seg), &
   &  i_seg, ix, iy, ll, io_seg, kxy, kpxy, nx_c, ny_c
   real dist, t_ref, t_max, t_min, delta, dip_s, stk_s
!
!     Input Fault position to memory
!
   open(12, file='Fault.time', status='old')
   read(12,*) nxs0, nys0, c_depth
   read(12,*) n_seg, dxs, dys, nx_p, ny_p, v_min, v_max, tbl, tbr
   read(12,*) ta0, dta, msou, v_ref, io_v_d
   do i_seg = 1, n_seg
      write(*,*) i_seg, stk_seg(1), stk_seg(2)
      read(12,*) io_seg, dip_seg(i_seg), stk_seg(i_seg)
      write(*,*) i_seg, io_seg, stk_seg(1), stk_seg(2)
      read(12,*) nxs_sub(i_seg), nys_sub(i_seg), delay_seg(i_seg)
      nxy = nxs_sub(i_seg) * nys_sub(i_seg)
      do ll = 1, nxy
         read(12,*) slip0(ll, i_seg), rake0(ll, i_seg), rupt_time0(ll, i_seg), tl0(ll, i_seg), tr0(ll, i_seg)
!  magic
         slip0(ll, i_seg) = int(slip0(ll, i_seg))
         rake0(ll, i_seg) = int(rake0(ll, i_seg))
!
      end do
   end do
   close(12)
!
!     Input Fault model to memory
!
   open(12, file='Fault.pos', status='old')
   do i_seg = 1, n_seg
      read(12,*) io_seg, dip_s, stk_s
      if ((abs(dip_s-dip_seg(i_seg)) .gt. 1.e-2).or. &
      &  (abs(stk_s-stk_seg(i_seg)) .gt. 1.e-2)) then
         write(*,*)'the value in Fault.pos is not matched with that in Fault.das'
         write(*,*) i_seg
         write(*,*) dip_s, dip_seg(i_seg)
         write(*,*) stk_s, stk_seg(i_seg)
      end if

      kxy = 0
      do IYS = 1, NYS_sub(i_seg)
         do IXS = 1, NXS_sub(i_seg)
            kxy = kxy+1
            kpxy = 0
            do IY = 1, ny_p
               do IX = 1, nx_p
                  kpxy = kpxy+1
                  read(12,*)(point_sources(k, ix, iy, ixs, iys, i_seg), k = 1, 7)
               end do
            end do
         end do
      end do
   end do
   close(12)
!
!  minimum and maximum subfault rupture arrival time
!
   nx_c = max(int(nx_p/2.0+0.51), 1)
   ny_c = max(int(ny_p/2.0+0.51), 1)
   kpxy = (ny_c-1)*nx_p+nx_c
   do i_seg = 1, n_seg
      do iys = 1, nys_sub(i_seg)
         do ixs = 1, nxs_sub(i_seg)
            kxy = (iys-1)*nxs_sub(i_seg)+ixs
            dist = point_sources(4, nx_c, ny_c, ixs, iys, i_seg)
            t_ref = dist/v_ref
            t_max = dist/v_min
            if (t_max .gt. t_latest) t_max = t_latest
            t_min = dist/v_max
            delta = t_min-t_ref
            if (tbl .lt. delta) then
               time_min(kxy, i_seg) = delta
            else
               time_min(kxy, i_seg) = tbl
            end if
            delta = t_max-t_ref
            if (tbr .gt. delta) then
               time_max(kxy, i_seg) = delta
            else
               time_max(kxy, i_seg) = tbr
            end if
         end do
      end do
   end do
!  
!  Input shear modulous Model
!
   open(12, file='Niu_model', status='old')
   read(12,*) n_s
   if (n_s.ne.n_seg) then
      write(*,*)'Amount of fault segments is different between mu file and Fault file'
      stop
   end if
   do i_s = 1, n_s
      read(12,*) io_seg, io_x, io_y
      nxys(i_s) = io_x*io_y
      read(12,*)(cniu(k, i_s), k = 1, nxys(i_s))
   end do 
   close(12)
   end subroutine get_faults_data


   subroutine write_model(dd, aa, tt, tl, tr)
   real :: dd(nnxy, max_seg), aa(nnxy, max_seg), tt(nnxy, max_seg)
   real :: tl(nnxy, max_seg), tr(nnxy, max_seg)
   real :: latitude_ep, longitude_ep, t_ref, moment_sol
   integer :: i_seg, iys, ixs, iy, ix, i_g, kp
   latitude_ep = 0.0
   longitude_ep = 0.0
   do i_seg = 1, n_seg
      do IYS = 1, NYS_sub(i_seg)
         do IXS = 1, NXS_sub(i_seg)
            do IY = 1, ny_p
               do IX = 1, nx_p
!
! we have found the epicenter
!
                  if (abs(point_sources(4, ix, iy, ixs, iys, i_seg)) .le. 1e-3) then 
                     latitude_ep = point_sources(1, ix, iy, ixs, iys, i_seg)
                     longitude_ep = point_sources(2, ix, iy, ixs, iys, i_seg)
                  end if
               end do
            end do
         end do
      end do
   end do


   open(13, file='Solucion.txt')
   write(13,*)'#Total number of fault_segments=', n_seg
   do i_g = 1, n_seg
      write(13,131)'#Fault_segment =', i_g, ' nx(Along-strike)=', &
     &  nxs_sub(i_g), ' Dx = ', dxs, 'km ny(downdip)= ', nys_sub(i_g), &
     & ' Dy = ', dys, 'km'
131           format(a, i4, a, i3, a, f5.2, a, i3, a, f5.2, a)
      write(13,132)'#Boundary of Fault_segment ', i_g, &
     & '. EQ in cell (', nxs0, ',', nys0, '). Lon:', longitude_ep, &
     & '   Lat:', latitude_ep
132           format(a, i4, a, i2, a, i2, a, f10.4, a, f10.4)
      write(13,*)'#Lon.  Lat.  Depth'
      write(13,*) point_sources(2, 1, 1, 1, 1, i_g), &! ix, iy, ixs, iys, i_g),
     &  point_sources(1, 1, 1, 1, 1, i_g), point_sources(3, 1, 1, 1, 1, i_g)
      write(13,*) point_sources(2, 1, ny_p, 1, nys_sub(i_g), i_g),  &
     &  point_sources(1, 1, ny_p, 1, nys_sub(i_g), i_g),  &
     &  point_sources(3, 1, ny_p, 1, nys_sub(i_g), i_g) 
      write(13,*) &
     &  point_sources(2, nx_p, ny_p, nxs_sub(i_g), nys_sub(i_g), i_g), &
     &  point_sources(1, nx_p, ny_p, nxs_sub(i_g), nys_sub(i_g), i_g), &
     &  point_sources(3, nx_p, ny_p, nxs_sub(i_g), nys_sub(i_g), i_g)
      write(13,*) point_sources(2, nx_p, 1, nxs_sub(i_g), 1, i_g), &
     &  point_sources(1, nx_p, 1, nxs_sub(i_g), 1, i_g), &
     &  point_sources(3, nx_p, 1, nxs_sub(i_g), 1, i_g)
      write(13,*) point_sources(2, 1, 1, 1, 1, i_g),  &
     &  point_sources(1, 1, 1, 1, 1, i_g), point_sources(3, 1, 1, 1, 1, i_g)
      write(13,*)'#Lat. Lon. depth slip rake strike dip t_rup t_ris t_fal mo'
      kp = 0
      ix = int(nx_p / 2) + 1
      iy = int(ny_p / 2) + 1
      do iys = 1, nys_sub(i_g)
         do ixs = 1, nxs_sub(i_g)
            kp = kp + 1
            t_ref = point_sources(5, ix, iy, ixs, iys, i_g)
            moment_sol = dd(kp, i_g) * cniu(kp, i_g) * dxs * dys * (10.0 ** 10.0)
            write(13, 133) point_sources(1, ix, iy, ixs, iys, i_g), &
         &  point_sources(2, ix, iy, ixs, iys, i_g), point_sources(3, ix, iy, ixs, iys, i_g), &
         &  dd(kp, i_g), aa(kp, i_g), stk_seg(i_g), dip_seg(i_g), &
         &  tt(kp, i_g) + t_ref, tl(kp, i_g), tr(kp, i_g), moment_sol
133  format(f14.6, f14.6, f14.6, f14.6, f14.6, f14.6, f14.6, f14.6, f14.6, f14.6, e14.6)
         end do
      end do
   end do
   close(13)
   end subroutine write_model


   subroutine get_model_space()
   implicit none
!
!  boundary conditions? or what?
!
   integer :: nx0, ny0
   parameter(nx0 = nnxs*2, ny0 = nnys*2)
   integer :: io_surf, nblock, npa, k, i_seg, io_right, io_left, io_up, io_down 
   real :: delt_x, delt_y, zmed_max, zmed_min, zleft_max, zleft_min
   real :: zright_max, zright_min, zup_max, zup_min, zdown_max, zdown_min
   real :: angle_max, angle_min, vel_max, vel_min
   real :: ddx1, ddx2, ddy1, ddy2
   integer :: io_seg, nmed, nleft2, nright2, nup2, ndown2, nangle, npv, nb, nsour
   integer :: nx, ny, i, j, k_s, i_ss, i_x, i_y
   real :: surface(1000, 4), xyb(nx0, ny0, 3), xr(5), u0(nt1, 4)
   open(17, file='bound.in', status='old')
   read(17,*) io_surf
   if (io_surf .eq. 1) then
      open(16, file='surface.constrain', status='old')
      read(16,*) nblock, io_surf
      do i = 1, nblock
         read(16,*)(surface(i, k), k = 1, 4)
      end do
      close(16)
   end if
   npa = 0
   k = 0
   do i_seg = 1, n_seg
      read(17,*) io_seg
      read(17,*) delt_x, delt_y
      read(17,*) io_left, io_right, io_up, io_down
      read(17,*) zmed_max, zmed_min, nmed
      read(17,*) zleft_max, zleft_min, nleft2
      read(17,*) zright_max, zright_min, nright2
      read(17,*) zup_max, zup_min, nup2
      read(17,*) zdown_max, zdown_min, ndown2
      read(17,*) angle_max, angle_min, nangle
      read(17,*) vel_max, vel_min, npv
      read(17,*) nb, nsour
      nx = nxs_sub(i_seg)
      ny = nys_sub(i_seg)
      rake_min = angle_min
!
!  check whether the contrains is frighting eath other
!
      if (io_right .eq. 1 .or. io_left .eq. 1) then
         ddx1 = (zmed_max-zleft_max)/(int((nx+1)/2))
         ddx2 = (zmed_max-zright_max)/(int((nx+1)/2))
         if (delt_x .lt. ddx1 .or. delt_x .lt. ddx2) then
            write(*,*)'the constrain of left or right boundary is error'
            stop
         end if
      end if
      if (io_down .eq. 1 .or. io_up .eq. 1) then
         ddy1 = (zmed_max-zup_max)/(int((ny+1)/2))
         ddy2 = (zmed_max-zdown_max)/(int((ny+1)/2))
         if (delt_y .lt. ddy1 .or. delt_y .lt. ddy2) then
            write(*,*)'the constrain of left or right boundary is error'
            stop
         end if
      end if

!  First:  Give the range of mediate part of space

      do i = 2, nxs_sub(i_seg)-1
         do j = 2, nys_sub(i_seg)-1
            xyb(i, j, 1) = zmed_min
            xyb(i, j, 2) = zmed_max
            xyb(i, j, 3) = nmed
         end do
      end do
!
!     Second: Give the value range of left and right part
!
      do j = 1, nys_sub(i_seg)
         if (io_left .eq. 1) then
            xyb(1, j, 1) = zleft_min
            xyb(1, j, 2) = zleft_max
            xyb(1, j, 3) = nleft2
         else
            xyb(1, j, 1) = zmed_min
            xyb(1, j, 2) = zmed_max
            xyb(1, j, 3) = nmed
         end if
         if (io_right .eq. 1) then
            xyb(nxs_sub(i_seg), j, 1) = zright_min
            xyb(nxs_sub(i_seg), j, 2) = zright_max
            xyb(nxs_sub(i_seg), j, 3) = nright2
         else
            xyb(nxs_sub(i_seg), j, 1) = zmed_min
            xyb(nxs_sub(i_seg), j, 2) = zmed_max
            xyb(nxs_sub(i_seg), j, 3) = nmed
         end if
      end do
!
!  finally: Give the value range of up and down part
!
      do i = 2, nxs_sub(i_seg)-1
         if (io_up .eq. 1) then
            xyb(i, 1, 1) = zup_min
            xyb(i, 1, 2) = zup_max
            xyb(i, 1, 3) = nup2
         else
            xyb(i, 1, 1) = zmed_min
            xyb(i, 1, 2) = zmed_max
            xyb(i, 1, 3) = nmed
         end if
         if (io_down .eq. 1) then
            xyb(i, nys_sub(i_seg), 1) = zdown_min
            xyb(i, nys_sub(i_seg), 2) = zdown_max
            xyb(i, nys_sub(i_seg), 3) = ndown2
         else
            xyb(i, nys_sub(i_seg), 1) = zmed_min
            xyb(i, nys_sub(i_seg), 2) = zmed_max
            xyb(i, nys_sub(i_seg), 3) = nmed
         end if
      end do
!
!  Recheck the range of mediate part
!
      do i = 2, nxs_sub(i_seg)-1
         do j = 2, nys_sub(i_seg)-1
            xr(1) = (i-1)*delt_x+zleft_max
            xr(2) = (nxs_sub(i_seg)-i)*delt_x+zright_max
            xr(3) = (j-1)*delt_y+zup_max
            xr(4) = (nys_sub(i_seg)-j)*delt_y+zdown_max
            xr(5) = xyb(i, j, 2)
            call bbsort(xr, 1, 5)
            xyb(i, j, 2) = xr(1)
         end do
      end do
!  
!  check the surface constrain
!  
      if (io_surf .eq. 1) then
         do k_s = 1, nblock
            i_ss = int(surface(k_s, 1)+0.1)
            if (i_ss .eq. i_seg) then
               i_y = 1
               i_x = int(surface(k_s, 2)+0.1)
               xyb(i_x, i_y, 1) = surface(k_s, 3)
               xyb(i_x, i_y, 2) = surface(k_s, 4)
            end if
         end do
      end if
!
!  Change xy_range to xyb    
!
      npa = npa+4*nxs_sub(i_seg)*nys_sub(i_seg)
      do j = 1, nys_sub(i_seg)
         do i = 1, nxs_sub(i_seg)
            k = k+1
            u0(k, 1) = xyb(i, j, 1)
            u0(k, 2) = xyb(i, j, 2)
            u0(k, 3) = xyb(i, j, 3)
            k = k+1
            u0(k, 1) = angle_min
            u0(k, 2) = angle_max
            u0(k, 3) = nangle
            k = k+1
            u0(k, 1) = vel_min
            u0(k, 2) = vel_max
            u0(k, 3) = npv
            k = k+1
            u0(k, 1) = 1
            u0(k, 2) = MSOU
            u0(k, 3) = MSOU
         end do
      end do
   end do
   close(17)

   do k = 1, npa
      beg(k) = u0(k, 1)
      np(k) = int(u0(k, 3) + 0.1)
      if (np(k) .gt. 1) then
         dp(k) = (u0(k, 2)-u0(k, 1))/(np(k)-1)
      else
         dp(k) = 0
      end if
   end do
   end subroutine get_model_space

   
   subroutine bbsort(a, mm, nn)
   implicit none
   real :: a(*), d
   integer :: mm, nn, m, j, i
   m = nn-mm+1
   do while (m .gt. 0)
! 10      if (m .gt. 0) then
      j = m+mm-2
      m = 0
      do i = mm, j
         if (a(i) .gt. a(i+1)) then
            d = a(i)
            a(i) = a(i+1)
            a(i+1) = d
            m = i-mm+1
         end if
      end do
   end do
!        go to 10
!        end if
   end subroutine bbsort


   subroutine get_special_boundaries()
   implicit none
   integer :: nm, nn, iss, ixs, iys, ll, i_g, kxy
   real :: dd(nnxy, max_seg), aa(nnxy, max_seg)
!
!  special boundary
!
   open(12, file='bound.special', status='old')
   read(12,*) nm
   do nn = 1, nm
      read(12,*) iss, ixs, iys
      ll = 0
      do i_g = 1, iss-1
         ll = ll+nxs_sub(i_g)*nys_sub(i_g)
      end do
      kxy = ixs+(iys-1)*nxs_sub(iss)
      ll = ll+ixs+(iys-1)*nxs_sub(iss)
      np(4*(ll-1)+1) = 2
      dp(4*(ll-1)+1) = 10
      beg(4*(11-1)+1) = 1
      dd(kxy, iss) = beg(4*(ll-1)+1)
      np(4*(ll-1)+2) = 2
      aa(kxy, iss) = beg(4*(ll-1)+2)
      np(4*(ll-1)+3) = 2
      np(4*(ll-1)+4) = 2
   end do
   close(12)
   end subroutine get_special_boundaries

   
   subroutine subfault_positions()
   implicit none
   integer k, nn_use
   integer n_is
   integer ix, iy, nxx, nyy, io_change
   character(len=80) aaaa
   integer i_s
!
!  We detect adjacent subfaults and their relative location
!  
   do i_s = 1, n_seg
      do ix = 1, nxs_sub(i_s)
         do iy = 1, nys_sub(i_s)
            nup(1, iy, ix, i_s) = i_s
            nup(2, iy, ix, i_s) = ix
            nup(3, iy, ix, i_s) = iy-1
            ndown(1, iy, ix, i_s) = i_s
            ndown(2, iy, ix, i_s) = ix
            ndown(3, iy, ix, i_s) = iy+1
            nleft(1, iy, ix, i_s) = i_s
            nleft(2, iy, ix, i_s) = ix-1
            nleft(3, iy, ix, i_s) = iy
            nright(1, iy, ix, i_s) = i_s
            nright(2, iy, ix, i_s) = ix+1
            nright(3, iy, ix, i_s) = iy
         end do
      end do
   end do
   open(22, file='continue', status='old')
   do I_s = 1, n_seg 
!       
      read(22,*)
      read(22,*) n_is, nyy
      do ix = 1, nxs_sub(I_s)
         nup(1, 1, ix, I_s) = n_is
         nup(2, 1, ix, I_s) = ix
         nup(3, 1, ix, I_s) = nyy
      end do
!       n
      read(22,*) n_is, nyy
      nn_use = nys_sub(i_s)
      do ix = 1, nxs_sub(I_s)
         ndown(1, nn_use, ix, I_s) = n_is
         ndown(2, nn_use, ix, i_s) = ix
         ndown(3, nn_use, ix, i_s) = nyy
      end do
!       t
      read(22,*) n_is, nxx
      do iy = 1, nys_sub(I_s)
         nleft(1, iy, 1, I_s) = n_is
         nleft(2, iy, 1, i_s) = nxx
         nleft(3, iy, 1, i_s) = iy
      end do
!       ht
      read(22,*) n_is, nxx
      nn_use = nxs_sub(i_s)
      do iy = 1, nys_sub(I_s)
         nright(1, iy, nn_use, i_s) = n_is
         nright(2, iy, nn_use, i_s) = nxx
         nright(3, iy, nn_use, i_s) = iy
      end do
   end do
   close(22)   
   open(22, file='continue.special', status='old')
   read(22,*) nn_use
   do k = 1, nn_use
      read(22,'(a)')aaaa
! 1  format(a)
      write(*,'(a)')aaaa
      read(22,*) i_s, ix, iy
      read(22,*) io_change, n_is, nxx, nyy
      if (io_change .eq. 1) then
         nup(1, iy, ix, i_s) = n_is
         nup(2, iy, ix, i_s) = nxx
         nup(3, iy, ix, i_s) = nyy
      end if
      read(22,*) io_change, n_is, nxx, nyy
      if (io_change .eq. 1) then
         ndown(1, iy, ix, i_s) = n_is
         ndown(2, iy, ix, i_s) = nxx
         ndown(3, iy, ix, i_s) = nyy
      end if
      read(22,*) io_change, n_is, nxx, nyy
      if (io_change .eq. 1) then
         nleft(1, iy, ix, i_s) = n_is
         nleft(2, iy, ix, i_s) = nxx
         nleft(3, iy, ix, i_s) = nyy
      end if
      read(22,*) io_change, n_is, nxx, nyy
      if (io_change .eq. 1) then
         nright(1, iy, ix, i_s) = n_is
         nright(2, iy, ix, i_s) = nxx
         nright(3, iy, ix, i_s) = nyy
      end if
   end do
   close(22)
   end subroutine subfault_positions


end module model_parameters
