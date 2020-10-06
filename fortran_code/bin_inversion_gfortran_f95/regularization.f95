module regularization
   

   use constants, only : max_seg, nnxy, dpi
   use model_parameters, only : nxs_sub, nys_sub, nleft, nright, nup, ndown, rake_min, n_seg 
   implicit none
   real :: slip_field(nnxy, max_seg, 2)
   integer :: nnn, nxys(max_seg)


contains
   
        
   subroutine define_slip_field(slip, rake)
!
!   Laplacian regularization of slip vector field
!
   implicit none
   real, intent(in) :: slip(nnxy, max_seg), rake(nnxy, max_seg)
   real angle
   integer i_s, i
!  
!  we define the slip vector field
!  
   nnn = 0
   do i_s = 1, n_seg
      nxys(i_s) = nxs_sub(i_s)*nys_sub(i_s)
      do i = 1, nxys(i_s)
         angle = (rake(i, i_s)-rake_min)*dpi
         slip_field(i, i_s, 1) = slip(i, i_s)*cos(angle)
         slip_field(i, i_s, 2) = slip(i, i_s)*sin(angle)
      end do
      nnn = nnn+nxys(i_s)
   end do
   end subroutine define_slip_field

   
   subroutine modify_slip_field(nn_sub, d_sub, a_sub)
!
!   Laplacian regularization of slip vector field
!
   implicit none
   integer, intent(in) :: nn_sub
   real, intent(in) :: d_sub, a_sub
   integer n_total, ll_s
   real angle
   integer i_s, i_ss
   n_total = 0
   do i_ss = 1, n_seg
      n_total = nxys(i_ss)+n_total
      if (nn_sub .le. n_total) then
         i_s = i_ss
         ll_s = nn_sub
         exit
      end if
   end do
   do i_ss = 1, i_s-1
      ll_s = ll_s-nxys(i_ss)
   end do
   angle = (a_sub-rake_min)*dpi
   slip_field(ll_s, i_s, 1) = d_sub*cos(angle)
   slip_field(ll_s, i_s, 2) = d_sub*sin(angle)
   end subroutine modify_slip_field


   subroutine lap(err)
!
!   Laplacian regularization of slip vector field
!
   implicit none
   real, intent(out) :: err
   integer n_sub
   integer n_is, ll
   integer nxx, nyy, jj, nx, ny
   real d1, d2, d3, d4, error!, kahan_y, kahan_t, kahan_c
   real(8) :: err2
   integer i_s
   
   err2 = 0.d0
!   kahan_y = 0.0
!   kahan_t = 0.0
!   kahan_c = 0.0
   do jj = 1, 2
      do i_s = 1, n_seg
         do n_sub = 1, nxys(i_s)
            ny = int(n_sub/nxs_sub(i_s))+1
            nx = n_sub-(ny-1)*nxs_sub(i_s)
            if (nx .eq. 0) then
               nx = nxs_sub(i_s)
               ny = ny-1
            end if
!                 write(*,*)"i_s", i_s, nx, ny, n_sub, jj
!       left
            n_is = nleft(1, ny, nx, i_s)
            nxx = nleft(2, ny, nx, i_s)
            nyy = nleft(3, ny, nx, i_s)
            if (n_is .eq. 0) then
               d1 = 0.0
            else
               ll = nxx+(nyy-1)*nxs_sub(n_is)
               d1 = slip_field(ll, n_is, jj)
            end if
!        write(*,*) n_is, nxx, nyy, ll,"left"
!       right   
            n_is = nright(1, ny, nx, i_s)
            nxx = nright(2, ny, nx, i_s)
            nyy = nright(3, ny, nx, i_s)
            if (n_is .eq. 0) then
               d3 = 0.0
            else
               ll = nxx+(nyy-1)*nxs_sub(n_is)
               d3 = slip_field(ll, n_is, jj)
            end if
!        write(*,*) n_is, nxx, nyy, ll,"right"
!       up    
            n_is = nup(1, ny, nx, i_s)
            nxx = nup(2, ny, nx, i_s)
            nyy = nup(3, ny, nx, i_s)
            if (n_is .eq. 0) then
               d2 = 0.0
            else
               ll = nxx+(nyy-1)*nxs_sub(n_is)
               d2 = slip_field(ll, n_is, jj)
            end if
!                    write(*,*) n_is, nxx, nyy, ll,"up"
!       down
            n_is = ndown(1, ny, nx, i_s)
            nxx = ndown(2, ny, nx, i_s)
            nyy = ndown(3, ny, nx, i_s)
            if (n_is .eq. 0) then
               d4 = 0.0
            else
               ll = nxx+(nyy-1)*nxs_sub(n_is)
               d4 = slip_field(ll, n_is, jj)
            end if
!                   write(*,*) n_is, nxx, nyy, ll,"down"
!
            error = slip_field(n_sub, i_s, jj)-0.25*(d1+d2+d3+d4)
            error = error*error
            err2 = err2+error
!            kahan_y = error-kahan_c
!            kahan_t = err+kahan_y
!            kahan_c = (kahan_t-err)-kahan_y
!            err = kahan_t
         end do
      end do
   end do
   err2 = sqrt(err2/nnn)
   err = real(err2)
   end subroutine lap


   pure subroutine tlap(tt, err)
!
!   Laplacian regularization of rupture initiation time
!
   implicit none
   real, intent(in) :: tt(nnxy, max_seg) 
   real, intent(out) :: err
   integer n_sub
   integer n_is
   integer nxx, nyy, nx, ny
   real d1, d2, d3, d4, error!, kahan_y, kahan_t, kahan_c
   real(8) :: err2
   integer i_s, ll

   err2 = 0.d0
!   kahan_y = 0.0
!   kahan_c = 0.0
!   kahan_t = 0.0
!   err = 0.0
   do i_s = 1, n_seg
      do n_sub = 1, nxys(i_s)
         ny = int(n_sub/nxs_sub(i_s))+1
         nx = n_sub-(ny-1)*nxs_sub(i_s)       
         if (nx .eq. 0) then
            nx = nxs_sub(i_s)
            ny = ny-1
         end if

!       left
         n_is = nleft(1, ny, nx, i_s)
         nxx = nleft(2, ny, nx, i_s)
         nyy = nleft(3, ny, nx, i_s)
         if (n_is .eq. 0) then
            d1 = tt(n_sub, i_s)
         else
            ll = nxx+(nyy-1)*nxs_sub(n_is)
            d1 = tt(ll, n_is)
         end if
!       right   
         n_is = nright(1, ny, nx, i_s)
         nxx = nright(2, ny, nx, i_s)
         nyy = nright(3, ny, nx, i_s)
         if (n_is .eq. 0) then
            d3 = tt(n_sub, i_s)
         else
            ll = nxx+(nyy-1)*nxs_sub(n_is)
            d3 = tt(ll, n_is)
         end if
!       up    
         n_is = nup(1, ny, nx, i_s)
         nxx = nup(2, ny, nx, i_s)
         nyy = nup(3, ny, nx, i_s)
         if (n_is .eq. 0) then
            d2 = tt(n_sub, i_s)
         else
            ll = nxx+(nyy-1)*nxs_sub(n_is)
            d2 = tt(ll, n_is)
         end if
!       down
         n_is = ndown(1, ny, nx, i_s)
         nxx = ndown(2, ny, nx, i_s)
         nyy = ndown(3, ny, nx, i_s)
         if (n_is .eq. 0) then
            d4 = tt(n_sub, i_s)
         else
            ll = nxx+(nyy-1)*nxs_sub(n_is)
            d4 = tt(ll, n_is)
         end if
         error = tt(n_sub, i_s)-0.25*(d1+d2+d3+d4)
         error = error*error
         err2 = err2+error
!         kahan_y = error-kahan_c
!         kahan_t = err+kahan_y
!         kahan_c = (kahan_t-err)-kahan_y
!         err = kahan_t
      end do
   end do
   err2 = sqrt(err2/nnn)
   err = real(err2)
   end subroutine tlap


end module regularization
   
