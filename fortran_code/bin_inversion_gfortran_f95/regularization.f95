module regularization
   

   use constants, only : max_seg, max_subf, dpi
   use model_parameters, only : nxs_sub, nys_sub, nleft, nright, nup, ndown, rake_min, segments 
   implicit none
   real :: slip_field(max_subf, max_seg, 2)
   integer :: subfaults, subfaults_segment(max_seg)


contains
   
        
   subroutine define_slip_field(slip, rake)
!
!   Laplacian regularization of slip vector field
!
   implicit none
   real, intent(in) :: slip(max_subf, max_seg), rake(max_subf, max_seg)
   real angle
   integer segment, i
!  
!  we define the slip vector field
!  
   subfaults = 0
   do segment = 1, segments
      subfaults_segment(segment) = nxs_sub(segment)*nys_sub(segment)
      do i = 1, subfaults_segment(segment)
         angle = (rake(i, segment)-rake_min)*dpi
         slip_field(i, segment, 1) = slip(i, segment)*cos(angle)
         slip_field(i, segment, 2) = slip(i, segment)*sin(angle)
      end do
      subfaults = subfaults+subfaults_segment(segment)
   end do
   end subroutine define_slip_field

   
   subroutine modify_slip_field(subfault, slip_subf, rake_subf)
!
!   Laplacian regularization of slip vector field
!
   implicit none
   integer, intent(in) :: subfault
   real, intent(in) :: slip_subf, rake_subf
   integer n_total, subfault_seg
   real angle
   integer segment, i
   n_total = 0
   do i = 1, segments
      n_total = subfaults_segment(i)+n_total
      if (subfault .le. n_total) then
         segment = i
         subfault_seg = subfault
         exit
      end if
   end do
   do i = 1, segment-1
      subfault_seg = subfault_seg-subfaults_segment(i)
   end do
   angle = (rake_subf-rake_min)*dpi
   slip_field(subfault_seg, segment, 1) = slip_subf*cos(angle)
   slip_field(subfault_seg, segment, 2) = slip_subf*sin(angle)
   end subroutine modify_slip_field


   subroutine slip_laplace(err)
!
!   Laplacian regularization of slip vector field
!
   implicit none
   real, intent(out) :: err
   integer n_sub
   integer n_is, ll
   integer nxx, nyy, jj, nx, ny
   real d1, d2, d3, d4!, error
   real*8 :: err2, error
   integer segment
   
   err2 = 0.d0
   do jj = 1, 2
      do segment = 1, segments
         do n_sub = 1, subfaults_segment(segment)
            ny = int(n_sub/nxs_sub(segment))+1
            nx = n_sub-(ny-1)*nxs_sub(segment)
            if (nx .eq. 0) then
               nx = nxs_sub(segment)
               ny = ny-1
            end if
!                 write(*,*)"segment", segment, nx, ny, n_sub, jj
!       left
            n_is = nleft(1, ny, nx, segment)
            nxx = nleft(2, ny, nx, segment)
            nyy = nleft(3, ny, nx, segment)
            if (n_is .eq. 0) then
               d1 = 0.0
            else
               ll = nxx+(nyy-1)*nxs_sub(n_is)
               d1 = slip_field(ll, n_is, jj)
            end if
!        write(*,*) n_is, nxx, nyy, ll,"left"
!       right   
            n_is = nright(1, ny, nx, segment)
            nxx = nright(2, ny, nx, segment)
            nyy = nright(3, ny, nx, segment)
            if (n_is .eq. 0) then
               d3 = 0.0
            else
               ll = nxx+(nyy-1)*nxs_sub(n_is)
               d3 = slip_field(ll, n_is, jj)
            end if
!        write(*,*) n_is, nxx, nyy, ll,"right"
!       up    
            n_is = nup(1, ny, nx, segment)
            nxx = nup(2, ny, nx, segment)
            nyy = nup(3, ny, nx, segment)
            if (n_is .eq. 0) then
               d2 = 0.0
            else
               ll = nxx+(nyy-1)*nxs_sub(n_is)
               d2 = slip_field(ll, n_is, jj)
            end if
!                    write(*,*) n_is, nxx, nyy, ll,"up"
!       down
            n_is = ndown(1, ny, nx, segment)
            nxx = ndown(2, ny, nx, segment)
            nyy = ndown(3, ny, nx, segment)
            if (n_is .eq. 0) then
               d4 = 0.0
            else
               ll = nxx+(nyy-1)*nxs_sub(n_is)
               d4 = slip_field(ll, n_is, jj)
            end if
!                   write(*,*) n_is, nxx, nyy, ll,"down"
!
            error = slip_field(n_sub, segment, jj)-0.25*(d1+d2+d3+d4)
            error = error*error
            err2 = err2+error
!            kahan_y = error-kahan_c
!            kahan_t = err+kahan_y
!            kahan_c = (kahan_t-err)-kahan_y
!            err = kahan_t
         end do
      end do
   end do
   err2 = sqrt(err2/subfaults)
   err = real(err2)
   end subroutine slip_laplace


   pure subroutine time_laplace(tt, err)
!
!   Laplacian regularization of rupture initiation time
!
   implicit none
   real, intent(in) :: tt(max_subf, max_seg) 
   real, intent(out) :: err
   integer n_sub
   integer n_is
   integer nxx, nyy, nx, ny
   real d1, d2, d3, d4!, error
   real*8 :: err2, error
   integer segment, ll

   err2 = 0.d0
!   kahan_y = 0.0
!   kahan_c = 0.0
!   kahan_t = 0.0
!   err = 0.0
   do segment = 1, segments
      do n_sub = 1, subfaults_segment(segment)
         ny = int(n_sub/nxs_sub(segment))+1
         nx = n_sub-(ny-1)*nxs_sub(segment)       
         if (nx .eq. 0) then
            nx = nxs_sub(segment)
            ny = ny-1
         end if

!       left
         n_is = nleft(1, ny, nx, segment)
         nxx = nleft(2, ny, nx, segment)
         nyy = nleft(3, ny, nx, segment)
         if (n_is .eq. 0) then
            d1 = tt(n_sub, segment)
         else
            ll = nxx+(nyy-1)*nxs_sub(n_is)
            d1 = tt(ll, n_is)
         end if
!       right   
         n_is = nright(1, ny, nx, segment)
         nxx = nright(2, ny, nx, segment)
         nyy = nright(3, ny, nx, segment)
         if (n_is .eq. 0) then
            d3 = tt(n_sub, segment)
         else
            ll = nxx+(nyy-1)*nxs_sub(n_is)
            d3 = tt(ll, n_is)
         end if
!       up    
         n_is = nup(1, ny, nx, segment)
         nxx = nup(2, ny, nx, segment)
         nyy = nup(3, ny, nx, segment)
         if (n_is .eq. 0) then
            d2 = tt(n_sub, segment)
         else
            ll = nxx+(nyy-1)*nxs_sub(n_is)
            d2 = tt(ll, n_is)
         end if
!       down
         n_is = ndown(1, ny, nx, segment)
         nxx = ndown(2, ny, nx, segment)
         nyy = ndown(3, ny, nx, segment)
         if (n_is .eq. 0) then
            d4 = tt(n_sub, segment)
         else
            ll = nxx+(nyy-1)*nxs_sub(n_is)
            d4 = tt(ll, n_is)
         end if
         error = tt(n_sub, segment)-0.25*(d1+d2+d3+d4)
         error = error*error
         err2 = err2+error
!         kahan_y = error-kahan_c
!         kahan_t = err+kahan_y
!         kahan_c = (kahan_t-err)-kahan_y
!         err = kahan_t
      end do
   end do
   err2 = sqrt(err2/subfaults)
   err = real(err2)
   end subroutine time_laplace


end module regularization
   
