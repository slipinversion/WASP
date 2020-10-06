module modelling_inputs


   implicit none
   integer :: idum, io_data, n_iter
   real :: t0, cooling_rate, t_stop, emin, smooth_moment, smooth_slip, smooth_time
   real :: t_mid, t_latest, cm_point
   integer :: io_re, io_func, io


contains


   subroutine get_annealing_param()
   implicit none
   open(10,file='HEAT.IN')
   read(10,*) n_iter, idum, io_data, cm_point
   read(10,*) t0, cooling_rate, t_stop, emin, smooth_moment, smooth_slip, smooth_time
   read(10,*) io_re, t_mid, io_func, t_latest
   read(10,*) io
   close(10)
   end subroutine get_annealing_param
   
   
end module modelling_inputs 
