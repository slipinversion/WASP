module modelling_inputs


   implicit none
   integer :: idum, io_data, n_iter
   real :: t0, cooling_rate, t_stop, emin0, smooth_moment, smooth_slip, smooth_time
   real :: t_mid, t_latest, moment_input
   integer :: io_re, io_func, io, events
   real :: moment_event(10)


contains


   subroutine get_annealing_param()
   implicit none
   write(*,*)'Store annealing and modelling parameters...'
   open(10, file='annealing.txt', status='old')
   read(10,*) n_iter, idum, io_data, moment_input
   read(10,*) t0, cooling_rate, t_stop, emin0, smooth_moment, smooth_slip, smooth_time
   read(10,*) io_re, t_mid, io_func, t_latest
   read(10,*) io
   close(10)
   end subroutine get_annealing_param
   
   
   subroutine moment_events()
   implicit none
   integer :: io, i
   open(10, file='moment_events.txt', status='old')
   events = 0
   do
      read(10,*,iostat=io)
      if (io/=0) exit
      events = events + 1
   enddo
   rewind (10)
   do i=1, events
      read(10,*)moment_event(i)
   enddo
   close(10)
   end subroutine moment_events
   

end module modelling_inputs 
