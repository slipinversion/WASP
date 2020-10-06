program run_modelling


   use constants, only : max_seg, nt1, nnxy
   use model_parameters, only : get_faults_data, get_model_space, get_special_boundaries, subfault_positions, &
                            &   write_model
   use modelling_inputs, only : get_annealing_param, n_iter, io_re, cooling_rate, t_stop, t_mid, t0
   use get_stations_data, only : get_data
   use retrieve_gf, only : get_gf
   use save_forward, only : write_forward
   use rise_time, only : get_source_fun
   use annealing, only : initial_model, initial_regularization, annealing_iter
   implicit none
   integer :: io_data2, i
   real :: dd(nnxy, max_seg), aa(nnxy, max_seg), tt(nnxy, max_seg), tl(nnxy, max_seg)
   real :: tr(nnxy, max_seg)
   real :: er, t

   io_data2 = 3
   call get_annealing_param()
   call get_faults_data()
   call get_model_space()
   call get_special_boundaries()
   call subfault_positions()
   call get_gf(io_data2)
   call get_data(io_data2)
   call get_source_fun()
   call initial_model(dd, aa, tt, tl, tr)
   call initial_regularization(dd, aa, tt, tl, tr)
   t = t_mid
   if (io_re.eq.0) t=t0
   do i=1, n_iter
      call annealing_iter(dd, aa, tt, tl, tr, er, t)
      write(*,*)'iter: ', i
      if (t.lt.t_stop) then
         t=t*0.995
      else
         t=t*cooling_rate
      endif
   enddo
   call write_forward(dd, aa, tt, tl, tr, io_data2)
   call write_model(dd, aa, tt, tl, tr)


end program run_modelling

