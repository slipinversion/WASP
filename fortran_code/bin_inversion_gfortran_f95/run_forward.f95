program run_forward


   use constants, only : max_seg, nt1, nnxy
   use model_parameters, only : get_faults_data, slip0, rake0, rupt_time0, tl0, tr0, write_model
   use modelling_inputs, only : get_annealing_param
   use get_stations_data, only : get_data
   use retrieve_gf, only : get_gf
   use save_forward, only : write_forward
   use static_data, only : initial_gps
   implicit none
   integer i
   character(len=10) :: input
   logical :: static, strong, cgps, body, surf, dart

   static = .False.
   strong = .False.
   cgps = .False.
   body = .False.
   surf = .False.
   dart = .False.
   do i = 1, iargc()
      call getarg(i, input)
      input = trim(input)
      if (input .eq.'gps') static = .True.
      if (input .eq.'strong') strong = .True.
      if (input .eq.'cgps') cgps = .True.
      if (input .eq.'body') body = .True.
      if (input .eq.'surf') surf = .True.
      if (input .eq.'dart') dart = .True.
   end do
   call get_annealing_param()
   call get_faults_data()
   call get_gf(strong, cgps, body, surf, dart)
   call get_data(strong, cgps, body, surf, dart)
   call write_forward(slip0, rake0, rupt_time0, tl0, tr0, strong, cgps, body, surf, dart)
   if (static) call initial_gps(slip0, rake0)
   call write_model(slip0, rake0, rupt_time0, tl0, tr0)


end program run_forward
