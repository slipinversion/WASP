program run_modelling


   use constants, only : max_seg, max_subfaults2, max_subf
   use model_parameters, only : get_faults_data, get_model_space, get_special_boundaries, subfault_positions, &
                            &   write_model, deallocate_ps
   use modelling_inputs, only : get_annealing_param, n_iter, io_re, cooling_rate, t_stop, t_mid, t0, idum
   use get_stations_data, only : get_data
   use retrieve_gf, only : get_gf, deallocate_gf
   use save_forward, only : write_forward
   use rise_time, only : get_source_fun, deallocate_source
   use static_data, only : initial_gps
   use insar_data, only : initial_insar, get_insar_gf, deallocate_insar_gf, &
                    &   get_insar_data, is_ramp, initial_ramp
   use random_gen, only : start_seed
   use annealing, only : initial_model, print_summary, &
                     &   annealing_iter3, annealing_iter4, n_threads
   implicit none
   integer :: i
   real :: slip(max_subf, max_seg), rake(max_subf, max_seg), rupt_time(max_subf, max_seg)
   real :: t_rise(max_subf, max_seg), t_fall(max_subf, max_seg)
   real :: t
   real*8 :: ramp(18)
   logical :: static, strong, cgps, dart, body, surf, auto
   logical :: get_coeff, insar, ramp_gf_file
   character(len=10) :: input

   write(*,'(/A/)')"CHEN-JI'S WAVELET KINEMATIC MODELLING METHOD"
   static = .False.
   insar = .False.
   get_coeff = .True.
   strong = .False.
   cgps = .False.
   body = .False.
   surf = .False.
   dart = .False.
   auto = .False.
   do i = 1, iargc()
      call getarg(i, input)
      input = trim(input)
      if (input .eq.'gps') static = .True.
      if (input .eq.'insar') insar = .True.
      if (input .eq.'strong') strong = .True.
      if (input .eq.'cgps') cgps = .True.
      if (input .eq.'body') body = .True.
      if (input .eq.'surf') surf = .True.
      if (input .eq.'dart') dart = .True.
      if (input .eq.'auto') auto = .True.
   end do
   call n_threads(auto)
   call get_annealing_param()
   call start_seed(idum)
   call get_faults_data()
   call get_model_space()
   call get_special_boundaries()
   call subfault_positions()
   call get_data(strong, cgps, body, surf, dart)
   call get_source_fun()
   call get_gf(strong, cgps, body, surf, dart)
   call initial_model(slip, rake, rupt_time, t_rise, t_fall)
   t = t_mid
   if (io_re .eq. 0) t = t0
   if (static) call initial_gps(slip, rake)
   if (insar .eqv. .False.) then
      call print_summary(slip, rake, rupt_time, t_rise, t_fall, static, &
           &  insar, get_coeff)
      write(*,*)'Start simmulated annealing...'
      do i = 1, n_iter
         if (static) then
            call annealing_iter4(slip, rake, rupt_time, t_rise, t_fall, &
               &  t, static, insar)
         else
            call annealing_iter3(slip, rake, rupt_time, t_rise, t_fall, t)
         endif   
         write(*,*)'iter: ', i
         if (t .lt. t_stop) then
            t = t*0.995
         else
            t = t*cooling_rate
         end if
      end do
      get_coeff = .False.
      call print_summary(slip, rake, rupt_time, t_rise, t_fall, static, &
              &  insar, get_coeff)
   else
      call get_insar_gf()
      call get_insar_data()
      call is_ramp(ramp_gf_file)
      if (ramp_gf_file .eqv. .False.) then
         call initial_insar(slip, rake)
         call print_summary(slip, rake, rupt_time, t_rise, t_fall, static, &
              &  insar, get_coeff)
         write(*,*)'Start simmulated annealing...'
         do i = 1, n_iter
            call annealing_iter4(slip, rake, rupt_time, t_rise, t_fall, &
               &  t, static, insar)
            write(*,*)'iter: ', i
            if (t .lt. t_stop) then
               t = t*0.995
            else
               t = t*cooling_rate
            end if
         end do
         get_coeff = .False.
         call print_summary(slip, rake, rupt_time, t_rise, t_fall, static, &
              &  insar, get_coeff)
         call initial_insar(slip, rake)
      else
         call initial_ramp(ramp) 
         call initial_insar(slip, rake, ramp)
         call print_summary(slip, rake, rupt_time, t_rise, t_fall, static, &
              &  insar, get_coeff, ramp)
         write(*,*)'Start simmulated annealing...'
         do i = 1, n_iter
            call annealing_iter4(slip, rake, rupt_time, t_rise, t_fall, &
               &  t, static, insar, ramp)
            write(*,*)'iter: ', i
            if (t .lt. t_stop) then
               t = t*0.995
            else
               t = t*cooling_rate
            end if
         end do
         get_coeff = .False.
         call print_summary(slip, rake, rupt_time, t_rise, t_fall, static, &
              &  insar, get_coeff, ramp)
         call initial_insar(slip, rake, ramp)
      endif
   end if
   call write_forward(slip, rake, rupt_time, t_rise, t_fall, strong, cgps, body, surf)
   if (static) call initial_gps(slip, rake)
   call write_model(slip, rake, rupt_time, t_rise, t_fall)
   write(*,'(/A/)')"END CHEN-JI'S WAVELET KINEMATIC MODELLING METHOD"
   call deallocate_source()
   call deallocate_gf()
   call deallocate_ps()
   if (insar) call deallocate_insar_gf()


end program run_modelling

