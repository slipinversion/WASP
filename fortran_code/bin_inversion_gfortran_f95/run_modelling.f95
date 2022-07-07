program run_modelling


   use constants, only : max_seg, max_subfaults, max_subf
   use model_parameters, only : get_faults_data, get_model_space, get_special_boundaries, &
         & subfault_positions, deallocate_ps
   use modelling_inputs, only : get_annealing_param, n_iter, io_re, cooling_rate, t_stop, t_mid, t0, idum
   use save_forward, only : write_forward
   use random_gen, only : start_seed
   use annealing, only : n_threads
   use ffm_methods, only: check_waveforms, check_static, waveform_ffm, &
                     &  mixed_ffm, static_ffm 
   implicit none
   integer :: i
   real :: slip(max_subfaults), rake(max_subfaults), rupt_time(max_subfaults)
   real :: t_rise(max_subfaults), t_fall(max_subfaults)
   real :: t
   real*8 :: ramp(18)
   logical :: static, strong, cgps, dart, body, surf, auto
   logical :: insar, ramp_gf_file
   logical :: use_waveforms, use_static
   character(len=10) :: input

   write(*,'(/A/)')"CHEN-JI'S WAVELET KINEMATIC MODELLING METHOD"
   static = .False.
   insar = .False.
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
   call check_waveforms(strong, cgps, body, surf, dart, use_waveforms)
   call check_static(static, insar, use_static)
   call n_threads(auto)
   call get_annealing_param()
   call start_seed(idum)
   call get_faults_data()
   call get_model_space()
   call get_special_boundaries()
   call subfault_positions()
   if ((use_waveforms) .and. (use_static .eqv. .False.)) then
      call waveform_ffm(strong, cgps, body, surf, dart, &
       & slip, rake, rupt_time, t_rise, t_fall)
   elseif ((use_static) .and. (use_waveforms .eqv. .False.)) then
      call static_ffm(slip, rake, static, insar)
   elseif ((use_static) .and. (use_waveforms)) then
      call mixed_ffm(strong, cgps, body, surf, dart, &
       & static, insar, slip, rake, rupt_time, t_rise, t_fall)
   endif
   call deallocate_ps()
   write(*,'(/A/)')"END CHEN-JI'S WAVELET KINEMATIC MODELLING METHOD"


end program run_modelling

