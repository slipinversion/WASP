module ffm_methods


   use constants, only : max_seg, max_subfaults2, max_subf
   use model_parameters, only : write_model, deallocate_ps
   use modelling_inputs, only : n_iter, io_re, cooling_rate, t_stop, t_mid, t0
   use get_stations_data, only : get_data
   use retrieve_gf, only : get_gf, deallocate_gf
   use save_forward, only : write_forward
   use rise_time, only : get_source_fun, deallocate_source
   use static_data, only : initial_gps
   use insar_data, only : initial_insar, get_insar_gf, deallocate_insar_gf, &
                    &   get_insar_data, is_ramp, initial_ramp
   use annealing, only : initial_model, print_summary, &
                     &   annealing_iter3, annealing_iter4, n_threads
   use annealing_static, only : print_static_summary, annealing_iter
   implicit none
   real :: rupt_time0(max_subf, max_seg)
   real :: t_rise0(max_subf, max_seg), t_fall0(max_subf, max_seg)
   real :: t
   real*8 :: ramp(18)
   integer :: i
   character(len=10) :: input


contains


   subroutine check_waveforms(strong, cgps, body, surf, dart, use_waveforms)
   implicit none
   logical :: strong, cgps, dart, body, surf, use_waveforms
   use_waveforms = .False.
   if (strong) use_waveforms = .True.
   if (cgps) use_waveforms = .True.
   if (body) use_waveforms = .True.
   if (surf) use_waveforms = .True.
   if (dart) use_waveforms = .True.
   end subroutine check_waveforms


   subroutine check_static(static, insar, use_static)
   implicit none
   logical :: static, insar, use_static
   use_static = .False.
   if (static) use_static = .True.
   if (insar) use_static = .True.
   end subroutine check_static


   subroutine waveform_ffm(strong, cgps, body, surf, dart, &
       & slip, rake, rupt_time, t_rise, t_fall)
   implicit none
   real :: slip(:, :), rake(:, :), rupt_time(:, :)
   real :: t_rise(:, :), t_fall(:, :)
   logical :: strong, cgps, dart, body, surf
   logical :: get_coeff, static, insar
   get_coeff = .True.
   static = .False.
   insar = .False.
   call get_data(strong, cgps, body, surf, dart)
   call get_source_fun()
   call get_gf(strong, cgps, body, surf, dart)
   call initial_model(slip, rake, rupt_time, t_rise, t_fall)
   t = t_mid
   if (io_re .eq. 0) t = t0
   call print_summary(slip, rake, rupt_time, t_rise, t_fall, static, &
        &  insar, get_coeff)
   write(*,*)'Start simmulated annealing...'
   do i = 1, n_iter
      call annealing_iter3(slip, rake, rupt_time, t_rise, t_fall, t)
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
   call write_forward(slip, rake, rupt_time, t_rise, t_fall, strong, cgps, body, surf)
   call write_model(slip, rake, rupt_time, t_rise, t_fall)
   call deallocate_source()
   call deallocate_gf()
   end subroutine waveform_ffm
   

   subroutine mixed_ffm(strong, cgps, body, surf, dart, &
       & static, insar, slip, rake, rupt_time, t_rise, t_fall)
   implicit none
   real :: slip(:, :), rake(:, :), rupt_time(:, :)
   real :: t_rise(:, :), t_fall(:, :)
   logical :: static, strong, cgps, dart, body, surf
   logical :: get_coeff, insar, ramp_gf_file
   get_coeff = .True.
   call get_data(strong, cgps, body, surf, dart)
   call get_source_fun()
   call get_gf(strong, cgps, body, surf, dart)
   call initial_model(slip, rake, rupt_time, t_rise, t_fall)
   t = t_mid
   if (io_re .eq. 0) t = t0
   if (static) call initial_gps(slip, rake)
   if (insar) call get_insar_gf()
   if (insar) call get_insar_data()
   if (insar) call is_ramp(ramp_gf_file)
   if (ramp_gf_file .eqv. .False.) then
      if (insar) call initial_insar(slip, rake)
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
      if (insar) call initial_insar(slip, rake)
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
   end if
   call write_forward(slip, rake, rupt_time, t_rise, t_fall, strong, cgps, body, surf)
   if (static) call initial_gps(slip, rake)
   call write_model(slip, rake, rupt_time, t_rise, t_fall)
   call deallocate_source()
   call deallocate_gf()
   if (insar) call deallocate_insar_gf()
   end subroutine mixed_ffm


   subroutine static_ffm(slip, rake, static, insar)
   implicit none
   real :: slip(:, :), rake(:, :)
   logical :: static, get_coeff, insar, ramp_gf_file
   get_coeff = .True.
   call initial_model(slip, rake, rupt_time0, t_rise0, t_fall0)
   t = t_mid
   if (io_re .eq. 0) t = t0
   if (static) call initial_gps(slip, rake)
   if (insar) call get_insar_gf()
   if (insar) call get_insar_data()
   if (insar) call is_ramp(ramp_gf_file)
   if (ramp_gf_file .eqv. .False.) then
      if (insar) call initial_insar(slip, rake)
      call print_static_summary(slip, rake, static, insar, get_coeff)
      write(*,*)'Start simmulated annealing...'
      do i = 1, n_iter
         call annealing_iter(slip, rake, t, static, insar)
         write(*,*)'iter: ', i
         if (t .lt. t_stop) then
            t = t*0.995
         else
            t = t*cooling_rate
         end if
      end do
      get_coeff = .False.
      call print_static_summary(slip, rake, static, insar, get_coeff)
      if (insar) call initial_insar(slip, rake)
   else
      call initial_ramp(ramp) 
      call initial_insar(slip, rake, ramp)
      call print_static_summary(slip, rake, static, insar, get_coeff, ramp)
      write(*,*)'Start simmulated annealing...'
      do i = 1, n_iter
         call annealing_iter(slip, rake, t, static, insar, ramp)
         write(*,*)'iter: ', i
         if (t .lt. t_stop) then
            t = t*0.995
         else
            t = t*cooling_rate
         end if
      end do
      get_coeff = .False.
      call print_static_summary(slip, rake, static, insar, get_coeff, ramp)
      call initial_insar(slip, rake, ramp)
   end if
   if (static) call initial_gps(slip, rake)
   call write_model(slip, rake, rupt_time0, t_rise0, t_fall0)
   if (insar) call deallocate_insar_gf()
   end subroutine static_ffm


end module ffm_methods

