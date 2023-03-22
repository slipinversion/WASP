module constants

   implicit none
   integer, parameter :: max_dip_psources = 25, max_stk_psources = 25
   integer, parameter :: max_stk_subfaults = 50, max_dip_subfaults = 20
   integer, parameter :: max_rise_time_range = 25
!   integer, parameter :: mpx = stk_psources * stk_subfaults, mpy = dip_psources * dip_subfaults
   integer, parameter :: max_subf = max_stk_subfaults * max_dip_subfaults
   integer, parameter :: max_psources = max_dip_psources * max_stk_psources, log2_pts = 10
   integer, parameter :: wave_pts = 2 ** log2_pts, wave_pts2 = 2 * wave_pts
   integer, parameter :: n_data = 10000, max_stations = 500
   integer, parameter :: max_seg = 30, max_subfaults = 1000, max_subfaults2 = 5000
   integer, parameter :: npuse = 1 + wave_pts/2, block_stg = npuse * 130
   integer, parameter :: block_far = 256 * 1024, ltde = 320000
   real*8, parameter :: pi = 4.d0*atan(1.d0), twopi = 2.d0 * pi
   real*8, parameter :: dpi = pi / 180.d0
   real, parameter :: pi_0 = 4.0*atan(1.0), twopi_0 = 2.0 * pi_0, dpi_0 = pi_0 / 180.0
!       maximum layer, distances, and time samples
!       maximum size of green's function (nnx*nny)
   integer, parameter :: nlay = 50, ndis = 1100, nt = 1025
   integer, parameter :: nny = 310

end module constants
