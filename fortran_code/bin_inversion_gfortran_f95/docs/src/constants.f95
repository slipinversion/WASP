module constants

   implicit none
   integer, parameter :: nnpy = 25, nnpx = 25, nnxs = 50, nnys = 20, mmsou = 25
   integer, parameter :: mpx = nnpx * nnxs, mpy = nnpy * nnys, nnxy = nnxs * nnys
   integer, parameter :: nnpxy = nnpx * nnpy, l2 = 10
   integer, parameter :: npth = 2 ** l2, inptd = 2 * npth, n_data = 10000
   integer, parameter :: nnsta = 300
   integer, parameter :: max_seg = 5, nnxy_m = 660, nt1 = max_seg * nnxy_m
   integer, parameter :: npuse = 513, block_stg = npuse * 130
   integer, parameter :: block_far = 256 * 1024, ltde = 320000
   real*8, parameter :: pi = 4.d0*atan(1.d0), twopi = 2.d0 * pi
   real*8, parameter :: dpi = pi / 180.d0
   real, parameter :: pi_0 = 4.0*atan(1.0), twopi_0 = 2.0 * pi_0, dpi_0 = pi_0 / 180.0
!       maximum layer, distances, and time samples
!       maximum size of green's function (nnx*nny)
   integer, parameter :: nlay = 50, ndis = 1100, nt = 1025
   integer, parameter :: nny = 310

end module constants
