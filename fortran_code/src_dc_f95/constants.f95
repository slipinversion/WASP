module constants

   implicit none
   integer, parameter :: nnpy = 25, nnpx = 25, nnxs = 50, nnys = 20
   integer, parameter :: mpx = nnpx*nnxs, mpy = nnpy * nnys, nnxy = nnxs * nnys
   integer, parameter :: l2 = 10
   integer, parameter :: npth = 2 ** l2, inptd = 2 * npth
   integer, parameter :: nnsta = 200
   integer, parameter :: max_seg = 5
   integer, parameter :: npuse = 513, block_stg = npuse * 130
   real*8, parameter :: zero=0.d0, one=1.d0, two=2.d0
   real*8, parameter :: pi = 4.d0*atan(one), twopi = two*pi, epsilon=0.0001d0
   real*8, parameter :: dpi = pi/180.d0
!       maximum layer, distances, and time samples
!       maximum size of green's function (nnx*nny)
   integer, parameter :: nlay = 50, ndis = 2000, nt = 2049
   integer, parameter :: nny = 310

end module constants
