Wavelet and simulated Annealing SliP inversion (WASP)

This code uses a wavelet based method to ivnert for slip on a fault using both regional and tele-seismic data. The code is absed on the original tele-seismic approach of Ji et al. (2002) and modified to include regional data (GNSS and strong motion) and for real-time operations by the Centro Sismologico Nacional of Chile. Details of the implementation can be found in Koch et al. (2019) and Goldberg et al. (2022).

Users of this code should consider citing the following relevant publications:
- Ji, C., D. J. Wald, and D. V. Helmberger (2002). Source description of the 1999 Hector Mine, California, earthquake, Part I: Wavelet domain inversion theory and resolution analysis, Bulletin of the Seismological Society of America, 92, no. 4, 1192–1207, https://doi.org/10.1785/0120000916.
- Koch, P., F. Bravo, S. Riquelme, and J. G. F. Crempien (2019). Near-real-time finite-fault inversions for large earthquakes in Chile using strong-motion data, Seismological Research Letters, 90, no. 5, 1971–1986, https://doi.org/10.1785/0220180294.
-Goldberg, D. E., P. Koch, D. Melgar, S. Riquelme, and W. L. Yeck (2022). Beyond the Teleseism: Introducing Regional Seismic and Geodetic Data into Routine USGS FiniteFault Modeling, Seismological Research Letters, 93, 3308–3323, https://doi.org/10.1785/0220220047.
-Zhu, L., & Rivera, L. A. (2002). A note on the dynamic and static displacements from a point source in multilayered media: A note on the dynamic and static displacements from a point source. Geophysical Journal International, 148(3), 619–627. https://doi.org/10.1046/j.1365-246X.2002.01610.x.

The surface wave GF bank (fd_bank) is now available in Zenodo, where the user needs to search for "Surface wave Green's functions for WASP slip inversion code" or go directly to: https://zenodo.org/record/7236739#.Y9q4BOzMKDV
Put the fd_bank file in fortran_code/gfs_nm/long/.

The LITHO1.0.nc can be downloaded here: https://ds.iris.edu/ds/products/emc-litho10/
Put the LITHO1.0.nc file in fortran_code/info/.

Example json files for the modelling of multi-segment faults are available now in multi_segment_example. With these files, the user can re-create the geometry for the solution of the 2011 Tohoku event published by the USGS.
