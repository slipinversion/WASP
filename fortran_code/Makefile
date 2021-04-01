all: finite_fault get_gf compute_gf clean


finite_fault: 
	$(MAKE) -C bin_inversion_gfortran_f95


get_gf: 
	$(MAKE) -C bin_str_f95


compute_gf: 
	$(MAKE) -C src_dc_f95


clean:
	\rm -f bin_inversion_gfortran_f95/*.o
	\rm -f bin_inversion_gfortran_f95/*.mod
	\rm -f bin_str_f95/*.o
	\rm -f bin_str_f95/*.mod
	\rm -f src_dc_f95/*.o
	\rm -f src_dc_f95/*.mod
