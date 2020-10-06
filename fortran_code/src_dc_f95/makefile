##  If you want to change the size of GF bank, please just modify the size.h
FFLAGS=-O
MY_BIN=.
F77=gfortran -ffixed-line-length-132 -mcmodel=large -O2 -Wall -g#-fbacktrace -ffpe-trap=underflow,denormal
# linux or Mac
#F77= ifort -132 -O2 -assume byterecl
#F77= gfortran -fno-rerun-loop-opt -fstrength-reduce -ffixed-line-length-132 -fno-range-check -mcmodel=large
#

modules = constants.o retrieve_gf.o bessel2.o vel_model_data.o wave_travel.o fk_source.o haskell.o prop.o fk_kernel.o fk.o bessel.o $(SUBS)
modules2 = constants.o retrieve_gf.o bessel2.o vel_model_data.o wave_travel.o fk_source.o haskell.o prop.o fk_kernel.o fk_openmp.o bessel.o $(SUBS)
modules3 = constants.o retrieve_gf.o bessel2.o vel_model_data.o wave_travel.o fk_source.o st_haskell.o prop.o fk_st_kernel.o fk_static.o bessel.o $(SUBS)
SUBS = fft.o Complex.o 
FKSUBS = green_bank_dc_3n.o cmodel_fk.o trav_fk.o fk_dc_4.o kernel.o prop.o haskell.o source.o bessel.o $(SUBS) # fk_dc_3.o
FKSUBS2 = green_bank_dc_31n.o cmodel_fk.o trav_fk.o fk_dc_41.o kernel.o prop.o haskell.o source.o bessel.o $(SUBS) # fk_dc_3.o
#FKSUBS_STA = gf_static.o ddis.o cmodel_fk.o trav_fk.o fk_dc_3.o kernel.o prop.o st_haskell.o source.o bessel.o $(SUBS) 
all: green_bank_openmp gf_static clean 


green_bank_fk: $(modules) green_bank_fk.o
	$(F77) -o $(MY_BIN)/green_bank_f95 $(modules) green_bank_fk.o

green_bank_openmp: $(modules2) green_bank_fk_openmp.o
	$(F77) -fopenmp -o $(MY_BIN)/green_bank_openmp_f95 $(modules2) green_bank_fk_openmp.o

gf_static: $(modules3) geodesics.o rad_pattern.o gf_static.o
	$(F77) -o $(MY_BIN)/gf_static_f95 $(modules3) geodesics.o rad_pattern.o gf_static.o


fk_openmp.mod: fk_openmp.o fk_openmp.f95
	$(F77) -fopenmp -c fk_openmp.f95

fk.mod: fk.o fk.f95
	$(F77) -c fk.f95

fk_static.mod: fk_static.o fk_static.f95
	$(F77) -c fk_static.f95

rad_pattern.mod: rad_pattern.o rad_pattern.f95
	$(F77) -c rad_pattern.f95

geodesics.mod: geodesics.o geodesics.f95
	$(F77) -c geodesics.f95

prop.mod: prop.o prop.f95
	$(F77) -c prop.f95

fk_kernel.mod: fk_kernel.o fk_kernel.f95
	$(F77) -c fk_kernel.f95

fk_st_kernel.mod: fk_st_kernel.o fk_st_kernel.f95
	$(F77) -c fk_st_kernel.f95

st_haskell.mod: st_haskell.o st_haskell.f95
	$(F77) -c st_haskell.f95

haskell.mod: haskell.o haskell.f95
	$(F77) -c haskell.f95

fk_source.mod: fk_source.o fk_source.f95
	$(F77) -c fk_source.f95

wave_travel.mod: wave_travel.o wave_travel.f95
	$(F77) -c wave_travel.f95

vel_model_data.mod: vel_model_data.o vel_model_data.f95
	$(F77) -c vel_model_data.f95

bessel2.mod: bessel2.o bessel2.f95
	$(F77) -c bessel2.f95

retrieve_gf.mod: retrieve_gf.o retrieve_gf.f95
	$(F77) -c retrieve:gf.f95

constants.mod: constants.o constants.f95
	$(F77) -c constants.f95


fft.o: fft.c 
	cc -c fft.c

Complex.o: Complex.c Complex.h
	cc -c Complex.c

rad_pattern.o: rad_pattern.f95
	$(F77) -c rad_pattern.f95

geodesics.o: geodesics.f95
	$(F77) -c geodesics.f95

fk_openmp.o: fk_openmp.f95
	$(F77) -fopenmp -c fk_openmp.f95

fk.o: fk.f95
	$(F77) -c fk.f95

fk_static.o: fk_static.f95
	$(F77) -c fk_static.f95

prop.o: prop.f95
	$(F77) -c prop.f95

fk_kernel.o: fk_kernel.f95
	$(F77) -c fk_kernel.f95

fk_st_kernel.o: fk_st_kernel.f95
	$(F77) -c fk_st_kernel.f95

st_haskell.o: st_haskell.f95
	$(F77) -c st_haskell.f95

haskell.o: haskell.f95
	$(F77) -c haskell.f95

fk_source.o: fk_source.f95
	$(F77) -c fk_source.f95

wave_travel.o: wave_travel.f95
	$(F77) -c wave_travel.f95

vel_model_data.o: vel_model_data.f95
	$(F77) -c vel_model_data.f95

bessel2.o: bessel2.f95
	$(F77) -c bessel2.f95

retrieve_gf.o: retrieve_gf.f95
	$(F77) -c retrieve_gf.f95

constants.o: constants.f95
	$(F77) -c constants.f95

bessel.o: bessel.f
	$(F77) -c bessel.f

bessel.f: bessel.FF
	cpp -traditional-cpp $< > $@

green_bank_fk.o: green_bank_fk.f95
	$(F77) -c green_bank_fk.f95

green_bank_fk_openmp.o: green_bank_fk_openmp.f95
	$(F77) -c green_bank_fk_openmp.f95

gf_static.o: gf_static.f95
	$(F77) -c gf_static.f95

clean:
	rm -f *.mod
	rm -f *.o
