MY_BIN=.

modules = constants.o modelling_inputs.o model_parameters.o wavelet_param.o wavelets.o get_stations_data.o retrieve_surf_gf.o rad_pattern.o geodesics.o retrieve_gf.o rise_time.o save_forward.o
modules2 = $(modules) static_data.o misfit_eval.o regularization.o random_gen.o annealing.o

#Linux
#F95= ifort -shared-intel -mcmodel=large -O2 -assume byterecl -p#-p -132 -shared-intel
#F77= ifort -O3 -132  -assume byterecl
F95= gfortran -O3 -Wall -fno-range-check -mcmodel=large -march=x86-64 -g# -ftree-loop-vectorize -ftree-loop-distribution#-fno-omit-frame-pointer -ffpe-trap=overflow,invalid#-ffpe-trap=overflow -funroll-loops -ffast-math -fcheck=all
F77= gfortran -O -fno-range-check -mcmodel=large -g -fbacktrace


all: run_forward run_modelling green_tele clean


#test: test.o
#	gfortran -o test test.o

run_forward: $(modules) static_data.o run_forward.o
	$(F95) -o run_forward $(modules) static_data.o run_forward.o

run_modelling: $(modules2) run_modelling.o
	$(F95) -fopenmp -o run_modelling $(modules2) run_modelling.o #-fopenmp

green_tele: bpfilter.o ddis.o green_tele.o
	$(F77) -o green_tele bpfilter.o ddis.o green_tele.o


annealing.mod: annealing.o annealing.f95
	$(F95) -fopenmp -c annealing.f95 #-fopenmp

%.mod: %.o %.f95
	$(F95) -c $<

#regularization.mod: regularization.o regularization.f95
#	$(F95) -c regularization.f95
#
#misfit_eval.mod: misfit_eval.o misfit_eval.f95
#	$(F95) -c misfit_eval.f95
#
#static_data.mod: static_data.o static_data.f95
#	$(F95) -c static_data.f95
#
#random_gen.mod: random_gen.o random_gen.f95
#	$(F95) -c random_gen.f95
#
#rise_time.mod: rise_time.o rise_time.f95
#	$(F95) -c rise_time.f95
#
#save_forward.mod: save_forward.o save_forward.f95
#	$(F95) -c save_forward.f95
#
#retrieve_gf.mod: retrieve_gf.o retrieve_gf.f95
#	$(F95) -c retrieve_gf.f95
#
#retrieve_surf_gf.mod: retrieve_surf_gf.o retrieve_surf_gf.f95
#	$(F95) -c retrieve_surf_gf.f95
#
#wavelet_param.mod: wavelet_param.o wavelets_param.f95
#	$(F95) -c wavelet_param.f95
#
#wavelets.mod: wavelets.o wavelets.f95
#	$(F95) -c wavelets.f95
#
#geodesics.mod: geodesics.o geodesics.f95
#	$(F95) -c geodesics.f95
#
#rad_pattern.mod: rad_pattern.o rad_pattern.f95
#	$(F95) -c rad_pattern.f95
#
#get_stations_data.mod: get_stations_data.o get_stations_data.f95
#	$(F95) -c get_stations_data.f95
#
#model_parameters.mod: model_parameters.o model_parameters.f95
#	$(F95) -c model_parameters.f95
#
#modelling_inputs.mod: modelling_inputs.o modelling_inputs.f95
#	$(F95) -c modelling_inputs.f95
#
#constants.mod: constants.o constants.f95
#	$(F95) -c constants.f95


bpfilter.o : bpfilter.f
	$(F77) -c bpfilter.f

green_tele.o : green_tele.f
	$(F77) -c green_tele.f

ddis.o : ddis.f
	$(F77) -c ddis.f

#test.o: test.f95
#	gfortran -c test.f95

annealing.o: annealing.f95
	$(F95) -fopenmp -c annealing.f95 #-fopenmp

%.o: %.f95
	$(F95) -c $<

#run_forward.o: run_forward.f95
#	$(F95) -c run_forward.f95
#
#run_modelling.o: run_modelling.f95
#	$(F95) -c run_modelling.f95
#
#annealing.o: annealing.f95
#	$(F95) -fopenmp -c annealing.f95 #-fopenmp
#
#regularization.o: regularization.f95
#	$(F95) -c regularization.f95
#
#misfit_eval.o: misfit_eval.f95
#	$(F95) -c misfit_eval.f95
#
#static_data.o: static_data.f95
#	$(F95) -c static_data.f95
#
#random_gen.o: random_gen.f95
#	$(F95) -c random_gen.f95 #-fopenmp
#
#rise_time.o: rise_time.f95
#	$(F95) -c rise_time.f95
#
#save_forward.o: save_forward.f95
#	$(F95) -c save_forward.f95
#
#retrieve_gf.o: retrieve_gf.f95
#	$(F95) -c retrieve_gf.f95
#
#retrieve_surf_gf.o: retrieve_surf_gf.f95
#	$(F95) -c retrieve_surf_gf.f95
#
#wavelet_param.o: wavelet_param.f95
#	$(F95) -c wavelet_param.f95
#
#wavelets.o: wavelets.f95
#	$(F95) -c wavelets.f95
#
#geodesics.o: geodesics.f95
#	$(F95) -c geodesics.f95
#
#rad_pattern.o: rad_pattern.f95
#	$(F95) -c rad_pattern.f95
#
#get_stations_data.o: get_stations_data.f95
#	$(F95) -c get_stations_data.f95
#
#model_parameters.o: model_parameters.f95
#	$(F95) -c model_parameters.f95
#
#modelling_inputs.o: modelling_inputs.f95
#	$(F95) -c modelling_inputs.f95 #-fopenmp
#
#constants.o: constants.f95
#	$(F95) -c constants.f95

 
clean:
	\rm -f *.o
	\rm -f *.mod
