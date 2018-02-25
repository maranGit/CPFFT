#
#               makefile for FFT_finite_3d program
#
mkllib=$INTEL/mkl/lib/em64t
mklinc=$INTEL/mkl/include
FCCFLAG= -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_lapack95_lp64 -liomp5 -lpthread
#
#                all compiled .mod and .o files
#
AllModule = main_data.mod fft.mod elem_block_data.mod mkl_dfti.mod mkl_dft_type.mod
AllDotO = mod_main.o mod_fft.o mod_eleblocks.o mkl_dfti.o FFT_finite_3d.o rstgp1.o mm01.o recstr_allocate.o init.o
#
#                          link
#
FFT_finite_3d: $(AllDotO)
	ifort -o FFT_finite_3d.exe -O3 $(AllDotO) -I$(mklinc) -L$(mkllib) $(FCCFLAG) -g -traceback -gen-interfaces -warn interfaces -check -O0 -qopenmp
#
#                        compile
#
main_data.mod: mod_main.o mod_main.f
	ifort -c mod_main.f
mod_main.o: mod_main.f
	ifort -c mod_main.f
fft.mod: mod_fft.o mod_fft.f
	ifort -c mod_fft.f
mod_fft.o: mod_fft.f
	ifort -c mod_fft.f
elem_block_data.mod: mod_eleblocks.o mod_eleblocks.f
	ifort -c mod_eleblocks.f
mod_eleblocks.o: mod_eleblocks.f
	ifort -c mod_eleblocks.f
mkl_dfti.mod: mkl_dfti.o mkl_dfti.f
	ifort -c mkl_dfti.f
mkl_dft_type.mod: mkl_dfti.o mkl_dfti.f
	ifort -c mkl_dfti.f
mkl_dfti.o: mkl_dfti.f
	ifort -c mkl_dfti.f
FFT_finite_3d.o: $(AllModule) mm01.o rstgp1.o FFT_finite_3d.f
	ifort -c FFT_finite_3d.f -qopenmp
mm01.o: elem_block_data.mod mm01.f
	ifort -c mm01.f
rstgp1.o: elem_block_data.mod mm01.o rstgp1.f
	ifort -c rstgp1.f
recstr_allocate.o: elem_block_data.mod mm01.o recstr_allocate.f
	ifort -c recstr_allocate.f
init.o: init.f
	ifort -c init.f
#
#                        remove
#
clean:
	rm *.mod *.o FFT_finite_3d
