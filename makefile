#
#               makefile for FFT_finite_3d program
#
mkllib=$INTEL/mkl/lib/em64t
mklinc=$INTEL/mkl/include
FCCFLAG= -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_lapack95_lp64
OPENMP     = -qopenmp
#
#
#
#                all compiled .mod and .o files
#
AllModule = fft.mod elem_block_data.mod mkl_dfti.mod mkl_dft_type.mod
AllDotO = mod_fft.o mod_eleblocks.o mkl_dfti.o FFT_finite_3d.o\
				 	rstgp1.o mm01.o recstr_allocate.o init.o rplstr.o update.o\
					polar.o qmply1.o dupstr.o cep2A.o gptns1.o polarDecom.o
#
#                          link
#
FFT_finite_3d: $(AllDotO)
	ifort -o FFT_finite_3d.exe -O0 $(AllDotO) -I$(mklinc) -L$(mkllib) $(FCCFLAG) -g -traceback -gen-interfaces -warn interfaces $(OPENMP)
#
#                        compile
#
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
FFT_finite_3d.o: $(AllModule) mm01.o rstgp1.o cep2A.o gptns1.o FFT_finite_3d.f
	ifort -c FFT_finite_3d.f
mm01.o: elem_block_data.mod mm01.f
	ifort -c mm01.f
rstgp1.o: elem_block_data.mod mm01.o rstgp1.f
	ifort -c rstgp1.f
recstr_allocate.o: elem_block_data.mod mm01.o recstr_allocate.f
	ifort -c recstr_allocate.f
init.o: elem_block_data.mod mm01.o init.f
	ifort -c init.f
rplstr.o: elem_block_data.mod rplstr.f
	ifort -c rplstr.f
update.o: elem_block_data.mod update.f
	ifort -c update.f
polar.o: polar.f
	ifort -c polar.f
qmply1.o: qmply1.f
	ifort -c qmply1.f
dupstr.o: dupstr.f
	ifort -c dupstr.f
cep2A.o: cep2A.f
	ifort -c cep2A.f
gptns1.o: gptns1.f
	ifort -c gptns1.f
polarDecom.o: polarDecom.f
	ifort -c polarDecom.f
#
#                        remove
#
clean:
	rm *.mod *.o FFT_finite_3d.exe
