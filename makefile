#
#               makefile for FFT_finite_3d program
#
mkllib=$INTEL/mkl/lib/em64t
mklinc=$INTEL/mkl/include
FCCFLAG= -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_lapack95_lp64
OPENMP     = -qopenmp
foption= -xHOST -ip -fpconstant -traceback -align array64byte \
				 -qopenmp -reentrancy -no-inline-factor
#
#
#
#                all compiled .mod and .o files
#
AllModule = fft.mod elem_block_data.mod mkl_dfti.mod mkl_dft_type.mod
AllDotO = mod_fft.o mod_eleblocks.o mkl_dfti.o FFT_finite_3d.o\
				 	rstgp1.o mm01.o recstr_allocate.o init.o rplstr.o update.o\
					polar.o qmply1.o dupstr.o cep2A.o gptns1.o
#
#                          link
#
FFT_finite_3d: $(AllDotO)
	ifort -o FFT_finite_3d.exe -O3 $(AllDotO) -I$(mklinc) -L$(mkllib) $(FCCFLAG) -g -traceback -gen-interfaces -warn interfaces $(OPENMP)
#
#                        compile
#
fft.mod: mod_fft.o mod_fft.f
	ifort $(foption) -c mod_fft.f
mod_fft.o: mod_fft.f
	ifort $(foption) -c mod_fft.f
elem_block_data.mod: mod_eleblocks.o mod_eleblocks.f
	ifort $(foption) -c mod_eleblocks.f
mod_eleblocks.o: mod_eleblocks.f
	ifort $(foption) -c mod_eleblocks.f
mkl_dfti.mod: mkl_dfti.o mkl_dfti.f
	ifort $(foption) -c mkl_dfti.f
mkl_dft_type.mod: mkl_dfti.o mkl_dfti.f
	ifort $(foption) -c mkl_dfti.f
mkl_dfti.o: mkl_dfti.f
	ifort $(foption) -c mkl_dfti.f
FFT_finite_3d.o: $(AllModule) mm01.o rstgp1.o cep2A.o gptns1.o FFT_finite_3d.f
	ifort $(foption) -c FFT_finite_3d.f
mm01.o: elem_block_data.mod mm01.f
	ifort $(foption) -c mm01.f
rstgp1.o: elem_block_data.mod mm01.o rstgp1.f
	ifort $(foption) -c rstgp1.f
recstr_allocate.o: elem_block_data.mod mm01.o recstr_allocate.f
	ifort $(foption) -c recstr_allocate.f
init.o: elem_block_data.mod mm01.o init.f
	ifort $(foption) -c init.f
rplstr.o: elem_block_data.mod rplstr.f
	ifort $(foption) -c rplstr.f
update.o: elem_block_data.mod update.f
	ifort $(foption) -c update.f
polar.o: polar.f
	ifort $(foption) -c polar.f
qmply1.o: qmply1.f
	ifort $(foption) -c qmply1.f
dupstr.o: dupstr.f
	ifort $(foption) -c dupstr.f
cep2A.o: cep2A.f
	ifort $(foption) -c cep2A.f
gptns1.o: gptns1.f
	ifort $(foption) -c gptns1.f
#
#                        remove
#
clean:
	rm *.mod *.o FFT_finite_3d.exe
