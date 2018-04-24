#
#               makefile for FFT_finite_3d program
#
mkllib=$INTEL/mkl/lib/em64t
mklinc=$INTEL/mkl/include
LIBMKL= -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_lapack95_lp64
OPENMP     = -qopenmp
COMPILER = ifort
RMC      = /bin/rm -f
foption= -xHOST -ip -fpconstant -traceback -align array64byte \
				 -qopenmp -reentrancy -no-inline-factor
#
INCLUDES = common.main param_def include_sig_up
MODS = mod_fft.o mod_eleblocks.o mkl_dfti.o mod_file_info.o
OBJ = $(MODS) FFT_finite_3d.o scan.o splunj.o star_com.o errmsg.o\
	rstgp1.o mm01.o recstr_allocate.o eblk_init.o rplstr.o update.o\
	polar.o qmply1.o dupstr.o cep2A.o gptns1.o drive_eps_sig.o\
	G_K_dF.o mm00.o FFT_init.o FFT_nr3.o mpi_code.o inmat.o\
	user_list.o inelem.o
#
#                          link
#
FFT_finite_3d: $(OBJ)
	$(COMPILER) -o FFT_finite_3d.exe -O3 $(OBJ) -I$(mklinc) -L$(mkllib) \
	$(LIBMKL) -g -traceback -gen-interfaces -warn interfaces $(OPENMP)
#
#       ---  default rule to make .o from .f ------
#            ===============================
#
#            => for simplicity we assume all the *.f depend on all
#               modules, includes, common.main, param_def
#
%.o: %.f $(MODS) $(INCLUDES)	 
	$(RMC)  $@
	$(COMPILER) $(foption) $< -c -o $@
#
#    ------   system level module files ----------------
#             =========================
#
mod_eleblocks.o : mod_eleblocks.f
	$(RMC)  $@
	$(COMPILER) $(foption) -c $<
mod_fft.o: mod_fft.f
	$(RMC)  $@
	$(COMPILER) $(foption) -c $<
mkl_dfti.o: mkl_dfti.f
	$(RMC)  $@
	$(COMPILER) $(foption) -c $<
mod_file_info.o: mod_file_info.f
	$(RMC)  $@
	$(COMPILER) $(foption) -c $<
#
#                        remove
#
clean:
	rm *.mod *.o FFT_finite_3d.exe
