#
#               makefile for FFT_finite_3d program
#
O        = .o
OD       = /home/rma8/ti6264/FFT/RanFFT/testAll/CPFFT/obj_linux
RMC      = /bin/rm -f
MVC      = /bin/mv -f
#
#               mkl library
#
mkllib=$INTEL/mkl/lib/em64t
mklinc=$INTEL/mkl/include
LIBMKL= -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_lapack95_lp64
#
#               openmp option
#
OPENMP     = -qopenmp
COMPILER = ifort
foption= -xHOST -ip -fpconstant -traceback -align array64byte \
				 -qopenmp -reentrancy -no-inline-factor
#
#               all target files
#
INCLUDES = common.main param_def include_sig_up
MODS = $(OD)/mod_fft.o \
	$(OD)/mod_eleblocks.o \
	$(OD)/mkl_dfti.o \
	$(OD)/mod_file_info.o
OBJ = $(MODS) \
	$(OD)/FFT_finite_3d.o \
	$(OD)/scan.o \
	$(OD)/splunj.o \
	$(OD)/star_com.o \
	$(OD)/errmsg.o\
	$(OD)/rstgp1.o \
	$(OD)/mm01.o \
	$(OD)/recstr_allocate.o \
	$(OD)/eblk_init.o \
	$(OD)/rplstr.o \
	$(OD)/update.o\
	$(OD)/polar.o \
	$(OD)/qmply1.o \
	$(OD)/dupstr.o \
	$(OD)/cep2A.o \
	$(OD)/gptns1.o \
	$(OD)/drive_eps_sig.o\
	$(OD)/G_K_dF.o \
	$(OD)/mm00.o \
	$(OD)/FFT_init.o \
	$(OD)/FFT_nr3.o \
	$(OD)/mpi_code.o \
	$(OD)/inmat.o\
	$(OD)/user_list.o \
	$(OD)/inelem.o \
	$(OD)/inlod.o \
	$(OD)/indypm.o \
	$(OD)/oudriv.o \
	$(OD)/ouneut.o \
	$(OD)/oumodel.o
#
#                          link
#
FFT_finite_3d.exe: $(OBJ)
	$(COMPILER) -o $@ -O3 $(OBJ) -I$(mklinc) -L$(mkllib) \
	$(LIBMKL) -g -traceback -gen-interfaces -warn interfaces -module $(OD) $(OPENMP)
	$(MVC) $@ ./../obj_linux/$@
	chmod ugo+rx  ./../obj_linux/$@
#
#       ---  default rule to make .o from .f ------
#            ===============================
#
#            => for simplicity we assume all the *.f depend on all
#               modules, includes, common.main, param_def
#
$(OD)/%.o: %.f $(MODS) $(INCLUDES)	 
	$(RMC)  $@
	$(COMPILER) $(foption) $< -c -o $@
#
#    ------   system level module files ----------------
#             =========================
#
$(OD)/mod_eleblocks$O : mod_eleblocks.f
	$(RMC)  $@
	$(COMPILER) $(foption) $< -c -o $@
$(OD)/mod_fft$O: mod_fft.f
	$(RMC)  $@
	$(COMPILER) $(foption) $< -c -o $@
$(OD)/mkl_dfti$O: mkl_dfti.f
	$(RMC)  $@
	$(COMPILER) $(foption) $< -c -o $@
$(OD)/mod_file_info$O: mod_file_info.f
	$(RMC)  $@
	$(COMPILER) $(foption) $< -c -o $@
#
#                        remove
#
clean:
	rm $(OD)/*