#
# Makefile to run SiB2 as a subroutine
#
# Last revision: 20180208 by nelsonvn
#
# Usage:    make
#

# This line replaces the OBJDIR and MODDIR default
OBJDIR = ./obj
MODDIR = ./obj

# Default options
FC = gfortran
FFLAGS = -ffixed-line-length-none -freal-4-real-8 -O2 -J$(MODDIR)

MODEL_SRC = \
	comsibc_h.f90 \
	pardif_h.f90 \
	cal2jul.f90 \
	derive_trans_new.f \
	aero.f90 \
	momopt_stockli.f90 \
	SIB2sub.f90 \
	sib2_offline_lib.f90
#	sib2_offline_lib_test.f90

MODEL_OBJ = $(patsubst %.f90,%.o, $(patsubst %.f,%.o,$(MODEL_SRC)))

%.o : %.f90
	$(FC) -fpic $(FFLAGS) -c $< -o $(OBJDIR)/$(@F)

%.o : %.f
# there was missing -fpic
#	$(FC) $(FFLAGS) -c $< -o $(OBJDIR)/$(@F)
	$(FC) -fpic $(FFLAGS) -c $< -o $(OBJDIR)/$(@F)


all:	$(OBJDIR) $(MODEL_OBJ)
	$(FC) -shared -Wl,-z,relro -o sib2r.so $(OBJDIR)/*.o
#	$(FC) $(FFLAGS) sib2_offline_lib_test.f90 -L. -lsib2_my -o sib2_offline_lib_test.out

rlib:	
	MAKEFLAGS="FCFLAGS=-ffixed-line-length-none\ -freal-4-real-8\ -O2 FFLAGS=-ffixed-line-length-none\ -freal-4-real-8\ -O2 FLIBS= LIBR0= LIBR1= LIBR="  R CMD SHLIB -o sib2r.so $(MODEL_SRC)

$(OBJDIR):
	mkdir -p $(OBJDIR)

clean:
	rm -f $(MODDIR)/* $(OBJDIR)/* *.o *.mod

show:
#	The following command show the model source code list
	@echo $(MODEL_SRC)

tar:
	rm -f sib2_library.tar.bz2
	tar -jcf sib2_library.tar.bz2 $(MODEL_SRC) sib2_offline_lib_test.f90 \
	data2 sib2_library.par sib2_offline_lib.r sib2_offline_lib_test.r \
	Makefile_libsib2 00README.TXT
