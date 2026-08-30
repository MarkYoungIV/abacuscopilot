.SUFFIXES: .f90

#FC = ifort
FC = gfortran

FFLAGS = -O2

LINK = -static

OBJS = kind_mod.o       \
       matrix_mod.o     \
       ions_mod.o       \
       options_mod.o    \
       charge_mod.o     \
       chgcar_mod.o     \
       cube_mod.o       \
       io_mod.o         \
       bader_mod.o      \
       voronoi_mod.o    \
       multipole_mod.o  \
       weight_mod.o     \
       dsyevj3_mod.o    \
       jacobi_mod.o     \
       critpoint_mod.o  \
       automatedcp_mod.o



%.o %.mod : %.f90
	$(FC) $(FFLAGS) -c -Wunused  $*.f90

bader: $(OBJS) main.o
	rm -f bader 
	#$(FC) $(LINK) main.o -o $@ $(OBJS) 
	$(FC)  main.o -o $@ $(OBJS) 

dist: bader
	tar -cf bader_lnx_64.tar bader
	gzip -9 bader_lnx_64.tar

clean:
	rm -f *.o *.mod bader bader_lnx_64.tar.gz

