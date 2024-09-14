FC=gfortran
OBJ=engines.f90 Generate_Rates.f90
compute_rate.exe: $(OBJ)
	$(FC) -o $@ $(OBJ)

clean:
	rm -f *.o
	rm -f *.mod
	rm -f *.exe
	rm -f *.out
	rm -f MQCT_TACS.dat

