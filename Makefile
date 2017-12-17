# COMPILE='ifort -O2 -r8 -check all'

COMPILE=gfortran -O2

dft.x: obj/main.o obj/diskOperations.o obj/dftOperations.o obj/mathLib.o
	${COMPILE}   obj/main.o obj/diskOperations.o obj/dftOperations.o obj/mathLib.o -o dft.x

obj/main.o: src/main.f90 obj/dftOperations.o obj/diskOperations.o obj
	${COMPILE}  -c -o obj/main.o src/main.f90

obj/diskOperations.o: src/diskOperations.f90 obj
	${COMPILE}  -c -o obj/diskOperations.o src/diskOperations.f90

obj/dftOperations.o: obj/mathLib.o obj obj/diskOperations.o src/dftOperations.f90
	${COMPILE}   -c -o obj/dftOperations.o src/dftOperations.f90

obj/mathLib.o: src/mathLib.f90 obj
	${COMPILE}  -c -o obj/mathLib.o src/mathLib.f90

obj:
	mkdir obj

clean:
	rm -rf obj
	rm -f *.mod
	rm -f dft.x

