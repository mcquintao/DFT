
dft.x: obj/main.o obj/diskOperations.o obj/dftOperations.o obj/mathLib.o
	gfortran obj/main.o obj/diskOperations.o obj/dftOperations.o obj/mathLib.o -o dft.x

obj/main.o: src/main.f90 obj
	gfortran -c -o obj/main.o src/main.f90

obj/diskOperations.o: src/diskOperations.f90 obj
	gfortran -c -o obj/diskOperations.o src/diskOperations.f90

obj/dftOperations.o: obj/mathLib.o obj src/dftOperations.f90
	gfortran -c -o obj/dftOperations.o src/dftOperations.f90

obj/mathLib.o: src/mathLib.f90 obj
	gfortran -c -o obj/mathLib.o src/mathLib.f90

obj:
	mkdir obj

clean:
	rm -f obj/*.o

clean-all:
	rm -rf obj
	rm -f dft.x

