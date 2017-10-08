
dft.x: obj/main.o obj/subHello.o obj/diskOperations.o
	gfortran obj/main.o obj/subHello.o obj/diskOperations.o -o dft.x

obj/main.o: src/main.f90 obj
	gfortran -c -o obj/main.o src/main.f90

obj/subHello.o: src/subHello.f90 obj
	gfortran -c -o obj/subHello.o src/subHello.f90

obj/diskOperations.o: src/diskOperations.f90 obj
	gfortran -c -o obj/diskOperations.o src/diskOperations.f90

obj:
	mkdir obj

clean:
	rm -f obj/*.o

clean-all:
	rm -rf obj
	rm -f dft.x

