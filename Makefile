
dft.x: obj/main.o obj/subHello.o
	gfortran obj/main.o obj/subHello.o -o dft.x

obj/main.o: src/main.f90 obj
	gfortran -c -o obj/main.o src/main.f90

obj/subHello.o: src/subHello.f90 obj
	gfortran -c -o obj/subHello.o src/subHello.f90

obj:
	mkdir obj

clean:
	rm -f obj/*.o

clean-all:
	rm -rf obj
	rm -f dft.x

