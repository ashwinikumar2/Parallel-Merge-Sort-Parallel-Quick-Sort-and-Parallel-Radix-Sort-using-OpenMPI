
output: tester.o sort.o
	mpiCC tester.o sort.o -o output

tester.o: tester.cpp
	mpiCC -c tester.cpp

sort.o: sort.cpp sort.h
	mpiCC -c sort.cpp

run:
	mpirun -np 6 ./output
clean:
	rm *.o output