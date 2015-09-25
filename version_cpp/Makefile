CC=mpicxx
CFLAGS=-O3 -Wall -std=c++0x
OBJ=file.o force.o thermostat.o \
	 main.o my_malloc.o cpff.o bicg_stab.o

%.o: %.cpp
	$(CC) -c -o $@ $< $(CFLAGS)

capacmd: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

clean:  
	rm -f mpi-main main capacmd *.o
