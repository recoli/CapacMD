CC=mpicxx
CFLAGS=-O2 -std=c++0x -Wall -lm 
OBJ=file.o force.o thermostat.o \
	 mpi-main.o my_malloc.o cpff.o bicg_stab.o

%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)

mpi-main: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

clean:  
	rm -f mpi-main *.o
