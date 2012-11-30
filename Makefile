CC = gcc
CFLAGS = -g -O3 -march=corei7
F90 = gfortran
FFLAGS = -g -O3 -march=corei7
LIBS =

OBJ       = wigner3j.o wigner6j.o wigner9j.o mod_wigner.o

OBJTEST   = test_wigner.o test_wignerf.o

default: lib

lib: $(OBJ)
	ar -r libwignercoupling.a $(OBJ)

test: $(OBJ) $(OBJTEST)
	$(CC) $(CFLAGS) $(LIBS) $(OBJ) test_wigner.o -o test_wigner -L./ -lwignercoupling -lgfortran
	$(F90) $(FFLAGS) $(LIBS) $(OBJ) test_wignerf.o -o test_wignerf -L./ -lwignercoupling

%.o: %.c
	$(CC) $(CFLAGS) $(LIBS) -c $<

%.o: %.f90
	$(F90) $(FFLAGS) $(LIBS) -c $<
clean:
	-rm -f *.o *.a core *.mod 

