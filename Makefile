CC = gcc
CFLAGS = -g -O3 -march=native
LIBS =

OBJ       = wigner3j.o wigner6j.o wigner9j.o

OBJTEST   = test_wigner.o

default: lib

lib: $(OBJ)
	ar -r libwignercoupling.a $(OBJ)

test: $(OBJ) $(OBJTEST)
	$(CC) $(CFLAGS) $(LIBS) $(OBJ) $(OBJTEST) -o test_wigner -L./ -lwignercoupling

%.o: %.c
	$(CC) $(CFLAGS) $(LIBS) -c $<

clean:
	-rm -f *.o *.a core 

