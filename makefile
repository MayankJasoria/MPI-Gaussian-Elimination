CC=mpicc

CFLAGS=-o

DEPS_PROG=gaussian_parallel.c

program:
	$(CC) $(DEPS_PROG) $(CFLAGS) a.out

generate:
	gcc input_generate.c -o generate.out

clean:
	rm -rf *.o *.out