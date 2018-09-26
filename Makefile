SRCS = domain_decomposition.f90
EXECS = $(SRCS:.f90=)
MPIFC = mpiifort
FLAGS =

all: $(EXECS)

%: %.f90
	$(MPIFC) $(FLAGS) -o $@ $<

run:
	mpirun -n 2 ./$(EXECS)

clean:
	rm -fv $(EXECS) *.mod
