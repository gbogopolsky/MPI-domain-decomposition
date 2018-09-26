SRCS = domain_decomposition.f90
EXECS = $(SRCS:.f90=)
MPIFC = mpiifort

all: $(EXECS)

%: %.f90
	$(MPIFC) -o $@ $<

clean:
	rm -fv $(EXECS) *.mod
