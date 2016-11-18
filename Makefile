integrators = euler.o
equations = diffusion.o
datatypes = grid.o
objects_mpi = diffusion_solve_mpi.o $(integrators) $(equations) $(datatypes)

CXX=mpic++
CXXFLAGS = -g -Wall -fopenmp

all: heat_mpi

heat_mpi : $(objects_mpi)
	$(MPICXX) -o $@ $^

clean:
	$(RM) *.o
	$(RM) .depend

depend:
	$(CXX) -MM $(CXXFLAGS) *.cc > .depend

-include .depend
