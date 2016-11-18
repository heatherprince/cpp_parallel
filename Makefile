integrators = euler.o
equations = diffusion.o
datatypes = grid.o
datatypes_omp = grid_omp.o
objects = diffusion_solve.o $(integrators) $(equations) $(datatypes)
objects_omp = diffusion_solve.o $(integrators) $(equations) $(datatypes_omp)
objects_mpi = diffusion_solve_mpi.o $(integrators) $(equations) $(datatypes)
test_grid_objects = test_grid.o $(datatypes)

CXX=icpc
CXXFLAGS = -g -Wall
heat_omp: CXXFLAGS = -g -Wall -qopenmp

all: heat_serial heat_omp test_grid

heat_serial : $(objects)
	$(CXX) -o $@ $^

heat_omp : $(objects_omp)
	$(CXX) -qopenmp -o $@ $^

test_grid : $(test_grid_objects)
		$(CXX) -o $@ $^

clean:
	$(RM) *.o
	$(RM) .depend

depend:
	$(CXX) -MM $(CXXFLAGS) *.cc > .depend

-include .depend
