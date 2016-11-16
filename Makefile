integrators = euler.o
equations = diffusion.o
datatypes_serial = grid.o
datatypes_omp = grid_omp.o
objects = diffusion_solve.o $(integrators) $(equations) $(datatypes_serial)
objects_omp = diffusion_solve.o $(integrators) $(equations) $(datatypes_omp)
test_grid_objects = test_grid.o $(datatypes_serial)

CXX=g++-6
CXXFLAGS = -g -Wall -fopenmp

all: heat_serial heat_omp test_grid

heat_serial : $(objects)
	$(CXX) -o $@ $^

heat_omp : $(objects_omp)
	$(CXX) -fopenmp -o $@ $^

test_grid : $(test_grid_objects)
		$(CXX) -o $@ $^

clean:
	$(RM) *.o
	$(RM) .depend

depend:
	$(CXX) -MM $(CXXFLAGS) *.cc > .depend

-include .depend
