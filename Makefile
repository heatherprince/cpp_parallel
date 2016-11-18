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
export SHELL:=/bin/bash

heat_omp: CXXFLAGS = -g -Wall -qopenmp

heat_mpi: CXX=mpic++

all: heat_serial heat_omp heat_mpi

heat_serial : $(objects)
	. /usr/share/Modules/init/bash; module load intel/16.0; $(CXX) -o $@ $^

heat_omp : $(objects_omp)
	. /usr/share/Modules/init/bash; module load intel/16.0; \
	$(CXX) -qopenmp -o $@ $^


heat_mpi : $(objects_mpi)
	. /usr/share/Modules/init/bash; module load openmpi/intel-16.0 intel/16.0; \
	$(CXX) -o $@ $^

test_grid : $(test_grid_objects)
		$(CXX) -o $@ $^

clean:
	$(RM) *.o
	$(RM) .depend

depend:
	. /usr/share/Modules/init/bash; module load openmpi/intel-16.0 intel/16.0; \
	$(CXX) -MM $(CXXFLAGS) *.cc > .depend

-include .depend
