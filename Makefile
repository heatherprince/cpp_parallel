integrators = euler.o
equations = diffusion.o
datatypes = grid.o
objects = diffusion_solve.o $(integrators) $(equations) $(datatypes)
test_grid_objects = test_grid.o $(datatypes)

CXXFLAGS = -g -Wall

all: heat_serial test_grid

heat_serial : $(objects)
	$(CXX) -o $@ $^

test_grid : $(test_grid_objects)
		$(CXX) -o $@ $^

clean:
	$(RM) *.o
	$(RM) .depend

depend:
	$(CXX) -MM $(CXXFLAGS) *.cc > .depend

-include .depend
