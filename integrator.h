#ifndef INTEGRATOR_H_
#define INTEGRATOR_H_
//for a grid
class Grid;
class Integrator {
 public:
  virtual ~Integrator() {}
  virtual int Step(double t, Grid &x) = 0;
};

#endif  // INTEGRATOR_H_
