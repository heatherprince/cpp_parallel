#ifndef EULER_H_
#define EULER_H_

#include "integrator.h"
class Model;
//euler on 2D grid
class Euler : public Integrator {
 public:
  Euler(double dt, const Model &model);
  ~Euler();
  int Step(double t, Grid &T);
 private:
  const int dimen_;                     // dimension of state x
  const int dimen_y_;                     // second dimension of state x
  const double dt_;                     // timestep
  const Model &model_;                  // functor to evaluate f(x,t)
  Grid *fx_;                          // temporary space to hold f(x,t)
};

#endif  // EULER_H
