#include <stdio.h>
#include <math.h>
#include "euler.h"
#include "model.h"
#include "grid.h"

//implementation of Euler that works on a square grid
Euler::Euler(double dt, const Model &model)
    : dimen_(model.dimen_()),
      dimen_y_(model.dimen_y_()),
      dt_(dt),
      model_(model) {
  fx_ = new Grid(dimen_, dimen_y_, M_PI, M_PI); // don't actually need xmax, ymax, so the last two arguments are unimportant
}

Euler::~Euler() {
  delete fx_;
}

int Euler::Step(double t, Grid &x) {
  model_.rhs(t, x, *fx_);
  x.AddGridTimesConst(dt_, *fx_);
  return 0;
}
