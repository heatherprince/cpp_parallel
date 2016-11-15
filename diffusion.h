#ifndef DIFFUSION_H_
#define DIFFUSION_H_

#include "model.h"
// Heat diffusion equation:
//   \dot T = kappa * grad^2 T
// implement on 2D grid
class Diffusion : public Model {
 public:
  Diffusion(double kappa, int nside, int xmax); //assumes that x and y both have same maximum xmax and same number of divisions nside
  ~Diffusion();
  int rhs(double t, const Grid &T, Grid &fx) const;   //T should be const, figure out how
  int dimen() const { return nside_; }
 private:
  const double kappa_;
  const double nside_;    //grid size is nside^2
  const double xmax_;
};

#endif  // DIFFUSION_H_
