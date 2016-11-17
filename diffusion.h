#ifndef DIFFUSION_H_
#define DIFFUSION_H_

#include "model.h"
// Heat diffusion equation:
//   \dot T = kappa * grad^2 T
// implement on 2D grid
class Diffusion : public Model {
 public:
  Diffusion(double kappa, int nside_x, int nside_y); //assumes that x and y both have same maximum xmax and same number of divisions nside
  ~Diffusion();
  int rhs(double t, const Grid &T, Grid &fx) const;   //T should be const, figure out how
  int dimen_() const { return nside_x_; }
  int dimen_y_() const { return nside_y_; }
 private:
  const double kappa_;
  const int nside_x_;
  const int nside_y_;
};

#endif  // DIFFUSION_H_
