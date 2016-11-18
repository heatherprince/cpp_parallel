#include "diffusion.h"
#include "grid.h"
#include <math.h>

//xmax and ymax to get dx, dy
Diffusion::Diffusion(double kappa, int nside_x, int nside_y)
    : kappa_(kappa),
      nside_x_(nside_x),
      nside_y_(nside_y){
}

Diffusion::~Diffusion(){
}


int Diffusion::rhs(double t, const Grid &T, Grid &fx) const {   //T should be const, figure out how
  T.GradSq(fx); 
  fx.MultiplyByConstant(kappa_); // set fx to kappa*grad^2 T
  return 0;
}
