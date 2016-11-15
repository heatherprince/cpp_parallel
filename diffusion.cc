#include "diffusion.h"
#include "grid.h"
#include <math.h>

//xmax and ymax to get dx, dy
Diffusion::Diffusion(double kappa, int nside, int xmax)
    : kappa_(kappa),
      nside_(nside),
      xmax_(xmax){
}

Diffusion::~Diffusion(){
}

int Diffusion::rhs(double t, const Grid &T, Grid &fx) const {   //T should be const, figure out how
  T.GradSq(fx); //check usage, should result in fx being grad sq T
  fx.MultiplyByConstant(kappa_); // set fx to kappa*grad^2 T
  return 0;
}
