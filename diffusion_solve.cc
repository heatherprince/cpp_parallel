#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
//differential equation
#include "diffusion.h"
//solver
#include "euler.h"
//grid
#include "grid.h"




int main(int argc, char *argv[]) {
  if (argc != 2) {
    printf("USAGE: %s <nx> \n", argv[0]);
    exit(1);
  }

  const int nside = atoi(argv[1]);  //check that it is an int

  double x_max=M_PI;
  //double dx=x_max/nside;
  double kappa=1.;  //what should this be? does it matter?
  double t_max=0.5*M_PI*M_PI/kappa;


  //choose dt less than dx^2/4kappa, but must divide into tmax exactly
  //double dt_max = dx*dx/(4*kappa)
  //int N_min=floor(t_max/dt_max)
  //int nsteps=N_min*2
  int nsteps=2*2*nside*nside; //minumum nsteps is 2N^2 from def of t_max, put in more steps to be safe

  double dt=t_max/nsteps;

  Model *model=new Diffusion(kappa, nside, x_max);

  //const int dimen = model->dimen();

  //initial condition on T --> where to put this?? make T_grid class?
  //T = new double[nside][nside]; //can't do this because nside only given at runtime

  Grid *T=new Grid(nside, nside, x_max, x_max); //initializes grid to zero
  T->InitializeTEdges();                         //boundary conditions: starts with cos^2(x) and sin^2(x) at opposite edges

  Integrator *integrator = new Euler(dt, *model);

  double t = 0;
  for (int i = 0; i < nsteps; ++i) {
    integrator->Step(t, *T);
    T->InitializeTEdges(); //maintain BC
    t = (i+1) * dt;
  }

  //output to file
  T->WriteToFile("T_out.txt");
  delete T;
  delete integrator;
  delete model;
  return 0;
}
