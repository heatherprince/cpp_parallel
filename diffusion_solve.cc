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
#include <ctime>



int main(int argc, char *argv[]) {
  if (argc != 2) {
    printf("USAGE: %s <nx> \n", argv[0]);
    exit(1);
  }
  time_t start = time(NULL);

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

  Model *model=new Diffusion(kappa, nside, nside);

  //const int dimen = model->dimen();

  //initial condition on T --> where to put this?? make T_grid class?
  //T = new double[nside][nside]; //can't do this because nside only given at runtime

  Grid *T=new Grid(nside, nside, x_max, x_max); //initializes grid to zero
  T->InitializeTEdges();                         //boundary conditions: starts with cos^2(x) and sin^2(x) at opposite edges

  Integrator *integrator = new Euler(dt, *model);

  time_t start_integration = time(NULL);
  double t = 0;
  for (int i = 0; i < nsteps; ++i) {
    integrator->Step(t, *T);
    T->InitializeTEdges(); //maintain BC
    t = (i+1) * dt;
  }
  time_t end_integration = time(NULL);

  //output to file
  char filename[] = "T_out.txt";
  T->WriteToFile(filename);
  double mean_temp=T->GetMean();
  printf("The mean temperature for nside=%4d is: %8.4f \n", nside, mean_temp);

  delete T;
  delete integrator;
  delete model;

  time_t end = time(NULL);

  double time_s=difftime(end,start);

  double time_integration_s=difftime(end_integration,start_integration);

  printf("Time elapsed running main of diffusion_solve is: %8.2f seconds \n", time_s);
  printf("Time elapsed integrating the diffusion equation is: %8.2f seconds \n", time_integration_s);

  return 0;
}
