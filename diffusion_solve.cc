#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <ctime>
#include <sstream>
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
  time_t start = time(NULL);

  const int nside = atoi(argv[1]);  //check that it is an int

  double x_max=M_PI;
  double kappa=1.;  //what should this be? does it matter?

  double t_max=0.5*M_PI*M_PI/kappa;
  //choose dt less than dx^2/4kappa, but must divide into tmax exactly
  int nsteps=2*2*nside*nside; //minumum nsteps is ~2 N^2 from def of t_max, put in more steps to be safe
  double dt=t_max/nsteps;

  Model *model=new Diffusion(kappa, nside, nside);
  Integrator *integrator = new Euler(dt, *model);

  Grid *T=new Grid(nside, nside, 0, 0, x_max, x_max);    //initializes grid to zero
  T->InitializeTEdges();                                 //boundary conditions: starts with cos^2(x) and sin^2(x) at opposite edges

  double t = 0;
  for (int i = 0; i < nsteps; ++i) {
    integrator->Step(t, *T);
    T->InitializeTEdges(); //maintain BC
    t = (i+1) * dt;
  }

  //output to file
  std::stringstream filename;
  filename << "OutputDatafiles/T_out_nside"<<nside<<"_serial.txt";
  T->WriteToFile(filename.str().c_str());


  double mean_temp=T->GetMean();
  printf("The mean temperature for nside=%4d is: %8.4f \n", nside, mean_temp);

  delete T;
  delete integrator;
  delete model;

  time_t end = time(NULL);
  double time_s=difftime(end,start);
  printf("Time elapsed running main of diffusion_solve is: %8.2f seconds \n", time_s);

  return 0;
}
