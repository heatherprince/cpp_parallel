#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <ctime>
#include <mpi.h>
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
  time_t start = time(nullptr);

  const int nside = atoi(argv[1]);  //check that it is an int
  double x_max=M_PI;
  double kappa=1.;  //what should this be? does it matter?
  double t_max=0.5*M_PI*M_PI/kappa;
  //choose dt less than dx^2/4kappa, but must divide into tmax exactly
  int nsteps=2*2*nside*nside; //minumum nsteps is 2N^2 from def of t_max, put in more steps to be safe
  double dt=t_max/nsteps;

  int my_rank, size, prev, next, tag1=1, tag2=2;


  MPI_Request reqs[4];
  MPI_Status stats[4];

  MPI_Init (&argc, &argv); //initialize MPI library
  MPI_Comm_size(MPI_COMM_WORLD, &size); //get number of processes
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); //get my process id

  prev = my_rank-1;
  next = my_rank+1;
  if (my_rank == 0) prev = size - 1;
  if (my_rank == (size - 1)) next = 0;

  printf("I am process %3d with neighbours %3d and %3d", my_rank, previous, next)



  Model *model=new Diffusion(kappa, nside, x_max);

  int nside_x=nside;
  int nside_y=nside/size; //check if it goes in exactly!
  y_max=x_max/size;

  Grid *T=new Grid(nside_x, nside_y, x_max, y_max); //initializes grid to zero
  T->InitializeTEdges();                         //boundary conditions: starts with cos^2(x) and sin^2(x) at opposite edges
  //make MPI versions of Grid functions!


  Integrator *integrator = new Euler(dt, *model);

  time_t start_integration = time(nullptr);
  double t = 0;
  for (int i = 0; i < nsteps; ++i) {
    integrator->Step(t, *T);
    //pass edge columns in ring formation
    T->InitializeTEdges(); //maintain BC
    t = (i+1) * dt;
  }
  time_t end_integration = time(nullptr);

  //output to file
  char filename[] = "T_out.txt";
  T->WriteToFile(filename);
  double mean_temp=T->GetMean();
  printf("The mean temperature for nside=%4d is: %8.4f \n", nside, mean_temp);

  delete T;
  delete integrator;
  delete model;

  MPI_Finalize(); //MPI cleanup

  time_t end = time(nullptr);

  double time_s=difftime(end,start);

  double time_integration_s=difftime(end_integration,start_integration);

  printf("Time elapsed running main of diffusion_solve is: %8.2f seconds \n", time_s);
  printf("Time elapsed integrating the diffusion equation is: %8.2f seconds \n", time_integration_s);

  return 0;
}
