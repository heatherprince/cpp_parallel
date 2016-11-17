#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
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


  const int nside = atoi(argv[1]);  //check that it is an int
  const double x_max=M_PI;
  const double kappa=1.;  //what should this be? does it matter?
  const double t_max=0.5*M_PI*M_PI/kappa;
  //choose dt less than dx^2/4kappa, but must divide into tmax exactly
  const int nsteps=2*2*nside*nside; //minumum nsteps is 2N^2 from def of t_max, put in more steps to be safe
  const double dt=t_max/nsteps;
  const double dx=xmax/(nside-1);

  int my_rank, size, prev, next, tag1=1, tag2=2, root_process=0;
  double mean_temp;



  MPI_Request reqs[4];
  MPI_Status stats[4];

  MPI_Init (&argc, &argv); //initialize MPI library
  MPI_Comm_size(MPI_COMM_WORLD, &size); //get number of processes
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); //get my process id

  double t1, t2;
  t1 = MPI_Wtime();



  prev = my_rank-1;
  next = my_rank+1;
  if (my_rank == 0) prev = size - 1;
  if (my_rank == (size - 1)) next = 0;

  printf("I am process %3d with neighbours %3d and %3d", my_rank, previous, next)

  Model *model=new Diffusion(kappa, nside, x_max);

  //divide up x's between ranks
  int nside_x=nside/size+2; //check if it goes in exactly! //divide equally then give an extra column on each side
  double my_x_min=rank*(xmax/size)-dx;
  double my_x_max=(rank+1)*(xmax/size);

  int nside_y=nside;
  double y_min=0.;
  double y_max=x_max;

  Grid *T=new Grid(nside_x, nside_y, my_x_min, my_x_max, ymin, y_max); //initializes grid to zero
  T->InitializeTEdges();                         //boundary conditions: cos^2(x) and sin^2(x) at opposite edges

  Integrator *integrator = new Euler(dt, *model);

  double t = 0;
  for (int i = 0; i < nsteps; ++i) {
    integrator->Step(t, *T);
    //pass edge columns in ring formation
    T->InitializeTEdges(); //maintain BC : cos^2(x) and sin^2(x) at opposite edges

    //pass end columns in a ring topology

    //done passing end columns
    t = (i+1) * dt;
  }

  if (my_rank==root_process){
    //T_full=new full T grid with data from all nodes
    char filename[] = "T_out.txt";
    //T_full->WriteToFile(filename);
    T_full->WriteToFile(filename);
  }
  double my_mean_temp=T->GetMean();
  //sum across all processes
  MPI_Reduce(&my_mean_temp, &mean_temp, 1, MPI_DOUBLE, MPI_SUM, root_process, MPI_COMM_WORLD);
  //reduction --> get mean of all different nodes and get master to write it
  if (my_rank==root_process){
    printf("The mean temperature for nside=%4d is: %8.4f \n", nside, mean_temp);
  }

  delete T;
  delete integrator;
  delete model;
  t2 = MPI_Wtime();
  if (my_rank==root_process){
    printf( "Elapsed time is %f\n", t2 - t1 );
  }
  MPI_Finalize(); //MPI cleanup


  return 0;
}
