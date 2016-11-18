#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <mpi.h>
#include <sstream>
//differential equation
#include "diffusion.h"
//solver
#include "euler.h"
//grid
#include "grid_mpi.h"




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
  const double dx=x_max/(nside-1);

  int my_rank, size, prev, next, tag1=1, tag2=2, root_process=0;
  double mean_temp;


  double *col_for_prev=new double[nside];
  double *col_for_next=new double[nside];
  double *col_from_prev=new double[nside];
  double *col_from_next=new double[nside];




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

  printf("I am process %3d with neighbours %3d and %3d\n", my_rank, prev, next);


  //divide up x's between ranks
  int nside_x=nside/size+2; //check if it goes in exactly! //divide equally then give an extra column on each side
  double my_x_min=dx*(my_rank*nside/size - 1);  //first rank will have xmin negative, ok because BC cos^2, sin^2 periodic
  double my_x_max=dx*((my_rank+1)*nside/size);  //last rank will have xmax larger than pi, ok because BC cos^2, sin^2 periodic

  printf("I am process %3d with xmin %5.2f, xmax %5.2f, and nside_x %d\n", my_rank, my_x_min, my_x_max, nside_x);

  int nside_y=nside;
  double my_y_min=0.;
  double my_y_max=x_max;
  Model *model=new Diffusion(kappa, nside_x, nside_y);
  Grid *T=new Grid(nside_x, nside_y, my_x_min, my_y_min,  my_x_max, my_y_max); //initializes grid to zero
  T->InitializeTEdges();                         //boundary conditions: cos^2(x) and sin^2(x) at opposite edges
  //if (my_rank==root_process){
  //  for(int i=0; i<nside_x; i++){
  //    for(int j=0; j<nside_y; j++){
  //    printf("T i: %4d, j: %4d, T: %5.2f \n",i,j,T->Get(i,j));
  //    }
  //  }
  //}
  Integrator *integrator = new Euler(dt, *model);

  double t = 0;
  for (int i = 0; i < nsteps; ++i) {
    integrator->Step(t, *T);
    //pass edge columns in ring formation
    T->InitializeTEdges(); //maintain BC : cos^2(x) and sin^2(x) at opposite edges

    //get end columns to pass to previous and next processes
    T->GetYColumn(1, col_for_prev); //first column apart from border
    T->GetYColumn(nside_x-2, col_for_next); //last column before border

    //pass end columns in a ring topology: col_for_prev is sent to col_from_next of previous processor, col_for_next to col_from_prev of next processor
    MPI_Irecv(col_from_prev, nside_y, MPI_DOUBLE, prev, tag1, MPI_COMM_WORLD, &reqs[0]);
    MPI_Irecv(col_from_next, nside_y, MPI_DOUBLE, next, tag2, MPI_COMM_WORLD, &reqs[1]);
    MPI_Isend(col_for_prev, nside_y, MPI_DOUBLE, prev, tag2, MPI_COMM_WORLD, &reqs[2]);
    MPI_Isend(col_for_next, nside_y, MPI_DOUBLE, next, tag1, MPI_COMM_WORLD, &reqs[3]);
    MPI_Waitall(4, reqs, stats);

    //put columns received in border
    T->SetYColumn(0, col_from_prev);
    T->SetYColumn(nside_x-1, col_from_next);

    //printf("Time iteration %d out of %d: process %3d sent and updated end columns \n", i, nsteps, my_rank);
    //done passing end columns
    t = (i+1) * dt;
  }
  //if (my_rank==root_process){
  //  for(int i=0; i<nside_x; i++){
  //    for(int j=0; j<nside_y; j++){
  //    printf("T i: %4d, j: %4d, T: %5.2f \n",i,j,T->Get(i,j));
  //    }
  //  }
  //}
  //if (my_rank==root_process){
    //T_full=new full T grid with data from all nodes
    //char filename[] = "T_out.txt";
    //T_full->WriteToFile(filename);
    //T->WriteToFile(filename);
  //}
  //I decided to write my files separately so that I could check that the border column swop was working
  //and also because I know that I can combine them really easily in Python using numpy before plotting
  //but if I had more time I would use MPI_IO to write to the same file from each process
  std::stringstream filename;

  filename << "OutputDatafiles/T_out_nside"<<nside<<"_process" << my_rank  <<  "_numproc" << size << ".txt";
  T->WriteToFile(filename.str().c_str());
  double my_mean_temp=T->GetMeanExcludeBorders();
  if (my_rank==root_process){
    printf("The mean temperature for the root process for nside=%4d is: %8.4f \n", nside, my_mean_temp);
  }
  //sum across all processes
  MPI_Reduce(&my_mean_temp, &mean_temp, 1, MPI_DOUBLE, MPI_SUM, root_process, MPI_COMM_WORLD);
  mean_temp/=size;
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
