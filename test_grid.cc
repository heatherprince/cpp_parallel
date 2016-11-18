#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctime>
#include "grid_mpi.h"

int main(int argc, char *argv[]) {
  if (argc != 2) {
    printf("USAGE: %s <nx> \n", argv[0]);
    exit(1);
  }

  time_t start = time(NULL);

  const int nside = atoi(argv[1]);  //check that it is an int

  double x_max=M_PI;
  double kappa=100.;



  Grid *T=new Grid(nside, nside, 0,0,x_max, x_max); //initializes grid to zero



  for(int j=0; j<nside; j++){
    for(int i=0; i<nside; i++){
      printf("T i: %4d, j: %4d, T: %5.2f \n",i,j,T->Get(i,j));
    }
  }

  T->InitializeTEdges();

  char filename[] = "file.txt";
  T->WriteToFile(filename);

  for(int j=0; j<nside; j++){
    for(int i=0; i<nside; i++){
      printf("T bound i: %4d, j: %4d, T: %5.2f \n",i,j,T->Get(i,j));
    }
  }

  Grid *gradsq=new Grid(nside, nside, 0,0, x_max, x_max);

  T->GradSq(*gradsq);
  for(int j=0; j<nside; j++){
    for(int i=0; i<nside; i++){
      printf("GradsqT i: %4d, j: %4d, T: %5.2f \n",i,j,gradsq->Get(i,j));
    }
  }

  gradsq->MultiplyByConstant(kappa);
  for(int j=0; j<nside; j++){
    for(int i=0; i<nside; i++){
      printf("kappa*GradsqT i: %4d, j: %4d, T: %5.2f \n",i,j,gradsq->Get(i,j));
    }
  }

  T->AddGridTimesConst(1., *gradsq);
  for(int j=0; j<nside; j++){
    for(int i=0; i<nside; i++){
      printf("T+kappa*GradsqT i: %4d, j: %4d, T: %5.2f \n",i,j,T->Get(i,j));
    }
  }

  delete T;
  delete gradsq;


  time_t end = time(NULL);

  double time_s=difftime(end,start);

  printf("Time elapsed running main of test_grid is: %8.2f seconds \n", time_s);

  return 0;
}
