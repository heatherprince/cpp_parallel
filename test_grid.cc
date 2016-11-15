#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>

#include "grid.h"

int main(int argc, char *argv[]) {
  if (argc != 2) {
    printf("USAGE: %s <nx> \n", argv[0]);
    exit(1);
  }

  const int nside = atoi(argv[1]);  //check that it is an int

  double x_max=M_PI;
  double kappa=100.;

  Grid *T=new Grid(nside, nside, x_max, x_max); //initializes grid to zero

  for(int j=0; j<nside; j++){
    for(int i=0; i<nside; i++){
      printf("T i: %4d, j: %4d, T: %5.2f \n",i,j,T->Get(i,j));
    }
  }

  T->InitializeTEdges();

  T->WriteToFile("file.txt");

  for(int j=0; j<nside; j++){
    for(int i=0; i<nside; i++){
      printf("T bound i: %4d, j: %4d, T: %5.2f \n",i,j,T->Get(i,j));
    }
  }

  Grid *gradsq=new Grid(nside, nside, x_max, x_max);

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

  return 0;
}
