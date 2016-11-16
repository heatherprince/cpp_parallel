#include <stdio.h>
#include <math.h>
#include "grid.h"

//implementation of Euler that works on a square grid
Grid::Grid(int nside_x, int nside_y, double x_max, double y_max)
    : dimen_x_(nside_x),
      dimen_y_(nside_y),
      max_x_(x_max),
      max_y_(y_max),
      dx_(x_max/(nside_x-1.)),
      dy_(y_max/(nside_y-1.)) {
  grid_ = new double*[dimen_y_];
  grid_[0] = new double[dimen_x_*dimen_y_];
  for (int j = 1; j < dimen_y_; j++){
    grid_[j] = grid_[j-1] + dimen_x_;
  } //I want x strips to be continuous in memory so I can divide up y in domain decomposition because of how boundary conditions are set

  //initialize all to zero by default, initialize edges for T in a separate function
  #pragma omp parallel for default(none)
  for(int j=0; j<dimen_y_; j++){
    for(int i=0; i<dimen_x_; i++){
      grid_[j][i]=0.;
    }
  }


}

Grid::~Grid() {
  delete[] grid_[0];
  delete[] grid_;
}

//i is x index, j is y index
double Grid::Get(int i, int j)const {
  //check i, i in correct range
  return grid_[j][i];
}

int Grid::Set(int i, int j, double val){
  //check i, i in correct range
  grid_[j][i]=val;
  return 0;
}

double Grid::GetMean()const {
  //check i, i in correct range
  //calculate mean temp
  double sum=0;
  #pragma omp parallel for default(none) reduction(+:sum)
  for(int j=0; j<dimen_y_; j++){
    for(int i=0; i<dimen_x_; i++){
      sum+=grid_[j][i];
    }
  }
  double mean=sum/(dimen_x_*dimen_y_);
  return mean;
}

int Grid::InitializeTEdges(){
  double x;
  #pragma omp parallel for default(none) private(x)
  for(int i=0; i<dimen_x_; i++){
    x=i*dx_;
    grid_[0][i]=cos(x)*cos(x);
    grid_[dimen_x_-1][i]=sin(x)*sin(x);
  }
  return 0;
}

int Grid::MultiplyByConstant(double c){
  #pragma omp parallel for default(none) shared(c)
  for(int j=0; j<dimen_y_; j++){
    for(int i=0; i<dimen_x_; i++){
      grid_[j][i]*=c;
    }
  }
  return 0;
}

int Grid::AddGridTimesConst(double c, Grid &grid2){
  #pragma omp parallel for default(none) shared(c, grid2)
  for(int j=0; j<dimen_y_; j++){
    for(int i=0; i<dimen_x_; i++){
      grid_[j][i]+=c*grid2.Get(i,j);
    }
  }
  return 0;
}

int Grid::GradSq(Grid &grad_sq_T)const  {
  double f;
  #pragma omp parallel for default(none) shared(grad_sq_T) private(f)
  for(int j=0; j<dimen_y_; j++){
    for(int i=0; i<dimen_x_; i++){
      //grad^2 T using finite differences for spatial gradient, 2 or 3 sided differences at corners/edges of grid
      if(i>0 && j>0 && i<dimen_x_-1 && j<dimen_y_-1){
        f= (grid_[j][i-1]+grid_[j][i+1]+grid_[j-1][i]+grid_[j+1][i]-4*grid_[j][i])/(dx_*dx_);
      } else if(i==0 and j==0){ //should I take periodicity into account?
        f= (grid_[j][i+1]+grid_[j+1][i]-2*grid_[j][i])/(dx_*dx_);
      } else if(i==0 and j==dimen_y_-1){
        f= (grid_[j][i+1]+grid_[j-1][i]-2*grid_[j][i])/(dx_*dx_);
      } else if(i==dimen_x_-1 and j==0){
        f= (grid_[j][i-1]+grid_[j+1][i]-2*grid_[j][i])/(dx_*dx_);
      } else if (i==dimen_x_-1 and j==dimen_y_-1){
        f= (grid_[j][i-1]+grid_[j-1][i]-2*grid_[j][i])/(dx_*dx_);
      } else if(i==0){
        f= (grid_[j][i+1]+grid_[j-1][i]+grid_[j+1][i]-3*grid_[j][i])/(dx_*dx_);
      } else if(i==dimen_x_-1){
        f= (grid_[j][i-1]+grid_[j-1][i]+grid_[j+1][i]-3*grid_[j][i])/(dx_*dx_);
      } else if(j==0){
        f= (grid_[j][i-1]+grid_[j][i+1]+grid_[j+1][i]-3*grid_[j][i])/(dx_*dx_);
      } else if(j==dimen_y_-1){
        f= (grid_[j][i-1]+grid_[j][i+1]+grid_[j-1][i]-3*grid_[j][i])/(dx_*dx_);
      }
      grad_sq_T.Set(i,j,f);
    }
  }
  return 0;
}


int Grid::WriteToFile(char *fname) const{
  FILE *f_out=fopen(fname,"w");
  for(int j=0; j<dimen_y_; j++){
    for(int i=0; i<dimen_x_; i++){
      fprintf(f_out, "%15.8f", grid_[j][i]);
    }
    fprintf(f_out, "\n");
  }
  return 0;
}
