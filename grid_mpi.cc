#include <stdio.h>
#include <math.h>
#include "grid_mpi.h"

//implementation of Euler that works on a square grid
Grid::Grid(int nside_x, int nside_y, double x_min, double y_min, double x_max, double y_max)
    : dimen_x_(nside_x),
      dimen_y_(nside_y),
      min_x_(x_min),
      min_y_(y_min),
      max_x_(x_max),
      max_y_(y_max),
      dx_((x_max-x_min)/(nside_x-1.)),
      dy_((y_max-y_min)/(nside_y-1.)) {
  grid_ = new double*[dimen_x_];
  grid_[0] = new double[dimen_x_*dimen_y_];
  for (int i = 1; i < dimen_x_; i++){
    grid_[i] = grid_[i-1] + dimen_y_;
  }

  //initialize all by default, initialize edges for T in a separate function
  for(int i=0; i<dimen_x_; i++){
    for(int j=0; j<dimen_y_; j++){
      grid_[i][j]=0.;
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
  return grid_[i][j];
}

int Grid::SetYColumn(int i, double *val){
  //check i, i in correct range
  grid_[i]=val;
  return 0;
}

//i is x index, j is y index
int Grid::GetYColumn(int i, double *val)const {  //sets val to column
  //check i, i in correct range
  val=grid_[i];
  return 0;
}

int Grid::Set(int i, int j, double val){
  //check i, i in correct range
  grid_[i][j]=val;
  return 0;
}

double Grid::GetMean()const {
  //check i, i in correct range
  //calculate mean temp
  double sum=0;
  for(int i=0; i<dimen_x_; i++){
    for(int j=0; j<dimen_y_; j++){
      sum+=grid_[i][j];
    }
  }
  double mean=sum/(dimen_x_*dimen_y_);
  return mean;
}

int Grid::InitializeTEdges(){
  for(int i=0; i<dimen_x_; i++){
    double x=i*dx_+min_x_;
    grid_[i][0]=cos(x)*cos(x);
    grid_[i][dimen_y_-1]=sin(x)*sin(x);
  }
  return 0;
}

int Grid::MultiplyByConstant(double c){
  for(int i=0; i<dimen_x_; i++){
    for(int j=0; j<dimen_y_; j++){
      grid_[i][j]*=c;
    }
  }
  return 0;
}

int Grid::AddGridTimesConst(double c, Grid &grid2){
  for(int i=0; i<dimen_x_; i++){
    for(int j=0; j<dimen_y_; j++){
      grid_[i][j]+=c*grid2.Get(i,j);
    }
  }
  return 0;
}

int Grid::GradSq(Grid &grad_sq_T)const  {
  double f;
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
