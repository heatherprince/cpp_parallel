#include <stdio.h>
#include <math.h>
#include "grid.h"

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
  //I want y strips to be continuous in memory so I can divide up x in domain decomposition
  //I made this choice so it would be easier to deal with boundary conditions in ring topology -> no special cases for ends
  grid_ = new double*[dimen_x_];
  grid_[0] = new double[dimen_x_*dimen_y_];
  for (int i = 1; i < dimen_x_; i++){
    grid_[i] = grid_[i-1] + dimen_y_;
  }

  //initialize all to zero by default, initialize edges for T in a separate function
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
  for(int j=0; j<dimen_y_; j++){
    grid_[i][j]=val[j];
  }
  return 0;
}

//i is x index, j is y index
int Grid::GetYColumn(int i, double *val)const {  //sets val to column
  //check i, i in correct range
  for(int j=0; j<dimen_y_; j++){
    val[j]=grid_[i][j];
  }
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

double Grid::GetMeanExcludeBorders()const {
  //check i, i in correct range
  //calculate mean temp
  double sum=0;
  for(int i=1; i<dimen_x_-1; i++){  //exclude first and last y column to get correct mean for each process
    for(int j=0; j<dimen_y_; j++){
      sum+=grid_[i][j];
    }
  }
  double mean=sum/((dimen_x_-2)*dimen_y_);
  return mean;
}

int Grid::InitializeTEdges(){
  double x;
  for(int i=0; i<dimen_x_; i++){
    x=i*dx_+min_x_;
    grid_[i][0]=cos(x)*cos(x);
    grid_[i][dimen_y_-1]=sin(x)*sin(x);
    //printf("i: %4d, dimen_x: %4d, dx: %5.2f, x_min: %5.2f, x_max: %5.2f, x: %5.2f, (x_max-x_min)/(n_x-1):%5.2f \n",i,dimen_x_,dx_, min_x_, max_x_, x, (max_x_-min_x_)/(dimen_x_-1));
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
  for(int i=0; i<dimen_x_; i++){
    for(int j=0; j<dimen_y_; j++){
      //grad^2 T using finite differences for spatial gradient, 2 or 3 sided differences at corners/edges of grid
      if(i>0 && j>0 && i<dimen_x_-1 && j<dimen_y_-1){
        f= (grid_[i-1][j]+grid_[i+1][j]+grid_[i][j-1]+grid_[i][j+1]-4*grid_[i][j])/(dx_*dx_);
      } else if(i==0 and j==0){ //should I take periodicity into account?
        f= (grid_[i+1][j]+grid_[i][j+1]-2*grid_[i][j])/(dx_*dx_);
      } else if(i==0 and j==dimen_y_-1){
        f= (grid_[i+1][j]+grid_[i][j-1]-2*grid_[i][j])/(dx_*dx_);
      } else if(i==dimen_x_-1 and j==0){
        f= (grid_[i-1][j]+grid_[i][j+1]-2*grid_[i][j])/(dx_*dx_);
      } else if (i==dimen_x_-1 and j==dimen_y_-1){
        f= (grid_[i-1][j]+grid_[i][j-1]-2*grid_[i][j])/(dx_*dx_);
      } else if(i==0){
        f= (grid_[i+1][j]+grid_[i][j-1]+grid_[i][j+1]-3*grid_[i][j])/(dx_*dx_);
      } else if(i==dimen_x_-1){
        f= (grid_[i-1][j]+grid_[i][j-1]+grid_[i][j+1]-3*grid_[i][j])/(dx_*dx_);
      } else if(j==0){
        f= (grid_[i-1][j]+grid_[i+1][j]+grid_[i][j+1]-3*grid_[i][j])/(dx_*dx_);
      } else if(j==dimen_y_-1){
        f= (grid_[i-1][j]+grid_[i+1][j]+grid_[i][j-1]-3*grid_[i][j])/(dx_*dx_);
      }
      grad_sq_T.Set(i,j,f);
    }
  }
  return 0;
}


int Grid::WriteToFile(const char *fname) const{
  FILE *f_out=fopen(fname,"w");
  for(int i=0; i<dimen_x_; i++){
    for(int j=0; j<dimen_y_; j++){
      fprintf(f_out, "%15.8f", grid_[i][j]);
    }
    fprintf(f_out, "\n");
  }
  return 0;
}
