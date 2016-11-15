#ifndef GRID_H_
#define GRID_H_

//2D grid
class Grid{
 public:
  Grid(int nside_x, int nside_y, double x_max, double y_max);
  ~Grid();
  double Get(int i, int j) const;
  int Set(int i, int j, double val);
  int InitializeTEdges();
  int MultiplyByConstant(double c);
  int AddGridTimesConst(double c, Grid &grid2);
  int GradSq(Grid &grad_sq_T) const;
  double GetMean() const;
  int WriteToFile(char *fname) const;
 private:
  const int dimen_x_;                     // dimension of x axis
  const int dimen_y_;                     // dimension of y axis
  const double max_x_;                    // maximum on x axis
  const double max_y_;                    // maximum on x axis
  const double dx_;                    // step on x axis
  const double dy_;                    // step on x axis
  double **grid_;                         // grid=2D array of doubles - public for testing
};


#endif  // EULER_H
