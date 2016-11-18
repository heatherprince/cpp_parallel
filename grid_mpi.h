#ifndef GRID_H_
#define GRID_H_


//2D grid
class Grid{
 public:
  Grid(int nside_x, int nside_y, double x_min, double y_min,  double x_max, double y_max);
  ~Grid();
  double Get(int i, int j) const;                  //returns value at x index i, y index j
  int GetYColumn(int i, double *val)const;   //sets val to point to y column at x index i
  int Set(int i, int j, double val);               //sets value at x index i, y index j to val
  int SetYColumn(int i, double *val);  //sets y column at x index i to val
  int InitializeTEdges();                          //initializes y=0 to cos^2(x), y=pi to sin^2(x)
  int MultiplyByConstant(double c);
  int AddGridTimesConst(double c, Grid &grid2);
  int GradSq(Grid &grad_sq_T) const;
  double GetMean() const;
  double GetMeanExcludeBorders() const;
  int WriteToFile(char *fname) const;
 private:
  const int dimen_x_;                     // dimension of x axis
  const int dimen_y_;                     // dimension of y axis
  const double min_x_;                    // minimum on x axis
  const double min_y_;                    // minimum on y axis
  const double max_x_;                    // maximum on x axis
  const double max_y_;                    // maximum on y axis
  const double dx_;                    // step on x axis
  const double dy_;                    // step on x axis
  double **grid_;                         // grid=2D array of doubles - public for testing
};


#endif  // EULER_H
