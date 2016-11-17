#ifndef MODEL_H_
#define MODEL_H_
class Grid;
class Model {
 public:
  virtual ~Model() {}

  // fx = f(x,t)
  virtual int rhs(double t, const Grid &x, Grid &fx) const = 0; //x should be const, figure out how

  // number of states (size of x)
  virtual int dimen_() const = 0;
  //if 2D num states in second dimension
  virtual int dimen_y_() const = 0;
};

#endif  // MODEL_H_
