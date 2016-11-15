#ifndef MODEL_H_
#define MODEL_H_
class Grid;
class Model {
 public:
  virtual ~Model() {}

  // fx = f(x,t)
  virtual int rhs(double t, const Grid &x, Grid &fx) const = 0; //x should be const, figure out how

  // number of states (size of x)
  virtual int dimen() const = 0;
};

#endif  // MODEL_H_
