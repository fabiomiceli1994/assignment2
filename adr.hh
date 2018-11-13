#ifndef CLASS_ADR
#define ADR

#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>    // std::find

#include "sparsematrix.hh"

class Adr
{
public:
  Adr(); //default construct
  Adr( unsigned int const M, double const a, double const b, double const c, double l, double b_0, double b_L ); //sets value of private data
  Adr( const Adr& source ); //copy constructor
  ~Adr(); //destructor

  unsigned int getJ () const; //returns the number of points on the interior of the set
  double getAlpha () const; //returns diffusion parameter
  double getBeta () const; //returns advection parameter
  double getGamma () const; //returns reaction parameter
  double getL () const; //returns the length of the interval
  double getu0 () const; //returns the boundary value at 0
  double getuL () const; //returns the boundary value at L
  SparseMatrix MatrixBuild (); //builds the matrix for the complete ADR equation
  SparseMatrix LaplacianMatrixBuild (); //builds the matrix for the Laplacian
  SparseMatrix GradientMatrixBuild (); //builds the matrix for the gradient
  std::vector<double> Ad_sol (); //returns the solution for advection diffusion each point of the grid
  std::vector<double> Dr_sol (); //returns the solution for diffusion reaction each point of the grid

private:
  unsigned int J_; //number of points in the interior of the set. Boundary are excluded
  double alpha_; //diffusion parameter
  double beta_; //advection parameter
  double gamma_; //reaction parameter
  double L_; //length of the interval
  double u0_; //boundary condition at 0
  double uL_; //boundary condition at L
};


#endif
