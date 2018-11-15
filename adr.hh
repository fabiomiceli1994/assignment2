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


private:
  unsigned int J_; //number of points in the interior of the set. Boundary are excluded
  double alpha_; //diffusion parameter
  double beta_; //advection parameter
  double gamma_; //reaction parameter
  double L_; //length of the interval
  double u0_; //boundary condition at 0
  double uL_; //boundary condition at L
};

std::vector<double> An_sol (unsigned int J, double const alpha, double const beta, double const gamma, double const L); //returns the analitical solution of the equation on the grid
std::vector<double> Solver ( unsigned int const J, double const alpha, double const beta, double const gamma, double const L, double const u0, double const uL, double const tol, double const itCheck, int const MaxIter); //solves the ADR equation
void ADR_Test (std::vector<unsigned int>& Jvec, double const alpha, double const beta, double const gamma, double const L, double const u0, double const uL, double const tol, double const itCheck, int const MaxIter, std::string ErrorFileName ); //error

#endif
