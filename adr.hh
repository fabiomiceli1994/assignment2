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
  Adr( unsigned int const M, double const a, double const b, double const c ); //sets value of private data
  Adr( const Adr& source ); //copy constructor
  ~Adr(); //destructor

  unsigned int getN () const; //returns the number of points on the interior of the set
  double getAlpha () const; //returns diffusion parameter
  double getBeta () const; //returns advection parameter
  double getGamma () const; //returns response parameter




private:
  SparseMatrix A_; //class of matrix needed to implement the solver
  unsigned int N_; //number of points in the interior of the set. Boundary are excluded
  double alpha_; //diffusion parameter
  double beta_; //advection parameter
  double gamma_; //reaction parameter

};

#endif
