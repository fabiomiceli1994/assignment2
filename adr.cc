#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <chrono>

//#include "sparsematrix.hh"
#include "adr.hh"

//default constructor
Adr::Adr ()
{
  A_ = SparseMatrix (102, 102); // I initalise it to take into account also the boundary points, so N+2
  N_ = 100; //A hundred points by default
  alpha_ = 1.; //by default the parameters are set to one
  beta_ = 1.;
  gamma_ = 1.;
}

//sets the value of the private data
Adr::Adr ( unsigned int const M, double const a, double const b, double const c ) // M being the # of points in the interior of the interval
{
  if( M>0 )
  {
    A_ = SparseMatrix (M+2,M+2); // I initalise it to take into account also the boundary points, so M+2
    N_ = M;
    alpha_ = a;
    beta_ = b;
    gamma_ = c;
  }else
  {
    std::cout << "Error. The number of points inserted by input is not strictly positive." << std::endl;
    exit(EXIT_FAILURE);
  }
}

// Copy constructor
Adr::Adr(const Adr& source )
{
  A_ = source.A_;
  N_ = source.N_;
  alpha_ = source.alpha_;
  beta_ = source.beta_;
  gamma_ = source.gamma_;
}

//destructor
Adr::~Adr()
{
  A_.~SparseMatrix(); //destroys the SparseMatrix
}

//returns the number of points in the interior of the set
unsigned int Adr::getN () const
{
  return N_;
}

//returns the diffusion parameter
double Adr::getAlpha () const
{
  return alpha_;
}

//returns the advection parameter
double Adr::getBeta () const
{
  return beta_;
}


//returns the response parameter
double Adr::getGamma () const
{
  return gamma_;
}
