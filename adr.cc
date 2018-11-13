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
  J_ = 100; //A hundred points by default
  alpha_ = 1.; //by default the parameters are set to one
  beta_ = 1.;
  gamma_ = 1.;
  L_ = 1.;
  u0_ =  0.;
  uL_ = 1.;
}

//sets the value of the private data
Adr::Adr ( unsigned int const M, double const a, double const b, double const c, double l, double b_0, double b_L ) // M being the # of points in the interior of the interval
{
  if( M>0 )
  {
    J_ = M;
    alpha_ = a;
    beta_ = b;
    gamma_ = c;
    L_ = l;
    u0_ = b_0;
    uL_ = b_L;
  }else
  {
    std::cout << "Error. The number of points inserted by input is not strictly positive." << std::endl;
    exit(EXIT_FAILURE);
  }
}

// Copy constructor
Adr::Adr(const Adr& source )
{
  J_ = source.J_;
  alpha_ = source.alpha_;
  beta_ = source.beta_;
  gamma_ = source.gamma_;
  L_ = source.L_;
  u0_ = source.u0_;
  uL_ = source.uL_;
}

//destructor
Adr::~Adr()
{

}

//returns the number of points in the interior of the set
unsigned int Adr::getJ () const
{
  return J_;
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

//returns the reaction parameter
double Adr::getGamma () const
{
  return gamma_;
}

//returns the length of the interval
double Adr::getL () const
{
  return L_;
}

//returns the boundary condition at 0
double Adr::getu0 () const
{
  return u0_;
}

//returns the boundary condition at L
double Adr::getuL () const
{
  return uL_;
}

//builds the matrix for the Laplacian
SparseMatrix Adr::MatrixBuild ()
{
  SparseMatrix A = SparseMatrix( J_, J_);
  double h = L_/J_; //mesh size
  double c_lap = -alpha_/(h*h); //multiplicative scalar for laplacian
  double c_grad = beta_/(2.*h); //multiplivative scalar for gradient
  for(unsigned int i=0; i<J_; ++i)
  {
    for(unsigned int j=0; j<J_; ++j)
    {
      if( (i-1) == j )
      {
        A.addEntry( i, j, c_lap-c_grad);
      }else if( i == j )
      {
        A.addEntry( i, j, -2*c_lap+gamma_);
      }else if( (i+1) == j )
      {
        A.addEntry( i, j, c_lap+c_grad);
      }else
      {
        continue;
      }
    }
  }
  return A;
}

//builds the matrix for the Laplacian
SparseMatrix Adr::LaplacianMatrixBuild ()
{
  SparseMatrix A = SparseMatrix( J_, J_);
  double h = L_/J_; //mesh size
  double c = -alpha_/(h*h); //multiplicative scalar
  for(unsigned int i=0; i<J_; ++i)
  {
    for(unsigned int j=0; j<J_; ++j)
    {
      if( (i-1) == j )
      {
        A.addEntry( i, j, c );
      }else if( i == j )
      {
        A.addEntry( i, j, -2*c);
      }else if( (i+1) == j )
      {
        A.addEntry( i, j, c);
      }else
      {
        continue;
      }
    }
  }
  return A;
}

//builds the matrix for the gradient
SparseMatrix Adr::GradientMatrixBuild ()
{
  SparseMatrix A = SparseMatrix( J_, J_);
  double h = L_/J_; //mesh size
  double c = beta_/(2*h); //multiplicative scalar
  for(unsigned int i=0; i<J_; ++i)
  {
    for(unsigned int j=0; j<J_; ++j)
    {
      if( (i-1) == j )
      {
        A.addEntry( i, j, c );
      }else if( (i+1) == j )
      {
        A.addEntry( i, j, -c );
      }else
      {
        continue;
      }
    }
  }
  return A;
}

//returns the solution of the equation for the advection diffusion type
std::vector<double> Adr::Ad_sol ()
{
  std::vector<double> sol ( J_+2 );
  sol[0]=u0_;
  sol[J_+1]=uL_;
  double c = beta_/alpha_; //exponent
  double h = L_/J_; //mesh size
  for(unsigned int i=1; i<J_+1; ++i)
  {
    sol[i] = (1.-exp(c*h*i))/(1.-exp(c*L_));
  }
  return sol;
}

//returns the solution of the equation for the diffusion reaction type
std::vector<double> Adr::Dr_sol ()
{
  std::vector<double> sol ( J_+2 );
  sol[0]=u0_;
  sol[J_+1]=uL_;
  double c = gamma_/alpha_; //exponent
  double h = L_/J_; //mesh size
  for(unsigned int i=1; i<J_+1; ++i)
  {
    sol[i] = (sinh(c*h*i))/(sinh(c*L_));
  }
  return sol;
}
