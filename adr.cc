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
  double h = L_/(J_+1); //mesh size
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

//returns the analitical solution of the equation on the grid
std::vector<double> An_sol ( unsigned int J, double const alpha, double const beta, double const gamma, double const L )
{
  std::vector<double> sol ( J );
  if(alpha!=0 && gamma==0) //checking non degeneracy and type of equation selected
  {
    double c = beta/alpha; //exponent
    double h = L/(J+1); //mesh size
    for(unsigned int i=0; i<J; ++i)
    {
      sol[i] = (1.-exp(c*h*(i+1)))/(1.-exp(c*L));
    }
  }else if(alpha!=0 && beta==0) //checking non degeneracy and type of equation selected
  {
    double c = sqrt(gamma/alpha); //exponent
    double h = L/(J+1); //mesh size
    for(unsigned int i=0; i<J; ++i)
    {
      sol[i] = (sinh(c*h*(i+1)))/(sinh(c*L));
    }
  }else
  {
    std::cout << "Error. The values of alpha, beta and gamma provided neither define an advection-diffusion nor a diffusion reaction one." << std::endl;
    exit(EXIT_FAILURE);
  }
  return sol;
}


//solving the ADR equation
void Solver ( std::vector<double>& u_x, unsigned int const J, double const alpha, double const beta, double const gamma, double const L, double const u0, double const uL, double const tol, double const itCheck, int const MaxIter)
{
  if( u_x.size() != J ) //checks consistency
  {
    std::cout << "Error. Size of vector solution and number of points J do not match." << std::endl;
    exit(EXIT_FAILURE);
  }
  Adr eq = Adr(J, alpha, beta, gamma, L, u0, uL);
  double h=L/(J+1); //mesh size
  double P_number = ( fabs(beta)*L )/(2*alpha); //P number
  double D_number = (gamma/alpha); //D number

  SparseMatrix A = eq.MatrixBuild(); //building the matrix to invert
  std::vector<double> b(J); //vector to invert against
  b[0]=u0*( ( alpha/(h*h) ) + ( beta/(2*h) ) ); //first boundary condition
  b[J-1]=uL*( ( alpha/(h*h) ) - ( beta/(2*h) ) ); //second boundary condition


  if(alpha!=0 && gamma==0) //checking non degeneracy and which type of equation has been chosen
  {

    std::string Filename = "P_numb_" + std::to_string (P_number) + "_J_" + std::to_string (J) + "_Residual_"; //file for residual analysis
    std::string Filename2 = "P_numb_" + std::to_string (P_number) + "_J_" + std::to_string (J) + "_Solution_"; //file for solution printing
    A.Gauss_Seidel(u_x, b, tol, itCheck, Filename, Filename2, MaxIter); //solving
    //outputting the parameters on the terminal every time a simulation is concluded
    std::cout << "Simulation for j=" << J << ", alpha=" << alpha << ", beta=" << beta << ", gamma=" << gamma << ", Pe=" << P_number << ", performed." << std::endl;

  }else if(alpha!=0 && beta==0)//checking non degeneracy and which type of equation has been chosen

  {
    std::string Filename = "D_numb_" + std::to_string (D_number) + "_J_" + std::to_string (J) + "_Residual_";
    std::string Filename2 = "D_numb_" + std::to_string (D_number) + "_J_" + std::to_string (J) + "_Solution_";
    A.Gauss_Seidel(u_x, b, tol, itCheck, Filename, Filename2, MaxIter);
    std::cout << "Simulation for j=" << J << ", alpha=" << alpha << ", beta=" << beta << ", gamma=" << gamma << ", Da=" << D_number << ", performed." << std::endl;

  }else //if the values inserted do not define AD or DR equation, exit
  {
    std::cout << "Error. The values of alpha, beta and gamma provided neither define an advection-diffusion nor a diffusion reaction one." << std::endl;
    exit(EXIT_FAILURE);
  }
}

//comparing with analytical solution
void ADR_Test (std::vector<unsigned int>& Jvec, double const alpha, double const beta, double const gamma, double const L, double const u0, double const uL, double const tol, double const itCheck, int const MaxIter, std::string ErrorFileName )
{

  double P_number = ( fabs(beta)*L )/(2*alpha); //P number
  double D_number = (gamma/alpha); //D number
  std::string FileVarName;
  if(alpha!=0 && gamma==0) //calling the outfile differently wrt which of equations has been chosen
  {
    FileVarName = "_P_numb_" + std::to_string(P_number) + ".txt";
  }else if (alpha!=0 && beta==0)
  {
    FileVarName = "_D_numb_" + std::to_string(D_number) + ".txt";
  }

  std::ofstream myOutFile (ErrorFileName + FileVarName ); //opening the file
  if ( !myOutFile.good() )
  {
    std::cout << "Failed to open the file." <<std::endl;
  }
  //making the file human readable
  myOutFile << "#Details of numerical method ADR equation through Gauss-Seidel algorithm for alpha=" << alpha << ", beta=" << beta << ", gamma=" << gamma << "." <<std::endl;
  myOutFile.width(25);
  myOutFile << std::left << "# 1-Number of points" ;
  myOutFile.width(25);
  myOutFile << std::left << "2-Error (LinfNorm)" ;
  myOutFile.width(25);
  myOutFile << std::left << "3-Scaling order of covergence" << std::endl;

  for(int j: Jvec) //lopping over the chosen values of j
  {
    double h = L/(j+1);
    std::vector<double> u_x(j);
    Solver(u_x, j, alpha, beta, gamma, L, u0, uL, tol, itCheck, MaxIter);
    double scaling = 12/(h*h)*(fabs(1.-exp(beta*L/alpha)))/(exp(beta*j*h/alpha))*(alpha/beta)*(alpha/beta)*(alpha/beta)*(alpha/beta);
    double error = LinfNorm( vectorSub( u_x, An_sol(j, alpha, beta, gamma, L) ) ); //Linf is defined in SparseMatrix
    myOutFile.width(25);
    myOutFile << std::left << j;
    myOutFile.width(25);
    myOutFile << std::left << error;
    myOutFile.width(25);
    myOutFile << std::left << error*scaling << std::endl;
  }

  myOutFile.close();
}
