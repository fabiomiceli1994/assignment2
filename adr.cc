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
  u0_ =  0.; //BCs given in the assignment
  uL_ = 1.;
}

//sets the value of the private data
Adr::Adr ( unsigned int const M, double const a, double const b, double const c, double l, double b_0, double b_L ) // M being the # of points in the interior of the interval
{
  if( M>0 ) //number of points in the interior of the interval has to be strictly positive
  {
    J_ = M;
    alpha_ = a;
    beta_ = b;
    gamma_ = c;
    L_ = l;
    u0_ = b_0;
    uL_ = b_L;
  }else //exits with warning
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
  SparseMatrix A = SparseMatrix( J_, J_); //constructor of SparseMatrix
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
std::vector<double> Adr::An_sol ( )
{
  std::vector<double> sol ( J_ );
  if(alpha_!=0 && beta_!=0 && gamma_==0) //checking non degeneracy and type of equation selected
  {
    double c = beta_/alpha_; //exponent
    double h = L_/(J_+1); //mesh size
    for(unsigned int i=0; i<J_; ++i) //assigning values in the interval
    {
      sol[i] = (1.-exp(c*h*(i+1)))/(1.-exp(c*L_));
    }
  }else if(alpha_!=0 && beta_==0 && gamma_!=0) //checking non degeneracy and type of equation selected
  {
    double c = sqrt(gamma_/alpha_); //exponent
    double h = L_/(J_+1); //mesh size
    for(unsigned int i=0; i<J_; ++i)
    {
      sol[i] = (sinh(c*h*(i+1)))/(sinh(c*L_));
    }
  }else if( alpha_!=0 && beta_==0 && gamma_==0 ) //straight line case
  {
    for(unsigned int i=0; i<J_; ++i)
    {
      double h = L_/(J_+1); //mesh size
      sol[i] = (i+1)*h;
    }
  }else
  {
    std::cout << "Error. The values of alpha, beta and gamma provided neither define an advection-diffusion nor a diffusion reaction one." << std::endl;
    exit(EXIT_FAILURE);
  }

  //inserting the BCs
  sol.insert(sol.begin(), u0_);
  sol.push_back(uL_);
  return sol;
}


//solving the ADR equation
void Adr::Solver ( std::vector<double>& u_x, double const tol, double const itCheck, int const MaxIter)
{
  if( u_x.size() != J_ ) //checks consistency
  {
    std::cout << "Error. Size of vector solution and number of points J do not match." << std::endl;
    exit(EXIT_FAILURE);
  }
  double h=L_/(J_+1); //mesh size
  double P_number = ( fabs(beta_)*L_ )/(2*alpha_); //P number
  double D_number = (gamma_/alpha_); //D number

  SparseMatrix A = (*this).MatrixBuild(); //building the matrix to invert
  std::vector<double> f(J_); //vector to invert against
  f[0] = u0_*( ( alpha_/(h*h) ) + ( beta_/(2*h) ) ); //first boundary condition
  f[J_-1] = uL_*( ( alpha_/(h*h) ) - ( beta_/(2*h) ) ); //second boundary condition

  std::string Filename, Filename2; //needed to have a variable output if the parameters vary

  if(alpha_!=0 && gamma_==0) //checking non degeneracy and which type of equation has been chosen
  {

    Filename = "P_numb_" + std::to_string (P_number) + "_Residual_" + "J_" ;
    Filename2 = "P_numb_" + std::to_string (P_number) + "_Solution_" + "J_" + std::to_string(J_) + ".txt";
    A.Gauss_Seidel(u_x, f, tol, itCheck, Filename, MaxIter); //solving
    //outputting the parameters on the terminal every time a simulation is concluded
    std::cout << "Solved for J=" << J_ << ", alpha=" << alpha_ << ", beta=" << beta_ << ", gamma=" << gamma_ << ", Pe=" << P_number << "." << std::endl;

  }else if(alpha_!=0 && beta_==0)//checking non degeneracy and which type of equation has been chosen

  {
    Filename = "D_numb_" + std::to_string (D_number) + "_Residual_" + "J_" ;
    Filename2 = "D_numb_" + std::to_string (D_number) + "_Solution_" + "J_" + std::to_string(J_) + ".txt";
    A.Gauss_Seidel(u_x, f, tol, itCheck, Filename, MaxIter);
    std::cout << "Solved for J=" << J_ << ", alpha=" << alpha_ << ", beta=" << beta_ << ", gamma=" << gamma_ << ", Da=" << D_number << "." << std::endl;
  }else //if the values inserted do not define AD or DR equation, exit
  {
    std::cout << "Error. The values of alpha, beta and gamma provided neither define an advection-diffusion nor a diffusion reaction one." << std::endl;
    exit(EXIT_FAILURE);
  }

  std::ofstream myFile (Filename2); //opening the file
  if ( !myFile.good() )
  {
    std::cout << "Failed to open the file." << std::endl;
  }

  myFile << "#Details of numerical method ADR equation through Gauss-Seidel algorithm for alpha=" << alpha_ << ", beta=" << beta_ << ", gamma=" << gamma_ << "." <<std::endl;
  myFile.width(25);
  myFile << std::left << "#1-x_j/L" ;
  myFile.width(25);
  myFile << std::left << "2-Value" << std::endl;

  //Gauss-Seidel solves in the interior of the interval. I put the boundary conditions in now
  u_x.insert(u_x.begin(), u0_);
  u_x.push_back(uL_);

  //printing the solution at each point of the grid
  for(unsigned int i=0; i<u_x.size(); ++i)
  {
    myFile.width(25);
    myFile << std::left << ((double)i)*L_/(J_+1) ;
    myFile.width(25);
    myFile << std::left << u_x[i] << std::endl;
  }

  myFile.close();

}

//comparing with analytical solution
void ADR_Test (std::vector<unsigned int>& Jvec, double const alpha, double const beta, double const gamma, double const L, double const u0, double const uL, double const tol, double const itCheck, int const MaxIter, std::string ErrorFileName )
{

  double P_number = ( fabs(beta)*L )/(2*alpha); //P number
  double D_number = (gamma/alpha); //D number
  double error, err_ratio, err_old=1; // variables to compute errors ratios.Err_old is set to 1 by default since at the first iterations no ratio will be printed.
  std::string FileVarName;
  if(alpha!=0 && gamma==0) //calling the outfile differently wrt which of equations has been chosen
  {
    FileVarName = "_P_numb_" + std::to_string(P_number) + ".txt";
  }else if (alpha!=0 && beta==0)
  {
    FileVarName = "_D_numb_" + std::to_string(D_number) + ".txt";
  }

  std::ofstream myOutFile (ErrorFileName + FileVarName); //opening the file
  if ( !myOutFile.good() )
  {
    std::cout << "Failed to open the file." <<std::endl;
  }
  //making the file human readable
  myOutFile << "#Details of numerical method for ADR equation through Gauss-Seidel algorithm for alpha=" << alpha << ", beta=" << beta << ", gamma=" << gamma << "." <<std::endl;
  myOutFile.width(25);
  myOutFile << std::left << "#1-Mesh size" ;
  myOutFile.width(25);
  myOutFile << std::left << "2-Error (LinfNorm)" ;
  myOutFile.width(25);
  myOutFile << std::left << "3-log2(Error Ratios)" << std::endl;
  for(unsigned int j=0; j<Jvec.size(); ++j) //looping over the different number of points chosen
  {
    Adr eq = Adr(Jvec[j], alpha, beta, gamma, L, u0, uL);
    std::vector<double> u_x(Jvec[j]); //initial guess
    eq.Solver(u_x, tol, itCheck, MaxIter); //solves the equation
    error = LinfNorm( vectorSub( u_x, eq.An_sol( ) ) ); //Linf is defined in SparseMatrix
    err_ratio = err_old/error;
    err_old = error;
    myOutFile.width(25);
    myOutFile << std::left << 1./(Jvec[j]+1);
    myOutFile.width(25);
    myOutFile << std::left << error;
    myOutFile.width(25);
    if( j == 0 || error == 0 ) //skips the first call since no ratio can be calculated and avoids printing "nan" if error is qual to zero
    {
      myOutFile << std::left << "\t" << std::endl;
    }else
    {
      myOutFile << std::left << log2(err_ratio) << std::endl;
    }
  }

  myOutFile.close();
}
