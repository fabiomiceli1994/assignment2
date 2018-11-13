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

int main()

{
  auto start = std::chrono::high_resolution_clock::now();

  Adr a = Adr (10, 1., 1., 0., 1., 0., 1.);
  SparseMatrix A = a.MatrixBuild();
  std::vector<unsigned int> Jvec = {10, 50, 100, 200, 300, 400, 500};
  std::vector<double> alpha = {1000, 0.25, 1.};
  std::vector<double> beta = {1., 0.5, 21};
  double gamma = 0;
  double L = 1;
  double u0 = 0;
  double uL = 1.;

  for(unsigned int i=0; i<alpha.size(); ++i)
  {
      ErrorAnalysis (Jvec, alpha[i], beta[i], gamma, L, u0, uL, 1e-6, 10000, 1e9, "Errors");
  }

  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  std::cout<< "Elapsed time: " << elapsed.count() << "s\n";



  return 0;
}
