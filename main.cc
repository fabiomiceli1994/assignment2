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

  std::vector<unsigned int> Jvec = {9, 19, 39, 79, 159, 319, 639}; //different values of J for testing
  std::vector<double> alpha = {1., 1000., 5., 0.25, 1.}; //values of alpha for testing
  std::vector<double> beta = {0., 1., 0.1, 0.5, 21.}; // values of beta for testing
  std::vector<double> gamma = {0., 0., 0., 0., 0.}; //gamma is fixed
  double L = 1; //size of the domain
  double u0 = 0; //BC at 0
  double uL = 1.; // BC at L

  for(unsigned int i=0; i<alpha.size(); ++i)
  {
      ADR_Test (Jvec, alpha[i], beta[i], gamma[i], L, u0, uL, 1e-6, 10000, 1e9, "Errors");
  }

  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  std::cout<< "Elapsed time: " << elapsed.count() << "s\n";



  return 0;
}
