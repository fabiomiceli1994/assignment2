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
  std::vector<unsigned int> Jvec = {10, 50, 100};


  //A.printMatrix();

  ErrorAnalysis (Jvec, 1., 1., 0., 1., 0., 1., 1e-6, 10000, 1e9, "Error.dat");



  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  std::cout<< "Elapsed time: " << elapsed.count() << "s\n";



  return 0;
  //std::cout << "gamma=" << a.getGamma() << std::endl;
}
