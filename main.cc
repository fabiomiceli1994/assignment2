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

  Adr a = Adr (10, 1., 1., 1., 1., 1., 1.);

  std::cout << "J=" << a.getJ() << std::endl;
  std::cout << "alpha=" << a.getAlpha() << std::endl;
  std::cout << "beta=" << a.getBeta() << std::endl;
  std::cout << "gamma=" << a.getGamma() << std::endl;
  std::cout << "L=" << a.getL() << std::endl;
  std::cout << "U0=" << a.getAlpha() << std::endl;
  std::cout << "UL=" << a.getBeta() << std::endl;

  SparseMatrix A = a.MatrixBuild();

  A.printMatrix();



  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  std::cout<< "Elapsed time: " << elapsed.count() << "s\n";



  return 0;
  //std::cout << "gamma=" << a.getGamma() << std::endl;
}
