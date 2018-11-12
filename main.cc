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

  Adr a = Adr (100, 1, 1, 1);

  // std::cout << "N=" << a.getN() << std::endl;
  // std::cout << "alpha=" << a.getAlpha() << std::endl;
  // std::cout << "beta=" << a.getBeta() << std::endl;
  // std::cout << "gamma=" << a.getGamma() << std::endl;

  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  std::cout<< "Elapsed time: " << elapsed.count() << "s\n";



  return 0;
  //std::cout << "gamma=" << a.getGamma() << std::endl;
}
