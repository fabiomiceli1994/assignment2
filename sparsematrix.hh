#ifndef CLASS_SPARSEMATRIX
#define SPARSEMATRIX

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

class SparseMatrix
{
public:
  SparseMatrix(); //default construct
  SparseMatrix( int rowSize, int colSize ); //sets value of private data
  SparseMatrix( const SparseMatrix& source ); //copy constructor
  ~SparseMatrix(); //destructor

  //operators within the class
  std::vector<double> operator*( const std::vector<double>& input ) const;

  unsigned int getRowSize () const; //gets the number of rows
  unsigned int getColSize () const; //gets the number of columns
  void addEntry ( unsigned int rowNumb, unsigned int colNumb, double newValue); //adds an entry to the matrix in the location [rowNumb][colNumb]
  double getValue (unsigned int x, unsigned int y) const; //return the element (x, y) of the real matrix
  void printMatrix (); //prints the matrix
  void Gauss_Seidel ( std::vector<double>& x_0, const std::vector<double>& b, const double tol, const int itCheck, std::string fileName, const int MaxIter ); //Gauss_Seidel


private:
  unsigned int rowSize_; //number of rows of the matrix
  unsigned int colSize_; //number of columns of the matrix
  std::vector<std::vector<double>*>* rows_; //vector of vector: first vector is matrix, its slots are the rows, whose slots are the entries
  std::vector<std::vector<int>*>* colsInd_; //vector of vector: first vector is matrix, its slots are rows, whose entries are the column indexes

};

std::vector<double> vectorSum ( std::vector<double> v1, std::vector<double> v2 ); //sums two vectors
std::vector<double> vectorSub ( std::vector<double> v1, std::vector<double> v2 ); //subtracts two vectors
double LinfNorm ( std::vector<double> v ); //returns the maximum of a vector. In this case I use it to construct the LinfNorm


#endif
