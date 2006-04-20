#include <iostream>
#include <list>
#include <string>
#include <algebra/block_matrix.h>
#include <algebra/matrix.h>

using std::cout;
using std::endl;
using namespace MathTL;

int main()
{
  cout << "Testing the class BlockMatrix ..." << endl;

  BlockMatrix<double> defaultM;
  cout << "- a default block matrix:" << endl << defaultM;

  BlockMatrix<double> N(1);
  cout << "- a BlockMatrix(1):" << endl << N;

  BlockMatrix<double> M(1,2);
  cout << "- a BlockMatrix(1,2):" << endl << M;

  Matrix<double>* Mblock = new Matrix<double>(2, 3);
  M.set_block(0, 0, Mblock);
  M.resize_block_row(0, Mblock->row_dimension());
  M.resize_block_column(0, Mblock->column_dimension());
  cout << "- the same block matrix after set_block(0,0,*):" << endl << M;

  return 0;
}
