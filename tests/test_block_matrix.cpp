#include <iostream>
#include <list>
#include <string>
#include <algebra/block_matrix.h>
#include <algebra/matrix.h>
#include <algebra/vector.h>

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
  cout << "- a BlockMatrix(1,2) M:" << endl << M;

  Matrix<double>* Mblock = new Matrix<double>(2, 3, "1 2 3 4 5 6");
  M.set_block(0, 0, Mblock);
  M.resize_block_row(0, Mblock->row_dimension());
  M.resize_block_column(0, Mblock->column_dimension());
  Mblock = new Matrix<double>(2, 2, "7 8 9 10");
  M.set_block(0, 1, Mblock);
  M.resize_block_column(1, Mblock->column_dimension());
  cout << "- the same block matrix M after 2 * set_block(0,*,*):" << endl << M;

  Vector<double> x(5, "1 2 3 4 5"), y(2);
  M.apply(x, y);
  cout << "- a vector x=" << x << ", Mx=" << y << endl;

  x.resize(2);
  x[0] = 1;
  x[1] = 2;
  y.resize(5);
  M.apply_transposed(x, y);
  cout << "- a vector x=" << x << ", (M^T)x=" << y << endl;

  return 0;
}
