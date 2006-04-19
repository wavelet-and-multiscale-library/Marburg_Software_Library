#include <iostream>
#include <list>
#include <string>
#include <algebra/block_matrix.h>

using std::cout;
using std::endl;
using namespace MathTL;

int main()
{
  cout << "Testing the class BlockMatrix ..." << endl;

  BlockMatrix<double> defaultM;
  cout << "- a default block matrix:" << endl << defaultM;

  BlockMatrix<double> M(1);
  cout << "- a BlockMatrix(1):" << endl << M;

  return 0;
}
