/*
 * Test precompute_stiffness_matrices.cpp
 * Generate some matrices. Then use compose_matrix (and check the output).
 * 
 * Check comment in precompute_stiffness_matrices.
 * I stored the computed matrices in the shared folder 
 * (slower IO, but useful if you intend to use several PCs in parallel later)
 */
#include <galerkin/precompute_stiffness_matrices.cpp>

#include <iostream>
#include <map>
#include <time.h>

#include <algebra/symmetric_matrix.h>
#include <algebra/sparse_matrix.h>
#include <numerics/sturm_bvp.h>
#include <numerics/iteratsolv.h>

#include <interval/i_index.h>
//#include <interval/ds_basis.h>
#include <interval/p_basis.h>
#include <adaptive/cdd1.h>
#include <adaptive/cdd2.h>
#include <geometry/sampled_mapping.h>

#include <cube/tbasis.h>
#include <galerkin/tbasis_equation.h>
#include <galerkin/cached_tproblem.h>
#include <cube/tbasis_evaluate.h>
#include <cube/tbasis_indexplot.h>



using namespace std;
using namespace MathTL;
using namespace WaveletTL;

using MathTL::SimpleSturmBVP;
using MathTL::CG;


int main()
{
    // test compose_matrix

    // store some simple matrices

    //sprintf(matrix_filename, "%s%s%s%d_%d",a_namestart,a_bcstring,a_nameend,ai,aj);
    SparseMatrix<double> A00(1,1);
    //const char* stiffnessMatrixStorageFolder = {"/import/shared/friedrich/source/precomputed/stiffness_matrices"};
    const char* A00_filename = "/import/shared/friedrich/source/precomputed/tmp/A_0_0";
    A00.set_entry(0,0,1);
    cout << "A00 = " << endl;
    A00.print(cout,2,4);
    cout << "storing matrix to file " << A00_filename << endl;
    A00.matlab_output(A00_filename,"A",1);

    SparseMatrix<double> A01(1,1);
    const char* A01_filename = "/import/shared/friedrich/source/precomputed/tmp/A_0_1";
    A01.set_entry(0,0,2);
    cout << "A01 = " << endl;
    A01.print(cout,2,4);
    cout << "storing matrix to file " << A01_filename << endl;
    A01.matlab_output(A01_filename,"A",1);

    SparseMatrix<double> A02(1,2);
    const char* A02_filename = "/import/shared/friedrich/source/precomputed/tmp/A_0_2";
    A02.set_entry(0,0,3);
    A02.set_entry(0,1,3);
    cout << "A02 = " << endl;
    A02.print(cout,2,4);
    cout << "storing matrix to file " << A02_filename << endl;
    A02.matlab_output(A02_filename,"A",1);

    SparseMatrix<double> A10(1,1);
    const char* A10_filename = "/import/shared/friedrich/source/precomputed/tmp/A_1_0";
    A10.set_entry(0,0,4);
    cout << "A10 = " << endl;
    A10.print(cout,2,4);
    cout << "storing matrix to file " << A10_filename << endl;
    A10.matlab_output(A10_filename,"A",1);

    SparseMatrix<double> A11(1,1);
    const char* A11_filename = "/import/shared/friedrich/source/precomputed/tmp/A_1_1";
    A11.set_entry(0,0,5);
    cout << "A11 = " << endl;
    A11.print(cout,2,4);
    cout << "storing matrix to file " << A11_filename << endl;
    A11.matlab_output(A11_filename,"A",1);

    SparseMatrix<double> A12(1,2);
    const char* A12_filename = "/import/shared/friedrich/source/precomputed/tmp/A_1_2";
    A12.set_entry(0,0,6);
    A12.set_entry(0,1,6);
    cout << "A12 = " << endl;
    A12.print(cout,2,4);
    cout << "storing matrix to file " << A12_filename << endl;
    A12.matlab_output(A12_filename,"A",1);

    SparseMatrix<double> A20(2,1);
    const char* A20_filename = "/import/shared/friedrich/source/precomputed/tmp/A_2_0";
    A20.set_entry(0,0,7);
    A20.set_entry(1,0,7);
    cout << "A20 = " << endl;
    A20.print(cout,2,4);
    cout << "storing matrix to file " << A20_filename << endl;
    A20.matlab_output(A20_filename,"A",1);

    SparseMatrix<double> A21(2,1);
    const char* A21_filename = "/import/shared/friedrich/source/precomputed/tmp/A_2_1";
    A21.set_entry(0,0,8);
    A21.set_entry(1,0,8);
    cout << "A21 = " << endl;
    A21.print(cout,2,4);
    cout << "storing matrix to file " << A21_filename << endl;
    A21.matlab_output(A21_filename,"A",1);

    SparseMatrix<double> A22(2,2);
    const char* A22_filename = "/import/shared/friedrich/source/precomputed/tmp/A_2_2";
    A22.set_entry(0,0,9);
    A22.set_entry(0,1,9);
    A22.set_entry(1,0,9);
    A22.set_entry(1,1,9);
    cout << "A22 = " << endl;
    A22.print(cout,2,4);
    cout << "storing matrix to file " << A22_filename << endl;
    A22.matlab_output(A22_filename,"A",1);

    SparseMatrix<double> B00(1,1);
    const char* B00_filename = "/import/shared/friedrich/source/precomputed/tmp/B_0_0";
    B00.set_entry(0,0,10);
    cout << "B00 = " << endl;
    B00.print(cout,2,4);
    cout << "storing matrix to file " << B00_filename << endl;
    B00.matlab_output(B00_filename,"A",1);

    SparseMatrix<double> B01(1,1);
    const char* B01_filename = "/import/shared/friedrich/source/precomputed/tmp/B_0_1";
    B01.set_entry(0,0,11);
    cout << "B01 = " << endl;
    B01.print(cout,2,4);
    cout << "storing matrix to file " << B01_filename << endl;
    B01.matlab_output(B01_filename,"A",1);

    SparseMatrix<double> B02(1,2);
    const char* B02_filename = "/import/shared/friedrich/source/precomputed/tmp/B_0_2";
    B02.set_entry(0,0,12);
    B02.set_entry(0,1,12);
    cout << "B02 = " << endl;
    B02.print(cout,2,4);
    cout << "storing matrix to file " << B02_filename << endl;
    B02.matlab_output(B02_filename,"A",1);

    SparseMatrix<double> B10(1,1);
    const char* B10_filename = "/import/shared/friedrich/source/precomputed/tmp/B_1_0";
    B10.set_entry(0,0,13);
    cout << "B10 = " << endl;
    B10.print(cout,2,4);
    cout << "storing matrix to file " << B10_filename << endl;
    B10.matlab_output(B10_filename,"A",1);

    SparseMatrix<double> B11(1,1);
    const char* B11_filename = "/import/shared/friedrich/source/precomputed/tmp/B_1_1";
    B11.set_entry(0,0,14);
    cout << "B11 = " << endl;
    B11.print(cout,2,4);
    cout << "storing matrix to file " << B11_filename << endl;
    B11.matlab_output(B11_filename,"A",1);

    SparseMatrix<double> B12(1,2);
    const char* B12_filename = "/import/shared/friedrich/source/precomputed/tmp/B_1_2";
    B12.set_entry(0,0,15);
    B12.set_entry(0,1,15);
    cout << "B12 = " << endl;
    B12.print(cout,2,4);
    cout << "storing matrix to file " << B12_filename << endl;
    B12.matlab_output(B12_filename,"A",1);

    SparseMatrix<double> B20(2,1);
    const char* B20_filename = "/import/shared/friedrich/source/precomputed/tmp/B_2_0";
    B20.set_entry(0,0,16);
    B20.set_entry(1,0,16);
    cout << "B20 = " << endl;
    B20.print(cout,2,4);
    cout << "storing matrix to file " << B20_filename << endl;
    B20.matlab_output(B20_filename,"A",1);

    SparseMatrix<double> B21(2,1);
    const char* B21_filename = "/import/shared/friedrich/source/precomputed/tmp/B_2_1";
    B21.set_entry(0,0,17);
    B21.set_entry(1,0,17);
    cout << "B21 = " << endl;
    B21.print(cout,2,4);
    cout << "storing matrix to file " << B21_filename << endl;
    B21.matlab_output(B21_filename,"A",1);

    SparseMatrix<double> B22(2,2);
    const char* B22_filename = "/import/shared/friedrich/source/precomputed/tmp/B_2_2";
    B22.set_entry(0,0,18);
    B22.set_entry(0,1,18);
    B22.set_entry(1,0,18);
    B22.set_entry(1,1,18);
    cout << "B22 = " << endl;
    B22.print(cout,2,4);
    cout << "storing matrix to file " << B22_filename << endl;
    B22.matlab_output(B22_filename,"A",1);


    // Compose a matrix

    //void compose_matrix(SparseMatrix<double> & A,
    //    const char* a_namestart, const char* a_bcstring, const char* a_nameend,
    //    const char* b_namestart, const char* b_bcstring, const char* b_nameend,
    //    //const char* out_namestart, const char* out_nameend,
    //    const unsigned int offset,
    //    const unsigned int number_of_gens))
    SparseMatrix<double> C (8,8); // o=0 => s =1, o=1=>s=3, o=2=> s=8
    compose_matrix(C,"/import/shared/friedrich/source/precomputed/tmp/A_","","",
                     "/import/shared/friedrich/source/precomputed/tmp/B_","","",
                     //"/import/shared/friedrich/source/precomputed/tmp/C")
                     2,
                     8); // achtung, dieser Wert macht erstmal nicht viel sinn
    cout << "C = " << endl;
    C.print(cout,4,4);
    C.add(1,transpose(C));
    cout << "C+C' = " << endl;
    C.print(cout,4,4);

}
