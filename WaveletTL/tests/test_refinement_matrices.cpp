/*
 * Test refinement.cpp
 * Generate some refinement matrices. Then use compose_matrix (and check the output).
 * 
 
 */


#include <iostream>
#include <map>
#include <time.h>

#include <algebra/symmetric_matrix.h>
#include <algebra/sparse_matrix.h>
#include <interval/i_index.h>
#include <interval/i_indexplot.h>
#include <interval/i_q_index.h>
#include <interval/i_q_indexplot.h>
#include <interval/ds_basis.h>
#include <interval/p_basis.h>
#include <interval/pq_frame.h>
#include <geometry/sampled_mapping.h>

#include <sstream>
#include <string>







using namespace std;
using namespace MathTL;
using namespace WaveletTL;




int main()
{
    const int d  = 3;
    const int dT = 3;
    bool bc = false;
    typedef PQFrame<d,dT> Frame;
    Frame frame(bc, bc);
    typedef Frame::Index Index;
    SparseMatrix<double> Mj0, Mj1;
    Mj0=frame.get_Mj0();
    Mj1=frame.get_Mj1();
    cout << "Mj0: " << endl << Mj0 << endl;
    cout << "Mj1: " << endl << Mj1 << endl;
    
    
    if(bc){
        const char* Mj0_filename = "refinement_matrices/bc/Mj0_3";
        const char* Mj1_filename = "refinement_matrices/bc/Mj1_3_3";
        cout << "storing matrix to file " << Mj0_filename << endl;
        Mj0.matlab_output(Mj0_filename,"Mj0",0);
        cout << "storing matrix to file " << Mj1_filename << endl;
        Mj1.matlab_output(Mj1_filename,"Mj1",0);
    }
    else{
        const char* Mj0_filename = "refinement_matrices/nobc/Mj0_3";  
        const char* Mj1_filename = "refinement_matrices/nobc/Mj1_3_3";
        cout << "storing matrix to file " << Mj0_filename << endl;
        Mj0.matlab_output(Mj0_filename,"Mj0",0);
        cout << "storing matrix to file " << Mj1_filename << endl;
        Mj1.matlab_output(Mj1_filename,"Mj1",0);
    }
    

}
