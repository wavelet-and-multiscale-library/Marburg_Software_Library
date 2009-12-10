#include <iostream>


#include <utils/fixed_array1d.h>
#include <cube/tbasis.h>
#include <utils/multiindex.h>
#include <cube/tbasis_evaluate.h>
#include <interval/p_basis.h>
#include <interval/p_evaluate.h>
#include <interval/ds_basis.h>
#include <interval/ds_evaluate.h>
using namespace std;
using namespace WaveletTL;

using MathTL::FixedArray1D;

int main()
{
	cout << "Testing tbasis." << endl;

	//typedef DSBasis<3,5> Basis1d;
        typedef PBasis<3,5> Basis1d;
	typedef TensorBasis<Basis1d,2> Basis;
	typedef Basis::Index Index;
        typedef Basis::Support Support;


	cout << "default constructor (with homogeneous b.c.)" << endl;
        Basis basisH;
	cout << "done" << endl;

	
	cout << "Constructor with specified boundary condition orders"<<endl;
	cout << "Initialize FixedArray<int,4> with b.c.'s" << endl;

	FixedArray1D<int,4> s,st;
	s[0]=s[1]=s[2]=s[3]=0;
	st[0]=st[1]=st[2]=st[3]=0;

	cout << "Call constructor with s" << endl;
	Basis basisS(s);
	cout << "done" << endl;
/*
	cout << "Call constructor with s,st" << endl;
	Basis basisSST(s,st);
	cout << "done" << endl;
*/
	cout << "Constructor with specified Dirichlet boundary conditions for the primal functions" << endl;
	cout << "Initialize FixedArray<bool,4> with b.c.'s" << endl;
	FixedArray1D<bool,4> bc;
	bc[0]=bc[1]=bc[2]=bc[3]=true;
	cout << "Call constructor with bc" << endl;
	Basis basisBC(bc);
	cout << "done" << endl;

	cout << "Constructor with precomputed instances of the 1D bases" << endl;
	cout << "Initialize 1d bases" << endl;
	//Basis1d bas0();
	Basis1d bas11(false,false);
	bas11.set_jmax(5);
	FixedArray1D<Basis1d*,2> basesArray;
	basesArray[0] = basesArray[1] = &bas11;
	cout << "Call constructor with bas11" << endl;
	Basis basisBas(basesArray);
	cout << "done" << endl;


        // Use this basis for testing:
        //Basis basis;
        //Basis basis(s);
        //Basis basis(s,st);
        //Basis basis(bc);
        Basis basis(basesArray);


        cout << "Testing j0()" << endl
             << "basisH " << basisH.j0() << endl
             << "basisS " << basisS.j0() << endl
//             << "basisSST "<< basisSST.j0() << endl
             << "basisBC "<< basisBC.j0() << endl
             << "basisBas "<< basisBas.j0() << endl;
 
        MultiIndex<int,2u> jmax;
        jmax[0]=jmax[1]=5;
        cout << "Testing set_jmax" << endl;
        cout << "basisH ...";
        basisH.set_jmax(jmax);
        cout << "done" << endl;
       
        cout << "basisS ...";
        basisS.set_jmax(jmax);
        cout << "done" << endl;

        /*
        cout << "basisSST ...";
        basisSST.set_jmax(jmax);
        cout << "done" << endl;
        */

        cout << "basisBC ...";
        basisBC.set_jmax(jmax);
        cout << "done" << endl;
        
        cout << "basisBas ...";
        basisBas.set_jmax(jmax);
        cout << "done" << endl;

        cout << "Testing support" << endl;
        // void support(const Index& lambda, Support& supp) const;
        Support supp;
        
        //cout << "DeltaLmax()  " << basis.bases()[0]->DeltaLmax() << " " << basis.bases()[1]->DeltaLmax() << endl;
        //cout << basisH.bases()[0]->DeltaLmax() << endl;
        //cout << "basis.first_generator();"  << basis.first_generator()<< endl;
        //basisH.first_generator();
        int i=0;
        for (Index it(basis.first_generator());it < basis.last_wavelet(jmax);++it,i++)
        {
            basis.support(it,supp);
            //cout << "supp = " << supp << endl;
            cout << "it= " << it << " has support 2^{-"
                 << "j= "<< it.j() << " + "<< " e= " << it.e()
                 << "}["
                 << supp.a[0] << ", "<<supp.b[0]
                 << "]x["
                 << supp.a[1] << ", "<<supp.b[1]
                 << "]"
                 << endl;
            if (i>10) break;
        }


#if 0

        // static double primal_regularity() { return IBASIS::primal_regularity(); }
    	// static unsigned int primal_polynomial_degree() { return IBASIS::primal_polynomial_degree(); }
    	// static unsigned int primal_vanishing_moments() { return IBASIS::primal_vanishing_moments(); }
        cout << "Testing primal_ functions"<<endl;
        cout << "primal_regularity "<< basis.primal_regularity() << endl
             << "primal_polynomial_degree " << basis.primal_polynomial_degree() << endl
             << "primal_vanishing_moments " << basis.primal_vanishing_moments() << endl;
#endif

        cout << "HIER" << endl;
        typedef Index::level_type level_type;
        typedef Index::type_type type_type;
        typedef Index::translation_type translation_type;

        level_type my_j;
        type_type my_e;
        translation_type my_k;
        my_j[0]=basis.bases()[0]->j0(); // == 4
        my_j[1]=basis.bases()[1]->j0();
        my_e[0]=0;
        my_e[1]=0;
        my_k[0]=basis.bases()[0]->DeltaLmin()+0; // DeltaLmin ==2, Deltasize(4) == 12
        my_k[1]=basis.bases()[0]->DeltaLmin()+1;

        //TensorIndex(const level_type& j, const type_type& e, const translation_type& k, const TENSORBASIS* basis);
        Index myindex (my_j,my_e,my_k,&basis);
        int res(7);
        
        typedef Basis1d::Index IIndex;
        InfiniteVector<double, IIndex> coeff1d;
        IIndex index1d(first_generator(basis.bases()[0], basis.bases()[0]->j0()));
        coeff1d.set_coefficient(index1d,1.0);
        int dil=8;
        Array1D<double> grid1d;
        Array1D<double> funval1d;
        grid1d.resize((1<<dil)+1);
        for(unsigned int k=0;k<grid1d.size();k++)
            grid1d[k]=k*(1.0/(1<<dil));
        Grid<1>gitter(grid1d);
        //evaluate(bas11,0,index1d,grid1d,funval1d);
        //evaluate(*(basesArray[0]),0,index1d,grid1d,funval1d);
        evaluate(*(basis.bases()[0]),0,index1d,grid1d,funval1d);
        SampledMapping<1> my_wavelet1d (gitter, funval1d);
        ofstream fileA("my_wavelet1d.m");
        my_wavelet1d.matlab_output(fileA);
        fileA.close();
        
cout << "HIER2" << endl;

cout << "j0 " <<basis.bases()[0]->j0() << endl;
cout << "deltasize " << basis.bases()[0]->Deltasize(basis.bases()[0]->j0())<<endl;
cout << "deltaLmin " << basis.bases()[0]->DeltaLmin() << endl;
cout << "firstgen " << first_generator(basis.bases()[0],basis.bases()[0]->j0()) << endl;
cout << Basis1d::Index(my_j[0],my_e[0],my_k[0],basis.bases()[0]) << endl;
        SampledMapping<1> my_wavelet1d00(evaluate(*(basis.bases()[0]), Basis1d::Index(my_j[0],my_e[0],my_k[0],basis.bases()[0]), false,res ));

//        SampledMapping<1> my_wavelet1d00(evaluate(*(basis.bases()[0]), first_generator(basis.bases()[0],4), false,res ));
        SampledMapping<1> my_wavelet1d01(evaluate(*(basis.bases()[0]), Basis1d::Index(my_j[0],my_e[0],my_k[0],basis.bases()[0]), true ,res ));
        SampledMapping<1> my_wavelet1d10(evaluate(*(basis.bases()[1]), Basis1d::Index(my_j[1],my_e[1],my_k[1],basis.bases()[1]), false,res ));
        SampledMapping<1> my_wavelet1d11(evaluate(*(basis.bases()[1]), Basis1d::Index(my_j[1],my_e[1],my_k[1],basis.bases()[1]), true ,res ));

        ofstream file1d00("my_wavelet1d00.m");
        ofstream file1d01("my_wavelet1d01.m");
        ofstream file1d10("my_wavelet1d10.m");
        ofstream file1d11("my_wavelet1d11.m");

        my_wavelet1d00.matlab_output(file1d00);
        my_wavelet1d01.matlab_output(file1d01);
        my_wavelet1d10.matlab_output(file1d10);
        my_wavelet1d11.matlab_output(file1d11);

        file1d00.close();
        file1d01.close();
        file1d10.close();
        file1d11.close();

        cout << "HIER3" << endl;

        SampledMapping<2> my_wavelet0(evaluate(basis, myindex,false,res));
        SampledMapping<2> my_wavelet1(evaluate(basis, myindex,true,res));
        
        ofstream file0("my_wavelet0.m");
        ofstream file1("my_wavelet1.m");
        my_wavelet0.matlab_output(file0);
        my_wavelet1.matlab_output(file1);
        file0.close();
        file1.close();
	return 0;
}