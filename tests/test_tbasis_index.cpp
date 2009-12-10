#include <iostream>

#include <interval/ds_basis.h>
#include <utils/fixed_array1d.h>
#include <cube/tbasis.h>
#include <utils/multiindex.h>

using namespace std;
using namespace WaveletTL;

using MathTL::FixedArray1D;

int main()
{
    cout << "Testing tbasis_index" << endl;
    typedef DSBasis<2,2> Basis1d;
    typedef TensorBasis<Basis1d,2> Basis;
    //typedef CubeBasis<Basis1d,2> Basis;
    typedef Basis::Index Index;
    Basis basisH;
    FixedArray1D<int,4> s; //,st;
    s[0]=s[1]=s[2]=s[3]=0;
    // st[0]=st[1]=st[2]=st[3]=0;
    Basis basisS(s);
    FixedArray1D<bool,4> bc;
    bc[0]=bc[1]=bc[2]=bc[3]=true;
    Basis basisBC(bc);
    Basis1d bas11(true,true);
    Basis1d bas00(false,false);
    bas11.set_jmax(5);
    bas00.set_jmax(5);
    FixedArray1D<Basis1d*,2> basesArray;
    basesArray[0] = basesArray[1] = &bas11;
    Basis basisBas(basesArray);

    // Run tests with this basis:
    //Basis basis;
    //Basis basis(s);
    //Basis basis(bc);
    Basis basis(basesArray);

    cout << "1d bases for basisBas: bas11 in each dimension" << endl;
    cout << "DeltaLmin = " << bas11.DeltaLmin()<< endl
         << "DeltaLmax = " << bas11.DeltaLmax() << endl
         << "Delta0min = " << bas11.Delta0min() << endl
         << "Delta0max = " << bas11.Delta0max(3) << endl
         << "DeltaRmin = " << bas11.DeltaRmin(3) << endl
         << "DeltaRmax = " << bas11.DeltaRmax(3) << endl
         << "Deltasize = " << bas11.Deltasize(3)<< endl
         << "Nablamin =  " << bas11.Nablamin()<< endl
         << "Nablamax =  " << bas11.Nablamax(3)<< endl
         << "Nablasize3= " << bas11.Nablasize(3)<< endl
         << "Nablasize4= " << bas11.Nablasize(4)<< endl
         << "Nablasize5= " << bas11.Nablasize(5)<< endl
         << "Nablasize6= " << bas11.Nablasize(6)<< endl;
    cout     << TensorIndex<Basis1d,2> (1600,&basisH) << " "
             << TensorIndex<Basis1d,2> (1600,&basisS) << " "
             << TensorIndex<Basis1d,2> (1600,&basisBC) << " "
             << TensorIndex<Basis1d,2> (1600,&basisBas) << endl;
    
    cout << "DeltaLmax()  " << basis.bases()[0]->DeltaLmax() << " " << basis.bases()[1]->DeltaLmax() << endl
         << "DeltLmin()   " << basis.bases()[0]->DeltaLmin() << " " << basis.bases()[1]->DeltaLmin() << endl
         << "Deltasize(3) " << basis.bases()[0]->Deltasize(3)<< " " << basis.bases()[1]->Deltasize(3)<< endl
         << "Nablamin()   " << basis.bases()[0]->Nablamin()  << " " << basis.bases()[1]->Nablamin()  << endl
         << "Nablamax(3)  " << basis.bases()[0]->Nablamax(3) << " " << basis.bases()[1]->Nablamax(3) << endl
         << "Nablasize(3) " << basis.bases()[0]->Nablasize(3)<< " " << basis.bases()[1]->Nablasize(3)<< endl
         << "Nablamax(4)  " << basis.bases()[0]->Nablamax(4) << " " << basis.bases()[1]->Nablamax(4) << endl
         << "Nablasize(4) " << basis.bases()[0]->Nablasize(4)<< " " << basis.bases()[1]->Nablasize(4)<< endl
         << "Nablamax(5)  " << basis.bases()[0]->Nablamax(5) << " " << basis.bases()[1]->Nablamax(5) << endl
         << "Nablasize(5) " << basis.bases()[0]->Nablasize(5)<< " " << basis.bases()[1]->Nablasize(5)<< endl
         << "Nablamax(6)  " << basis.bases()[0]->Nablamax(6) << " " << basis.bases()[1]->Nablamax(6) << endl
         << "Nablasize(6) " << basis.bases()[0]->Nablasize(6)<< " " << basis.bases()[1]->Nablasize(6)<< endl;
     

#if 1

/*
	first_generator(const TENSORBASIS* basis);
	last_generator(const TENSORBASIS* basis);
	first_wavelet(const TENSORBASIS* basis, const typename TensorIndex<IBASIS,DIM,TENSORBASIS>::level_type j);
	last_wavelet(const TENSORBASIS* basis, const typename TensorIndex<IBASIS,DIM,TENSORBASIS>::level_type j);
	first_generator_num(const TENSORBASIS* basis);
	last_generator_num(const TENSORBASIS* basis);
	first_wavelet_num(const TENSORBASIS* basis, const typename TensorIndex<IBASIS,DIM,TENSORBASIS>::level_type j);
	last_wavelet_num(const TENSORBASIS* basis, const typename TensorIndex<IBASIS,DIM,TENSORBASIS>::level_type j);
*/
    cout << "Testing first/last routines" << endl;


    cout << "first_generator = " << first_generator<Basis1d,2,Basis> (&basis) << " .number()= "<< first_generator<Basis1d,2,Basis> (&basis).number() << endl;
    cout << "last_generator = " << last_generator<Basis1d,2,Basis> (&basis) << " .number()= "<< last_generator<Basis1d,2,Basis> (&basis).number() << endl;

		TensorIndex<Basis1d,2,Basis>::level_type j;
		TensorIndex<Basis1d,2,Basis>::type_type e;
		TensorIndex<Basis1d,2,Basis>::translation_type k;
/*
                j = basis.j0();
		for (unsigned int i = 0; i < 2; i++)
		{
			e[i] = 0;
			k[i] = basis.bases()[i]->DeltaRmax(j[i]);
		}
                TensorIndex<Basis1d,2,Basis> temp_lg (j, e, k, &basis) ;
                cout << "last generator?? "<<  temp_lg << " .number()= "<< temp_lg.number() << endl;
*/

    //typedef TensorIndex<Basis1d,2,Basis>::level_type level_type;
    //level_type j;
    j[0]=3;j[1]=3;
    TensorIndex<Basis1d,2,Basis> temp_fw;
    //first_wavelet<Basis1d,2,Basis> (&basis,j);

    temp_fw = first_wavelet<Basis1d,2,Basis> (&basis,j);
    cout << "first_wavelet("<<j<<") = " << temp_fw;
    cout << " .number()= "<< temp_fw.number() << endl;

    cout << "last_wavelet("<<j<<") = " << last_wavelet<Basis1d,2,Basis> (&basis,j) << " .number()= "<< last_wavelet<Basis1d,2,Basis> (&basis,j).number() << endl;
    j[0]=6;j[1]=4;
    cout << "first_wavelet("<<j<<") = " << first_wavelet<Basis1d,2,Basis> (&basis,j) << " .number()= "<< first_wavelet<Basis1d,2,Basis> (&basis,j).number() << endl;
    cout << "last_wavelet("<<j<<") = " << last_wavelet<Basis1d,2,Basis> (&basis,j) << " .number()= "<< last_wavelet<Basis1d,2,Basis> (&basis,j).number() << endl;

    cout << "first_generator_num = " << first_generator_num<Basis1d,2,Basis> (&basis) << endl;
    cout << "last_generator_num = " << last_generator_num<Basis1d,2,Basis> (&basis) << endl;
    j[0]=3;j[1]=3;
    cout << "first_wavelet_num("<<j<<") = " << first_wavelet_num<Basis1d,2,Basis> (&basis,j) << endl;
    cout << "last_wavelet_num("<<j<<") = " << last_wavelet_num<Basis1d,2,Basis> (&basis,j) << endl;
    j[0]=6;j[1]=4;
    cout << "first_wavelet_num("<<j<<") = " << first_wavelet_num<Basis1d,2,Basis> (&basis,j) << endl;
    cout << "last_wavelet_num("<<j<<") = " << last_wavelet_num<Basis1d,2,Basis> (&basis,j) << endl;

#endif


#if 1
    cout << "Testing TensorIndex constructors" << endl;
      	// TensorIndex(const TENSORBASIS* basis = 0);
    	// TensorIndex(const level_type& j, const type_type& e, const translation_type& k, const TENSORBASIS* basis);
    	// TensorIndex(const TensorIndex& lambda);
    	// TensorIndex(const TensorIndex* lambda);
    	// TensorIndex(const int number, const TENSORBASIS* basis);
    int i=0, imax(4000);
    cout << "Testing ++ and construction by number" << endl;    
    
    /*
        cout << "by number" << endl;
    for (unsigned int i = 499; i < 1000; i++)
        cout << "i= "<<i<<" index = "<< TensorIndex<Basis1d,2> (i,&basis) << endl;
    */
/*
    for (TensorIndex<Basis1d,2> it(first_generator<Basis1d,2,Basis> (&basis));i<imax;++i, ++it)
    {
        cout << " it= "<< it << " index= "<< TensorIndex<Basis1d,2> (i, &basis) <<"i/it.n/ind.n="<<i<<" "<< it.number()<<" "<<TensorIndex<Basis1d,2> (i, &basisH).number()<<endl;
    }
    */
    for (TensorIndex<Basis1d,2> it(first_generator<Basis1d,2,Basis> (&basis));i<imax;++i, ++it)
        if (it != TensorIndex<Basis1d,2> (i, &basis))
        {
            cout << "Problem entdeckt fuer Index Nummer "<<i<<endl;
            cout << " it= "<< it << " index= "<< TensorIndex<Basis1d,2> (i, &basis) <<"i/it.n/ind.n="<<i<<" "<< it.number()<<" "<<TensorIndex<Basis1d,2> (i, &basisH).number()<<endl;
            break;
        }

    
    cout << "Testing copy constructors" << endl;
    for (TensorIndex<Basis1d,2> it(first_generator<Basis1d,2,Basis> (&basis));it.number()<imax;++it)
    {
        TensorIndex<Basis1d,2> temp_index1(it);
        TensorIndex<Basis1d,2> temp_index2(&it);
        if (temp_index1 != it || temp_index1.number() != it.number())
        {
            cout << "1::Problem entdeckt fuer Index Nummer it.n/t_it.n"<<it.number()<< " "<<temp_index1.number() << endl;
            cout << " it= "<< it << " t_ind1= "<< temp_index1<<endl;
            break;
        }
        if (temp_index2 != it || temp_index2.number() != it.number())
        {
            cout << "2::Problem entdeckt fuer Index Nummer it.n(t_it.n"<<it.number()<< " "<<temp_index2.number() << endl;
            cout << " it= "<< it << " t_ind2= "<< temp_index2<<endl;
            break;
        }
    }

    cout << "Testing constructor by j,e,k" << endl;
    // TensorIndex(const level_type& j, const type_type& e, const translation_type& k, const TENSORBASIS* basis);
    for (TensorIndex<Basis1d,2> it(first_generator<Basis1d,2,Basis> (&basis));it.number()<imax;++it)
    {
        TensorIndex<Basis1d,2> temp_index(it.j(),it.e(),it.k(),&basis);
        if (temp_index != it || temp_index.number() != it.number())
        {
            cout << "Problem entdeckt fuer Index Nummer it.n/t_it.n"<<it.number()<< " "<<temp_index.number() << endl;
            cout << " it= "<< it << " temp_index= "<< temp_index<< "it.n, t_i.n = " << it.number() << ", " << temp_index.number() << endl;
            break;
        }
    }
#endif

#if 0
    bool error = false;
    cout << "Testing operator <=" << endl;
    for (int i=0;i<200;i++)
    {
        cout << i << endl;
        TensorIndex<Basis1d,2> outer(i*20,&basis);
        for (int j=0;j<4000;j++)
        {
            TensorIndex<Basis1d,2> inner(j,&basis);
            if ((outer<=inner) != (i*20 <= j))
            {
                cout << "Problem entdeckt fuer (i*10,j) = (" << (i*20) <<", "<<j<<")"<<endl;
                cout << " outer <= inner "<< outer << " <= "<< inner << " = " << (outer <= inner) << endl;
                cout << " i*10 <= j "<< (i*20) << " <= "<< j << " = " << (i*20 <= j) << endl;
                cout << "((outer<=inner) != (i*20 <= j)) = " << ((outer<=inner) != (i*20 <= j)) << endl;
                error = true;
                break;
            }
        }
        if (error == true) break;
    }
#endif

    return 0;
}
