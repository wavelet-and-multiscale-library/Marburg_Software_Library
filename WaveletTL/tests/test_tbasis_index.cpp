#include <iostream>

#include <interval/p_basis.h>
#include <interval/ds_basis.h>
#include <utils/fixed_array1d.h>
#include <cube/tbasis.h>
#include <cube/tbasis_index.h>
#include <utils/multiindex.h>

using namespace std;
using namespace WaveletTL;

using MathTL::FixedArray1D;

int main()
{
    cout << "Testing tbasis_index" << endl;
#if 0
    //typedef DSBasis<2,2> Basis1d;
    typedef PBasis<2,2> Basis1d;
    typedef TensorBasis<Basis1d,2> Basis;
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
    const unsigned int dim = 2;
#else

    const int d  = 2;
    const int dT = 2;
    //const unsigned int dim = 2; const int levelrange(4);
    const unsigned int dim = 2; const int levelrange(1);
    typedef PBasis<d,dT> Basis1d;
    //typedef DSBasis<d,dT> Basis1d;
    typedef TensorBasis<Basis1d,dim> Basis;
//    typedef Basis::Index Index;
//    typedef TensorIndex<Basis1d,dim>::level_type level_type;
//    typedef TensorIndex<Basis1d,dim>::type_type type_type;
//    typedef TensorIndex<Basis1d,dim>::translation_type translation_type;
    
    Basis basisH;
    FixedArray1D<int,2*dim> s; //,st;
    for (unsigned int i(0);i<(2*dim);i++) //s[0]=s[1]=s[2]=s[3]=s[4]=s[5]=0;
    {
        s[i]=0;
    }
    // st[0]=st[1]=st[2]=st[3]=0;
    Basis basisS(s);
    FixedArray1D<bool,2*dim> bc;
    for (unsigned int i=0;i<2*dim;i++)
    {
        bc[i]=false;
    }
    //bc[0]=bc[1]=bc[2]=bc[3]=bc[4]=bc[5]=true;
    Basis basisBC(bc);
    Basis1d bas00(false,false); Basis1d bas10(true,false); Basis1d bas01(false,true); Basis1d bas11(true,true);
    bas00.set_jmax(5); bas10.set_jmax(5); bas01.set_jmax(5); bas11.set_jmax(5);
    FixedArray1D<Basis1d*,dim> basesArray;
    for(unsigned int i=0;i<dim;i++)
    {
        basesArray[i] = &bas00;
    }
    
#endif

    // Run tests with this basis:
    //Basis basis;
    //Basis basis(s);
    //Basis basis(bc);
    Basis basis(basesArray);
    basis.set_jmax(multi_degree(basis.j0())+levelrange);
    
    Basis basisBas(basesArray); basisBas.set_jmax(multi_degree(basisBas.j0())+levelrange);
    basesArray[0] = &bas10;
    Basis basisBas2(basesArray);basisBas2.set_jmax(multi_degree(basisBas2.j0())+levelrange);
    basesArray[0] = &bas01;
    Basis basisBas3(basesArray);basisBas3.set_jmax(multi_degree(basisBas3.j0())+levelrange);
    basesArray[dim-1] = &bas10;
    Basis basisBas4(basesArray);basisBas4.set_jmax(multi_degree(basisBas4.j0())+levelrange);
    basesArray[1] = &bas01;
    Basis basisBas5(basesArray);basisBas5.set_jmax(multi_degree(basisBas5.j0())+levelrange);
    basesArray[0] = &bas11; basesArray[dim-1] = &bas01;
    Basis basisBas6(basesArray);basisBas6.set_jmax(multi_degree(basisBas6.j0())+levelrange);

    FixedArray1D<Basis*,6> TBasisArray;
    TBasisArray[0] = &basisBas;
    TBasisArray[1] = &basisBas2;
    TBasisArray[2] = &basisBas3;
    TBasisArray[3] = &basisBas4;
    TBasisArray[4] = &basisBas5;
    TBasisArray[5] = &basisBas6;
    
    cout << "1d bases bas11 in each dimension" << endl;
    cout << "DeltaLmin  = " << bas11.DeltaLmin()  << endl
         << "DeltaLmax  = " << bas11.DeltaLmax()  << endl
         << "Delta0min  = " << bas11.Delta0min()  << endl
         << "Delta0max3 = " << bas11.Delta0max(3) << endl
         << "DeltaRmin3 = " << bas11.DeltaRmin(3) << endl
         << "DeltaRmax3 = " << bas11.DeltaRmax(3) << endl
         << "Deltasize3 = " << bas11.Deltasize(3) << endl
         << "Nablamin   = " << bas11.Nablamin()   << endl
         << "Nablamax3  = " << bas11.Nablamax(3)  << endl
         << "Nablasize3 = " << bas11.Nablasize(3) << endl
         << "Nablasize4 = " << bas11.Nablasize(4) << endl
         << "Nablasize5 = " << bas11.Nablasize(5) << endl
         << "Nablasize6 = " << bas11.Nablasize(6) << endl;
    cout << TensorIndex<Basis1d,dim> (1600,&basisH) << " ";
    cout << TensorIndex<Basis1d,dim> (1600,&basisS) << " ";
    cout << TensorIndex<Basis1d,dim> (1600,&basisBC) << " ";
    cout << TensorIndex<Basis1d,dim> (1600,&basisBas) << endl;
    cout << "minimal levels. H = " << basisH.j0() << " S = " << basisS.j0() << " BC = " << basisBC.j0() << " Bas = " << basisBas.j0() << endl;
    
    cout << "display informations for the selected basis" << endl;
    cout << "DeltaLmax()  " << basis.bases()[0]->DeltaLmax();  if (dim > 1) cout << " " << basis.bases()[1]->DeltaLmax();  if (dim > 2) cout << " " << basis.bases()[2]->DeltaLmax();  cout << endl;
    cout << "DeltLmin()   " << basis.bases()[0]->DeltaLmin();  if (dim > 1) cout << " " << basis.bases()[1]->DeltaLmin();  if (dim > 2) cout << " " << basis.bases()[2]->DeltaLmin();  cout << endl;
    cout << "Deltasize(" << basis.j0() << ") " << basis.bases()[0]->Deltasize(basis.j0()[0]); 
    if (dim > 1) cout << " " << basis.bases()[1]->Deltasize(basis.j0()[1]); 
    if (dim > 2) cout << " " << basis.bases()[2]->Deltasize(basis.j0()[2]); cout << endl;
    
    cout << "Nablamin()   " << basis.bases()[0]->Nablamin();   if (dim > 1) cout << " " << basis.bases()[1]->Nablamin();   if (dim > 2) cout << " " << basis.bases()[2]->Nablamin();   cout << endl;
    
    cout << "Nablamax(" << basis.j0() << ") " << basis.bases()[0]->Nablamax(basis.j0()[0]); 
    if (dim > 1) cout << " " << basis.bases()[1]->Nablamax(basis.j0()[1]);  
    if (dim > 2) cout << " " << basis.bases()[2]->Nablamax(basis.j0()[2]);  cout << endl;
    cout << "Nablasize(" << basis.j0() << ") " << basis.bases()[0]->Nablasize(basis.j0()[0]); 
    if (dim > 1) cout << " " << basis.bases()[1]->Nablasize(basis.j0()[1]);  
    if (dim > 2) cout << " " << basis.bases()[2]->Nablasize(basis.j0()[2]);  cout << endl;
    
    cout << "Nablamax(" << basis.j0() << "+ 1 ) " << basis.bases()[0]->Nablamax(basis.j0()[0]+1); 
    if (dim > 1) cout << " " << basis.bases()[1]->Nablamax(basis.j0()[1]+1);  
    if (dim > 2) cout << " " << basis.bases()[2]->Nablamax(basis.j0()[2]+1);  cout << endl;
    cout << "Nablasize(" << basis.j0() << "+ 1 ) " << basis.bases()[0]->Nablasize(basis.j0()[0]+1); 
    if (dim > 1) cout << " " << basis.bases()[1]->Nablasize(basis.j0()[1]+1);  
    if (dim > 2) cout << " " << basis.bases()[2]->Nablasize(basis.j0()[2]+1);  cout << endl;
    
    cout << "Nablamax(" << basis.j0() << "+ 2 ) " << basis.bases()[0]->Nablamax(basis.j0()[0]+2); 
    if (dim > 1) cout << " " << basis.bases()[1]->Nablamax(basis.j0()[1]+2);
    if (dim > 2) cout << " " << basis.bases()[2]->Nablamax(basis.j0()[2]+2);  cout << endl;
    cout << "Nablasize(" << basis.j0() << "+ 2 ) " << basis.bases()[0]->Nablasize(basis.j0()[0]+2);
    if (dim > 1) cout << " " << basis.bases()[1]->Nablasize(basis.j0()[1]+2);
    if (dim > 2) cout << " " << basis.bases()[2]->Nablasize(basis.j0()[2]+2);  cout << endl;
    
    cout << "Nablamax(" << basis.j0() << "+ 3 ) " << basis.bases()[0]->Nablamax(basis.j0()[0]+3); 
    if (dim > 1) cout << " " << basis.bases()[1]->Nablamax(basis.j0()[1]+3);
    if (dim > 2) cout << " " << basis.bases()[2]->Nablamax(basis.j0()[2]+3);  cout << endl;
    cout << "Nablasize(" << basis.j0() << "+ 3 ) " << basis.bases()[0]->Nablasize(basis.j0()[0]+3);
    if (dim > 1) cout << " " << basis.bases()[1]->Nablasize(basis.j0()[1]+3);
    if (dim > 2) cout << " " << basis.bases()[2]->Nablasize(basis.j0()[2]+3);  cout << endl;
    
    
    
#if 1
    // some output for debugging
    for (unsigned int i=0; (int)i<TBasisArray[0]->degrees_of_freedom(); ++i)
    {
        cout << "i = " << i << "; lambda(i) = " << *TBasisArray[0]->get_wavelet(i) << endl;
    }
    abort();
#endif
#if 0
    cout << "Testing TensorIndex constructors" << endl;
    //TensorIndex(const TENSORBASIS* basis = 0);
    TensorIndex<Basis1d,dim> tindex1;
    cout << "TensorIndex() " << tindex1 << "; .number() = " << tindex1.number() << endl;
    
    //TensorIndex(const level_type& j, const type_type& e, const translation_type& k, const TENSORBASIS* basis);
    
    level_type temp_j;
    type_type temp_e;
    translation_type temp_k;
    for (unsigned int i=0; i< dim; i++)
    {
        temp_j[i] = basis.j0()[i] + 2;
        temp_e[i] = 1;
        temp_k[i] = basis.bases()[i]->Nablamin();
    }
    //cout << "j = " << temp_j << "; e = " << temp_e << "; k = " << temp_k << endl;
    TensorIndex<Basis1d,dim> tindex2(temp_j,temp_e,temp_k,&basis);
    cout << "TensorIndex("<< temp_j << ", " << temp_e << ", " << temp_k << ", basis) " << tindex2 << "; .number() = " << tindex2.number() << endl;
    
    //TensorIndex(const int& j, const int& e, const int& k, const TENSORBASIS* empty);
    
    typedef TensorBasis<Basis1d,1> TBasis1d;
    //typedef Basis::Index Index;
    
    FixedArray1D<Basis1d*,1> basesArray1d_dummy;
    basesArray1d_dummy[0] = &bas00;
    TBasis1d basisBas1d_dummy(basesArray1d_dummy);
    TensorIndex<Basis1d,1> tindex3(4,1,3,&basisBas1d_dummy);
    cout << "TensorIndex<Basis1d,1> tindex3(4,1,3,bas00) " << tindex3 << "; .number() = " << tindex3.number() << endl;

    TensorIndex<Basis1d,1> tindex3a(3,1,0,&basisBas1d_dummy);
    cout << "TensorIndex<Basis1d,1> tindex3(4,1,3,bas00) " << tindex3a << "; .number() = " << tindex3a.number() << endl;
    TensorIndex<Basis1d,1> tindex3b(3,1,1,&basisBas1d_dummy);
    cout << "TensorIndex<Basis1d,1> tindex3(4,1,3,bas00) " << tindex3b << "; .number() = " << tindex3b.number() << endl;
    TensorIndex<Basis1d,1> tindex3c(3,1,2,&basisBas1d_dummy);
    cout << "TensorIndex<Basis1d,1> tindex3(4,1,3,bas00) " << tindex3c << "; .number() = " << tindex3c.number() << endl;
    TensorIndex<Basis1d,1> tindex3d(3,1,3,&basisBas1d_dummy);
    cout << "TensorIndex<Basis1d,1> tindex3(4,1,3,bas00) " << tindex3d << "; .number() = " << tindex3d.number() << endl;
    
    //TensorIndex(const TensorIndex& lambda);
    TensorIndex<Basis1d,dim> tindex4(tindex2);
    cout << "TensorIndex(" << tindex2 << ") " << tindex4 << "; .number() = " << tindex4.number() << endl;
    
    //TensorIndex(const TensorIndex* lambda);
    TensorIndex<Basis1d,dim> tindex5(tindex2);
    cout << "TensorIndex(" << tindex2 << ") " << tindex5 << "; .number() = " << tindex5.number() << endl;
    
    //TensorIndex(const level_type& j, const type_type& e, const translation_type& k, const int number, const TENSORBASIS* basis);
    TensorIndex<Basis1d,dim> tindex6(temp_j,temp_e,temp_k,47, &basis);
    cout << "TensorIndex("<< temp_j << ", " << temp_e << ", " << temp_k << ", 47, basis) " << tindex6 << "; .number() = " << tindex6.number() << endl;
    
    //TensorIndex(const int number, const TENSORBASIS* basis);
    TensorIndex<Basis1d,dim> tindex7(tindex2.number(),&basis);
    cout << "TensorIndex(" << tindex2.number() << ",basis) " << tindex7 << "; .number() = " << tindex7.number() << endl;
    
    TensorIndex<Basis1d,dim> tindex7a(tindex2.number()-1,&basis);
    cout << "TensorIndex(" << (tindex2.number()-1) << ",basis) " << tindex7a << "; .number() = " << tindex7a.number() << endl;
    TensorIndex<Basis1d,dim> tindex7b(tindex2.number()+2,&basis);
    cout << "TensorIndex(" << (tindex2.number()+2) << ",basis) " << tindex7b << "; .number() = " << tindex7b.number() << endl;
    
    cout << "Testing operators" << endl;
    cout << "=,==,!=,++" << endl;
    cout << "tindex1 = " << tindex1 << "; .number() = " << tindex1.number() << "; tindex2 = " << tindex2 << "; .number() = " << tindex2.number() << endl;
    tindex1 = tindex2;
    cout << "tindex1 = tindex2 => tindex1 = " << tindex1 << "; .number() = " << tindex1.number() << endl;
    cout << "tindex6 = " << tindex6 << "; .number() = " << tindex6.number() << endl;
    cout << "tindex2 == tindex6 => " << (tindex2 == tindex6) << endl;
    cout << "tindex2 != tindex6 => " << (tindex2 != tindex6) << endl;
    cout << "++tindex1 => " << (++tindex1) << endl; 
    cout << "; tindex1 = " << tindex1 << "; .number() = " << tindex1.number() <<  endl;
            
#endif
    
#if 0
    cout << "Testing ++ and construction by number" << endl;    
/*
    cout << "by number" << endl;
    for (unsigned int i = 0; i < 20; i++)
    {
        cout << "i= "<<i<<" index = "<< TensorIndex<Basis1d,dim> (i,&basis) << endl;
    }
*/
    //TensorIndex<Basis1d,dim> (15, &basis);

    int minlambdanum(0000), maxlambdanum(500);
    assert(maxlambdanum < basis.degrees_of_freedom());
    Index lambda_run(basis.get_wavelet(minlambdanum));
    
    //TensorIndex<Basis1d,dim> (272, &basis);
#if 1
    for (unsigned int i=minlambdanum; i< maxlambdanum; i++)
    {
        cout << "N = " << i << "; " << (*basis.get_wavelet(i)).number() << "; " << *basis.get_wavelet(i) << "; " << TensorIndex<Basis1d,dim> (i, &basis).number() << "; " << TensorIndex<Basis1d,dim> (i, &basis) << endl;
    }
#endif
    
    cout << last_wavelet<Basis1d,dim,Basis>(&basis,multi_degree(basis.j0())) << "; .number() = " << last_wavelet<Basis1d,dim,Basis>(&basis,multi_degree(basis.j0())).number() << endl;
    Index temp_ind(last_wavelet<Basis1d,dim,Basis>(&basis,multi_degree(basis.j0())));
    Index temp_ind2(temp_ind.j(), temp_ind.e(), temp_ind.k(), temp_ind.basis());
    
    cout <<  temp_ind2 << "; .number() = " << temp_ind2.number() << endl;
    
    
    for (unsigned int num(minlambdanum); num < maxlambdanum; ++num)
    {
        Index tempi1(basis.get_wavelet(num));
        Index tempi2(TensorIndex<Basis1d,dim> (num, &basis));
        Index tempi3(tempi1.j(),tempi1.e(),tempi1.k(),tempi1.basis());
        if ((((tempi1 != tempi2) || (tempi2 != tempi3))
            || (tempi1.number() != tempi2.number()))
            || (tempi2.number() != tempi3.number()))
        {
            cout << "Problem entdeckt!" << endl;
            cout << "N = "<< num <<" tempi1 = ("<< tempi1.number() << ")=" << tempi1 
                                 << "; tempi2 = ("<< tempi2.number() << ")=" << tempi2 
                                 << "; tempi3 = ("<< tempi3.number() << ")=" << tempi3 << endl; 
            abort();
            //break;
        }
    }
    for (TensorIndex<Basis1d,dim> it(basis.get_wavelet(minlambdanum));minlambdanum<maxlambdanum;++minlambdanum, ++it)
    {
        if ( ( (it != TensorIndex<Basis1d,dim> (minlambdanum, &basis)) || (it.number() != TensorIndex<Basis1d,dim> (minlambdanum, &basis).number()) ) || (it.number() != minlambdanum) )
        {
            cout << "Problem entdeckt!" << endl;
            cout << "i = "<< minlambdanum <<" it= ("<<it.number()<<") = "<< it << " index= ("<< TensorIndex<Basis1d,dim> (minlambdanum, &basis).number() << ") = "<<TensorIndex<Basis1d,dim> (minlambdanum, &basis) <<endl;
            abort();
            //break;
        }
    }
    
    cout << "Testing copy constructors" << endl;
    for (TensorIndex<Basis1d,dim> it(basis.get_wavelet(minlambdanum));it.number()<maxlambdanum;++it)
    {
        TensorIndex<Basis1d,dim> temp_index1(it);
        TensorIndex<Basis1d,dim> temp_index2(&it);
        if (temp_index1 != it || temp_index1.number() != it.number())
        {
            cout << "1::Problem entdeckt fuer Index Nummer it.n/t_it.n"<<it.number()<< " "<<temp_index1.number() << endl;
            cout << " it= "<< it << " t_ind1= "<< temp_index1<<endl;
            abort();
            //break;
        }
        if (temp_index2 != it || temp_index2.number() != it.number())
        {
            cout << "2::Problem entdeckt fuer Index Nummer it.n(t_it.n"<<it.number()<< " "<<temp_index2.number() << endl;
            cout << " it= "<< it << " t_ind2= "<< temp_index2<<endl;
            abort();
            //break;
        }
    }

    cout << "Testing constructor by j,e,k" << endl;
    // TensorIndex(const level_type& j, const type_type& e, const translation_type& k, const TENSORBASIS* basis);
    for (TensorIndex<Basis1d,dim> it(basis.get_wavelet(minlambdanum));it.number()<maxlambdanum;++it)
    {
        TensorIndex<Basis1d,dim> temp_index(it.j(),it.e(),it.k(),&basis);
        //cout << temp_index << endl;
        if (temp_index != it || temp_index.number() != it.number())
        {
            cout << "Problem entdeckt!"<<endl;
            cout << " it= ("<<it.number()<<") = "<< it << " copy= ("<< temp_index.number() << ") = "<<temp_index <<endl;
            abort();
            //break;
        }
    }
#endif

#if 0
    const int step(1),maxnumber(500); // 1,4000 = 3 min

    bool error = false;
    cout << "Testing operator <=" << endl;
    for (int i=0;i<maxnumber/step;i++)
    {
        //cout << i << endl;
        TensorIndex<Basis1d,dim> outer(i*step,&basis);
        for (int j=0;j<maxnumber;j++)
        {
            TensorIndex<Basis1d,dim> inner(j,&basis);
            if ((outer<=inner) != (i*step <= j))
            {
                cout << "Problem entdeckt fuer (i*"<<step<<",j) = (" << (i*step) <<", "<<j<<")"<<endl;
                cout << " outer <= inner "<< outer << " <= "<< inner << " = " << (outer <= inner) << endl;
                cout << " i*"<<step<<" <= j "<< (i*step) << " <= "<< j << " = " << (i*step <= j) << endl;
                cout << "((outer<=inner) != (i*"<<step<<" <= j)) = " << ((outer<=inner) != (i*step <= j)) << endl;
                error = true;
                break;
            }
        }
        if (error == true) break;
    }
#else
    cout << "skipping test of <= operator" << endl;
#endif
    
#if 0
    cout << "Testing first/last _ generator/wavelet" << endl;
    cout << "first_generator(&basis) = " << first_generator<Basis1d,dim,Basis>(&basis) << "; .number() = " << first_generator<Basis1d,dim,Basis>(&basis).number() << endl;
    cout << "last_generator(&basis) = "  << last_generator<Basis1d,dim,Basis>(&basis)  << "; .number() = " << last_generator<Basis1d,dim,Basis>(&basis).number() << endl;

    //first_wavelet(const TENSORBASIS* basis, const typename TensorIndex<IBASIS,DIM,TENSORBASIS>::level_type j);
    cout << "first_wavelet(&basis, " << basis.j0() << ") = " << first_wavelet<Basis1d,dim,Basis>(&basis, basis.j0()) << "; .number() = " << first_wavelet<Basis1d,dim,Basis>(&basis, basis.j0()).number() << endl;
    level_type temp_j2(basis.j0());
    temp_j2[dim-1]=temp_j2[dim-1]+1;
    cout << "first_wavelet(&basis, " << temp_j2 << ") = " << first_wavelet<Basis1d,dim,Basis>(&basis, temp_j2) << "; .number() = " << first_wavelet<Basis1d,dim,Basis>(&basis, temp_j2).number() << endl;
    
    //first_wavelet(const int level);
    for (int i=multi_degree(basis.j0()); i<multi_degree(basis.j0())+5; ++i)
    {
        cout << "first_wavelet(&basis,"<< i <<") = " << first_wavelet<Basis1d,dim,Basis>(&basis,i) << "; .number() = " << first_wavelet<Basis1d,dim,Basis>(&basis,i).number() << endl;
    }
    
    //last_wavelet(const TENSORBASIS* basis, const typename TensorIndex<IBASIS,DIM,TENSORBASIS>::level_type j);
    cout << "last_wavelet(&basis, " << basis.j0() << ") = " << last_wavelet<Basis1d,dim,Basis>(&basis, basis.j0()) << "; .number() = " << last_wavelet<Basis1d,dim,Basis>(&basis, basis.j0()).number() << endl;
    cout << "last_wavelet(&basis, " << temp_j2 << ") = " << last_wavelet<Basis1d,dim,Basis>(&basis, temp_j2) << "; .number() = " << last_wavelet<Basis1d,dim,Basis>(&basis, temp_j2).number() << endl;
    
    //last_wavelet(const TENSORBASIS* basis, const unsigned int level);
    for (int i=multi_degree(basis.j0()); i<multi_degree(basis.j0())+5; ++i)
    {
        cout << "last_wavelet("<< i <<") = " << last_wavelet<Basis1d,dim,Basis>(&basis,i) << "; .number() = " << last_wavelet<Basis1d,dim,Basis>(&basis,i).number() << endl;
    }
    
    //first_generator_num(const TENSORBASIS* basis);
    cout << "first_generator_num(&basis) = " << first_generator_num<Basis1d,dim,Basis>(&basis) << endl;
    
    //last_generator_num(const TENSORBASIS* basis);
    cout << "last_generator_num(&basis) = " << last_generator_num<Basis1d,dim,Basis>(&basis) << endl;
    
    //first_wavelet_num(const TENSORBASIS* basis, const typename TensorIndex<IBASIS,DIM,TENSORBASIS>::level_type j);
    cout << "first_wavelet_num(&basis, " << basis.j0() << ") = " << first_wavelet_num<Basis1d,dim,Basis>(&basis, basis.j0()) << endl;
    cout << "first_wavelet_num(&basis, " << temp_j2 << ") = " << first_wavelet_num<Basis1d,dim,Basis>(&basis, temp_j2) << endl;
    
    //last_wavelet_num(const TENSORBASIS* basis, const typename TensorIndex<IBASIS,DIM,TENSORBASIS>::level_type j);
    cout << "last_wavelet_num(&basis, " << basis.j0() << ") = " << last_wavelet_num<Basis1d,dim,Basis>(&basis, basis.j0()) << endl;
    cout << "last_wavelet_num(&basis, " << temp_j2 << ") = " << last_wavelet_num<Basis1d,dim,Basis>(&basis, temp_j2) << endl;
    
    //last_wavelet_num(const TENSORBASIS* basis, const unsigned int level);
    for (int i = multi_degree(basis.j0()); i < multi_degree(basis.j0())+5; ++i)
    {
        cout << "last_wavelet_num(&basis,"<< i <<") = " << last_wavelet_num<Basis1d,dim,Basis>(&basis,i) <<  endl;
    }
    
    
#endif
    
    
#if 0

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

    Index temp(first_generator<Basis1d,dim,Basis> (&basis));
    typedef Index::level_type index_lt;
    index_lt temp2(temp.j());
    cout << "first gen: " << temp << endl << "first_gen.j()= " << temp2 << endl << "first_gen.j.number= "<< temp2.number() << endl;

    typedef std::map<index_lt,int> map_type;
    typedef map_type::iterator map_it;
    map_type mapp((indexmapping(temp2)));
    for (map_it mit(mapp.begin()), mit_end(mapp.end()); mit!=mit_end;mit++)
    {
        cout << (*mit).first << " => " << (*mit).second << endl;
    }

    //cout << "indexnumbers for temp="<<temp2<< " gives " << endl << indexnumbers(temp2) << endl;
    
    cout << "first_generator = " << first_generator<Basis1d,dim,Basis> (&basis) << " .number()= "<< first_generator<Basis1d,dim,Basis> (&basis).number() << endl;
    cout << "last_generator = " << last_generator<Basis1d,dim,Basis> (&basis) << " .number()= "<< last_generator<Basis1d,dim,Basis> (&basis).number() << endl;

		TensorIndex<Basis1d,dim,Basis>::level_type j;
		TensorIndex<Basis1d,dim,Basis>::type_type e;
		TensorIndex<Basis1d,dim,Basis>::translation_type k;
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

    //typedef TensorIndex<Basis1d,dim,Basis>::level_type level_type;
    //level_type j;
    for (unsigned int k=0;k<dim;k++) j[k]=3;
    // j[0]=3;j[1]=3;
    TensorIndex<Basis1d,dim,Basis> temp_fw;
    //first_wavelet<Basis1d,2,Basis> (&basis,j);

    temp_fw = first_wavelet<Basis1d,dim,Basis> (&basis,j);
    cout << "first_wavelet("<<j<<") = " << temp_fw;
    cout << " .number()= "<< temp_fw.number() << endl;

    cout << "last_wavelet("<<j<<") = " << last_wavelet<Basis1d,dim,Basis> (&basis,j) << " .number()= "<< last_wavelet<Basis1d,dim,Basis> (&basis,j).number() << endl;
    j[0]=6;
    if (dim > 1) j[1]=4;
    cout << "first_wavelet("<<j<<") = " << first_wavelet<Basis1d,dim,Basis> (&basis,j) << " .number()= "<< first_wavelet<Basis1d,dim,Basis> (&basis,j).number() << endl;
    cout << "last_wavelet("<<j<<") = " << last_wavelet<Basis1d,dim,Basis> (&basis,j) << " .number()= "<< last_wavelet<Basis1d,dim,Basis> (&basis,j).number() << endl;

    cout << "first_generator_num = " << first_generator_num<Basis1d,dim,Basis> (&basis) << endl;
    cout << "last_generator_num = " << last_generator_num<Basis1d,dim,Basis> (&basis) << endl;
    j[0]=3;
    if (dim > 1) j[1]=3;
    cout << "first_wavelet_num("<<j<<") = " << first_wavelet_num<Basis1d,dim,Basis> (&basis,j) << endl;
    cout << "last_wavelet_num("<<j<<") = " << last_wavelet_num<Basis1d,dim,Basis> (&basis,j) << endl;
    j[0]=6;
    if (dim > 1) j[1]=4;
    cout << "first_wavelet_num("<<j<<") = " << first_wavelet_num<Basis1d,dim,Basis> (&basis,j) << endl;
    cout << "last_wavelet_num("<<j<<") = " << last_wavelet_num<Basis1d,dim,Basis> (&basis,j) << endl;

#endif

#if 0
    cout << "Testing modified code for first_wavelet(basis,j)" << endl;
    level_type temp_mi, temp_mi2;
    temp_mi = basis.j0();
    int temp_upto (20);
    for (unsigned int i=0; i< temp_upto; ++i)
    {
        for (unsigned int j=0; j<dim; ++j)
        {
            temp_mi[j] = basis.j0()[j]+temp_mi2[j];
        }
        cout << "first_wavelet(&basis," << temp_mi << ") = " << first_wavelet<Basis1d,dim,Basis>(&basis,temp_mi) << "; .number() = " << first_wavelet<Basis1d,dim,Basis>(&basis,temp_mi).number() << endl;
        ++temp_mi2;
    }
#endif
    
#if 0
    cout << "Testing modified code for first_wavelet(basis,level)" << endl;
    
    for (unsigned int level=multi_degree(basis.j0()); level< multi_degree(basis.j0())+15;++level)
    {
        cout << "first_wavelet(&basis," << level << ") = " << first_wavelet<Basis1d,dim,Basis>(&basis,level) << endl;
        cout << "       ; .number() = " << first_wavelet<Basis1d,dim,Basis>(&basis,level).number() << endl;
    }
#endif
    
#if 0
    cout << "Testing modified code for last_wavelet(basis,level)" << endl;
    
    for (unsigned int level=multi_degree(basis.j0()); level< multi_degree(basis.j0())+15;++level)
    {
        assert ((last_wavelet<Basis1d,dim,Basis>(&basis,level) == last_wavelet2<Basis1d,dim,Basis>(&basis,level))
                && (last_wavelet<Basis1d,dim,Basis>(&basis,level).number() == last_wavelet2<Basis1d,dim,Basis>(&basis,level).number()));
    }
    
    int repetitions(500), minlevel(multi_degree(basis.j0()));
    int maxlevel = minlevel + 15;
    
    Index lambda_run2;
    clock_t tstart, tend;
    double time1(0), time2(0);
    
    tstart = clock();
    for (unsigned int r=0; r < repetitions; ++r)
    {
        for (unsigned int level=minlevel; level< maxlevel;++level)
        {
            lambda_run2 = last_wavelet<Basis1d,dim,Basis>(&basis,level);
        }
    }
    tend = clock();
    time1 = (double)(tend-tstart)/CLOCKS_PER_SEC;
    
    tstart = clock();
    for (unsigned int r=0; r < repetitions; ++r)
    {
        for (unsigned int level=minlevel; level< maxlevel;++level)
        {
            lambda_run2 = last_wavelet2<Basis1d,dim,Basis>(&basis,level);
        }
    }
    tend = clock();
    time2 = (double)(tend-tstart)/CLOCKS_PER_SEC;
    
    cout << "Repetitions = " << repetitions << "; minlevel = " << minlevel << "; maxlevel = " << maxlevel << "; "
            << "last_wavelet1 " << time1 << "sec; last_wavelet2 " << time2 << "sec" << endl;
#endif
    
    
#if 0
    cout << "Testing modified code for TensorIndex(num,basis)" << endl;
    
    int minnum(0);
    //int maxnum(100);
    int maxnum(basis.degrees_of_freedom());
    
    for (unsigned int b=0; b< TBasisArray.size(); ++b)
    //int b=1;
    {
        for (unsigned int n = minnum; n < maxnum;++n)
        {
            Index temp_ind1(n,TBasisArray[b]);
            Index temp_ind2(0,n,TBasisArray[b]);
            /*
            cout << "N = " << n << "; temp_ind1 = " << temp_ind1 << "; .number() = " << temp_ind1.number() << "; "
                    "temp_ind2 = " << temp_ind2 << "; .number() = " << temp_ind2.number();

            if (!((temp_ind1 == temp_ind2) && (temp_ind1.number() == temp_ind2.number())))
                cout << " (Problem)";
            cout << endl;
             * */
            assert ((temp_ind1 == temp_ind2) && (temp_ind1.number() == temp_ind2.number()));
        }
    }
    
    int repetitions2(000), minnum2(0), maxnum2(basis.degrees_of_freedom());

    Index lambda_run3;
    clock_t tstart2, tend2;
    double time3(0), time4(0);
    
    tstart2 = clock();
    for (unsigned int r=0; r < repetitions2; ++r)
    {
        for (unsigned int n=minnum2; n < maxnum2;++n)
        {
            Index lambda_run3(n,&basis);
        }
    }
    tend2 = clock();
    time3 = (double)(tend2-tstart2)/CLOCKS_PER_SEC;
    
    tstart2 = clock();
    for (unsigned int r=0; r < repetitions2; ++r)
    {
        for (unsigned int n=minnum2; n < maxnum2;++n)
        {
            Index lambda_run3(0,n,&basis);
        }
    }
    tend2 = clock();
    time4 = (double)(tend2-tstart2)/CLOCKS_PER_SEC;
    
    cout << "Repetitions = " << repetitions2 << "; min lambdanum = " << minnum2 << "; max lambdanum = " << maxnum2 << "; "
            << "new code " << time3 << "sec; old code " << time4 << "sec" << endl;
#endif
    
#if 0
    cout << "Testing modified code for TensorIndex(j,e,k)" << endl;
    
    Index temp_ind1;
    //for (unsigned int b=0; b<TBasisArray.size(); ++b)
    for (unsigned int b=3; b<4; ++b)
    {
        cout << "Basis b = " << b << endl;
        for (unsigned int  i=0; i<TBasisArray[b]->degrees_of_freedom(); ++i)
        //for (unsigned int  i=6656; i<6670; ++i)
        //for (unsigned int  i=0; i<2000; ++i)
        {
            //cout << "j0 = " << TBasisArray[b]->j0() << "; bases()[d]->j0() = [" << TBasisArray[b]->bases()[0]->j0();
            //for (unsigned int d = 1; d < dim; ++d)
            //    cout << ", " << TBasisArray[b]->bases()[d]->j0();
            //cout << "]" << endl;
            temp_ind1 = TBasisArray[b]->get_wavelet(i);
            //cout << "i = " << i << "; temp_ind1 = " << temp_ind1 << endl;
            Index temp_ind2(temp_ind1.j(), temp_ind1.e(), temp_ind1.k(),TBasisArray[b]);
            //cout << "temp_ind2 = " << temp_ind2 << endl;
            Index temp_ind3(temp_ind1.j(), temp_ind1.e(), temp_ind1.k(),TBasisArray[b],0);
            if (   ( !( (temp_ind1 == temp_ind2) && (temp_ind2 == temp_ind3)) )
                || (! (((temp_ind1.number() == i) && (temp_ind2.number() == i)) && (temp_ind3.number() == i)) )
               )
            {
                cout << "i = " << i << "; temp_ind1 = " << temp_ind1 << "; .number() = " << temp_ind1.number() << "; temp_ind2 = " << temp_ind2 << "; .number() = " << temp_ind2.number() << "; temp_ind3 = " << temp_ind3 << "; .number() = " << temp_ind3.number() << endl;
                abort();
            }
            //cout << "i = " << i << "; temp_ind1 = " << temp_ind1 << "; .number() = " << temp_ind1.number() << "; temp_ind2 = " << temp_ind2 << "; .number() = " << temp_ind2.number() << "; temp_ind3 = " << temp_ind3 << "; .number() = " << temp_ind3.number() << endl;
            assert ( (temp_ind1 == temp_ind2) && (temp_ind2 == temp_ind3));
            assert (((temp_ind1.number() == i) && (temp_ind2.number() == i)) && (temp_ind3.number() == i));
        }
    }

    cout << "no differences detected" << endl;
    cout << "checking speedup" << endl;
    int repetitions(10);
    
    Index lambda_run2;
    clock_t tstart, tend;
    double time1(0), time2(0);
    
    cout << "testing old code" << endl;
    tstart = clock();
    for (unsigned int r=0; r < repetitions; ++r)
    {
        for (unsigned int b=0; b<TBasisArray.size(); ++b)
        {
            //cout << "Basis b = " << b << endl;
            for (unsigned int  i=0; i<TBasisArray[b]->degrees_of_freedom(); ++i)
            {
                temp_ind1 = TBasisArray[b]->get_wavelet(i);
                Index temp_ind2(temp_ind1.j(), temp_ind1.e(), temp_ind1.k(),TBasisArray[b]);
            }
        }
    }
    tend = clock();
    time1 = (double)(tend-tstart)/CLOCKS_PER_SEC;
    
    cout << "testing new code" << endl;
    tstart = clock();
    for (unsigned int r=0; r < repetitions; ++r)
    {
        for (unsigned int b=0; b<TBasisArray.size(); ++b)
        {
            //cout << "Basis b = " << b << endl;
            for (unsigned int  i=0; i<TBasisArray[b]->degrees_of_freedom(); ++i)
            {
                temp_ind1 = TBasisArray[b]->get_wavelet(i);
                Index temp_ind2(temp_ind1.j(), temp_ind1.e(), temp_ind1.k(),TBasisArray[b],0);
            }
        }
    }
    tend = clock();
    time2 = (double)(tend-tstart)/CLOCKS_PER_SEC;
    
    cout << "Repetitions = " << repetitions << "; old code " << time1 << "sec; new code " << time2 << "sec" << "; time1/time2 = " << (time1/time2) << endl;
#endif
    return 0;
}
