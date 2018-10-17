#include <cmath>
#include <time.h>
#include <iostream>
#include "tframe_index.h"

// implementation for tframe.h

namespace WaveletTL
{

  	template <class IFRAME, unsigned int DIM>
  	TensorFrame<IFRAME,DIM>::TensorFrame(const FixedArray1D<IFRAME*,DIM> frames) {
    	for (unsigned int i = 0; i < DIM; i++) {
      		frames_[i] = frames[i];
      		j0_[i] = frames_[i]->j0();
    	}
    	delete_pointers = false;
        setup_full_collection_ = false;
        precomputed_supports_ = false;
  	}

  	template <class IFRAME, unsigned int DIM>
  	TensorFrame<IFRAME,DIM>::~TensorFrame()
  	{
	    if (delete_pointers) {
			for (typename list<IFRAME*>::const_iterator it(frames_infact.begin()); it != frames_infact.end(); ++it)
				delete *it;
            }
  	}
        
        template <class IFRAME, unsigned int DIM>
        const int
        TensorFrame<IFRAME,DIM>::Deltasize(const int j) const {
            const unsigned int Deltaj = frames_[0]->Deltasize(j);
            return pow((Deltaj),DIM);
        }
        
        template <class IFRAME, unsigned int DIM>
        typename TensorFrame<IFRAME,DIM>::Index
        TensorFrame<IFRAME,DIM>::first_generator(const level_type& j, const polynomial_type& p) const
        {
            assert(j >= j0_);
            typename Index::type_type e;
            typename Index::translation_type k;
            for (unsigned int i = 0; i < DIM; i++) {
                e[i]=0;
                k[i]=frames_[i]->DeltaLmin();   //todo: +1?
            }
            return Index(p, j0_, e, k, p.number() * Nablasize_, this);
                          
        }
        
         template <class IFRAME, unsigned int DIM>
        typename TensorFrame<IFRAME,DIM>::Index
        TensorFrame<IFRAME,DIM>::last_generator(const level_type& j, const polynomial_type& p) const
        {
            if(DIM==2){
                assert(j >= j0_);
                typename Index::type_type e;
                // setup highest translation index for e=(0,0), p=4
                typename Index::translation_type k(frames_[0]->DeltaRmax(j[0]), frames_[1]->DeltaRmax(j[1]));
                return Index(p,j, e, k, p.number()* Nablasize_+Deltasize(j[0])-1, this); 
            }
//            else{
//            assert(j >= j0_);
//            typename Index::type_type e;
//            typename Index::translation_type k;
//            for (unsigned int i = 0; i < DIM; i++)
//            {
//                e[i] = 0;
//                k[i] = frames_[i]->DeltaRmax(j0_[i]);
//            }
//            int res=1;
//            for (unsigned int i = 0; i < DIM; i++)
//                res *= frames_[i]->Deltasize(j0_[i]);
//            return Index(p, j0_, e, k, p.number() * (last_quarklet(jmax_).number()+1)+(res-1), this);
//            }
        }
         
         template <class IFRAME, unsigned int DIM>
        typename TensorFrame<IFRAME,DIM>::Index
        TensorFrame<IFRAME,DIM>::first_quarklet(const level_type& j, const polynomial_type& p, const int& number) const
        {
             if(DIM==2){
                 assert(j >= j0());

                typename Index::type_type e;
                typename Index::translation_type k;
    
                bool sofar_only_generators = true;
    
                if (j[0] == j0_[0])
                {
                    e[0] = 0;
                    k[0] = frames_[0]->DeltaLmin();
                } else
                {
                    e[0] = 1;
                    k[0] = frames_[0]->Nablamin();
                    sofar_only_generators = false;
                }
    
                if ( (sofar_only_generators == true) || (j[1] != j0_[0]) )
                {
                    e[1] = 1;
                    k[1] = frames_[1]->Nablamin();
                    //sofar_only_generators = false;
                } else
                {
                    e[1] = 0;
                    k[1] = frames_[1]->DeltaLmin();
                }
    
                if (number==-1)
                {
                    level_type jdiff;
                    jdiff[0]= j[0]-j0_[0], jdiff[1]= j[1]-j0_[1];
                    int level = jdiff.number();
                    int altnumber = p.number()* Nablasize_+first_wavelet_numbers[level];
                    return Index(p,j, e, k, altnumber, this);
                }
                else
                return Index(p,j, e, k, number, this);
             }
//            MultiIndex<int,DIM> temp_mi(j);
//            for (unsigned int i=0; i<DIM; ++i)
//                temp_mi[i] -= j0_[i];
//            if(p.number()==0)
//            return first_wavelets[temp_mi.number()];
//            else{
//                MultiIndex<int,DIM> tempj,tempe,tempk;
//                tempj = first_wavelets[temp_mi.number()].j();
//                tempe = first_wavelets[temp_mi.number()].e();
//                tempk = first_wavelets[temp_mi.number()].k();
//                return Index(p,tempj,tempe,tempk, p.number() * (last_quarklet(jmax_).number()+1) 
//                        + first_wavelets[temp_mi.number()].number(), this);
//            }
        }
        
//        template <class IFRAME, unsigned int DIM>
//        typename TensorFrame<IFRAME,DIM>::Index
//        TensorFrame<IFRAME,DIM>::first_quarklet(const int levelsum, const typename Index::polynomial_type p) const
//        {            
//            MultiIndex<int,DIM> temp_mi;
//            temp_mi[0] = levelsum - multi_degree(j0_);
//            if(p.number()==0)
//            return first_wavelets[temp_mi.number()];
//            else{
//                MultiIndex<int,DIM> tempj,tempe,tempk;
//                tempj = first_wavelets[temp_mi.number()].j();
//                tempe = first_wavelets[temp_mi.number()].e();
//                tempk = first_wavelets[temp_mi.number()].k();
//                return Index(p,tempj,tempe,tempk, p.number() * (last_quarklet(jmax_).number()+1) 
//                        + first_wavelets[temp_mi.number()].number()-1, this);
//            }
//        }
        
        template <class IFRAME, unsigned int DIM>
        typename TensorFrame<IFRAME,DIM>::Index
        TensorFrame<IFRAME,DIM>::last_quarklet(const level_type& j, const polynomial_type& p, const int& number) const
        {   
//            cout<<"number:"<<number<<endl;
            if(DIM==2){
               assert(j >= j0());
    
               typename Index::type_type e(1, 1);

                // setup highest translation index for e=(1,1), p=2
                typename Index::translation_type k(frames_[0]->Nablamax(j[0]), frames_[1]->Nablamax(j[1]));
//              cout << "in last_quarklet level_type" << endl;
                if (number==-1)
                {
                    level_type jdiff;
                    jdiff[0]= j[0]-j0_[0], jdiff[1]= j[1]-j0_[1];
                    int level = jdiff.number();
                    int altnumber = p.number()* Nablasize_+last_wavelet_numbers[level];
                    return Index(p,j, e, k, altnumber, this);
                }
                else
    
                return Index(p,j, e, k, number, this); 
            }
//            MultiIndex<int,DIM> temp_mi(j);
//            for (unsigned int i=0; i<DIM; ++i)
//                temp_mi[i] -= j0_[i];
//            if(p.number()==0)
//            return last_wavelets[temp_mi.number()];
//            else{
//                MultiIndex<int,DIM> tempj,tempe,tempk;
//                tempj = last_wavelets[temp_mi.number()].j();
//                tempe = last_wavelets[temp_mi.number()].e();
//                tempk = last_wavelets[temp_mi.number()].k();
//                return Index(p,tempj,tempe,tempk, (p.number()+1) * (last_quarklet(jmax_).number()+1) - 1,this);
//            }
        }
        
        template <class IFRAME, unsigned int DIM>
        typename TensorFrame<IFRAME,DIM>::Index
        TensorFrame<IFRAME,DIM>::last_quarklet(const int& levelsum, const polynomial_type& p, const int& number) const
        {
            if(DIM==2){
                assert(levelsum >= (int) multi_degree(j0_));
    
                typename Index::type_type e(1, 1);
                typename Index::level_type j(levelsum - j0_[0],  j0_[0]);

                // setup highest translation index for e=(1,1), p=2
                typename Index::translation_type k(frames_[0]->Nablamax(j[0]), frames_[1]->Nablamax(j[1]));
//                cout << "in last_quarklet" << endl;
                if (number==-1)
                {
//        cout << "altnumber" << endl;
//        cout << "j: " << j << ", j.number" << j.number() << endl;
//        cout << "j0_: " << j0_ <<", j0_.number" << j0_.number() << endl;
                    level_type jdiff;
                    jdiff[0]= j[0]-j0_[0], jdiff[1]= j[1]-j0_[1];
                    int level = jdiff.number();
//        cout << "level: " << level << endl;
                    int altnumber = p.number()* Nablasize_+last_wavelet_numbers[level];
                    return Index(p,j, e, k, altnumber, this);
                }
                else{
//        cout << "number" << endl;
                return Index(p,j, e, k, number, this); 
                }
            }
//            MultiIndex<int,DIM> temp_mi;
//            temp_mi[0] = levelsum - multi_degree(j0_);
//            if(p.number()==0)
//            return last_wavelets[temp_mi.number()];
//            else{
//                
//                MultiIndex<int,DIM> tempj,tempe,tempk;
//                tempj = last_wavelets[temp_mi.number()].j();
//                tempe = last_wavelets[temp_mi.number()].e();
//                tempk = last_wavelets[temp_mi.number()].k();
//                return Index(p,tempj,tempe,tempk, (p.number()+1) * (last_quarklet(jmax_).number()+1) - 1,this);
//            }
                
        }
        
#if 0               
        template <class IFRAME, unsigned int DIM>
        void
        TensorFrame<IFRAME,DIM>::precompute_firstlast_wavelets()
        {
            MultiIndex<int,DIM> level_it;
            level_it[0] = jmax_ - multi_degree(j0_);
            const int numoflevels = level_it.number() + 1;
            level_it[0] = 0;
            first_wavelets.resize(numoflevels);
            last_wavelets.resize(numoflevels);
            typename Index::polynomial_type temp_p;
            typename Index::level_type temp_j(j0_);
            typename Index::type_type temp_e;
            typename Index::translation_type temp_k;
            
            
            
            // first_level
            unsigned int l = 0;
            if (DIM == 1)
            {
                temp_p[0] = 0;
                temp_e[0] = 1;
                temp_k[0] = frames_[0]->Nablamin();
            }
            else
            {
                for (unsigned int i = 0; i < DIM-1; i++) 
                {
                    temp_p[i] = 0;
                    temp_e[i] = 0;
                    temp_k[i] = frames_[i]->DeltaLmin();
                }
                temp_p[DIM-1] = 0;
                temp_e[DIM-1] = 1;
                temp_k[DIM-1] = frames_[DIM-1]->Nablamin();
            }
            first_wavelets[0] = Index(temp_p, temp_j, temp_e, temp_k, this);
            for (unsigned int i = 0; i < DIM; i++) 
                {
                    temp_p[i] = 0;
                    temp_e[i] = 1;
                    temp_k[i] = frames_[i]->Nablamax(temp_j[i]);
                }
            last_wavelets[0] = Index(temp_p, temp_j, temp_e, temp_k, this);
            ++level_it;
            ++l;
            while ((int)l < numoflevels)
            {
                //cout << "level_it == " << level_it << "; l == " << l <<  endl;
                for (unsigned int i = 0; i < DIM; i++) 
                {
                    temp_j[i] = j0_[i] + level_it[i];
                    if (level_it[i] == 0)
                    {
                        temp_p[i] = 0;
                        temp_e[i] = 0;
                        temp_k[i] = frames_[i]->DeltaLmin();
                    } else
                    {
                        temp_p[i] = 0;
                        temp_e[i] = 1;
                        temp_k[i] = frames_[i]->Nablamin();
                    }
                }
//                cout << "temp_p == " << temp_p << "; temp_j == " << temp_j << "; temp_e == " << temp_e << "; temp_k == " << temp_k << endl;
                first_wavelets[l] = Index(temp_p, temp_j, temp_e, temp_k, this);
                for (unsigned int i = 0; i < DIM; i++) 
                {
                    temp_p[i] = 0;
                    temp_e[i] = 1;
                    temp_k[i] = frames_[i]->Nablamax(temp_j[i]);
                }
//                cout << "temp_p == " << temp_p << "; temp_j == " << temp_j << "; temp_e == " << temp_e << "; temp_k == " << temp_k << endl;
                last_wavelets[l] = Index(temp_p, temp_j, temp_e, temp_k, this);
                ++level_it;
                ++l;
            }
        }
#endif
        template <class IFRAME, unsigned int DIM>
        void
        TensorFrame<IFRAME,DIM>::support(const int& lambda_num, Support& supp) const
        {

            if (precomputed_supports_){
                supp=all_supports_[lambda_num];
            }
            else {
                cout << "TensorFrame<IFRAME,DIM>::support(): precompute supports first!" << endl;
                abort();      
            }
        }

	template <class IFRAME, unsigned int DIM>
  	void
  	TensorFrame<IFRAME,DIM>::support(const Index& lambda, Support& supp) const
  	{
        if (precomputed_supports_){
         supp=all_supports_[lambda.number()];
        }
        else{    
            // PERFORMANCE?
            typename Index::level_type temp_j = lambda.j();
            typename Index::type_type temp_e = lambda.e();
            //MultiIndex<int,DIM> temp_j=lambda.j(),temp_e=lambda.e();
            for (unsigned int i(0); i < DIM; i++) {
                supp.j[i] = temp_j[i]+temp_e[i];
                //supp.j[i] = lambda.j()[i] + lambda.e()[i];
                frames_[i]->support(temp_j[i],
                                    temp_e[i],
                                    lambda.k()[i],
                                    supp.a[i],
                                    supp.b[i]);
            }
        }    
  	}
        
        template <class IFRAME, unsigned int DIM>
        SampledMapping<DIM>
        TensorFrame<IFRAME,DIM>::evaluate
            (const typename TensorFrame<IFRAME,DIM>::Index& lambda,
            const int resolution) const
        {
        SampledMapping<DIM> r;

        typedef typename TensorFrame<IFRAME,DIM>::Index Index;

        typename Index::type_type zero;
     
        Point<DIM> x;
        Point<DIM> y;
        FixedArray1D<Array1D<double>,DIM> values;

        for(int i=0;i<DIM;i++){
            x[i]=0;
            y[i]=1;
            values[i].resize((1<<resolution)+1); 
            values[i] = frames_[i]->evaluate(typename IFRAME::Index(lambda.p()[i],
                                                              lambda.j()[i],
							      lambda.e()[i],
							      lambda.k()[i],
							      frames_[i]),
				       resolution).values();
        }
 	r = SampledMapping<DIM>(x, y, values);
        return r;
        }
        
        template <class IFRAME, unsigned int DIM>
        SampledMapping<DIM>
        TensorFrame<IFRAME,DIM>::evaluate
            (const InfiniteVector<double, typename TensorFrame<IFRAME,DIM>::Index>& coeffs,
            const int resolution) const
        {
            // first prepare a plot of the zero function
            Point<DIM> x;
            Point<DIM> y;
            FixedArray1D<Array1D<double>,DIM> values;
            for(int i=0;i<DIM;i++){
                x[i]=0;
                y[i]=1;
                values[i].resize((1<<resolution)+1);
                for (int j = 0; j <= 1<<resolution; j++) {
                    values[i][j]  = 0;
                }
            }
            SampledMapping<DIM>  result;
            result = SampledMapping<DIM>(x,y, values);
          
            // add all plots of the single functions
            typedef typename TensorFrame<IFRAME,DIM>::Index Index;
            for (typename InfiniteVector<double,Index>::const_iterator it(coeffs.begin()),
            itend(coeffs.end()); it != itend; ++it) {
            SampledMapping<DIM>  temp(evaluate(it.index(), resolution));
            result.add(*it, temp);
            }
    
            return result;
        }
        
        template <class IFRAME, unsigned int DIM>
        void
        TensorFrame<IFRAME,DIM>::setup_full_collection()
        {
        if(DIM==2){
            if (jmax_ == -1 || jmax_ < (int) multi_degree(j0())) {
                cout << "TensorFrame<IFRAME,DIM>::setup_full_collection(): specify a maximal level of resolution first!" << endl;
                abort();
            }   
    
            cout << "setting up collection of quarklet indices..." << endl;

            MultiIndex<int,DIM> level_it;
            level_it[0] = jmax_ - multi_degree(j0());
            const int numoflevels = level_it.number() + 1;
            level_it[0] = 0;
            first_wavelet_numbers.resize(numoflevels);
            last_wavelet_numbers.resize(numoflevels);
    
            typename Index::polynomial_type p(0,0), pmax(pmax_,0), pmin(0,0);
            set<Index> Lambda;
            int waveletlevel=0;
            Index ind = first_generator(j0_, p);
            for (int k = 0;; k++) {
//        full_collection[k] = ind;
                Lambda.insert(ind);
//                cout << ind << ", " << endl;
                if(ind==first_quarklet(ind.j(),pmin,k)){
                    first_wavelet_numbers[waveletlevel]=k;
                
//                cout << "first: " << k << endl;
                }
                if(ind==last_quarklet(ind.j(),pmin,k)){
                    last_wavelet_numbers[waveletlevel]=k;
//            cout << "last: " << k << endl;
                    waveletlevel++;
                }
                if(ind==last_quarklet(jmax_, p,k)){
                    if(multi_degree(p)==0){
//                cout << "Nablasize wird gesetzt auf " << k+1 << endl;
                        Nablasize_=k+1;
                    }
                   if(p==pmax){
//                       cout<<"bin hier"<<endl;
//                       cout<<ind<<endl;
                       break;
                   }    
                    ++p;
                    ind=first_generator(j0_, p);
                }
                    else
                        ++ind;
            }
    
            full_collection.resize(Lambda.size());
            all_supports_.resize(Lambda.size());
            int i = 0;
            for (typename set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it, i++){
                full_collection[i] = *it;
                Support supp;
                support(*it,supp);
                all_supports_[i] = supp;
        
//        cout << *it << ", " << (*it).number() << endl;
            }
            precomputed_supports_ = true;
    
    
    
//    cout << "Nablasize in setup: " << Nablasize_ << endl;
    
            cout << "done setting up collection of quarklet indices and precompute supports..." << endl;
            cout << "total degrees of freedom between j0_ = " << j0_ << " and (jmax_= " << jmax_ << ", pmax_= " << pmax_ << ") is " << full_collection.size() << endl;
            setup_full_collection_ = true;
        }    
        }
}
        
        
  
//        template <class IFRAME, unsigned int DIM>
//        typename TensorFrame<IFRAME,DIM>::Index
//        TensorFrame<IFRAME,DIM>::first_generator(const typename Index::polynomial_type p) const
//        {
//            typename Index::type_type e;
//            typename Index::translation_type k;
//            for (unsigned int i = 0; i < DIM; i++) {
//                e[i]=0;
//                k[i]=frames_[i]->DeltaLmin();
//            }
//             return Index(p, j0_, e, k, p.number() * (last_quarklet(jmax_).number()+1), this);
//        }

//        template <class IFRAME, unsigned int DIM>
//        typename TensorFrame<IFRAME,DIM>::Index
//        TensorFrame<IFRAME,DIM>::first_generator(const unsigned int j, const typename Index::polynomial_type p) const
//        {
//            typename Index::type_type e;
//            typename Index::translation_type k;
//            for (unsigned int i = 0; i < DIM; i++) {
//                e[i]=0;
//                k[i]=frames_[i]->DeltaLmin();
//            }
//            return Index(p, j0_, e, k, p.number() * (last_quarklet(jmax_).number()+1),  this);
//        }

        

       

//#if _PRECOMPUTE_FIRSTLAST_WAVELETS
               
        
#if 0
        
        template <class IFRAME, unsigned int DIM>
        typename TensorFrame<IFRAME,DIM>::Index
        TensorFrame<IFRAME,DIM>::first_wavelet(const MultiIndex<int,DIM> j) const
        {
#if _TFRAME_DEBUGLEVEL_ >= 1
            assert(multi_degree(j) > multi_degree(j0_)
                    || (multi_degree(j) == multi_degree(j0_) && j >= j0_)
                  );
#endif
            typename Index::type_type e;
            typename Index::translation_type k;
            bool first_level = true;
            for (unsigned int i = 0; i < DIM; i++) {
                if (j[i] == j0_[i])
                {
                    p[i] = 0;
                    e[i] = 0;
                    k[i] = frames_[i]->DeltaLmin();
                } else
                {
                    p[i] = 0;
                    e[i] = 1;
                    k[i] = frames_[i]->Nablamin();
                    first_level = false;
                }
            }
            if (first_level == true)
            {
                p[DIM-1] = 0;
                e[DIM-1] = 1;
                k[DIM-1] = frames_[DIM-1]->Nablamin();
            }
            return Index(p, j, e, k, this);
        }

        template <class IFRAME, unsigned int DIM>
        typename TensorFrame<IFRAME,DIM>::Index
        TensorFrame<IFRAME,DIM>::first_wavelet(const int levelsum) const
        {            
#if _TFRAME_DEBUGLEVEL_ >= 1
            assert(levelsum >= multi_degree(j0_) );
#endif
            typename Index::level_type j;
            typename Index::type_type e;
            typename Index::translation_type k;

            p[DIM-1] = 0;
            j[DIM-1] = j0_[DIM-1] + levelsum - multi_degree(j0_);
            e[DIM-1] = 1;
            k[DIM-1] = frames_[DIM-1]->Nablamin();

            for (unsigned int i = 0; i < DIM-1; i++)
            {
                p[i] = 0;
                j[i] = j0_[i];
                e[i] = 0;
                k[i] = frames_[i]->DeltaLmin();
            }
            return Index(p, j, e, k, this);
        }

        template <class IFRAME, unsigned int DIM>
        typename TensorFrame<IFRAME,DIM>::Index
        TensorFrame<IFRAME,DIM>::last_wavelet(const MultiIndex<int,DIM> j) const
        {            
#if _TFRAME_DEBUGLEVEL_ >= 1
            assert(multi_degree(j) > multi_degree(j0_)
                    || ((multi_degree(j) == multi_degree(j0_)) && (j0_ <= j))
                  );
#endif
            typename Index::type_type e;
            typename Index::translation_type k;
            for (unsigned int i = 0; i < DIM; i++) {
                p[i] = 0;
                e[i] = 1;
                k[i] = frames_[i]->Nablamax(j[i]);
            }
            return Index(p, j, e, k, this);
        }

        template <class IFRAME, unsigned int DIM>
        typename TensorFrame<IFRAME,DIM>::Index
        TensorFrame<IFRAME,DIM>::last_wavelet(const int levelsum) const
        {
#if _TFRAME_DEBUGLEVEL_ >= 1
            assert(levelsum >= multi_degree(j0_));
#endif
            typename Index::polynomial_type p;
            typename Index::level_type j;
            typename Index::type_type e;
            typename Index::translation_type k;
            p[0] = 0;
            j[0]= j0_[0]+levelsum-multi_degree(j0_);
            e[0]=1;
            k[0]= frames_[0]->Nablamax(j[0]);
            for (unsigned int i = 1; i < DIM; i++) {
                p[i] = 0;
                j[i] = j0_[i];
                e[i] = 1;
                k[i] = frames_[i]->Nablamax(j0_[i]);
            }
            return Index(p, j, e, k, this);
        }
        
#endif // if _PRECOMPUTE_FIRSTLAST_WAVELETS
#if 0      
	template <class IFRAME, unsigned int DIM>
	void
	TensorFrame<IFRAME,DIM>::expand(const Function<DIM>* f,
                                        const bool primal,
                                        const MultiIndex<int,DIM> jmax,
                                        InfiniteVector<double,Index>& coeffs) const
  	{
            assert(primal == false); // only integrate against primal wavelets and generators
            for (unsigned int i=0; (int)i<degrees_of_freedom(); i++)
            {
                const double coeff = integrate(f, full_collection[i]);
                if (fabs(coeff)>1e-15)
                    coeffs.set_coefficient(full_collection[i], coeff);
            }
            /*
            for (Index lambda = first_generator(), lambda_end(last_wavelet(jmax));;++lambda)
            {
                const double coeff = integrate(f, lambda);
                if (fabs(coeff)>1e-15)
                    coeffs.set_coefficient(lambda, coeff);
                if (lambda == lambda_end)
                    break;
            }
            */
  	}

	template <class IFRAME, unsigned int DIM>
	void
	TensorFrame<IFRAME,DIM>::expand(const Function<DIM>* f,
                                        const bool primal,
                                        const unsigned int jmax,
                                        InfiniteVector<double,Index>& coeffs) const
  	{
            assert(primal == false); // only integrate against primal wavelets and generators
            /*
            for (Index lambda = first_generator(), lambda_end(last_wavelet(jmax));;++lambda)
            {
                const double coeff = integrate(f, lambda);
                if (fabs(coeff)>1e-15)
                    coeffs.set_coefficient(lambda, coeff);
                if (lambda == lambda_end)
                    break;
            }
            */
            for (unsigned int i=0; (int)i<degrees_of_freedom(); i++)
            {
                const double coeff = integrate(f, full_collection[i]);
                if (fabs(coeff)>1e-15)
                    coeffs.set_coefficient(full_collection[i], coeff);
            }
            
  	}
        
	template <class IFRAME, unsigned int DIM>
	double
	TensorFrame<IFRAME,DIM>::integrate(const Function<DIM>* f,
                                           const Index& lambda) const
  	{
            // f(v) = \int_0^1 g(t)v(t) dt
            double r = 0;
            // first compute supp(psi_lambda)
            Support supp;
            support(lambda, supp);
            // setup Gauss points and weights for a composite quadrature formula:
            const int N_Gauss = 5;
            FixedArray1D<double,DIM> h;
            for (unsigned int i=0; i<DIM;i++)
                h[i]=ldexp(1.0, -supp.j[i]); // granularity for the quadrature
                //const double h = ldexp(1.0, -supp.j); // granularity for the quadrature
            FixedArray1D<Array1D<double>,DIM> gauss_points, gauss_weights, v_values;
            for (unsigned int i = 0; i < DIM; i++) {
                gauss_points[i].resize(N_Gauss*(supp.b[i]-supp.a[i]));
                gauss_weights[i].resize(N_Gauss*(supp.b[i]-supp.a[i]));
                for (int patch = supp.a[i]; patch < supp.b[i]; patch++)
                    for (int n = 0; n < N_Gauss; n++) {
                        gauss_points[i][(patch-supp.a[i])*N_Gauss+n]
                                = h[i]*(2*patch+1+GaussPoints[N_Gauss-1][n])/2.;
                        gauss_weights[i][(patch-supp.a[i])*N_Gauss+n]
                                = h[i]*GaussWeights[N_Gauss-1][n];
                    }
            }
            // compute the point values of the integrand (where we use that it is a tensor product)
            for (unsigned int i = 0; i < DIM; i++)
                frames()[i]->evaluate(0,
                                     lambda.p()[i],
                                     lambda.j()[i],
                                     lambda.e()[i],
                                     lambda.k()[i],
                                     gauss_points[i],
                                     v_values[i]);
            // iterate over all points and sum up the integral shares
            int index[DIM]; // current multiindex for the point values
            for (unsigned int i = 0; i < DIM; i++)
                index[i] = 0;
            Point<DIM> x;
            while (true) {
                for (unsigned int i = 0; i < DIM; i++)
                    x[i] = gauss_points[i][index[i]];
                double share = f->value(x);
                for (unsigned int i = 0; i < DIM; i++)
                    share *= gauss_weights[i][index[i]] * v_values[i][index[i]];
                r += share;
                // "++index"
                bool exit = false;
                for (unsigned int i = 0; i < DIM; i++) {
                    if (index[i] == N_Gauss*(supp.b[i]-supp.a[i])-1) {
                        index[i] = 0;
                        exit = (i == DIM-1);
                    } else {
                        index[i]++;
                        break;
                    }
                }
                if (exit) break;
            }
            return r;
        }
#endif
//        template <class IFRAME, unsigned int DIM>
//  	double
//  	TensorFrame<IFRAME,DIM>::evaluate(const unsigned int derivative, 
//                                          const Index& lambda,
//                                          const Point<DIM> x) const
//  	{
//            double value = 1.0;
//            for (unsigned int i = 0; i < DIM; i++) // loop through components of the tensor product
//                value *= frames_[i]->evaluate(derivative,
//                                             lambda.p()[i], lambda.j()[i], lambda.e()[i], lambda.k()[i],
//                                             x[i]);
//            return value;
//        }


        
//            if (jmax_ < multi_degree(j0_) ) {
//                cout << "TensorFrame<IFRAME,DIM>::setup_full_collection(): the specified maximal level jmax is invalid. Specify a higher maximal level jmax_" << endl;
//                cout << "jmax_ = " << jmax_ << "; j0_ = " << j0_ << endl;
//                abort();
//            }
//            typedef typename Index::polynomial_type polynomial_type;
//            typedef  std::set<polynomial_type> pol_set;
//            pol_set temp(degree_indices<int,DIM>(pmax_));
//            typename pol_set::const_iterator itend(temp.end());
//            --itend;
//            int degrees_of_freedom = ((*itend).number()+1)*(last_quarklet_num<IFRAME,DIM,TensorFrame<IFRAME,DIM> >(this, jmax_) +1); // +1 since numbering begins at 0
//            cout << "total degrees of freedom between j0_ = " << j0_ << " and (jmax_= " << jmax_ << ", pmax_= " << pmax_ << ") is " << degrees_of_freedom << endl;
//            cout << "setting up collection of quarklet indices..." << endl;
//            full_collection.resize(degrees_of_freedom);
//            all_supports_.resize(degrees_of_freedom);
//            typename Index::polynomial_type p;
//            Index ind = first_generator(j0_, p);
//            for (int k = 0; k < degrees_of_freedom; k++) {
////                cout << ind << ", " << endl;
//                full_collection[k] = ind;
//                Support supp;
//                support(ind,supp);
//                all_supports_[k] = supp;
//                
//                if(ind==last_quarklet(jmax_, p)){
//                    if(multi_degree(p)==0){
////                cout << "Nablasize wird gesetzt auf " << k+1 << endl;
//                    Nablasize_=k+1;
//                    }
////                    if(p==pmax)
////                        break;
//                    ++p;
//                    ind=first_generator(j0_, p);
//                }
//                    else
//                        ++ind;
//            }
//            setup_full_collection_ = true;
//            precomputed_supports_ = true;
//            cout << "done setting up collection of quarklet indices..." << endl;
        


