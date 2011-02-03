// implementation for tbasis_support.h

#include <utils/multiindex.h>
#include <utils/fixed_array1d.h>
#include <iostream>

using MathTL::multi_degree;
using MathTL::FixedArray1D;

namespace WaveletTL
{
    template <class IBASIS, unsigned int DIM>
    inline
    void
    support(const TensorBasis<IBASIS,DIM>& basis,
            const typename TensorBasis<IBASIS,DIM>::Index& lambda,
            typename TensorBasis<IBASIS,DIM>::Support& supp)
    {
        basis.support(lambda, supp);
    }

    template <class IBASIS, unsigned int DIM>
    bool
    intersect_supports(const TensorBasis<IBASIS,DIM>& basis,
                       const typename TensorBasis<IBASIS,DIM>::Index& lambda,
                       const typename TensorBasis<IBASIS,DIM>::Index& mu,
                       typename TensorBasis<IBASIS,DIM>::Support& supp)
    {
        typename TensorBasis<IBASIS,DIM>::Support supp_lambda;
        WaveletTL::support<IBASIS,DIM>(basis, lambda, supp_lambda);
        typename TensorBasis<IBASIS,DIM>::Support supp_mu;
        WaveletTL::support<IBASIS,DIM>(basis, mu, supp_mu);
        // determine support intersection granularity,
        // adjust single support granularities if necessary
        for (unsigned int i=0;i<DIM;i++)
        {
            supp.j[i] = std::max(supp_lambda.j[i], supp_mu.j[i]);
            if (supp_lambda.j[i] > supp_mu.j[i]) {
                const int adjust = 1<<(supp_lambda.j[i]-supp_mu.j[i]);
                supp_mu.a[i] *= adjust;
                supp_mu.b[i] *= adjust;
            } else {
                const int adjust = 1<<(supp_mu.j[i]-supp_lambda.j[i]);
                supp_lambda.a[i] *= adjust;
                supp_lambda.b[i] *= adjust;
            }
            supp.a[i] = std::max(supp_lambda.a[i],supp_mu.a[i]);
            supp.b[i] = std::min(supp_lambda.b[i],supp_mu.b[i]);
            if (supp.a[i] >= supp.b[i])
                return false;
        }
        return true;
    }


    template <class IBASIS, unsigned int DIM>
    void intersecting_wavelets(const TensorBasis<IBASIS,DIM>& basis,
                               const typename TensorBasis<IBASIS,DIM>::Index& lambda,
                               const MultiIndex<int,DIM> j, const bool generators,
                               std::list<typename TensorBasis<IBASIS,DIM>::Index>& intersecting)
    {
        #if 0
        if (generators) assert (j==basis.j0());

        typedef typename TensorBasis<IBASIS,DIM>::Index Index;
        intersecting.clear();

        // the set of intersecting wavelets is a cartesian product from d sets from the 1D case,
        // so we only have to compute the relevant 1D indices
        typedef typename IBASIS::Index Index1D;
        FixedArray1D<std::list<Index1D>,DIM> intersecting_1d_generators, intersecting_1d_wavelets;
        // prepare all intersecting wavelets and generators in the i-th coordinate direction
        for (unsigned int i = 0; i < DIM; i++) {
            if ((j[i]==basis.j0()[i]) && (generators || DIM > 1) )
                intersecting_wavelets(*basis.bases()[i],
                                      Index1D(lambda.j()[i],
                                              lambda.e()[i],
                                              lambda.k()[i],
                                              basis.bases()[i]),
                                      j[i], true, intersecting_1d_generators[i]);
            if (!generators)
                intersecting_wavelets(*basis.bases()[i],
                                      Index1D(lambda.j()[i],
                                              lambda.e()[i],
                                              lambda.k()[i],
                                              basis.bases()[i]),
                                      j[i], false, intersecting_1d_wavelets[i]);
        }
        // generate all relevant tensor product indices with either e=(0,...,0) or e!=(0,...,0)

        // Version 1: unsortierter Output
        // PERFORMANCE: Reihenfolge entspricht nicht der von tbasis_index. 
        // Problem? Ja! Reihenfolge wird in apply(window,y,res) ausgenutzt,
        // um heraus zu finden, ob der Zeileniterator erhöht werden muss
        typedef std::list<FixedArray1D<Index1D,DIM> > list_type;
        list_type indices;
        FixedArray1D<Index1D,DIM> helpindex;
        if ((j[0] == basis.j0()[0]) && (generators || DIM > 1) ) {
            for (typename std::list<Index1D>::const_iterator it(intersecting_1d_generators[0].begin()),
                    itend(intersecting_1d_generators[0].end());
                    it != itend; ++it) {
                        helpindex[0] = *it;
                        indices.push_back(helpindex);
                    }
        }
        if (!(generators)) {
            for (typename std::list<Index1D>::const_iterator it(intersecting_1d_wavelets[0].begin()),
                    itend(intersecting_1d_wavelets[0].end());
                    it != itend; ++it) {
                        helpindex[0] = *it;
                        indices.push_back(helpindex);
                    }
        }
        for (unsigned int i = 1; i < DIM; i++) {
            list_type sofar;
            sofar.swap(indices);
            for (typename list_type::const_iterator itready(sofar.begin()), itreadyend(sofar.end());
                    itready != itreadyend; ++itready) {
                        helpindex = *itready;
                        unsigned int esum = 0;
                        for (unsigned int k = 0; k < i; k++)
                            esum += helpindex[k].e();
                        if (generators || ( (j[i] == basis.j0()[i]) && ( i < DIM-1 || (i == (DIM-1) && esum > 0) ) ) ) {
                            for (typename std::list<Index1D>::const_iterator it(intersecting_1d_generators[i].begin()),
                                    itend(intersecting_1d_generators[i].end());
                                    it != itend; ++it) {
                                        helpindex[i] = *it;
                                        indices.push_back(helpindex);
                                    }
                        }
                        if (!(generators)) {
                            for (typename std::list<Index1D>::const_iterator it(intersecting_1d_wavelets[i].begin()),
                                    itend(intersecting_1d_wavelets[i].end());
                                    it != itend; ++it) {
                                        helpindex[i] = *it;
                                        indices.push_back(helpindex);
                                    }
                        }
                    }
        }
        // compose the results
        typename Index::type_type help_e;
        typename Index::translation_type help_k;
        for (typename list_type::const_iterator it(indices.begin()), itend(indices.end()); it != itend; ++it)
        {
            for (unsigned int i = 0; i < DIM; i++)
            {
                help_e[i] = (*it)[i].e();
                help_k[i] = (*it)[i].k();
            }
            intersecting.push_back(Index(j, help_e, help_k, &basis));
        }
#else
     if (generators) assert (j==basis.j0());

    typedef typename TensorBasis<IBASIS,DIM>::Index Index;
    intersecting.clear();
 
 //   if (! generators) {
      FixedArray1D<int,DIM>
      minkwavelet, maxkwavelet, minkgen, maxkgen;
      typedef typename IBASIS::Index Index1D;
      int minkkkk;
      int maxkkkk;

      // prepare all intersecting wavelets and generators in the i-th coordinate direction
      for (unsigned int i = 0; i < DIM; i++) {
        get_intersecting_wavelets_on_level(*(basis.bases()[i]),
	    	Index1D(lambda.j()[i],
		    lambda.e()[i],
		    lambda.k()[i],
		    basis.bases()[i]),
		    j[i], true, minkkkk,maxkkkk);
        minkgen[i] = minkkkk;
        maxkgen[i] = maxkkkk;
        if (!(generators))
	  get_intersecting_wavelets_on_level(*basis.bases()[i],
		      Index1D(lambda.j()[i],
			      lambda.e()[i],
			      lambda.k()[i],
			      basis.bases()[i]),
		      j[i], false, minkkkk,maxkkkk);
        minkwavelet[i] = minkkkk;
        maxkwavelet[i] = maxkkkk;
       } // end for

      unsigned int result = 0;
      int deltaresult = 0;
      int genfstlvl = 0;
      bool gen = 0;
      //const Array1D<Index>* full_collection = &basis.full_collection;

      MultiIndex<int,DIM> type;
      type[DIM-1] = 1;
      int tmp = 1;
      bool exit = 0;
      bool nureinmal = 0;
      bool doch = 0;

      // determine how many wavelets there are on all the levels
      // below the level of this index

      MathTL::FixedArray1D<std::map<int, int>,DIM> sizes; // store number of basis elements. Generators on level j0 (0), wavelets on level j0 (1), j0+1 (2), ...
                int uptothislevel(0); // number of basis function below the current level
		int oncurrentlevel(1); // number of base elements on current level j
		int range(0); // determines wether a new entry has to be computed in "sizes"
		MultiIndex<int,DIM> currentlevel, j0(basis.j0());
		MultiIndex<int,DIM> currenttype;		
		for (unsigned int i=0; i < DIM; i++) {
                    currentlevel[i]=j0[i];
                    currenttype[i]=0;
                    sizes[i][0] = basis.bases()[i]->Deltasize(j0[i]); // Number of generators on level j0
                    sizes[i][1] = basis.bases()[i]->Nablasize(j0[i]); // Number of Wavelets on level j0
		}
		// iterate over all level up to j_ and add up number of basis functions up to the level
                while(!exit){
                exit=1;
                doch = 0;

                if (!generators)
                while (true)
		{
                   // if(!nureinmal || doch){
                      // compute number of functions on this level
                      oncurrentlevel = 1;
                      for (unsigned int i = 0; i < DIM; i++)
                      {
                          oncurrentlevel *= sizes[i][currentlevel[i]+currenttype[i]-j0[i]];
                      }
                      uptothislevel += oncurrentlevel;
                //    }
                    // increase index = (currentlevel,currenttype)
                    // "small loop"
                    // iterate over all combinations of generators/wavelets for all dimensions with currentlevel[i]=j0[i]
                    // this looks like binary addition: (in 3 Dim:) gwg is followed by gww (g=0,w=1)
                    bool done = true;
                    for (int i(DIM-1); i >= 0; i--)
                    {
                        // find first position on level j0
                        if (currentlevel[i] == j0[i])
                        {
                            if (currenttype[i] == 1)
                            {
                                currenttype[i]=0;
                            } else
                            {
                                currenttype[i]=1;
                                done = false;
                    //            doch = 1;
                                break;
                            }
                        }
                    }
                    // done == true means that all components with currentlevel[i]=j0[i] were wavelets.
                    // the level j has to be increased.
                    // iterate "big loop", meaning: "currentlevel++"

                    if (done == true)
                    {
                    for (int i(DIM-1); i >= 0; i--)
                    {
                        if (i != 0)
                        {
                            if (currentlevel[i] != j0[i])
                            {
                                // increase left neighbor
                                currentlevel[i-1]=currentlevel[i-1]+1;
                                if (currentlevel[i-1]-j0[i] == range) sizes[i-1][range+1]=basis.bases()[i]->Nablasize(currentlevel[i-1]); // if needed compute and store new size information
                                currenttype[i-1]=1;
                                int temp = currentlevel[i]-j0[i];
                                currentlevel[i]=j0[i];
                                currenttype[i]=0;
                                currentlevel[DIM-1]=j0[DIM-1]+temp-1;
                                currenttype[DIM-1]= (temp == 1?0:1);
  				break;
                            }
  			} else // i == 0. "big loop" arrived at the last index. We have to increase the level!
  			{
                            range = range +1;
                            if (DIM == 1)
                            {
                                currenttype[i] = 1;
                                currentlevel[i]=currentlevel[i]+1;
                                sizes[0][range+1]=basis.bases()[0]->Nablasize(currentlevel[0]); // if needed compute and store new size information
                            }
                            else
                            {
                                //currentlevel[DIM-1]=j0[DIM-1]+currentlevel[0]-j0[0]+1; currenttype[DIM-1]=1;
                                currentlevel[DIM-1]=j0[DIM-1]+range; currenttype[DIM-1]=1;
				currenttype[0]=0; currentlevel[0]=j0[0];
                                sizes[DIM-1][range+1]=basis.bases()[DIM-1]->Nablasize(currentlevel[DIM-1]); // if needed compute and store new size information
                            }
  			break; // unnoetig, da i==0 gilt.
  			}
                    } // end of "big loop"
                    }  // end if done == true
                    if (currentlevel == j || nureinmal) break; 
		}  // end while
                //berechnet wie viele Wavelets es gibt mit zu kleinem Indices

      FixedArray1D<int,DIM> help1, help2;
      
      for(unsigned int i = 0; i<DIM; i++)
         help1[i]=0;

      // berechnet wie viele indices mit einem zu kleinem translationstyp es gibt, so dass sich die Wavelets nicht schneiden
      unsigned int result2 = 0;
      for (unsigned int i = 0; i < DIM; i++) {  // begin for1
        int tmp = 1;

        for (unsigned int l = i+1; l < DIM; l++) 
	    tmp *= sizes[l][currentlevel[l]+currenttype[l]-j0[l]];;
        

        help2[i] = tmp;
        if (currenttype[i] == 0) {
	  if (minkgen[i] == (basis.bases()[i]->DeltaLmin()))
	    continue;
        }
        else
	  if (minkwavelet[i] == (basis.bases()[i]->Nablamin()))
	    continue;
      
        
        if (currenttype[i] == 0) {
	tmp *= minkgen[i]-basis.bases()[i]->DeltaLmin();
        }
        else
	  tmp *= minkwavelet[i]-basis.bases()[i]->Nablamin();

        result2 += tmp;
      }  // end for1

      tmp = 0;

      if (currenttype[DIM-1] == 0) {
	tmp = maxkgen[DIM-1] - minkgen[DIM-1]+1;
      }
      else{
	tmp = maxkwavelet[DIM-1] - minkwavelet[DIM-1]+1; 
      }
     
      bool exit2 = 0;

      while(!exit2){

      // fügt die Indizes ein die sich überlappen
      for (unsigned int i = uptothislevel + result2; i < uptothislevel + result2 + tmp; i++) {
        const Index* ind = basis.get_wavelet(i);   //&((*full_collection)[i]);
	intersecting.push_back(*ind);
      }

      for (unsigned int i = DIM-2; i >= 0; i--) {
            if(currenttype[i]==0){
	      if ( help1[i] < maxkgen[i]-minkgen[i]) {
	        help1[i]++;
                result2 = result2 + help2[i];
               // result2 = result2 +4;
                for (unsigned int j = i+1; j<=DIM-2;j++){
                    if(currenttype[i] == 0){
                       result2 = result2 - help2[j]*(maxkgen[j] - minkgen[j]+1);
                    } 
                    else
                       result2 = result2 - help2[j]*(maxkwavelet[j] - minkwavelet[j]+1);
                }
                break;
              }
              else {
                 help1[i]=0;
                 exit2 = (i==0);
                 break;
              }
            }
            else {
              if ( help1[i] < maxkwavelet[i] - minkwavelet[i]) {
	        help1[i]++;
                result2 = result2 + help2[i];
                //result2 = result2 +5;
                for (unsigned int j = i+1; j<=DIM-2;j++){
                    if(currenttype[i] == 0){
                       result2 = result2 - help2[j]*(maxkgen[j] - minkgen[j]+1);
                    }
                    else
                       result2 = result2 - help2[j]*(maxkwavelet[j] - minkwavelet[j]+1);
                }
                break;
              }
              else {
                 help1[i]=0;
                 exit2 = (i==0);
                 break;
              }
	    }
	  } //end for
      } //end while 2
  

     // berechnet den nächsten Typ 
   if(!generators)
     for (unsigned int i = DIM-1; i >= 0; i--) {
	    if (currenttype[i] == 0 ){  // ()   //  && currentlevel[i] == j0[i] && 
	     // currenttype[i]++;
              nureinmal = 1;
              exit = 0;
	      break;
	    }
            else if (currentlevel[i] = j0[i]){
              exit = (i==0);
              if(exit)
                break;
            }
	  } //end for
       } // end while 1

#endif
//#else
#if 0
                // a brute force solution
                typedef typename TensorBasis<IBASIS,DIM>::Support Support;
                Support supp;
                if (generators) {
                    for (Index mu = first_generator<IBASIS,DIM>(&basis, j);; ++mu) {
                        if (intersect_supports(basis, lambda, mu, supp))
                            intersecting.push_back(mu);
                        if (mu == last_generator<IBASIS,DIM>(&basis, j)) break;
                    }
                } else {
                    for (Index mu = first_wavelet<IBASIS,DIM>(&basis, j);; ++mu) {
                        if (intersect_supports(basis, lambda, mu, supp))
                            intersecting.push_back(mu);
                        if (mu == last_wavelet<IBASIS,DIM>(&basis, j)) break;
                    }
                }
#endif
    }

    template <class IBASIS, unsigned int DIM>
    void intersecting_elements(const TensorBasis<IBASIS,DIM>& basis,
                               const typename TensorBasis<IBASIS,DIM>::Index& lambda,
                               const MultiIndex<int,DIM> j,
                               std::list<typename TensorBasis<IBASIS,DIM>::Index>& intersecting)
    {
        typedef typename TensorBasis<IBASIS,DIM>::Index Index;
        intersecting.clear();

        // the set of intersecting wavelets is a cartesian product from d sets from the 1D case,
        // so we only have to compute the relevant 1D indices
        typedef typename IBASIS::Index Index1D;
        FixedArray1D<std::list<Index1D>,DIM> intersecting_1d_generators, intersecting_1d_wavelets;
        // prepare all intersecting wavelets and generators in the i-th coordinate direction

        /*
        typedef typename Index::type_type type_type;       
        type_type min_type,max_type;
        for (unsigned int i=0; i<DIM;i++)
        {
            min_type[i]= (basis.j0()[i] == j[i])? 0:1;
            max_type[i]=1;
        }
        */
        for (unsigned int i = 0; i < DIM; i++)
        {

            //if (mintype[i]==0)
            if (j[i] == basis.j0()[i])
                intersecting_wavelets(*basis.bases()[i],
                                      Index1D(lambda.j()[i],
                                              lambda.e()[i],
                                              lambda.k()[i],
                                              basis.bases()[i]),
                                      j[i], true, intersecting_1d_generators[i]);

            intersecting_wavelets(*basis.bases()[i],
                                  Index1D(lambda.j()[i],
                                          lambda.e()[i],
                                          lambda.k()[i],
                                          basis.bases()[i]),
                                  j[i], false, intersecting_1d_wavelets[i]);
        }
        // generate all tensor product indices

        // Iteration over the type-vektor e

        typename Index::type_type min_type,current_type;
        for (unsigned int i=0;i<DIM;i++)
        {
            min_type[i]=(j[i]==basis.j0()[i]) ?0:1;
            current_type[i]=min_type[i];
        }
//        cout << "tbasis type " << current_type << cout;cout.flush();

        
        // --------------------------------------------
        // speichere zu jedem Typ was bisher berechnet wurde. key = nummer des typs
        typedef std::list<FixedArray1D<Index1D,DIM> > list_type;
        typedef std::list<list_type> type_storage;
        type_storage storage;

        list_type temp_indices;
        FixedArray1D<Index1D,DIM> helpindex;
        if (j[0] == basis.j0()[0])
        {
            for (typename std::list<Index1D>::const_iterator it(intersecting_1d_generators[0].begin()),
                    itend(intersecting_1d_generators[0].end());
                    it != itend; ++it)
            {
                helpindex[0] = *it;
                temp_indices.push_back(helpindex);
            }
            storage.push_back(temp_indices);
            temp_indices.clear();
        }

        for (typename std::list<Index1D>::const_iterator it(intersecting_1d_wavelets[0].begin()),
                itend(intersecting_1d_wavelets[0].end());
                it != itend; ++it)
        {
            helpindex[0] = *it;
            temp_indices.push_back(helpindex);
        }
        storage.push_back(temp_indices);
        temp_indices.clear();

        for (int i = 1; i < DIM; i++)
        {
            type_storage bisher;
            bisher.swap(storage);
            for (typename type_storage::const_iterator it(bisher.begin()), itend(bisher.end()); it != itend; ++it)
            {
                // combine every element in "it" with every generator and wavelet in direction i

                //list_type sofar1;
                //sofar1.swap(indices1);

                if (j[i]==basis.j0()[i])
                { // add all possible generators
                    for (typename list_type::const_iterator itready((*it).begin()), itreadyend((*it).end()); itready != itreadyend; ++itready)
                    {
                        helpindex = *itready;
                        //unsigned int esum = 0;
                        //for (unsigned int k = 0; k < i; k++)
                        //esum += helpindex[k].e();
                        for (typename std::list<Index1D>::const_iterator it(intersecting_1d_generators[i].begin()), itend(intersecting_1d_generators[i].end()); it != itend; ++it)
                        {
                            helpindex[i] = *it;
                            temp_indices.push_back(helpindex);
                        }
                    }
                    storage.push_back(temp_indices);
                    temp_indices.clear();
                }
                // now add all possible wavelets
                for (typename list_type::const_iterator itready((*it).begin()), itreadyend((*it).end()); itready != itreadyend; ++itready)
                {
                    helpindex = *itready;
                    for (typename std::list<Index1D>::const_iterator it(intersecting_1d_wavelets[i].begin()), itend(intersecting_1d_wavelets[i].end()); it != itend; ++it)
                    {
                        helpindex[i] = *it;
                        temp_indices.push_back(helpindex);
                    }
                }
                storage.push_back(temp_indices);
                temp_indices.clear();
            }
        }
        // compose the results
        typename Index::type_type help_e;
        typename Index::translation_type help_k;
        for (typename type_storage::const_iterator it_store(storage.begin()),it_store_end(storage.end());it_store != it_store_end; it_store++)
        {
            //list_type somelist(*it_store);
            //typename list_type::const_iterator it22(*it_store.begin()), itend22(*it_store.end());
            for (typename list_type::const_iterator it((*it_store).begin()), itend((*it_store).end()); it != itend; ++it)
            {
                for (unsigned int i = 0; i < DIM; i++)
                {
                    help_e[i] = (*it)[i].e();
                    help_k[i] = (*it)[i].k();
                }
                intersecting.push_back(Index(j, help_e, help_k, &basis));
            }
        }
    };

    template <class IBASIS, unsigned int DIM>
    bool intersect_singular_support(const TensorBasis<IBASIS,DIM>& basis,
                                    const typename TensorBasis<IBASIS,DIM>::Index& lambda,
                                    const typename TensorBasis<IBASIS,DIM>::Index& mu)
    {
        // we have intersection of the singular supports if and only if
        // (cube_support:)   one of the components has this property in one dimension
        // (tbasis_support:) all of the components have this property
        typedef typename IBASIS::Index Index1D;
        for (unsigned int i = 0; i < DIM; i++) {
            if (!(intersect_singular_support(*basis.bases()[i],
                                             Index1D(lambda.j()[i], lambda.e()[i], lambda.k()[i], basis.bases()[i]),
                                             Index1D(mu.j()[i], mu.e()[i], mu.k()[i], basis.bases()[i]))
                          ))
                return false;
        }
        return true;
    }
}

