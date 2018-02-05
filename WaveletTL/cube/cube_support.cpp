// implementation for cube_support.h

#include <utils/multiindex.h>
#include <utils/fixed_array1d.h>

using MathTL::multi_degree;
using MathTL::FixedArray1D;

namespace WaveletTL
{
  template <class IBASIS, unsigned int DIM>
  inline
  void
  support(const CubeBasis<IBASIS,DIM>& basis,
	  const typename CubeBasis<IBASIS,DIM>::Index& lambda,
	  typename CubeBasis<IBASIS,DIM>::Support& supp)
  {
    basis.support(lambda, supp);
  }

  template <class IBASIS, unsigned int DIM>
  bool
  intersect_supports(const CubeBasis<IBASIS,DIM>& basis,
		     const typename CubeBasis<IBASIS,DIM>::Index& lambda,
		     const typename CubeBasis<IBASIS,DIM>::Index& mu,
		     typename CubeBasis<IBASIS,DIM>::Support& supp)
  {
    typename CubeBasis<IBASIS,DIM>::Support supp_lambda;
    WaveletTL::support<IBASIS,DIM>(basis, lambda, supp_lambda);

    typename CubeBasis<IBASIS,DIM>::Support supp_mu;
    WaveletTL::support<IBASIS,DIM>(basis, mu, supp_mu);

    // determine support intersection granularity,
    // adjust single support granularities if necessary
    supp.j = std::max(supp_lambda.j, supp_mu.j);
    if (supp_lambda.j > supp_mu.j)
    {
      const int adjust = 1<<(supp_lambda.j-supp_mu.j);
      for (unsigned int i = 0; i < DIM; i++)
      {
	    supp_mu.a[i] *= adjust;
	    supp_mu.b[i] *= adjust;
      }
    }
    else
    {
      const int adjust = 1<<(supp_mu.j-supp_lambda.j);
      for (unsigned int i = 0; i < DIM; i++)
      {
	    supp_lambda.a[i] *= adjust;
	    supp_lambda.b[i] *= adjust;
      }
    }

    for (unsigned int i = 0; i < DIM; i++)
    {
      supp.a[i] = std::max(supp_lambda.a[i],supp_mu.a[i]);
      supp.b[i] = std::min(supp_lambda.b[i],supp_mu.b[i]);

      if (supp.a[i] >= supp.b[i])
      {
	    return false;
      }
    }

    return true;
  }

  template <class IBASIS, unsigned int DIM>
  void intersecting_wavelets(const CubeBasis<IBASIS,DIM>& basis,
			     const typename CubeBasis<IBASIS,DIM>::Index& lambda,
			     const int j, const bool generators,
			     std::list<typename CubeBasis<IBASIS,DIM>::Index>& intersecting)
  {
    typedef typename CubeBasis<IBASIS,DIM>::Index Index;

    intersecting.clear();

#if 1
    // the set of intersecting wavelets is a cartesian product from d sets from the 1D case,
    // so we only have to compute the relevant 1D indices
    typedef typename IBASIS::Index Index1D;
    FixedArray1D<std::list<Index1D>,DIM>
      intersecting_1d_generators, intersecting_1d_wavelets;

    // prepare all intersecting wavelets and generators in the i-th coordinate direction
    for (unsigned int i = 0; i < DIM; i++)
    {
      intersecting_wavelets(*basis.bases()[i],
			    Index1D(lambda.j(),
				    lambda.e()[i],
				    lambda.k()[i],
				    basis.bases()[i]),
			    j, true, intersecting_1d_generators[i]);
      if (!(generators))
      {
	    intersecting_wavelets(*basis.bases()[i],
			      Index1D(lambda.j(),
				      lambda.e()[i],
				      lambda.k()[i],
				      basis.bases()[i]),
			      j, false, intersecting_1d_wavelets[i]);
      }
    }

    // generate all relevant tensor product indices with either e=(0,...,0) or e!=(0,...,0)
    typedef std::list<FixedArray1D<Index1D,DIM> > list_type;
    list_type indices;
    FixedArray1D<Index1D,DIM> helpindex;

    if (DIM > 1 || (DIM == 1 && generators))
    {
      for (typename std::list<Index1D>::const_iterator it(intersecting_1d_generators[0].begin()),
	       itend(intersecting_1d_generators[0].end());
	       it != itend; ++it)
	   {
	     helpindex[0] = *it;
	     indices.push_back(helpindex);
       }
    }
    if (!(generators))
    {
      for (typename std::list<Index1D>::const_iterator it(intersecting_1d_wavelets[0].begin()),
	     itend(intersecting_1d_wavelets[0].end());
	   it != itend; ++it)
	   {
	     helpindex[0] = *it;
	     indices.push_back(helpindex);
       }
    }
    for (unsigned int i = 1; i < DIM; i++)
    {
      list_type sofar;
      sofar.swap(indices);
      for (typename list_type::const_iterator itready(sofar.begin()), itreadyend(sofar.end());
	       itready != itreadyend; ++itready)
      {
	    helpindex = *itready;
	    unsigned int esum = 0;
	    for (unsigned int k = 0; k < i; k++)
	    {
	      esum += helpindex[k].e();
        }
	    if (generators || (i < DIM-1 || (i == (DIM-1) && esum > 0)))
	    {
	      for (typename std::list<Index1D>::const_iterator it(intersecting_1d_generators[i].begin()),
		       itend(intersecting_1d_generators[i].end());
	           it != itend; ++it)
	       {
	         helpindex[i] = *it;
	         indices.push_back(helpindex);
	       }
	    }
	    if (!(generators))
	    {
	      for (typename std::list<Index1D>::const_iterator it(intersecting_1d_wavelets[i].begin()),
		       itend(intersecting_1d_wavelets[i].end());
	           it != itend; ++it)
          {
	        helpindex[i] = *it;
	        indices.push_back(helpindex);
	      }
	    }
      }
    }

    // compose the results
    typename Index::type_type help_e;
    typename Index::translation_type help_k;
    for (typename list_type::const_iterator it(indices.begin()), itend(indices.end());
	     it != itend; ++it)
    {
      for (unsigned int i = 0; i < DIM; i++)
      {
	    help_e[i] = (*it)[i].e();
	    help_k[i] = (*it)[i].k();
      }
      intersecting.push_back(Index(j, help_e, help_k, &basis));
    }


#else

    typedef typename CubeBasis<IBASIS,DIM>::Index Index;

    int k = -1;
	if ( generators ) {
      k=0;
    }
    else {
      k=j-basis.j0()+1;
    }
    //std::list<typename Frame::Index> intersect_diff;

//! generators
    if (true) {
      FixedArray1D<int,DIM>
      minkwavelet, maxkwavelet, minkgen, maxkgen;
      typedef typename IBASIS::Index Index1D;
      int minkkkk;
      int maxkkkk;

      // prepare all intersecting wavelets and generators in the i-th coordinate direction
      for (unsigned int i = 0; i < DIM; i++) {
        get_intersecting_wavelets_on_level(*basis.bases()[i],
	    	Index1D(lambda.j(),
		    lambda.e()[i],
		    lambda.k()[i],
		    basis.bases()[i]),
		    j, true, minkkkk,maxkkkk);
        minkgen[i]=minkkkk;
        maxkgen[i] = maxkkkk;
        if (!(generators))
	  get_intersecting_wavelets_on_level(*basis.bases()[i],
		      Index1D(lambda.j(),
			      lambda.e()[i],
			      lambda.k()[i],
			      basis.bases()[i]),
		      j, false, minkkkk,maxkkkk);
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
      unsigned int tmp = 1;
      bool exit = 0;

      // determine how many wavelets there are on all the levels
      // below the level of this index
      if (! gen) {
        result = 0;
        genfstlvl =1;
        //generators on level j0
        for (unsigned int i = 0; i< DIM; i++)
	  genfstlvl *= (basis.bases()[i])->Deltasize((basis.bases()[i])->j0());
        //additional wavelets on level j
        //            =(#Gen[1]+#Wav[1])*...*(#Gen[Dim-1]+#Wav[Dim-1])
        //             -#Gen[1]*...*#Gen[Dim-1]
        for (int lvl= 0 ;
	     lvl < (j -basis.j0());
	     lvl++){
	  int genCurLvl = 1;
	  int addWav = 1;
	  for (unsigned int i = 0; i< DIM; i++) {
	    unsigned int curJ = basis.bases()[i]->j0()+lvl;
	    int genCurDim = (basis.bases()[i])->Deltasize(curJ);
	    genCurLvl *= genCurDim;
	    addWav *= genCurDim+ (basis.bases()[i])->Nablasize(curJ);
	  }
	  result += addWav-genCurLvl;
        }
        result += genfstlvl;
      }



      while(!exit){
      FixedArray1D<int,DIM> help1, help2;

      for(unsigned int i = 0; i<DIM; i++)
         help1[i]=0;

      // berechnet wie viele indices mit einem zu kleinem translationstyp es gibt, so dass sich die Wavelets nicht schneiden
      unsigned int result2 = 0;
      for (unsigned int i = 0; i < DIM; i++) {  // begin for1
        int tmp = 1;

        for (unsigned int l = i+1; l < DIM; l++) {
	  if (type[l] == 0)
	    tmp *= (basis.bases())[l]->Deltasize(j);
	  else
	    tmp *= (basis.bases())[l]->Nablasize(j);
        }

        help2[i] = tmp;

        if (type[i] == 0) {
	  if (minkgen[i] == (basis.bases())[i]->DeltaLmin())
	    continue;
        }
        else
	  if (minkwavelet[i] == (basis.bases())[i]->Nablamin())
	    continue;


        if (type[i] == 0) {
	tmp *= minkgen[i]-(basis.bases())[i]->DeltaLmin();
        }
        else
	  tmp *= minkwavelet[i]-(basis.bases())[i]->Nablamin();

        result2 += tmp;
      }  // end for1

      int tmp = 0;

      if (type[DIM-1] == 0) {
	tmp = maxkgen[DIM-1] - minkgen[DIM-1]+1;
      }
      else{
	tmp = maxkwavelet[DIM-1] - minkwavelet[DIM-1]+1;
      }

      bool exit2 = 0;

      while(!exit2){

      // fügt die Indizes ein die sich überlappen
      for (unsigned int i = result + result2; i < result + result2 + tmp; i++) {
        const Index* ind = basis.get_wavelet(i);   //&((*full_collection)[i]);
	intersecting.push_back(*ind);
      }

      for (unsigned int i = DIM-2; i >= 0; i--) {
            if(type[i]==0){
	      if ( help1[i] < maxkgen[i]-minkgen[i]) {
	        help1[i]++;
                result2 = result2 + help2[i];
                for (unsigned int j = i+1; j<=DIM-2;j++){
                    if(type[i] == 0){
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
                for (unsigned int j = i+1; j<=DIM-2;j++){
                    if(type[i] == 0){
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


      // berechnet wie viele Indizes von dem jeweiligen Typ in Patches p liegen
      tmp = 1;
      for (unsigned int i = 0; i < DIM; i++) {
	      if (type[i] == 0)
	        tmp *= (basis.bases())[i]->Deltasize(j);
	      else
	        tmp *= (basis.bases())[i]->Nablasize(j);
	    }

      result += tmp;


     // berechnet den nächsten Typ
     for (unsigned int i = DIM-1; i >= 0; i--) {
	    if ( type[i] == 1 ) {
	      type[i] = 0;
              exit = (i == 0);
	      if(exit)
               break;
	    }
	    else {
	      type[i]++;
	      break;
	    }
	  } //end for
       } // end while 1
      } // end if
    // } // end if

    else { // if generators
        // a brute force solution
    typedef typename CubeBasis<IBASIS,DIM>::Support Support;
    Support supp;
    if (generators) {
      for (Index mu = basis.first_generator (j);; ++mu) {
	if (intersect_supports(basis, lambda, mu, supp))
	  intersecting.push_back(mu);
	if (mu == basis.last_generator(j)) break;
      }
    }
    }


//*/
#endif
//#else
#if 0
    // a brute force solution
    typedef typename CubeBasis<IBASIS,DIM>::Support Support;
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
  bool intersect_singular_support(const CubeBasis<IBASIS,DIM>& basis,
				  const typename CubeBasis<IBASIS,DIM>::Index& lambda,
				  const typename CubeBasis<IBASIS,DIM>::Index& mu)
  {
    // we have intersection of the singular supports if and only if
    // one of the components have this property in one dimension
    typedef typename IBASIS::Index Index1D;
    for (unsigned int i = 0; i < DIM; i++) {
      if (intersect_singular_support
	  (*basis.bases()[i],
	   Index1D(lambda.j(), lambda.e()[i], lambda.k()[i], basis.bases()[i]),
	   Index1D(mu.j(), mu.e()[i], mu.k()[i], basis.bases()[i])))
	return true;
    }
    return false;
  }
}
