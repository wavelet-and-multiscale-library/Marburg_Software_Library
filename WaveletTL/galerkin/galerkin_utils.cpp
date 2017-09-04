// implementation for galerkin_utils.h

namespace WaveletTL
{
  //  template <class PROBLEM>
//   void setup_stiffness_matrix(PROBLEM& P,
// 			      const std::set<typename PROBLEM::Index>& Lambda,
// 			      SparseMatrix<double>& A_Lambda,
// 			      bool preconditioned)
//   {
//     A_Lambda.resize(Lambda.size(), Lambda.size());
    
//     typedef typename SparseMatrix<double>::size_type size_type;

//     size_type row = 0;
//     typedef typename PROBLEM::Index Index;
//     typedef std::list<Index> IntersectingList;
//     for (typename std::set<Index>::const_iterator it1(Lambda.begin()), itend(Lambda.end());
// 	 it1 != itend; ++it1, ++row)
//       {
// 	const double d1 = preconditioned ? P.D(*it1) : 1.0;
// 	std::list<size_type> indices;
// 	std::list<double> entries;

// 	// determine list of thise indices which intersect with *it1
//  	IntersectingList Nu;
//  	intersecting_wavelets(P.basis(), *it1, (*it1).p(), Lambda, Nu);

// 	size_type column = 0;
// #if _WAVELETTL_GALERKINUTILS_VERBOSITY >= 1
// 	cout << "setup_stiffness_matrix(): doing row " << row << " of " << Lambda.size()
// 	     << " (wavelet index " << *it1 << ")" << endl;
// #endif
// 	typename std::list<Index>::const_iterator it3 = Nu.begin(); 
// 	//cout << "GROESSE VON NU= " << Nu.size() << endl;
// // 	int count = 0;
// // 	if (Nu.size() != Lambda.size()) {
// // // 	  cout << "+++++++++++++++++++++++++++++++++++++++" << endl;
// // // 	  cout << Nu.size() << " " << Lambda.size() << endl;
// // 	  //abort();
// // 	}
//  	for (typename std::set<Index>::const_iterator it2(Lambda.begin());
//  	     (it2 != itend) && (it3 != Nu.end()) ; ++it2, ++it3, ++column)
// 	  // 	for (typename std::list<Index>::const_iterator it2(Nu.begin()), itend2(Nu.end());
// 	  // 	     it2 != itend2; ++it2, ++column)
// 	  {
// 	    //cout << *it2 << " " << *it3 << endl;
// 	    //cout << (*it2).number() << " " << (*it3).number() << endl;
// 	    while ((*it2).number() != (*it3).number()) {
// 	      ++it2;
// 	      ++column;
// 	    }
// 	    // 	    if (intersect_singular_support(P.basis(), *it1, *it2)) {
// 	    double entry = P.a(*it2, *it1);
// 	    //count++;
// 	    //const double entry = 0;
// #if _WAVELETTL_GALERKINUTILS_VERBOSITY >= 2
//  	    if (fabs(entry) > 1e-15) {
//  	      cout << " column: " << *it2 <<  ", value " << entry << endl;
//  	    }
// #endif
// 	    if (fabs(entry) > 1e-15) {
// 		indices.push_back(column);
// 		entries.push_back(entry / (d1 * (preconditioned ? P.D(*it2) : 1.0)));
// 	    }
// 	    // 	    }
// 	  }
// // 	if (count != Lambda.size()) {
// // 	  cout << "++++++++++++++++++++++" << endl;
// // 	}
// 	A_Lambda.set_row(row, indices, entries);
//       }
//   }
 

  template <class PROBLEM>
  void setup_stiffness_matrix(PROBLEM& P,
			      const std::set<typename PROBLEM::Index>& Lambda,
			      SparseMatrix<double>& A_Lambda,
			      bool preconditioned)
  {
      
    A_Lambda.resize(Lambda.size(), Lambda.size());

    typedef typename SparseMatrix<double>::size_type size_type;

   
    typedef typename PROBLEM::Index Index;
#if PARALLEL==1
    cout<<"parallel computing stiffness matrix"<<endl;
#pragma omp parallel
    {
#pragma omp for
    for(size_type row=0;row<Lambda.size();row++){
        typename std::set<Index>::const_iterator it1(Lambda.begin());
        advance(it1, row);
#else
    cout<<"sequentiell computing stiffness matrix"<<endl;
    size_type row=0;
    for (typename std::set<Index>::const_iterator it1(Lambda.begin()), itend(Lambda.end());
	 it1 != itend; ++it1, ++row)
      {
#endif
//        cout << "Index nu: " << *it1 << endl;
	const double d1 = preconditioned ? P.D(*it1) : 1.0;
        
	std::list<size_type> indices;
	std::list<double> entries;
#if 0 //parallelization of inner for loop is not useful
#pragma omp parallel for
        for(int column=0;column<Lambda.size();column++){
            typename std::set<Index>::const_iterator it2(Lambda.begin());
            advance(it2,column);
#else
	size_type column = 0;
	for (typename std::set<Index>::const_iterator it2(Lambda.begin()), itend(Lambda.end());
	     it2 != itend; ++it2, ++column)
	  {
#endif
	    // 	    if (intersect_singular_support(P.basis(), *it1, *it2)) {
            
            double entry = P.a(*it2, *it1);
	    
	    //const double entry = 0;
#if _WAVELETTL_GALERKINUTILS_VERBOSITY >= 2
 	    if (fabs(entry) > 1e-15) {
 	      cout << " column: " << *it2 <<  ", value " << entry << endl;
 	    }
#endif
	    if (fabs(entry) > 1e-15) {
		indices.push_back(column);
		entries.push_back(entry / (preconditioned ? d1 * P.D(*it2) : 1.0));
	    }
	    // 	    }
	  }
	A_Lambda.set_row(row, indices, entries);
      }
#if PARALLEL==1
    }
#endif
  }

    template <class PROBLEM>
    void setup_stiffness_matrix(PROBLEM& P,
            const std::set<int>& Lambda,
            SparseMatrix<double>& A_Lambda,
            bool preconditioned)
    {
        A_Lambda.resize(Lambda.size(), Lambda.size());
        typedef typename SparseMatrix<double>::size_type size_type;
        size_type row = 0;
        //typedef typename PROBLEM::Index Index;
        if (preconditioned)
        {
            for (typename std::set<int>::const_iterator it1(Lambda.begin()), itend(Lambda.end());
                    it1 != itend; ++it1, ++row)
            {
                const double d1 = P.D(*it1);
                std::list<size_type> indices;
                std::list<double> entries;

                size_type column = 0;
    
                for (typename std::set<int>::const_iterator it2(Lambda.begin());
                        it2 != itend; ++it2, ++column)
                {
                    double entry = P.a(*it2, *it1);
    #if _WAVELETTL_GALERKINUTILS_VERBOSITY >= 2
                    if (fabs(entry) > 1e-15) 
                    {
                        cout << " column: " << *it2 <<  ", value " << entry << endl;
                    }
    #endif
                    if (fabs(entry) > 1e-15) 
                    {
                        indices.push_back(column);
                        entries.push_back(entry / (d1 * P.D(*it2)));
                    }
                }
                A_Lambda.set_row(row, indices, entries);
            }
        }
        else
        {
            for (typename std::set<int>::const_iterator it1(Lambda.begin()), itend(Lambda.end());
                    it1 != itend; ++it1, ++row)
            {
                std::list<size_type> indices;
                std::list<double> entries;
                size_type column = 0;
   
                for (typename std::set<int>::const_iterator it2(Lambda.begin());
                        it2 != itend; ++it2, ++column)
                {
                    double entry = P.a(*it2, *it1);
    #if _WAVELETTL_GALERKINUTILS_VERBOSITY >= 2
                    if (fabs(entry) > 1e-15) 
                    {
                        cout << " column: " << *it2 <<  ", value " << entry << endl;
                    }
    #endif
                    if (fabs(entry) > 1e-15) 
                    {
                        indices.push_back(column);
                        entries.push_back(entry);
                    }
                }
                A_Lambda.set_row(row, indices, entries);
            }
        }
        
    }
  

  template <class PROBLEM>
  void setup_stiffness_matrix(PROBLEM& P,
			      const std::set<typename PROBLEM::Index>& Lambda1,
			      const std::set<typename PROBLEM::Index>& Lambda2,
			      SparseMatrix<double>& A_Lambda,
			      bool preconditioned)
  {
      
    A_Lambda.resize(Lambda1.size(), Lambda2.size());
    
    typedef typename SparseMatrix<double>::size_type size_type;


    size_type row = 0;
    typedef typename PROBLEM::Index Index;
    typedef std::list<Index> IntersectingList;
      
    for (typename std::set<Index>::const_iterator it1(Lambda1.begin()), itend(Lambda1.end());
	 it1 != itend; ++it1, ++row)
      {

	// determine list of thise indices which intersect with *it1
	IntersectingList Nu;
	intersecting_wavelets(P.basis(), *it1.p(), *it1.p(), Lambda2, Nu);

	const double d1 = preconditioned ? P.D(*it1) : 1.0;
	std::list<size_type> indices;
	std::list<double> entries;

	size_type column = 0;

// 	for (typename std::set<Index>::const_iterator it2(Lambda2.begin()), itend2(Lambda2.end());
// 	     it2 != itend2; ++it2, ++column)
	for (typename std::list<Index>::const_iterator it2(Nu.begin()), itend2(Nu.end());
	     it2 != itend2; ++it2, ++column)

	  {
	    //if (intersect_singular_support(P.basis(), *it1, *it2))
	    //double entry = P.a(*it2, *it1);
	    double entry = P.a(*it1, *it2);
	      
	    //const double entry = 0;
#if _WAVELETTL_GALERKINUTILS_VERBOSITY >= 2
 	    if (fabs(entry) > 1e-15) {
 	      cout << " column: " << *it2 <<  ", value " << entry << endl;
 	    }
#endif
	    if (fabs(entry) > 1e-15) {
		indices.push_back(column);
		entries.push_back(entry / (preconditioned ? d1 * P.D(*it2) : 1.0));
	    }
	    // 	    }
	  }
	A_Lambda.set_row(row, indices, entries);
      }
        
  }
//   template <class PROBLEM>
//   void setup_stiffness_matrix(PROBLEM& P,
// 			      const std::set<typename PROBLEM::Index>& Lambda1,
// 			      const std::set<typename PROBLEM::Index>& Lambda2,
// 			      SparseMatrix<double>& A_Lambda,
// 			      bool preconditioned)
//   {
//     A_Lambda.resize(Lambda1.size(), Lambda2.size());
    
//     typedef typename SparseMatrix<double>::size_type size_type;

//     size_type row = 0;
//     typedef typename PROBLEM::Index Index;
//     for (typename std::set<Index>::const_iterator it1(Lambda1.begin()), itend(Lambda1.end());
// 	 it1 != itend; ++it1, ++row)
//       {
// 	const double d1 = preconditioned ? P.D(*it1) : 1.0;
// 	std::list<size_type> indices;
// 	std::list<double> entries;

// 	size_type column = 0;
// #if _WAVELETTL_GALERKINUTILS_VERBOSITY >= 1
// 	cout << "setup_stiffness_matrix(): doing row " << row << " of " << Lambda1.size()
// 	     << " (wavelet index " << *it1 << ")" << endl;
// #endif
// 	for (typename std::set<Index>::const_iterator it2(Lambda2.begin()), itend2(Lambda2.end());
// 	     it2 != itend2; ++it2, ++column)
// 	  {
// 	    // 	    if (intersect_singular_support(P.basis(), *it1, *it2)) {
// 	    //double entry = P.a(*it2, *it1);
// 	    double entry = P.a(*it1, *it2);
	    
// 	    //const double entry = 0;
// #if _WAVELETTL_GALERKINUTILS_VERBOSITY >= 2
//  	    if (fabs(entry) > 1e-15) {
//  	      cout << " column: " << *it2 <<  ", value " << entry << endl;
//  	    }
// #endif
// 	    if (fabs(entry) > 1e-15) {
// 		indices.push_back(column);
// 		entries.push_back(entry / (d1 * (preconditioned ? P.D(*it2) : 1.0)));
// 	    }
// 	    // 	    }
// 	  }
// 	A_Lambda.set_row(row, indices, entries);
//       }
//   }


  template <class PROBLEM>
  void setup_righthand_side(PROBLEM& P,
			    const std::set<typename PROBLEM::Index>& Lambda,
			    Vector<double>& F_Lambda)
  {
    F_Lambda.resize(Lambda.size());

    typedef typename SparseMatrix<double>::size_type size_type;

    size_type row = 0;
    typedef typename PROBLEM::Index Index;
    for (typename std::set<Index>::const_iterator it(Lambda.begin()), itend(Lambda.end());
	 it != itend; ++it, ++row) {
      F_Lambda[row] = P.f(*it)/P.D(*it);
    }
  }
  

    /*
     * int variant.
     * faster.
     * Difference to Index variant is the assumption that P.f(int) is already preconditioned!
     */
    template <class PROBLEM>
    void setup_righthand_side(PROBLEM& P,
            const std::set<int>& Lambda,
            Vector<double>& F_Lambda)

    {
        F_Lambda.resize(Lambda.size());
        typedef typename SparseMatrix<double>::size_type size_type;
        size_type row = 0;
        for (typename std::set<int>::const_iterator it(Lambda.begin()), itend(Lambda.end()); it != itend; ++it, ++row) 
        {
            F_Lambda[row] = P.f(*it);
        }
    }
}
