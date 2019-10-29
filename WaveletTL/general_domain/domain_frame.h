// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Philipp Keding, Alexander Sieber                                   |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_DOMAIN_FRAME_H
#define _WAVELETTL_DOMAIN_FRAME_H

#include <map>

#include <algebra/vector.h>
#include <algebra/infinite_vector.h>
#include <algebra/sparse_matrix.h>
#include <algebra/block_matrix.h>
#include <geometry/sampled_mapping.h>
#include <utils/fixed_array1d.h>
#include <utils/multiindex.h>

#include <general_domain/domain_frame_index.h>

// for convenience, include also some functionality
//#include <general_domain/domain_support.h>

using std::map;
using std::set;
using MathTL::FixedArray1D;
using MathTL::InfiniteVector;
using MathTL::SparseMatrix;
using MathTL::BlockMatrix;

namespace WaveletTL
{
  /*!
    Template class for composite quarklet frames over decomposable domains. big cubes are not yet supported
   
    such that the primal and the dual frame fulfill homogeneous Dirichlet b.c.'s
    on the outer domain boundary.

    References:
     [CDFS13] Piecewise Tensor Product Wavelet Bases by Extensions and Approximation Rates
              N. Chegini, S. Dahlke, U. Friedrich, R. Stevenson
              Mathematics of Computation 82 (2013), 2157-2190. 
     
  */
  template <class IFRAME, int NPATCHES> //NPATCHES has to be accessible to support struct
  class DomainFrame
  {
  public:
      
    //! type of the interval frame
    typedef IFRAME IntervalFrame;
    
    
    //! quarklet index class
    typedef DomainFrameIndex<IntervalFrame, NPATCHES> Index;
    typedef typename Index::level_type level_type;
    typedef typename Index::polynomial_type polynomial_type;
    typedef typename Index::type_type type_type;
    

    //! default constructor
    DomainFrame(const IntervalFrame* frame1d, const IntervalFrame* frame1d_11, 
          const IntervalFrame* frame1d_01, const IntervalFrame* frame1d_10, 
            const Array1D<Point<2, int> >& corners, const Array1D<FixedArray1D<int,2> >& extensions);


    //! virtual destructor
    virtual ~DomainFrame() {}

    //! coarsest possible level j0
    inline const level_type& j0() const { return j0_; }
    
    

    //! geometric type of the support sets 
    //Maybe a constructor is necessary for correct initialization. Put it on, if errors occur @PHK
    typedef struct Support{
      int j[2];       // granularity
      int xmin[NPATCHES];
      int xmax[NPATCHES];
      int ymin[NPATCHES];
      int ymax[NPATCHES];
      Support() //general struct constructor
      {
          j[0]=0,j[1]=0;
          for(int i=0;i<NPATCHES;i++){
              xmin[i]=0,xmax[i]=0,ymin[i]=0,ymax[i]=0;
          }
      };
      
      
    } Support; 
    
    //alternative support
//    typedef struct {
//      		int j[2];
//      		int xmin[];
//                int xmax[];
//                int ymin[];
//                int ymax[];
//    	} Support;
    
//    inline std::ostream& operator << (std::ostream& os,
//                                      const typename  Support& supp)
//    {
//      using namespace std;
//      
//      os << "("
//         << supp.j[0]
//         << ","
//         << supp.j[1]
//         << ", ["
//         << supp.xmin[0]
//         << ","
//         << supp.xmax[0]
//         << "])";
//      return os;
//    }
    
    

    //! compute the support of psi_lambda, using the internal cache
    void support(const Index& lambda, Support& supp) const;
    
    //! compute the support of psi_lambda, using the internal cache; only works for precomputed supports
    void support(const int& lambda_num, Support& supp) const;
    
    //! critical Sobolev regularity for the primal generators/quarklets
    static double primal_regularity() { return IntervalFrame::primal_regularity(); } // dirty, we should use max(1.5,~) instead
    
    //! degree of polynomial reproduction for the primal generators/quarklets
    static unsigned int primal_polynomial_degree() { return IntervalFrame::primal_polynomial_degree(); }

    //! number of vanishing moments for the primal quarklets
    static unsigned int primal_vanishing_moments() { return IntervalFrame::primal_vanishing_moments(); }
    
    //! read access to the number of patches
    const int num_real_patches() const {return corners_.size(); }
    
    //! read access to the number of patches
    const int num_logical_patches() const {return extensions_.size(); }
    
    //! read access to corners
    const Array1D<Point<2, int> > get_corners() const {return corners_;}
    
    //! read access to extensions
    const Array1D<FixedArray1D<int,2> > get_extensions() const{return extensions_;}

    //! read access to the underlying 1D frame
    const IntervalFrame* frame1d() const { return frame1d_; }
    
    //! read access to the underlying 1D frame
    const IntervalFrame* frame1d_11() const { return frame1d_11_; }
    
    //! read access to the underlying 1D frame
    const IntervalFrame* frame1d_01() const { return frame1d_01_; }
    
    //! read access to the underlying 1D frame
    const IntervalFrame* frame1d_10() const { return frame1d_10_; }
    
    //! read access to the underlying 1D frame
    const IntervalFrame* frames(const int patch, const int dir) const { 
           //methode nicht sinnvoll für logische patches
        //loop over non-logical patches
        int s0=1,s1=1;
        for(unsigned int ext1=0;ext1<extensions_.size();ext1++){
            if(extensions_[ext1][0]==patch){ //extension from patch to another one
                int target_patch=extensions_[ext1][1];
//                cout<<"target_patch="<<target_patch<<endl;
                if(corners_[target_patch][0]==corners_[patch][0] && dir==1){ //extension in y-direction
                    if(corners_[target_patch][1]<corners_[patch][1]){ //extension from north to south
//                        return frame1d_01_;
                        s0=0;
//                        cout<<"s0="<<s0<<endl;
                    }
                    else{
//                        return frame1d_10_; //extension from south to north
                        s1=0;
//                        cout<<"s0="<<s0<<endl;
                    }
                }
                if(corners_[target_patch][1]==corners_[patch][1] && dir==0){ //extension in x-direction
                    if(corners_[target_patch][0]<corners_[patch][0]){ //extension from east to west
//                        return frame1d_01_;
                        s0=0;
//                        cout<<"s0="<<s0<<endl;
                    }
                    else{
//                        return frame1d_10_; //extension from west to east
                        s1=0;
//                        cout<<"s0="<<s0<<endl;
                    }
                }
            }
        }
//        cout<<s0<<s1<<endl;
        if(s0==0 && s1==0){
            return frame1d_;
        }
        if(s0==0 && s1==1){
            return frame1d_01_;
        }
        if(s0==1 && s1==0){
            return frame1d_10_;
        }
        //to do: extensions for domains such as big cubes, free boundary conditions for interface frames
        
        //no extension
        return frame1d_11_;
    }
    
//    //! size of Delta_j
    const int Deltasize(const int j) const;
//
//    //! sizes of the different quarklet index sets
//    const int Nabla01size(const int j) const;
//    const int Nabla10size(const int j) const;
//    const int Nabla11size(const int j) const;

    /*!
      The following routines provide read access to the diverse refinement (sub)matrices
      on a level j >= j0. The row and column indices follow the less<Index> ordering.
      Those matrices will be collected in an internal cache to provide faster access.
    */
//    const BlockMatrix<double>&  get_Mj0      (const int j) const;
//    const BlockMatrix<double>&  get_Mj0T     (const int j) const;
//    const SparseMatrix<double>& get_Mj1c_1d  (const int j) const; // initial stable completion in 1D
//    const BlockMatrix<double>&  get_Mj1c_01  (const int j) const;
//    const BlockMatrix<double>&  get_Mj1c_10  (const int j) const;
//    const BlockMatrix<double>&  get_Mj1c_11  (const int j) const;

    //! index of first generator on level j >= j0
    Index first_generator(const level_type& j, const polynomial_type& p = polynomial_type()) const;
      
    //! index of last generator on level j >= j0
    Index last_generator(const level_type& j, const polynomial_type& p = polynomial_type()) const;
      
    //! index of first quarklet on level j >= j0
    Index first_quarklet(const level_type& j, const polynomial_type& p, const int& number=-1) const;
      
    //! index of first quarklet on level j >= j0 with type e
//    Index first_quarklet(const level_type& j, const type_type& e, const polynomial_type& p) const;

    //! index of last quarklet on level j >= j0
    Index last_quarklet(const level_type& j, const polynomial_type& p = polynomial_type(), const int& number=-1) const;
    
    //! index of last quarklet on level j >= j0
    Index last_quarklet(const int& levelsum, const polynomial_type& p = polynomial_type(), const int& number=-1) const;
    
//    //! Index of last quarklet on level j >= ||j0||_1
//    Index last_quarklet(const int levelsum, const polynomial_type& p = polynomial_type()) const;

    //! RECONSTRUCT routine, simple version
    /*!
      Constructs for a given single quarklet index lambda a coefficient set c,
      such that
      \psi_lambda = \sum_{\lambda'}c_{\lambda'}\psi_{\lambda'}
      where always |\lambda'|>=j
    */
//    void reconstruct_1(const Index& lambda, const int j,
//		       InfiniteVector<double, Index>& c) const;
//
//    //! RECONSTRUCT routine, full version
//    /*!
//      Constructs for a given coefficient set c another one v,
//      such that
//      \sum_{\lambda}c_\lambda\psi_lambda = \sum_{\lambda'}v_{\lambda'}\psi_{\lambda'}
//      where always |\lambda'|>=j
//    */
//    void reconstruct(const InfiniteVector<double, Index>& c, const int j,
//		     InfiniteVector<double, Index>& v) const;

    /*!
      Evaluate a single primal generator or quarklet \psi_\lambda
      on a dyadic subgrid of the L-shaped domain
    */
    Array1D<SampledMapping<2> >
    evaluate
    (const typename DomainFrame<IntervalFrame, NPATCHES>::Index& lambda,
     const int resolution) const;
    
    /*!
      Evaluate an arbitrary linear combination of primal/dual quarklets
      on a dyadic subgrid of the L-shaped domain
    */
    Array1D<SampledMapping<2> >
    evaluate
    (const InfiniteVector<double, typename DomainFrame<IntervalFrame, NPATCHES>::Index>& coeffs,
     const int resolution) const;

    void set_jpmax(const int jmax, const int pmax) {
      jmax_ = jmax;
      pmax_ = pmax;
      setup_full_collection();
    }
    
    inline const int get_jmax() const {
            return jmax_;
        }
        
    inline const int get_pmax() const {
        return pmax_;
    }
    
    inline const int get_Nablasize() const {
        return Nablasize_;
    }

    inline const Array1D<int> get_first_wavelet_numbers() const {
        return first_wavelet_numbers;
    }
    
    inline const Array1D<int> get_last_wavelet_numbers() const {
        return last_wavelet_numbers;
    }
    
    //! get the quarklet index corresponding to a specified number
    const inline Index* get_quarklet (const int number) const {
      return &full_collection[number];
    }
    
    inline const bool get_setup_full_collection() const {
            return setup_full_collection_;
        }
    
    const inline Support& get_support (const int number) const {
        return all_supports_[number];
    }
    
    

    //! number of quarklets between coarsest and finest level
    const int degrees_of_freedom() const { return full_collection.size(); };
//    typedef std::map<int,Support> SupportCache;
//    const SupportCache suppcache() const{ return supp_cache;};
    
    //alternative for SupportCache
    Array1D<Support> all_supports_;


  protected:
      
    /*
    *  Collection of first and last wavelet numbers on all levels up to jmax
    *  Precomputed for speedup
    */
    Array1D<int> first_wavelet_numbers, last_wavelet_numbers;
//    Array1D<Array1D<int> > first_quarklet_numbers, last_quarklet_numbers, first_quark_numbers;
      
    //! Coarsest possible level j0
    level_type j0_;

    //! finest possible level
    int jmax_;
    
    //! highest possible polynomial
    int pmax_;

    //! setup full collectin of quarklets between j0_ and jmax_ as long as a jmax_ has been specified
    void setup_full_collection();
    
    //! setup_full_collection_ is set to 1 after initialising full collection
    bool setup_full_collection_;

    //! collection of all quarklets between coarsest and finest level
    Array1D<Index> full_collection;

    //! the interval 1d quarklet frame
    const IntervalFrame* frame1d_;
    const IntervalFrame* frame1d_11_;
    const IntervalFrame* frame1d_01_;
    const IntervalFrame* frame1d_10_;
    
    
    const Array1D<Point<2, int> > corners_;
    const Array1D<FixedArray1D<int,2> > extensions_;
    
    //! Degrees of freedom on level p=(0,0) 
    int Nablasize_;

//    //! support cache
////    typedef std::map<Index,Support> SupportCache;
//    mutable SupportCache supp_cache;
    
    
    
    bool precomputed_supports_;

  };
  


}

#include <general_domain/domain_frame.cpp>


#endif
