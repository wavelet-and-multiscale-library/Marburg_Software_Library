// implementation for error_H_scale.h

namespace FrameTL
{
  template <class IBASIS>
  double
  error_H_scale_interval (const int order,
			  const AggregatedFrame<IBASIS,1,1>& frame,
			  const InfiniteVector<double, typename AggregatedFrame<IBASIS,1,1>::Index>& coeffs,
			  const Function<1>& f) {
    double res = 0.;


    typedef typename AggregatedFrame<IBASIS,1,1>::Index Index;
    //typedef typename AggregatedFrame<IBASIS,1,1>::Support SuppType;

    const int granularity = 100;
    const double dx = 1./granularity;

    double approx = 0.;
    for (int i = 0; i < granularity; i++) {
      Point<1> p1(i*dx + 0.5*dx);
      for(typename InfiniteVector<double,Index>::const_iterator it2(coeffs.begin()),
	    itend(coeffs.end()); it2 != itend; ++it2) {
	
	//if (in_support(frame, it2.index(), p)) {
	Point<1> p_0_1;
	frame.atlas()->charts()[it2.index().p()]->map_point_inv(p1,p_0_1);

	double wavVal = 0.;
	switch (order) {
	  // L_2 error
	case 0: {
	  if (in_support(frame, it2.index(), p1)) {
	    wavVal = WaveletTL::evaluate(*(frame.bases()[it2.index().p()]->bases()[0]), 0,
					 typename IBASIS::Index(it2.index().j(),
								it2.index().e()[0],
								it2.index().k()[0],
								frame.bases()[it2.index().p()]->bases()[0]),
					 
					 p_0_1[0]);
	  }
	  
	  approx += (*it2)*( wavVal / frame.atlas()->charts()[(it2.index()).p()]->Gram_factor(p_0_1) );
	  break;
	}// end case 0
	case 1: {
	  if (in_support(frame, it2.index(), p1)) {
	    wavVal = WaveletTL::evaluate(*(frame.bases()[it2.index().p()]->bases()[0]), 1,
					  typename IBASIS::Index(it2.index().j(),
								 it2.index().e()[0],
								 it2.index().k()[0],
								 frame.bases()[it2.index().p()]->bases()[0]),
					 
					  p_0_1[0]);

	    approx += (*it2)*(wavVal
			       /
			       (
				frame.atlas()->charts()[(it2.index()).p()]->Gram_factor(p_0_1)
				*
				frame.atlas()->charts()[it2.index().p()]->a_i(0)
				)
			       );
	  }
	  break;
	}// end case 1
	case 2: {
	  // not implemented
	}// end case 2
	}// end switch
      }
      approx -= f.value(p1);
      approx *= approx;
      res += approx*dx;
      
      approx = 0;
    }

    return sqrt(res);  

  }


//   template <class IBASIS, int DIM>
//   double
//   error_H_scale_interval (const int order,
// 			  const AggregatedFrame<IBASIS,DIM,DIM>& frame,
// 			  const InfiniteVector<double, typename AggregatedFrame<IBASIS,DIM,DIM>::Index>& coeffs,
// 			  const Function<1>& f) {
//     double res = 0.;
//     const int N_Gauss = 4;

//     typedef typename AggregatedFrame<IBASIS,1,1>::Index Index;
    
//     list<double> grid;
    
//     typedef typename AggregatedFrame<IBASIS,1,1>::Support SuppType;
    

//     // setup grid for composite quadrature
//      for(typename InfiniteVector<double,Index>::const_iterator it(coeffs.begin()),
// 	   itend(coeffs.end()); it != itend; ++it) {
//        const SuppType* supp = &(frame.all_patch_supports[it.index().number()]);

//        // ?
//        int res = it.index().j()==frame.j0() ? it.index().j() : it.index().j();

//        double dx = 1./(1 << res);
//        int k = 0;
//        while (! (supp->a[0]+k*dx > supp->b[0])) {
// 	 grid.push_back(supp->a[0]+k*dx);
// 	 k++;
//        }
//        // to be on the save side
//        grid.push_back(supp->b[0]);
//      }

//      grid.push_back(0.5);
//      grid.sort();
//      grid.unique();



//      double u = 0.;
//      double v = 0.;
//      double approx = 0.;

//      for (typename list<double>::const_iterator it(grid.begin()),
// 	   itend(grid.end()); it != itend; ++it) {
//        u = *it;
//        ++it;
//        if (it != itend)
// 	 v = *it; 
//        else
// 	 break;

//        Array1D<double> gauss_points, gauss_weights;
//        gauss_points.resize(N_Gauss);
//        gauss_weights.resize(N_Gauss);
       
//        for (int n = 0; n < N_Gauss ; n++) {
// 	 gauss_points[n] = 0.5*((v-u)*GaussPoints[N_Gauss-1][n]+u+v);
// 	 gauss_weights[n] = (v-u)*GaussWeights[N_Gauss-1][n];
//        }


//        for (int n = 0; n < N_Gauss ; n++) {
// 	 Point<1> p(gauss_points[n]);
// 	 for(typename InfiniteVector<double,Index>::const_iterator it2(coeffs.begin()),
// 	       itend(coeffs.end()); it2 != itend; ++it2) {
	   
// 	   if (in_support(frame, it2.index(), p)) {
// 	     Point<1> p_0_1;
// 	     frame.atlas()->charts()[it2.index().p()]->map_point_inv(p,p_0_1);
// 	     double wavVal = 0.;
// 	     switch (order) {
// 	       // L_2 error
// 	     case 0: {
// 	       wavVal = WaveletTL::evaluate(*(frame.bases()[it2.index().p()]->bases()[0]), 0,
// 						   typename IBASIS::Index(it2.index().j(),
// 									  it2.index().e()[0],
// 									  it2.index().k()[0],
// 									  frame.bases()[it2.index().p()]->bases()[0]),
						   
// 						   p_0_1[0]);
	       
// 	       approx += (*it2)*( wavVal / frame.atlas()->charts()[(it2.index()).p()]->Gram_factor(p_0_1) );
// 	       break;
// 	     }
// 	       // H^1 error
// 	     case 1: {
// 	       wavVal = WaveletTL::evaluate(*(frame.bases()[it2.index().p()]->bases()[0]), 1,
// 						   typename IBASIS::Index(it2.index().j(),
// 									  it2.index().e()[0],
// 									  it2.index().k()[0],
// 									  frame.bases()[it2.index().p()]->bases()[0]),
						   
// 						   p_0_1[0]);
// 	       approx += (*it2)*(wavVal
// 				 /
// 				 (
// 				  frame.atlas()->charts()[(it2.index()).p()]->Gram_factor(p_0_1)
// 				  *
// 				  frame.atlas()->charts()[it2.index().p()]->a_i(0)
// 				  )
// 				 );
// 	       break;
// 	     }
// 	       // H^2 error
// 	     case 2: {
// 	       wavVal = WaveletTL::evaluate(*(frame.bases()[it2.index().p()]->bases()[0]), 2,
// 						   typename IBASIS::Index(it2.index().j(),
// 									  it2.index().e()[0],
// 									  it2.index().k()[0],
// 									  frame.bases()[it2.index().p()]->bases()[0]),
						   
// 						   p_0_1[0]);

// 	       approx += (*it2)*(wavVal
// 				 /
// 				 (
// 				  frame.atlas()->charts()[(it2.index()).p()]->Gram_factor(p_0_1)
// 				  *
// 				  frame.atlas()->charts()[it2.index().p()]->a_i(0)*frame.atlas()->charts()[it2.index().p()]->a_i(0)
// 				  )
// 				 );
// 	     }
//        	     }
// 	   }// end if
// 	 }
// 	 approx -= f.value(p);
// 	 approx *= approx;
// 	 res += gauss_weights[n]*approx;
// 	 approx = 0;
//        }
//      }     
//      return sqrt(res);  

//   }


  template <class IBASIS>
  double
  error_H_scale_Lshaped(const int order,
			const AggregatedFrame<IBASIS,2,2>& frame,
			const InfiniteVector<double, typename AggregatedFrame<IBASIS,2,2>::Index>& coeffs,
			const Function<2>& f)
  {

    // split L shaped into three equally large cubes and apply a composite midpoint rule
    double res = 0.;
    Vector<double> values(2);

    typedef AggregatedFrame<IBASIS,2,2> Frame2D;
    typedef typename Frame2D::Index Index;
    //typedef typename Frame2D::Support SuppType;
    
    const int granularity = 10; // formerly 200 (for the tests in the PhD thesis)
    
    
    double approx = 0.;
    Vector<double> approx_der(2);
    approx_der *= 0.;

    Vector<double> values_grad_wav(2);
    values_grad_wav *= 0.;
    Vector<double> values_grad_f(2);
    values_grad_f *= 0.;


    // first cube: (-1,0)^2
    Point<2> a(-1,-1);
    Point<2> b(0,0);
    double dx = (b[0]-a[0]) / granularity;
    double dy = (b[1]-a[1]) / granularity;
    for (int i = 0; i < granularity; i++) {
      for (int j = 0; j < granularity; j++) {
	Point<2> p(a[0]+i*dx+dx*0.5, a[1]+j*dy+dy*0.5);
	 for(typename InfiniteVector<double,Index>::const_iterator it2(coeffs.begin()),
	    itend(coeffs.end()); it2 != itend; ++it2) {
	   if (in_support(frame, it2.index(), p)) {
	     switch (order) {
	       // L_2 error
	     case 0: {
	       approx += *it2 * frame.evaluate(it2.index(), p);
	       break;
	     }
	     case 1: {
	       frame.evaluate_gradient(it2.index(), p, values_grad_wav);
	       values_grad_wav *= *it2;
	       approx_der += values_grad_wav;
	       break;
	     }
	     }
	   }
	 }// end for coefficients
	 switch (order) {
	   // L_2 error
	 case 0: {
	   approx -= f.value(p);
	   approx *= approx * dx * dy;
	   res += approx;
	   approx = 0.;
	   break;
	 }
	 case 1: {
	   f.vector_value(p, values_grad_f);
	   res += l2_norm_sqr(approx_der-values_grad_f) * dx * dy;
	   approx_der *= 0.;
	   break;
	 }
	 }
      }
    }

    // second cube: (0,1)x(-1,0)
    a[0] = 0;
    a[1] = -1;
    b[0] = 1;
    b[1] = 0;
    dx = (b[0]-a[0]) / granularity;
    dy = (b[1]-a[1]) / granularity;
    for (int i = 0; i < granularity; i++) {
      for (int j = 0; j < granularity; j++) {
	Point<2> p(a[0]+i*dx+dx*0.5, a[1]+j*dy+dy*0.5);
	 for(typename InfiniteVector<double,Index>::const_iterator it2(coeffs.begin()),
	    itend(coeffs.end()); it2 != itend; ++it2) {
	   if (in_support(frame, it2.index(), p)) {
	     switch (order) {
	       // L_2 error
	     case 0: {
	       approx += *it2 * frame.evaluate(it2.index(), p);
	       break;
	     }
	     case 1: {
	       frame.evaluate_gradient(it2.index(), p, values_grad_wav);
	       values_grad_wav *= *it2;
	       approx_der += values_grad_wav;
	       break;
	     }
	     }
	   }
	 }// end for coefficients
	 switch (order) {
	   // L_2 error
	 case 0: {
	   approx -= f.value(p);
	   approx *= approx * dx * dy;
	   res += approx;
	   approx = 0.;
	   break;
	 }
	 case 1: {
	   f.vector_value(p, values_grad_f);
	   res += l2_norm_sqr(approx_der-values_grad_f) * dx * dy;
	   approx_der *= 0.;
	   break;
	 }
	 }
      }
   }

    // third cube: (-1,0)x(0,1)
    a[0] = -1;
    a[1] = 0;
    b[0] = 0;
    b[1] = 1;
    dx = (b[0]-a[0]) / granularity;
    dy = (b[1]-a[1]) / granularity;
    for (int i = 0; i < granularity; i++) {
      for (int j = 0; j < granularity; j++) {
	Point<2> p(a[0]+i*dx+dx*0.5, a[1]+j*dy+dy*0.5);
	for(typename InfiniteVector<double,Index>::const_iterator it2(coeffs.begin()),
	    itend(coeffs.end()); it2 != itend; ++it2) {
	   if (in_support(frame, it2.index(), p)) {
	     switch (order) {
	       // L_2 error
	     case 0: {
	       approx += *it2 * frame.evaluate(it2.index(), p);
	       break;
	     }
	     case 1: {
	       frame.evaluate_gradient(it2.index(), p, values_grad_wav);
	       values_grad_wav *= *it2;
	       approx_der += values_grad_wav;
	       break;
	     }
	     }
	   }
	 }// end for coefficients
	 switch (order) {
	   // L_2 error
	 case 0: {
	   approx -= f.value(p);
	   approx *= approx * dx * dy;
	   res += approx;
	   approx = 0.;
	   break;
	 }
	 case 1: {
	   f.vector_value(p, values_grad_f);
	   res += l2_norm_sqr(approx_der-values_grad_f) * dx * dy;
	   approx_der *= 0.;
	   break;
	 }
	 }
      }
    }


    return sqrt(res);
  }

  double
  H_1_semi_norm_Lshaped(const Function<2>& gradient)
  {
    double res = 0.;
    // split L shaped into three equally large cubes and apply a composite midpoint rule
    const int granularity = 1024;

    Vector<double> values(2);
    // first cube: (-1,0)^2
    Point<2> a(-1,-1);
    Point<2> b(0,0);
    //Grid grid(a,b,granularity);
    double dx = (b[0]-a[0]) / granularity;
    double dy = (b[1]-a[1]) / granularity;
    for (int i = 0; i < granularity; i++) {
      for (int j = 0; j < granularity; j++) {
	Point<2> x(a[0]+i*dx+dx*0.5, a[1]+j*dy+dy*0.5);
	gradient.vector_value(x,values);
	res += ((values[0]*values[0])+(values[1]*values[1]))*dx*dy;
      }
    }
    
    // second cube: (0,1)x(-1,0)
    a[0] = 0;
    a[1] = -1;
    b[0] = 1;
    b[1] = 0;
    dx = (b[0]-a[0]) / granularity;
    dy = (b[1]-a[1]) / granularity;
    for (int i = 0; i < granularity; i++) {
      for (int j = 0; j < granularity; j++) {
	Point<2> x(a[0]+i*dx+dx*0.5, a[1]+j*dy+dy*0.5);
	gradient.vector_value(x,values);
	res += ((values[0]*values[0])+(values[1]*values[1]))*dx*dy;
      }
    }

    // third cube: (-1,0)x(0,1)
    a[0] = -1;
    a[1] = 0;
    b[0] = 0;
    b[1] = 1;

    dx = (b[0]-a[0]) / granularity;
    dy = (b[1]-a[1]) / granularity;
    for (int i = 0; i < granularity; i++) {
      for (int j = 0; j < granularity; j++) {
	Point<2> x(a[0]+i*dx+dx*0.5, a[1]+j*dy+dy*0.5);
	gradient.vector_value(x,values);
	res += ((values[0]*values[0])+(values[1]*values[1]))*dx*dy;
      }
    }
    return sqrt(res);
  }


}
