// implementation for cube_evaluate.h

#include <utils/array1d.h>
#include <utils/fixed_array1d.h>
#include <geometry/point.h>
#include <geometry/grid.h>
#include <geometry/sampled_mapping.h>

using namespace MathTL;

namespace WaveletTL
{
  template <class IBASIS, unsigned int DIM>
  SampledMapping<DIM>
  evaluate(const CubeBasis<IBASIS,DIM>& basis,
	   const typename CubeBasis<IBASIS,DIM>::Index& lambda,
	   const bool primal,
	   const int resolution)
  {
    FixedArray1D<Array1D<double>,DIM> values; // point values of the factors within psi_lambda
    for (unsigned int i = 0; i < DIM; i++)
      values[i] = evaluate(*(basis.bases()[i]),
			   typename IBASIS::Index(lambda.j(),
						  lambda.e()[i],
						  lambda.k()[i],
						  basis.bases()[i]),
			   primal,
			   resolution).values();

    SampledMapping<DIM> result(Point<DIM>(0), Point<DIM>(1), values);
    return result; // gcc 2.95 does not like these two lines melted into one
  }


  template <class IBASIS, unsigned int DIM>
  SampledMapping<DIM>
  evaluate(const CubeBasis<IBASIS,DIM>& basis,
	   const InfiniteVector<double, typename CubeBasis<IBASIS,DIM>::Index>& coeffs,
	   const bool primal,
	   const int resolution)
  {
    Grid<DIM> grid(Point<DIM>(0), Point<DIM>(1), 1<<resolution);
    SampledMapping<DIM> result(grid); // zero

    typedef typename CubeBasis<IBASIS,DIM>::Index Index;
    for (typename InfiniteVector<double,Index>::const_iterator it(coeffs.begin()),
	   itend(coeffs.end()); it != itend; ++it)
      result.add(*it, evaluate(basis, it.index(), primal, resolution));

    return result;
  }


  template <class IBASIS, unsigned int DIM>
  double
  evaluate(const CubeBasis<IBASIS,DIM>& basis,
	   const typename CubeBasis<IBASIS,DIM>::Index& lambda,
	   const Point<DIM> p)
  {
    double r(1);
    for (unsigned int i = 0; i < DIM; i++)
    {
      r *= evaluate(*(basis.bases()[i]), 0,
			   typename IBASIS::Index(lambda.j(),
						  lambda.e()[i],
						  lambda.k()[i],
						  basis.bases()[i]),
			              p[i]);
    }
    return r;
  }


  template <class IBASIS, unsigned int DIM>
  double
  evaluate(const CubeBasis<IBASIS,DIM>& basis,
       const InfiniteVector<double, typename CubeBasis<IBASIS,DIM>::Index>& coeffs,
	   const Point<DIM> p)
  {
    double r(0);
    for (typename InfiniteVector<double, typename CubeBasis<IBASIS,DIM>::Index>::const_iterator it(coeffs.begin()),
         itend(coeffs.end()); it != itend; ++it)
    {
        r += ( (*it) * evaluate(basis, it.index(), p) );
    }
    return r;
  }


  template <class IBASIS, unsigned int DIM>
  void
  evaluate(const CubeBasis<IBASIS,DIM>& basis,
	   const typename CubeBasis<IBASIS,DIM>::Index& lambda,
	   const Point<DIM> p,
	   FixedArray1D<double,DIM>& deriv_values)
  {
    for (unsigned int i = 0; i < DIM; i++)
    {
      deriv_values[i] = 1.;
      for (unsigned int j = 0; j < DIM; j++)
      {
        if (i != j)
        {
          deriv_values[i] *= evaluate(*(basis.bases()[j]), 0,
			   typename IBASIS::Index(lambda.j(),
						  lambda.e()[j],
						  lambda.k()[j],
						  basis.bases()[j]),
			              p[j]);
        }
        else
        {
          deriv_values[i] *= evaluate(*(basis.bases()[j]), 1,
			   typename IBASIS::Index(lambda.j(),
						  lambda.e()[j],
						  lambda.k()[j],
						  basis.bases()[j]),
			              p[j]);
        }
      }
    }
  }


  template <class IBASIS, unsigned int DIM>
  void
  evaluate(const CubeBasis<IBASIS,DIM>& basis,
           const InfiniteVector<double, typename CubeBasis<IBASIS,DIM>::Index>& coeffs,
	       const Point<DIM> p,
	       FixedArray1D<double,DIM>& deriv_values)
  {
    // set deriv_values to zero
    for (int i = 0; i < DIM; i++)
    {
      deriv_values[i] = 0.;
    }
    // dummy variable
    FixedArray1D<double,DIM> share;

    // iterate over all coeffs and sum the share
    for (typename InfiniteVector<double, typename CubeBasis<IBASIS,DIM>::Index>::const_iterator it(coeffs.begin()),
         itend(coeffs.end()); it != itend; ++it)
    {
        evaluate(basis, it.index(), p, share);
        //evaluate_tuned(basis, it.index(), p, share);    // tuned Version (Christoph)
        for (int i = 0; i < DIM; i++)
        {
          deriv_values[i] += ( (*it) * share[i] );
        }
    }
  }


  /*! C.Hartmann: added on 25.09.15
    Evaluate a single primal generator or wavelet \psi_\lambda
    on an arbitrary grid of [0,1]^d. The grid is given by a tensor product
    of 1-dimensional grids on [0,1].
    ATTENTION: Only dimension = 2 and evaluation of primal generators/wavelets supported a.t.m.
  */
  template <class IBASIS, unsigned int DIM>
  Matrix<double>
  evaluate(const CubeBasis<IBASIS,DIM>& basis,
	   const typename CubeBasis<IBASIS,DIM>::Index& lambda,
	   const FixedArray1D<Array1D<double>,DIM>& grid)
  {
    Matrix<double> r(grid[0].size(), grid[1].size());

    // values of Psi_x on x-grid
    MathTL::Array1D<double> values_x(grid[0].size());
    for (unsigned int i(0); i < values_x.size(); i++)
 	{
      values_x[i] = evaluate_primbs_22_bc11_v3(lambda.j(), lambda.e()[0], lambda.k()[0], grid[0][i]);
    }

    // values of Psi_y on y-grid
    MathTL::Array1D<double> values_y(grid[1].size());
    for (unsigned int i(0); i < values_y.size(); i++)
 	{
      values_y[i] = evaluate_primbs_22_bc11_v3(lambda.j(), lambda.e()[1], lambda.k()[1], grid[1][i]);
    }

    // compute tensor product
    for (unsigned int m(0); m < r.row_dimension(); m++)
      for (unsigned int n(0); n < r.column_dimension(); n++)
	    r(m,n) = values_x[m] * values_y[n];

    return r;
  }


  /*! C.Hartmann: added on 25.09.15
    Evaluate an arbitrary linear combination of primal wavelets
    on an arbitrary grid of [0,1]^d. The grid is given by a tensor product
    of 1-dimensional grids on [0,1].
    ATTENTION: Only dimension = 2 and evaluation of primal generators/wavelets supported a.t.m.
  */
  template <class IBASIS, unsigned int DIM>
  Matrix<double>
  evaluate(const CubeBasis<IBASIS,DIM>& basis,
           const InfiniteVector<double, typename CubeBasis<IBASIS,DIM>::Index>& coeffs,
           const FixedArray1D<Array1D<double>,DIM>& grid)
  {
    Matrix<double> r(grid[0].size(), grid[1].size());
    Matrix<double> dummy(grid[0].size(), grid[1].size());

    typedef typename CubeBasis<IBASIS,DIM>::Index Index;
    for (typename InfiniteVector<double,Index>::const_iterator it(coeffs.begin()),
	   itend(coeffs.end()); it != itend; ++it)
    {
      dummy = evaluate(basis, it.index(), grid);

      for (unsigned int m(0); m < r.row_dimension(); m++)
      {
        for (unsigned int n(0); n < r.column_dimension(); n++)
        {
	      r(m,n) += ( (*it) * dummy(m,n) );
        }
      }

    }

    return r;
  }


  /*! C.Hartmann: added on 26.09.15
    Evaluate the first derivative of a single primal generator or wavelet \psi_\lambda
    on an arbitrary grid of [0,1]^d. The grid is given by a tensor product of 1-dimensional grids on [0,1].
    Function values are stored in the matrices deriv_x and deriv_y. Matrix dimensions must match grid dimension!
    ATTENTION: Only dimension = 2 and evaluation of primal generators/wavelets supported a.t.m.
  */
  template <class IBASIS, unsigned int DIM>
  void
  evaluate(const CubeBasis<IBASIS,DIM>& basis,
	   const typename CubeBasis<IBASIS,DIM>::Index& lambda,
	   const FixedArray1D<Array1D<double>,DIM>& grid,
	   Matrix<double>& deriv_x, Matrix<double>& deriv_y)
  {

    // values of Psi_x on x-grid
    MathTL::Array1D<double> Psi_x(grid[0].size());
    for (unsigned int i(0); i < Psi_x.size(); i++)
 	{
      Psi_x[i] = evaluate_primbs_22_bc11_v3(lambda.j(), lambda.e()[0], lambda.k()[0], grid[0][i]);
    }

    // values of d/dx Psi_x on x-grid
    MathTL::Array1D<double> Psi_x_deriv(grid[0].size());
    for (unsigned int i(0); i < Psi_x_deriv.size(); i++)
 	{
      Psi_x_deriv[i] = evaluate_primbs_22_bc11_deriv(lambda.j(), lambda.e()[0], lambda.k()[0], grid[0][i]);
    }

    // values of Psi_y on y-grid
    MathTL::Array1D<double> Psi_y(grid[1].size());
    for (unsigned int i(0); i < Psi_y.size(); i++)
 	{
      Psi_y[i] = evaluate_primbs_22_bc11_v3(lambda.j(), lambda.e()[1], lambda.k()[1], grid[1][i]);
    }

    // values of d/dy Psi_y on y-grid
    MathTL::Array1D<double> Psi_y_deriv(grid[1].size());
    for (unsigned int i(0); i < Psi_y_deriv.size(); i++)
 	{
      Psi_y_deriv[i] = evaluate_primbs_22_bc11_deriv(lambda.j(), lambda.e()[1], lambda.k()[1], grid[1][i]);
    }

    // compute tensor products
    for (unsigned int m(0); m < grid[0].size(); m++)
    {
      for (unsigned int n(0); n < grid[1].size(); n++)
      {
	    deriv_x(m,n) = Psi_x_deriv[m] * Psi_y[n];
	    deriv_y(m,n) = Psi_x[m] * Psi_y_deriv[n];
      }
    }

  }



  /*! C.Hartmann: added on 26.09.15
    Evaluate the first derivative of an arbitrary linear combination of primal wavelets
    on an arbitrary grid of [0,1]^d. The grid is given by a tensor product of 1-dimensional grids on [0,1].
    Function values are stored in the matrices deriv_x and deriv_y. Matrix dimensions must match grid dimension!
    ATTENTION: Only dimension = 2 and evaluation of primal generators/wavelets supported a.t.m.
  */
  template <class IBASIS, unsigned int DIM>
  void
  evaluate(const CubeBasis<IBASIS,DIM>& basis,
           const InfiniteVector<double, typename CubeBasis<IBASIS,DIM>::Index>& coeffs,
           const FixedArray1D<Array1D<double>,DIM>& grid,
           Matrix<double>& deriv_x, Matrix<double>& deriv_y)
  {
    Matrix<double> dummy_x(grid[0].size(), grid[1].size());
    Matrix<double> dummy_y(grid[0].size(), grid[1].size());

    deriv_x.resize(grid[0].size(), grid[1].size());
    deriv_y.resize(grid[0].size(), grid[1].size());

    typedef typename CubeBasis<IBASIS,DIM>::Index Index;
    for (typename InfiniteVector<double,Index>::const_iterator it(coeffs.begin()),
	   itend(coeffs.end()); it != itend; ++it)
    {
      evaluate(basis, it.index(), grid, dummy_x, dummy_y);

      for (unsigned int m(0); m < grid[0].size(); m++)
      {
        for (unsigned int n(0); n < grid[1].size(); n++)
        {
	      deriv_x(m,n) += ( (*it) * dummy_x(m,n) );
	      deriv_y(m,n) += ( (*it) * dummy_y(m,n) );
        }
      }

    }

  }


  //! tuned Versions (Christoph)

  template <class IBASIS, unsigned int DIM>
  void
  evaluate_tuned(const CubeBasis<IBASIS,DIM>& basis,
	   const typename CubeBasis<IBASIS,DIM>::Index& lambda,
	   const Point<DIM> p,
	   FixedArray1D<double,DIM>& deriv_values)
  {
    for (unsigned int i = 0; i < DIM; i++)
    {
      deriv_values[i] = 1.;
      for (unsigned int j = 0; j < DIM; j++)
      {
        if (i != j)
        {
          deriv_values[i] *= evaluate(*(basis.bases()[j]), 0,
                                      lambda.j(), lambda.e()[j], lambda.k()[j],
			                          p[j]);
        }
        else
        {
          deriv_values[i] *= evaluate(*(basis.bases()[j]), 1,
			                          lambda.j(), lambda.e()[j], lambda.k()[j],
			                          p[j]);
        }
      }
    }
  }


  template <class IBASIS, unsigned int DIM>
  void
  evaluate_tuned(const CubeBasis<IBASIS,DIM>& basis,
           const InfiniteVector<double, typename CubeBasis<IBASIS,DIM>::Index>& coeffs,
	       const Point<DIM> p,
	       FixedArray1D<double,DIM>& deriv_values)
  {

    FixedArray1D<double,DIM> share;
    typename CubeBasis<IBASIS,DIM>::Index lambda;
    unsigned int i, j;

    // set deriv_values to zero
    for (i = 0; i < DIM; i++)
    {
      deriv_values[i] = 0.;
    }


    // iterate over all coeffs and sum the share
    for (typename InfiniteVector<double, typename CubeBasis<IBASIS,DIM>::Index>::const_iterator it(coeffs.begin()),
         itend(coeffs.end()); it != itend; ++it)
    {
        lambda = it.index();

        for (i = 0; i < DIM; i++)
        {
          share[i] = 1.;
          for (j = 0; j < DIM; j++)
          {
            if (i != j)
            {
              share[i] *= evaluate(*(basis.bases()[j]), 0,
                                          lambda.j(), lambda.e()[j], lambda.k()[j],
			                              p[j]);
            }
            else
            {
              share[i] *= evaluate(*(basis.bases()[j]), 1,
			                              lambda.j(), lambda.e()[j], lambda.k()[j],
			                              p[j]);
            }
          }
        }

        for (i = 0; i < DIM; i++)
        {
          deriv_values[i] += ( (*it) * share[i] );
        }
    }
  }


  template <class IBASIS>
  void
  evaluate_tuned_dim2(const CubeBasis<IBASIS,2>& basis,
           const InfiniteVector<double, typename CubeBasis<IBASIS,2>::Index>& coeffs,
	       const Point<2> p,
	       FixedArray1D<double,2>& deriv_values)
  {
    typename CubeBasis<IBASIS,2>::Index lambda;
    int k1, k2, scale;
    // set deriv_values to zero
    deriv_values[0] = 0.;
    deriv_values[1] = 0.;

    clock_t begin_calc_deriv_outer = clock();

    // iterate over all coeffs and sum the share
    for (typename InfiniteVector<double, typename CubeBasis<IBASIS,2>::Index>::const_iterator it(coeffs.begin()),
         itend(coeffs.end()); it != itend; ++it)
    {
        lambda = it.index();

        // check if the support of Psi_x contains x-value of p
        (basis.bases()[0])->support(lambda.j(), lambda.e()[0], lambda.k()[0], k1, k2);
        scale = 1u << (lambda.j() + lambda.e()[0]);
        if ( ((scale*p[0]) < k1) || ((scale*p[0]) > k2) )
        {
          continue;
        }

        // check if the support of Psi_y contains y-value of p
        (basis.bases()[1])->support(lambda.j(), lambda.e()[1], lambda.k()[1], k1, k2);
        scale = 1u << (lambda.j() + lambda.e()[1]);
        if ( ((scale*p[1]) < k1) || ((scale*p[1]) > k2) )
        {
          continue;
        }

        // measure time consumption of function calls "evaluate"
        clock_t begin_calc_deriv = clock();

        // sum the share
        deriv_values[0] += ( (*it) * evaluate(*(basis.bases()[0]), 1, lambda.j(), lambda.e()[0], lambda.k()[0], p[0]) *
                                     evaluate(*(basis.bases()[1]), 0, lambda.j(), lambda.e()[1], lambda.k()[1], p[1]) );

        deriv_values[1] += ( (*it) * evaluate(*(basis.bases()[0]), 0, lambda.j(), lambda.e()[0], lambda.k()[0], p[0]) *
                                     evaluate(*(basis.bases()[1]), 1, lambda.j(), lambda.e()[1], lambda.k()[1], p[1]) );

        // measure time consumption of function calls "evaluate"
        clock_t end_calc_deriv = clock();
        double elapsed_secs = double(end_calc_deriv - begin_calc_deriv) / CLOCKS_PER_SEC;
        time_consumption_of_calc_deriv += elapsed_secs;

    }


    clock_t end_calc_deriv_outer = clock();
    double elapsed_secs_outer = double(end_calc_deriv_outer - begin_calc_deriv_outer) / CLOCKS_PER_SEC;
    time_consumption_of_calc_deriv_outer += elapsed_secs_outer;

  }


  template <class IBASIS>
  void
  evaluate_tuned_dim2_lin_primbs(const CubeBasis<IBASIS,2>& basis,
           const InfiniteVector<double, typename CubeBasis<IBASIS,2>::Index>& coeffs,
	       const Point<2> p,
	       FixedArray1D<double,2>& deriv_values)
  {
    typename CubeBasis<IBASIS,2>::Index lambda;

    int scale;
    int k1_h, k2_h;
    // set deriv_values to zero
    deriv_values[0] = 0.;
    deriv_values[1] = 0.;

#if 1
    clock_t begin_calc_deriv_outer = clock();
#endif

    // iterate over all coeffs and sum the share
    for (typename InfiniteVector<double, typename CubeBasis<IBASIS,2>::Index>::const_iterator it(coeffs.begin()),
         itend(coeffs.end()); it != itend; ++it)
    {
        lambda = it.index();
        scale = 1u << lambda.j();

        // check if the support of Psi_x contains x-value of p
        // (basis.bases()[0])->support(lambda.j(), lambda.e()[0], lambda.k()[0], k1, k2);

        // compute the support of Psi_lambda (HARD CODED, ONLY FOR LINEAR PRIMBS BASIS WITH HOMOGENEOUS BC's!!!)
        if (lambda.e()[0])   // Wavelet
        {
          if (lambda.k()[0] > 0)
          {
            if (lambda.k()[0] < (scale-1))  // inner wavelet
            {
              k1_h = lambda.k()[0] - 1;
              k2_h = lambda.k()[0] + 2;
            }
            else        // right boundary wavelet
            {
              k1_h = scale - 3;
              k2_h = scale;
            }
          }
          else          // left boundary wavelet
          {
            k1_h = 0;
            k2_h = 3;
          }
        }
        else                // Generator
        {
          k1_h = lambda.k()[0] - 1;
          k2_h = lambda.k()[0] + 1;
        }


        if ( ((scale*p[0]) < k1_h) || ((scale*p[0]) > k2_h) )
        {
          continue;
        }



        // check if the support of Psi_y contains y-value of p
        // (basis.bases()[1])->support(lambda.j(), lambda.e()[1], lambda.k()[1], k1, k2);

        // compute the support of Psi_lambda (HARD CODED, ONLY FOR LINEAR PRIMBS BASIS WITH HOMOGENEOUS BC's!!!)
        if (lambda.e()[1])   // Wavelet
        {
          if (lambda.k()[1] > 0)
          {
            if (lambda.k()[1] < (scale-1))  // inner wavelet
            {
              k1_h = lambda.k()[1] - 1;
              k2_h = lambda.k()[1] + 2;
            }
            else        // right boundary wavelet
            {
              k1_h = scale - 3;
              k2_h = scale;
            }
          }
          else          // left boundary wavelet
          {
            k1_h = 0;
            k2_h = 3;
          }
        }
        else                // Generator
        {
          k1_h = lambda.k()[1] - 1;
          k2_h = lambda.k()[1] + 1;
        }


        if ( ((scale*p[1]) < k1_h) || ((scale*p[1]) > k2_h) )
        {
          continue;
        }


#if 0
        // measure time consumption of function calls "evaluate"
        clock_t begin_calc_deriv = clock();
#endif

#if 0
        // sum the share
        deriv_values[0] += ( (*it) * evaluate(*(basis.bases()[0]), 1, lambda.j(), lambda.e()[0], lambda.k()[0], p[0]) *
                                     evaluate(*(basis.bases()[1]), 0, lambda.j(), lambda.e()[1], lambda.k()[1], p[1]) );

        deriv_values[1] += ( (*it) * evaluate(*(basis.bases()[0]), 0, lambda.j(), lambda.e()[0], lambda.k()[0], p[0]) *
                                     evaluate(*(basis.bases()[1]), 1, lambda.j(), lambda.e()[1], lambda.k()[1], p[1]) );
#endif

#if 1
// sum the share, USE optimized versions 'evaluate_primbs_22_bc11' and 'evaluate_primbs_22_bc11_deriv'
        deriv_values[0] += ( (*it) * evaluate_primbs_22_bc11_deriv(lambda.j(), lambda.e()[0], lambda.k()[0], p[0]) *
                                     evaluate_primbs_22_bc11_v3(lambda.j(), lambda.e()[1], lambda.k()[1], p[1]) );

        deriv_values[1] += ( (*it) * evaluate_primbs_22_bc11_v3(lambda.j(), lambda.e()[0], lambda.k()[0], p[0]) *
                                     evaluate_primbs_22_bc11_deriv(lambda.j(), lambda.e()[1], lambda.k()[1], p[1]) );
#endif


# if 0
//! sanity check for 'evaluate_primbs_22_bc11' :
double diff_x = evaluate(*(basis.bases()[0]), 0, lambda.j(), lambda.e()[0], lambda.k()[0], p[0]) - evaluate_primbs_22_bc11_v3(lambda.j(), lambda.e()[0], lambda.k()[0], p[0]);
double diff_y = evaluate(*(basis.bases()[1]), 0, lambda.j(), lambda.e()[1], lambda.k()[1], p[1]) - evaluate_primbs_22_bc11_v3(lambda.j(), lambda.e()[1], lambda.k()[1], p[1]);

if (abs(diff_x) > 1e-14)
{
  cout <<  std::setprecision(20);
  cout << "ERROR in evaluate primbs: diff_x = " << diff_x << endl;
  cout << "lambda: " << lambda << endl;
  cout <<  std::setprecision(6);
  exit(1);
}
if (abs(diff_y) > 1e-14)
{
  cout << std::setprecision(20);
  cout << "ERROR in evaluate primbs: diff_y = " << diff_y << endl;
  cout << "lambda: " << lambda << endl;
  cout <<  std::setprecision(6);
  exit(1);
}
#endif

# if 0
//! sanity check for 'evaluate_primbs_22_bc11_deriv' :
double diff_x_deriv = evaluate(*(basis.bases()[0]), 1, lambda.j(), lambda.e()[0], lambda.k()[0], p[0]) - evaluate_primbs_22_bc11_deriv(lambda.j(), lambda.e()[0], lambda.k()[0], p[0]);
double diff_y_deriv = evaluate(*(basis.bases()[1]), 1, lambda.j(), lambda.e()[1], lambda.k()[1], p[1]) - evaluate_primbs_22_bc11_deriv(lambda.j(), lambda.e()[1], lambda.k()[1], p[1]);

if (abs(diff_x_deriv) > 1e-13)
{
  cout <<  std::setprecision(20);
  cout << "ERROR in evaluate primbs derivative: diff_x_deriv = " << diff_x_deriv << endl;
  cout << "lambda: " << lambda << " point: (" << p[0] << ", " << p[1] << ")" << endl;
  cout << "evaluate:           " << evaluate(*(basis.bases()[0]), 1, lambda.j(), lambda.e()[0], lambda.k()[0], p[0]) << endl;
  cout << "evaluate_primbs_22: " << evaluate_primbs_22_bc11_deriv(lambda.j(), lambda.e()[0], lambda.k()[0], p[0]) << endl;
  cout <<  std::setprecision(6);
  exit(1);
}
if (abs(diff_y_deriv) > 1e-13)
{
  cout <<  std::setprecision(20);
  cout << "ERROR in evaluate primbs derivative: diff_y_deriv = " << diff_y_deriv << endl;
  cout << "lambda: " << lambda << " point: (" << p[0] << ", " << p[1] << ")" << endl;
  cout << "evaluate:           " << evaluate(*(basis.bases()[1]), 1, lambda.j(), lambda.e()[1], lambda.k()[1], p[1]) << endl;
  cout << "evaluate_primbs_22: " << evaluate_primbs_22_bc11_deriv(lambda.j(), lambda.e()[1], lambda.k()[1], p[1]) << endl;
  cout <<  std::setprecision(6);
  exit(1);
}
#endif

        // measure time consumption of function calls "evaluate"
#if 0
        clock_t end_calc_deriv = clock();
        double elapsed_secs = double(end_calc_deriv - begin_calc_deriv) / CLOCKS_PER_SEC;
        time_consumption_of_calc_deriv += elapsed_secs;
#endif

    }

#if 1
    clock_t end_calc_deriv_outer = clock();
    double elapsed_secs_outer = double(end_calc_deriv_outer - begin_calc_deriv_outer) / double(CLOCKS_PER_SEC);
    time_consumption_of_calc_deriv_outer += elapsed_secs_outer;
#endif

  }


  template <class IBASIS>
  double
  evaluate_tuned_dim2_lin_primbs(const CubeBasis<IBASIS,2>& basis,
           const InfiniteVector<double, typename CubeBasis<IBASIS,2>::Index>& coeffs,
	       const Point<2> p)
  {
    typename CubeBasis<IBASIS,2>::Index lambda;

    int scale;
    int k1_h, k2_h;
    double r = 0.0;


    // iterate over all coeffs and sum the share
    for (typename InfiniteVector<double, typename CubeBasis<IBASIS,2>::Index>::const_iterator it(coeffs.begin()),
         itend(coeffs.end()); it != itend; ++it)
    {
        lambda = it.index();
        scale = 1u << lambda.j();

        // check if the support of Psi_x contains x-value of p
        // (basis.bases()[0])->support(lambda.j(), lambda.e()[0], lambda.k()[0], k1, k2);

        // compute the support of Psi_lambda (HARD CODED, ONLY FOR LINEAR PRIMBS BASIS WITH HOMOGENEOUS BC's!!!)
        if (lambda.e()[0])   // Wavelet
        {
          if (lambda.k()[0] > 0)
          {
            if (lambda.k()[0] < (scale-1))  // inner wavelet
            {
              k1_h = lambda.k()[0] - 1;
              k2_h = lambda.k()[0] + 2;
            }
            else        // right boundary wavelet
            {
              k1_h = scale - 3;
              k2_h = scale;
            }
          }
          else          // left boundary wavelet
          {
            k1_h = 0;
            k2_h = 3;
          }
        }
        else                // Generator
        {
          k1_h = lambda.k()[0] - 1;
          k2_h = lambda.k()[0] + 1;
        }


        if ( ((scale*p[0]) < k1_h) || ((scale*p[0]) > k2_h) )
        {
          continue;
        }



        // check if the support of Psi_y contains y-value of p
        // (basis.bases()[1])->support(lambda.j(), lambda.e()[1], lambda.k()[1], k1, k2);

        // compute the support of Psi_lambda (HARD CODED, ONLY FOR LINEAR PRIMBS BASIS WITH HOMOGENEOUS BC's!!!)
        if (lambda.e()[1])   // Wavelet
        {
          if (lambda.k()[1] > 0)
          {
            if (lambda.k()[1] < (scale-1))  // inner wavelet
            {
              k1_h = lambda.k()[1] - 1;
              k2_h = lambda.k()[1] + 2;
            }
            else        // right boundary wavelet
            {
              k1_h = scale - 3;
              k2_h = scale;
            }
          }
          else          // left boundary wavelet
          {
            k1_h = 0;
            k2_h = 3;
          }
        }
        else                // Generator
        {
          k1_h = lambda.k()[1] - 1;
          k2_h = lambda.k()[1] + 1;
        }


        if ( ((scale*p[1]) < k1_h) || ((scale*p[1]) > k2_h) )
        {
          continue;
        }


        // sum the share, USE optimized version 'evaluate_primbs_22_bc11_v3'
        r += ( (*it) * evaluate_primbs_22_bc11_v3(lambda.j(), lambda.e()[0], lambda.k()[0], p[0]) *
                       evaluate_primbs_22_bc11_v3(lambda.j(), lambda.e()[1], lambda.k()[1], p[1]) );


    }

    return r;
  }



}

