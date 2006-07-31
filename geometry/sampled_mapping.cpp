// implementation for sampled_mapping.h

namespace MathTL
{
  SampledMapping<1>::SampledMapping()
    : Grid<1>(), values_()
  {
  }

  SampledMapping<1>::SampledMapping(const SampledMapping<1>& sm)
    : Grid<1>(sm), values_(sm.values_)
  {
  }

  SampledMapping<1>::SampledMapping(const Grid<1>& grid)
    : Grid<1>(grid), values_(grid.size())
  {
    for (unsigned int i = 0; i < grid.size(); i++)
      values_[i] = 0;
  }

  SampledMapping<1>::SampledMapping(const Grid<1>& grid, const Array1D<double>& values)
    : Grid<1>(grid), values_(values)
  {
  }

  SampledMapping<1>::SampledMapping(const Grid<1>& grid, const Function<1>& f)
    : Grid<1>(grid), values_(grid.size())
  {
    for (unsigned int n(0); n < grid.size(); n++)
      values_[n] = f.value(Point<1>(grid_[n]));
  }
  
  SampledMapping<1>::SampledMapping(const MultiIndex<int,1>& a,
				    const MultiIndex<int,1>& b,
				    const InfiniteVector<double, MultiIndex<int,1> >& values,
				    const int resolution)
    : Grid<1>(a[0], b[0], (1<<resolution)*(b[0]-a[0]))
  {
    values_.resize(Grid<1>::size());
    for (int k(a[0]<<resolution), n(0); k <= (b[0]<<resolution); k++, n++)
      values_[n] = values.get_coefficient(MultiIndex<int,1>(k));
  }

  SampledMapping<1>::SampledMapping(const Point<1>& a,
				    const Point<1>& b,
				    const FixedArray1D<Array1D<double>,1>& values)
    : Grid<1>(a[0], b[0], values[0].size()-1), values_(values[0])
  {
  }
  
  SampledMapping<1>::SampledMapping(const Chart<1>& ch,
				    const FixedArray1D<Array1D<double>,1>& values,
				    const unsigned int resolution)
  {
    const unsigned int n_points = (1 << resolution)+1;
    const double h = 1.0 / (n_points-1);
    
    Point<1> x;
    Point<1> x_patch;
    
    grid_.resize(n_points);

    values_.resize(n_points);

    for (unsigned int i = 0; i < n_points; i++) {
      x[0] = h*i;
      ch.map_point(x,x_patch);  
      grid_[i] = x_patch[0];
      values_[i] = values[0][i] / ch.Gram_factor(x);
    }
      
  }

  SampledMapping<1>::SampledMapping(const Chart<1>& ch,
				    const unsigned int resolution)

  {
    const unsigned int n_points = (1 << resolution)+1;
    const double h = 1.0 / (n_points-1);
    
    Point<1> x;
    Point<1> x_patch;
    
    grid_.resize(n_points);
    values_.resize(n_points);

    for (unsigned int i = 0; i < n_points; i++) {
      x[0] = h*i;
      ch.map_point(x,x_patch);  
      grid_[i] = x_patch[0];
      values_[i] = 0;
    } 
  }

  SampledMapping<1>& 
  SampledMapping<1>::operator = (const SampledMapping<1>& sm)
  {
    Grid<1>::operator = (sm);
    values_ = sm.values_;
    return *this;
  }
  
  void
  SampledMapping<1>::add(const SampledMapping<1>& s)
  {
    assert(values_.size() == s.values_.size());
    for (unsigned int i = 0; i < values_.size(); i++)
      values_[i] += s.values_[i];
  }

  void
  SampledMapping<1>::add(const double alpha, const SampledMapping<1>& s)
  {
    assert(values_.size() == s.values_.size());
    for (unsigned int i = 0; i < values_.size(); i++)
      values_[i] += alpha * s.values_[i];
  }

  void
  SampledMapping<1>::mult(const double alpha)
  {
    for (unsigned int i = 0; i < values_.size(); i++)
      values_[i] *= alpha;
  }


  void
  SampledMapping<1>::matlab_output(std::ostream& os) const
  {
    Grid<1>::matlab_output(os);
    os << "y = " // here we can take y
       << values_
       << ";"
       << std::endl;
  }

  SampledMapping<2>::SampledMapping()
    : Grid<2>(), values_()
  {
  }

  SampledMapping<2>::SampledMapping(const SampledMapping<2>& sm)
    : Grid<2>(sm), values_(sm.values_)
  {
  }

  SampledMapping<2>::SampledMapping(const Grid<2>& grid)
    : Grid<2>(grid)
  {
    values_.resize(gridx_.row_dimension(), gridx_.column_dimension());
  }

  SampledMapping<2>::SampledMapping(const Grid<2>& grid, const Matrix<double>& values)
    : Grid<2>(grid), values_(values)
  {
  }

  SampledMapping<2>::SampledMapping(const Grid<2>& grid, const Function<2>& f)
    : Grid<2>(grid), values_()
  {
    values_.resize(gridx_.row_dimension(), gridx_.column_dimension());
    for (unsigned int m(0); m < values_.row_dimension(); m++)
      for (unsigned int n(0); n < values_.column_dimension(); n++)
	values_(m,n) = f.value(Point<2>(gridx_(m,n), gridy_(m,n)));
  }
  
  SampledMapping<2>::SampledMapping(const MultiIndex<int,2>& a,
				    const MultiIndex<int,2>& b,
				    const InfiniteVector<double, MultiIndex<int,2> >& values,
				    const int resolution)
  {
    // TODO: implement this
  }

  SampledMapping<2>::SampledMapping(const Point<2>& a,
				    const Point<2>& b,
				    const FixedArray1D<Array1D<double>,2>& values)
    : Grid<2>(a, b, values[0].size()-1, values[1].size()-1),
      values_(values[1].size(), values[0].size())
  {
    for (unsigned int m(0); m < values_.row_dimension(); m++)
      for (unsigned int n(0); n < values_.column_dimension(); n++)
	values_(m,n) = values[1][m] * values[0][n];
  }

  SampledMapping<2>::SampledMapping(const Chart<2>& ch,
				    const FixedArray1D<Array1D<double>,2>& values,
				    const unsigned int resolution)
  {

    const unsigned int  n_points = (1 << resolution)+1;
    const double h = 1.0 / (n_points-1);
    
    Point<2> x;
    Point<2> x_patch;
    
    gridx_.resize(n_points,n_points);
    gridy_.resize(n_points,n_points);

    values_.resize(n_points,n_points);

    // setup grid
    for (unsigned int i = 0; i < n_points; i++) {
      x[0] = h*i;
     for (unsigned int j = 0; j < n_points; j++) {
       x[1] = h*j;
       ch.map_point(x,x_patch);
       gridx_.set_entry(i,j,x_patch[0]);
       gridy_.set_entry(i,j,x_patch[1]);
       values_.set_entry(i,j,(values[0][i] * values[1][j]) / ch.Gram_factor(x));
     }
    }

//      values_.resize(n_points,n_points);
//      // setup values
//      for (unsigned int i = 0; i < n_points; i++) {
//        x[0] = h*i;
//        for (unsigned int j = 0; j < n_points; j++) {
// 	 x[1] = h*j;
// 	 values_.set_entry(i,j,(values[0][i] * values[1][j]) / ch.Gram_factor(x));
//        }
//      }
  }

  SampledMapping<2>::SampledMapping(const Chart<2>& ch,
				    const unsigned int resolution)

  {
    
    const unsigned int  n_points = (1 << resolution)+1;
    const double h = 1.0 / (n_points-1);
    
    Point<2> x;
    Point<2> x_patch;
    
    gridx_.resize(n_points,n_points);
    gridy_.resize(n_points,n_points);
    
    // setup grid
    for (unsigned int i = 0; i < n_points; i++) {
      x[0] = h*i;
      for (unsigned int j = 0; j < n_points; j++) {
	x[1] = h*j;
	ch.map_point(x,x_patch);
	//cout << "i=" << x_patch[0] << " j= " << x_patch[1] << endl;
	gridx_.set_entry(i,j,x_patch[0]);
	gridy_.set_entry(i,j,x_patch[1]);
      }
    }

    values_.resize(gridx_.row_dimension(), gridx_.column_dimension());    
    
  }


  SampledMapping<2>&
  SampledMapping<2>::operator = (const SampledMapping<2>& sm)
  {
    Grid<2>::operator = (sm);
    values_ = sm.values_;
    return *this;
  }

  void
  SampledMapping<2>::add(const SampledMapping<2>& s)
  {
    assert(values_.row_dimension() == s.values_.row_dimension()
	   && values_.column_dimension() == s.values_.column_dimension());
    for (unsigned int m(0); m < values_.row_dimension(); m++)
      for (unsigned int n(0); n < values_.column_dimension(); n++)
	values_(m,n) += s.values_(m,n); // Matrix does not (yet) have an add() method
  }

  void
  SampledMapping<2>::add(const double alpha, const SampledMapping<2>& s)
  {
    assert(values_.row_dimension() == s.values_.row_dimension()
	   && values_.column_dimension() == s.values_.column_dimension());
    for (unsigned int m(0); m < values_.row_dimension(); m++)
      for (unsigned int n(0); n < values_.column_dimension(); n++)
	values_(m,n) += alpha * s.values_(m,n); // Matrix does not (yet) have an add() method
  }

  void
  SampledMapping<2>::mult(const double alpha)
  {
    for (unsigned int m(0); m < values_.row_dimension(); m++)
      for (unsigned int n(0); n < values_.column_dimension(); n++)
	values_(m,n) *= alpha; // Matrix does not (yet) have an add() method    
  }


  void
  SampledMapping<2>::matlab_output(std::ostream& os) const
  {
    Grid<2>::matlab_output(os);
    os << "z = ";
    print_matrix(values_, os);
    os << ";"
       << std::endl;
  }

  template <unsigned int DIM>
  void matlab_output(std::ostream& os,
		     const Array1D<SampledMapping<DIM> >& values)
  {
    switch (DIM) {
    case 1: {
      for (unsigned int i = 0; i < values.size(); i++) {
	values[i].matlab_output(os);
	os << "hold on" << std::endl
	   << "plot(x,y)" << std::endl;
	if (i == (values.size()-1))
	  os << "hold off" << std::endl;
      }
      break;
    }
    case 2: {
      for (unsigned int i = 0; i < values.size(); i++) {
	values[i].matlab_output(os);
	os << "hold on" << std::endl
	   << "surf(x,y,z)" << std::endl;
	if (i == (values.size()-1))
	  os << "hold off" << std::endl;
      }
      break;
    }      
    }// end switch

  }

  
}
