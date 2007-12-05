// implementation for sampled_mapping.h

namespace MathTL
{
  template <class C>
  SampledMapping<1,C>::SampledMapping()
    : Grid<1>(), values_()
  {
  }

  template <class C>
  SampledMapping<1,C>::SampledMapping(const SampledMapping<1,C>& sm)
    : Grid<1>(sm), values_(sm.values_)
  {
  }

  template <class C>
  SampledMapping<1,C>::SampledMapping(const Grid<1>& grid)
    : Grid<1>(grid), values_(grid.size())
  {
    for (unsigned int i = 0; i < grid.size(); i++)
      values_[i] = 0;
  }

  template <class C>
  SampledMapping<1,C>::SampledMapping(const Grid<1>& grid, const Array1D<C>& values)
    : Grid<1>(grid), values_(values)
  {
  }

  template <class C>
  SampledMapping<1,C>::SampledMapping(const Grid<1>& grid, const Function<1,C>& f)
    : Grid<1>(grid), values_(grid.size())
  {
    for (unsigned int n(0); n < grid.size(); n++)
      values_[n] = f.value(Point<1,C>(grid_[n]));
  }
  
  template <class C>
  SampledMapping<1,C>::SampledMapping(const MultiIndex<int,1>& a,
				      const MultiIndex<int,1>& b,
				      const InfiniteVector<C, MultiIndex<int,1> >& values,
				      const int resolution)
    : Grid<1>(a[0], b[0], (1<<resolution)*(b[0]-a[0]))
  {
    values_.resize(Grid<1>::size());
    for (int k(a[0]<<resolution), n(0); k <= (b[0]<<resolution); k++, n++)
      values_[n] = values.get_coefficient(MultiIndex<int,1>(k));
  }

  template <class C>
  SampledMapping<1,C>::SampledMapping(const int a,
				      const int b,
				      const InfiniteVector<C, int>& values,
				      const int resolution)
    : Grid<1>(a, b, (1<<resolution)*(b-a))
  {
    values_.resize(Grid<1>::size());
    for (int k(a<<resolution), n(0); k <= (b<<resolution); k++, n++)
      values_[n] = values.get_coefficient(k);
  }
  
  template <class C>
  SampledMapping<1,C>::SampledMapping(const Point<1,C>& a,
				      const Point<1,C>& b,
				      const FixedArray1D<Array1D<C>,1>& values)
    : Grid<1>(a[0], b[0], values[0].size()-1), values_(values[0])
  {
  }
  
  template <class C>
  SampledMapping<1,C>::SampledMapping(const Chart<1>& ch,
				      const FixedArray1D<Array1D<C>,1>& values,
				      const unsigned int resolution)
  {
    const unsigned int n_points = (1 << resolution)+1;
    const C h = 1.0 / (n_points-1);
    
    Point<1,C> x;
    Point<1,C> x_patch;
    
    grid_.resize(n_points);

    values_.resize(n_points);

    for (unsigned int i = 0; i < n_points; i++) {
      x[0] = h*i;
      ch.map_point(x,x_patch);  
      grid_[i] = x_patch[0];
      values_[i] = values[0][i] / ch.Gram_factor(x);
    }
      
  }

  template <class C>
  SampledMapping<1,C>::SampledMapping(const Chart<1>& ch,
				      const unsigned int resolution)
    
  {
    const unsigned int n_points = (1 << resolution)+1;
    const C h = 1.0 / (n_points-1);
    
    Point<1,C> x;
    Point<1,C> x_patch;
    
    grid_.resize(n_points);
    values_.resize(n_points);

    for (unsigned int i = 0; i < n_points; i++) {
      x[0] = h*i;
      ch.map_point(x,x_patch);  
      grid_[i] = x_patch[0];
      values_[i] = 0;
    } 
  }

  template <class C>
  SampledMapping<1,C>& 
  SampledMapping<1,C>::operator = (const SampledMapping<1,C>& sm)
  {
    Grid<1>::operator = (sm);
    values_ = sm.values_;
    return *this;
  }
  
  template <class C>
  void
  SampledMapping<1,C>::add(const SampledMapping<1,C>& s)
  {
    assert(values_.size() == s.values_.size());
    for (unsigned int i = 0; i < values_.size(); i++)
      values_[i] += s.values_[i];
  }

  template <class C>
  void
  SampledMapping<1,C>::add(const C alpha, const SampledMapping<1,C>& s)
  {
    assert(values_.size() == s.values_.size());
    for (unsigned int i = 0; i < values_.size(); i++)
      values_[i] += alpha * s.values_[i];
  }

  template <class C>
  void
  SampledMapping<1,C>::mult(const C alpha)
  {
    for (unsigned int i = 0; i < values_.size(); i++)
      values_[i] *= alpha;
  }

  template <class C>
  void
  SampledMapping<1,C>::matlab_output(std::ostream& os) const
  {
    Grid<1>::matlab_output(os);
    os << "y = " // here we can take y
       << values_
       << ";"
       << std::endl;
  }

  template <class C>
  void
  SampledMapping<1,C>::gnuplot_output(std::ostream& os) const
  {
    unsigned int i;
    const unsigned int size = grid_.size();
    assert(values_.size() == size); // sizes of arrays must match

    for (i = 0; i < size; i++) // loop through entries
      os << grid_[i] << "\t" << values_[i] << std::endl; // format: one value pair per row, in each row x y
  }

  template <class C>
  SampledMapping<2,C>::SampledMapping()
    : Grid<2>(), values_()
  {
  }

  template <class C>
  SampledMapping<2,C>::SampledMapping(const SampledMapping<2,C>& sm)
    : Grid<2>(sm), values_(sm.values_)
  {
  }

  template <class C>
  SampledMapping<2,C>::SampledMapping(const Grid<2>& grid)
    : Grid<2>(grid)
  {
    values_.resize(gridx_.row_dimension(), gridx_.column_dimension());
  }

  template <class C>
  SampledMapping<2,C>::SampledMapping(const Grid<2>& grid, const Matrix<C>& values)
    : Grid<2>(grid), values_(values)
  {
  }

  template <class C>
  SampledMapping<2,C>::SampledMapping(const Grid<2>& grid, const Function<2,C>& f)
    : Grid<2>(grid), values_()
  {
    values_.resize(gridx_.row_dimension(), gridx_.column_dimension());
    for (unsigned int m(0); m < values_.row_dimension(); m++)
      for (unsigned int n(0); n < values_.column_dimension(); n++)
	values_(m,n) = f.value(Point<2,C>(gridx_(m,n), gridy_(m,n)));
  }
  
  template <class C>
  SampledMapping<2,C>::SampledMapping(const MultiIndex<int,2>& a,
				      const MultiIndex<int,2>& b,
				      const InfiniteVector<C, MultiIndex<int,2> >& values,
				      const int resolution)
  {
    // TODO: implement this
  }

  template <class C>
  SampledMapping<2,C>::SampledMapping(const Point<2,C>& a,
				      const Point<2,C>& b,
				      const FixedArray1D<Array1D<C>,2>& values)
    : Grid<2>(a, b, values[0].size()-1, values[1].size()-1),
      values_(values[1].size(), values[0].size())
  {
    for (unsigned int m(0); m < values_.row_dimension(); m++)
      for (unsigned int n(0); n < values_.column_dimension(); n++)
	values_(m,n) = values[1][m] * values[0][n];
  }

  template <class C>
  SampledMapping<2,C>::SampledMapping(const Chart<2>& ch,
				      const FixedArray1D<Array1D<C>,2>& values,
				      const unsigned int resolution)
  {
    
    const unsigned int  n_points = (1 << resolution)+1;
    const C h = 1.0 / (n_points-1);
    
    Point<2,C> x;
    Point<2,C> x_patch;
    
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

  template <class C>
  SampledMapping<2,C>::SampledMapping(const Chart<2>& ch,
				      const unsigned int resolution)

  {
    
    const unsigned int  n_points = (1 << resolution)+1;
    const C h = 1.0 / (n_points-1);
    
    Point<2,C> x;
    Point<2,C> x_patch;
    
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

  template <class C>
  SampledMapping<2,C>&
  SampledMapping<2,C>::operator = (const SampledMapping<2,C>& sm)
  {
    Grid<2>::operator = (sm);
    values_ = sm.values_;
    return *this;
  }

  template <class C>
  void
  SampledMapping<2,C>::add(const SampledMapping<2,C>& s)
  {
    assert(values_.row_dimension() == s.values_.row_dimension()
	   && values_.column_dimension() == s.values_.column_dimension());
    for (unsigned int m(0); m < values_.row_dimension(); m++)
      for (unsigned int n(0); n < values_.column_dimension(); n++)
	values_(m,n) += s.values_(m,n); // Matrix does not (yet) have an add() method
  }

  template <class C>
  void
  SampledMapping<2,C>::add(const C alpha, const SampledMapping<2,C>& s)
  {
    assert(values_.row_dimension() == s.values_.row_dimension()
	   && values_.column_dimension() == s.values_.column_dimension());
    for (unsigned int m(0); m < values_.row_dimension(); m++)
      for (unsigned int n(0); n < values_.column_dimension(); n++)
	values_(m,n) += alpha * s.values_(m,n); // Matrix does not (yet) have an add() method
  }

  template <class C>
  void
  SampledMapping<2,C>::mult(const C alpha)
  {
    for (unsigned int m(0); m < values_.row_dimension(); m++)
      for (unsigned int n(0); n < values_.column_dimension(); n++)
	values_(m,n) *= alpha; // Matrix does not (yet) have an add() method    
  }


  template <class C>
  void
  SampledMapping<2,C>::matlab_output(std::ostream& os) const
  {
    Grid<2>::matlab_output(os);
    os << "z = ";
    print_matrix(values_, os);
    os << ";"
       << std::endl;
  }

  template <class C>
  void
  SampledMapping<2,C>::octave_output(std::ostream& os) const
  {
    Grid<2>::octave_output(os);
    os << "z = ";
    print_matrix(values_, os);
//     os << ";"
//        << std::endl;
  }
  
  template <unsigned int DIM, class C>
  void matlab_output(std::ostream& os,
		     const Array1D<SampledMapping<DIM,C> >& values)
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

  template <unsigned int DIM, class C>
  void octave_output(std::ostream& os,
		     const Array1D<SampledMapping<DIM,C> >& values)
  {
    assert(DIM==2);

    for (unsigned int i = 0; i < values.size(); i++) {
      values[i].octave_output(os);
      os << "hold on" << std::endl
	 << "mesh(x,y,z)" << std::endl;
      if (i == (values.size()-1))
	os << "hold off" << std::endl;
    }
  }
  
}
