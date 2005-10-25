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
      values_(values[0].size(), values[1].size())
  {
    for (unsigned int m(0); m < values_.row_dimension(); m++)
      for (unsigned int n(0); n < values_.column_dimension(); n++)
	values_(m,n) = values[0][m] * values[1][n];
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
  SampledMapping<2>::matlab_output(std::ostream& os) const
  {
    Grid<2>::matlab_output(os);
    os << "z = ";
    print_matrix(values_, os);
    os << ";"
       << std::endl;
  }
  
}
