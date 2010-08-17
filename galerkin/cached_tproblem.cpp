// implementation for cached_tproblem.h

#include <cmath>
#include <algebra/vector.h>
#include <numerics/eigenvalues.h>
#include <galerkin/galerkin_utils.h>

namespace WaveletTL
{
  template <class PROBLEM>
  CachedTProblem<PROBLEM>::CachedTProblem(const PROBLEM* P,
					  const double estnormA,
					  const double estnormAinv)
    : problem(P), normA(estnormA), normAinv(estnormAinv)
  {
  }

  template <class PROBLEM>
  double
  CachedTProblem<PROBLEM>::a(const Index& lambda,
			     const Index& nu) const
  {
    double r = 0;

    if (problem->local_operator()) {
      const int lambda_num = lambda.number();
      const int nu_num = nu.number();

      // Be careful, there is no generator level in the cache! The situation from the MRA setting:
      // KEY OF GENERATOR LEVEL IS j0-1 NOT j0 !!!!
      // does not hold in the tensor setting. Generators and wavelets on
      // the minimal level are thrown together in one index set (componentwise),
      // this also applies to tensors of generators on the lowest level

      index_lt lambda_key(lambda.j()),first_level(basis().j0());
      for (int k=0;k<space_dimension;k++)
      {
          lambda_key[k] = lambda_key[k]-first_level[k];
      }
//TODO (PERFORMANCE): store the numbers of all levels up to jmax, do not compute anything here:
      int blocknumber(lambda_key.number());

      // check wether entry has already been computed
      typedef std::list<Index> IntersectingList;

      // search for column 'mu'
      typename ColumnCache::iterator col_lb(entries_cache.lower_bound(nu_num));
      typename ColumnCache::iterator col_it(col_lb);

      if (col_lb == entries_cache.end() ||
	  entries_cache.key_comp()(nu_num, col_lb->first))
	{
	  // insert a new column
	  typedef typename ColumnCache::value_type value_type;
	  col_it = entries_cache.insert(col_lb, value_type(nu_num, Column()));
	}

      Column& col(col_it->second);

      // check wether the level 'lambda' belongs to has already been calculated
      typename Column::iterator lb(col.lower_bound(blocknumber));
      typename Column::iterator it(lb);
      
      if (lb == col.end() ||
	  col.key_comp()(blocknumber, lb->first))
      {
          // no entries have ever been computed for this column and this level
	  // compute whole level block
	  // insert a new level

          typedef typename Column::value_type value_type;
          it = col.insert(lb, value_type(blocknumber, Block()));
          Block& block(it->second);

          if ((blocknumber == 0) && (lambda.j() == first_level))
          {
              // Insert Generators and Wavelets

              // also add Generators to the Cache
              // this case has to be considered seperatly because
              // intersecting_wavelets divides the case of a basis
              // function made entirely of generators and the cache does not.

              // TODO : produce nicer looking code:
              IntersectingList nusG,nusW;
              intersecting_wavelets(basis(), nu,
				    lambda.j(),
                                    true, // this argument isn't that helpful for the way the method is used in a() .  It doesn't play a role for DIM>1 and miserably fails for DIM=1. Need to work on the routine in tbasis_support!
				    nusG);
              intersecting_wavelets(basis(), nu,
				    lambda.j(),
                                    false, // this argument isn't that helpful for the way the method is used in a() .  It doesn't play a role for DIM>1 and miserably fails for DIM=1. Need to work on the routine in tbasis_support!
				    nusW);
              // compute entries
              for (typename IntersectingList::const_iterator it(nusG.begin()), itend(nusG.end());it != itend; ++it)
              {
                  const double entry = problem->a(*it, nu);
                  typedef typename Block::value_type value_type_block;
                  if (fabs(entry) > 1e-16 ) //(entry != 0.)
                  {
                      block.insert(block.end(), value_type_block((*it).number(), entry));
                      if ((*it).number() == lambda_num)
                      {
                          r = entry;
                      }
                  }
              }

              for (typename IntersectingList::const_iterator it(nusW.begin()), itend(nusW.end());it != itend; ++it)
              {
                  const double entry = problem->a(*it, nu);
                  typedef typename Block::value_type value_type_block;
                  if (fabs(entry) > 1e-16 ) //(entry != 0.)
                  {
                      // Insertion should be efficient, since the wavelets coma after the generators
                      block.insert(block.end(), value_type_block((*it).number(), entry));
                      if ((*it).number() == lambda_num)
                      {
                          r = entry;
                      }
                  }
              }
          }
          else
          {
              // there are no Generators

              IntersectingList nus;
              intersecting_wavelets(basis(),
                                    nu,
                                    //std::max(j, basis().j0()),
                                    lambda.j(),
                                    false, // this argument isn't that helpful for the way the method is used in a() .  It doesn't play a role for DIM>1 and miserably fails for DIM=1. Need to work on the routine in tbasis_support!
              //                    == (basis().j0()-1),

                                    nus);
              // compute entries
              for (typename IntersectingList::const_iterator it(nus.begin()), itend(nus.end());it != itend; ++it)
              {
                  const double entry = problem->a(*it, nu);
                  typedef typename Block::value_type value_type_block;
                  if (fabs(entry) > 1e-16 ) //(entry != 0.)
                  {
                      block.insert(block.end(), value_type_block((*it).number(), entry));
                      if ((*it).number() == lambda_num)
                      {
                          r = entry;
                      }
                  }
              }
          }
      }
      // level already exists --> extract row corresponding to 'lambda'
      else {
	Block& block(it->second);

 	//typename Block::iterator block_lb(block.lower_bound(lambda));
	typename Block::iterator block_lb(block.lower_bound(lambda_num));
 	typename Block::iterator block_it(block_lb);
	// level exists, but in row 'lambda' no entry is available ==> entry must be zero
	if (block_lb == block.end() ||
	    block.key_comp()(lambda_num, block_lb->first))
	  {
	    r = 0;
	  }
	else {
	  r = block_it->second;
	}
      }
    } // end of local operator
    else {
      // for nonlocal operators, we put full level blocks into the cache, regardless of support intersections

      const int lambda_num = lambda.number();
      const int nu_num = nu.number();

      // comment about the key applies here too.
      //typedef typename Index::type_type generator_type;
      //int j = (lambda.e() == generator_type()) ? (lambda.j()-1) : lambda.j();
      
      index_lt lambda_key(lambda.j()), first_level(basis().j0());
      for (int k=0;k<space_dimension;k++)
      {
          lambda_key[k] = lambda_key[k]-first_level[k];
      }
//TODO PERFORMANCE: store the numbrs of all levels up to jmax, do not compute anything here:
      int blocknumber(lambda_key.number());

      // search for column 'mu'
      typename ColumnCache::iterator col_lb(entries_cache.lower_bound(nu_num));
      typename ColumnCache::iterator col_it(col_lb);

      if (col_lb == entries_cache.end() ||
	  entries_cache.key_comp()(nu_num, col_lb->first))
	{
	  // insert a new column
	  typedef typename ColumnCache::value_type value_type;
	  col_it = entries_cache.insert(col_lb, value_type(nu_num, Column()));
	}

      Column& col(col_it->second);

      // check wether the level block which 'lambda' belongs to has already been calculated
      typename Column::iterator lb(col.lower_bound(blocknumber));
      typename Column::iterator it(lb);

      if (lb == col.end() ||
	  col.key_comp()(blocknumber, lb->first))
      {
	  // no entries have ever been computed for this column and this level
	  // compute whole level block

	  // insert a new level
	  typedef typename Column::value_type value_type;
	  it = col.insert(lb, value_type(blocknumber, Block()));

	  Block& block(it->second);

	  // collect all indices in the level block
	  typedef std::list<Index> IndexList;
	  IndexList nus;

	  if (lambda.j() == first_level)
          {
              // generators & wavelets on level j0
              for (Index lambda_it(basis().first_generator(first_level)), lambda_end(basis().last_wavelet(first_level));lambda_it != lambda_end; ++lambda_it)
              {
                  nus.push_back(lambda_it);
              }
          } else
          {
              // wavelets on level j > j0
              for (Index lambda_it(basis().first_wavelet(lambda.j())), lambda_end(basis().last_wavelet(lambda.j()));lambda_it != lambda_end; ++lambda_it)
              {
                  nus.push_back(lambda_it);
              }
          }
	  // compute entries
	  for (typename IndexList::const_iterator it(nus.begin()), itend(nus.end()); it != itend; ++it)
          {
              const double entry = problem->a(*it, nu);
              typedef typename Block::value_type value_type_block;
              if (fabs(entry) > 1e-16 ) {
                  block.insert(block.end(), value_type_block((*it).number(), entry));
                  if ((*it).number() == lambda_num)
                  {
                      r = entry;
                  }
              }
          }
      }
      // level already exists --> extract row corresponding to 'lambda'
      else
      {
          Block& block(it->second);

          typename Block::iterator block_lb(block.lower_bound(lambda_num));
          typename Block::iterator block_it(block_lb);
          if (block_lb == block.end() || block.key_comp()(lambda_num, block_lb->first))
          {
              // level exists, but in row 'lambda' no entry is available ==> entry must be zero
              r = 0;
          }
          else
          {
              r = block_it->second;
          }
      }
    }
    return r;
  }

  template <class PROBLEM>
  void
  CachedTProblem<PROBLEM>::add_ball(const Index& lambda,
				    //InfiniteVector<double, Index>& w,
				    Vector<double>& w,
				    const int radius,
				    const double factor,
				    const int maxlevel,
                                    const CompressionStrategy strategy) const
  {
      /*
       * For dim=1,2 it is easy to describe all levels in the (1-norm) ball of radius around lambda.
       * For higher dimensions this is done recursivly.
       */
      double d1 = D(lambda);
      if (space_dimension == 1)
      {
          int j0 = this->basis().j0()[0];
          Index mu;
          for (int level = max (j0, lambda.j()[0]-radius); level  < min(lambda.j()[0]+radius,maxlevel)+1;level++)
          {
              mu = this->basis().first_wavelet(level);
              // the result of the following call is that the cache holds the whole level block corresponing to all wavelets with level = |mu|
              // ideally one should replace this mehtod call with a slightly modified version of the code from a(,)
              // using the first generator on level j0 would meand no difference, since block 0 contains wavelets and generators
              a(mu,lambda);
              // add the level
              Block& block (entries_cache[lambda.number()][level-j0]);
#if 1
              for (typename Block::const_iterator it(block.begin()), itend(block.end()); it != itend; ++it)
              {
                  // high caching strategy:
                  // The call of D fills the diagonal block of each column in block
                  w[it->first]=w[it->first]+factor*(it->second)/d1/D(this->basis().get_wavelet(it->first));
                  // low caching strategy:
                  // by calling the uncached Version of D no clumn but the one of lambda is altered. This Version is slower than the high caching strategy.
                  //w[it->first]=w[it->first]+factor*(it->second)/d1/(problem->D(this->basis().get_wavelet(it->first)));
              }
#else
              // mid caching strategy:
              // maybe it is too much effort to compute the whole diagonal block (as in the high caching strategy).
              // For the columns in block other than the one of lambda we use D_2 to only compute and cache a(mu,mu) (a single entry)
              // Do do this a bool variable has to be introduced marking that only one entry in mus column is stored.
              // this caching strategy is useful for sparse signals with active coefficients on high levels
              // ... needs to be tested though ...
              typename Block::const_iterator it(block.begin()), itend(block.end());
              w[it->first]=w[it->first]+factor*(it->second)/d1/D(this->basis().get_wavelet(it->first));
              ++it;
              for (; it != itend; ++it)
              {
                  w[it->first]=w[it->first]+factor*(it->second)/d1/D_already(this->basis().get_wavelet(it->first));
              }
#endif
          }
      }
      else if (space_dimension == 2)
      {

          // The ball can be described of levellines consisting of levels with the same multidegree
          // The first is determined with the distance of lambda.j and j0. The last with maxlevel
          // The first level in a levelline is determined with minx = min(j0[0], lambda.j[0]-radius)
          // The last level in a levelline is determined with miny = min(j0[1], lambda.j[1]-radius)

          Index mu;
          index_lt j0(this->basis().j0());
          int lambdaline = lambda.j()[0]+lambda.j()[1];
          int lowestline = j0[0]+j0[1];
          int dist2j0=lambdaline-lowestline;
          int dist2maxlevel=maxlevel-lambdaline;
          MultiIndex<int,space_dimension> currentlevel;
          int xstart,xend,ystart,blocknumber;
          // iterate the levellines. offset relative to lambdas levelline
          for (int offset = -std::min(dist2j0,radius); offset < std::min(dist2maxlevel,radius)+1; offset++)
          {
              // iterate over the levels on the levelline
              xstart = lambda.j()[0]-radius+ceil((radius+offset)/2.0); //x coordinate of the first level, ignoring restrictions by j0[0]
              xend = lambda.j()[0]+floor((radius+offset)/2.0); // same for the last
              ystart = lambda.j()[1]+floor((radius+offset)/2.0); // and for the second dimension
              // the first level in the ball is denoted with steps=0.
              // The first level on the current levelline in the ball may have steps >0.
              for (int steps = max(0,j0[0]-xstart); steps <= min(xend-xstart,ystart-j0[1]) ; steps++ )
              {
                  currentlevel[0]= xstart+steps;
                  currentlevel[1]= ystart-steps;
                  mu = this->basis().first_wavelet(currentlevel);

                  a(mu,lambda); // computes & stores the block in the cache
                  blocknumber = currentlevel[0]-j0[0] + ((lambdaline-lowestline+offset)*(lambdaline-lowestline+offset+1))/2;

                  // add the level
                  Block& block (entries_cache[lambda.number()][blocknumber]);
#if 1
                  for (typename Block::const_iterator it(block.begin()), itend(block.end()); it != itend; ++it)
                  {
                      // high caching strategy
                      w[it->first]=w[it->first]+factor*(it->second)/d1/D(this->basis().get_wavelet(it->first));
                      // low caching strategy
                      //w[it->first]=w[it->first]+factor*(it->second)/d1/(problem->D(this->basis().get_wavelet(it->first)));
                  }
#else
                  // mid caching strategy
                  typename Block::const_iterator it(block.begin()), itend(block.end());
                  w[it->first]=w[it->first]+factor*(it->second)/d1/D(this->basis().get_wavelet(it->first));
                  ++it;
                  for (; it != itend; ++it)
                  {
                      w[it->first]=w[it->first]+factor*(it->second)/d1/D_already(this->basis().get_wavelet(it->first));
                  }
#endif
              }
          }
      }
      else // dim > 2. iteration over all levels in 'range' is done recursivly for arbitrary dimensions
      {
          add_level_recurse(lambda,w,radius,factor,lambda.j(),0,maxlevel,true,d1);
      }
  }

  

  template <class PROBLEM>
  void
  CachedTProblem<PROBLEM>::add_level_recurse(const Index& lambda,
                                            Vector<double>& w,
                                            const int radius,
                                            const double factor,
                                            const index_lt & current_level,
                                            const int current_dim,
                                            const int maxlevel,
                                            const bool legal,
                                            const double precond) const
  {

      if (current_dim < space_dimension)
      {
          /* we call a (multiindex) level legal, if its norm is <=maxlevel.
           * let lmin and lmax be the lowest and largest (integer) values
           * that result in a level with norm <= jmax if written at [curdimension]
           * in the current_level. then we have j0[curdim]<= lmin[curdim] (left boundary)
           * For the right boundary we observe that curlevel does not need to be a legal
           * level in order for the disc with radius around it to contain legal levels.
           * In such a case let r2 be the minimal radius such that the r2-ball around
           * curlevel contains a legal level. Denote the level resulting from curlevel, where
           * the entry [curdim] is increased by k by lvl_k. We observe that the minimal radius
           * around lvl_k that gives a ball that contains legal levels has a radius of r2+2*k.
           * r2 can be computed with lmax.
           */
          int current_j0 (this->basis().j0()[current_dim]);
          index_lt next_level(current_level);
          int leftrange = - min (next_level[current_dim]-current_j0,radius);
          // compiler complains about the line:
          // int rightrange = std::min(radius, (maxlevel - multi_degree(current_level)));

          int rightrange = (maxlevel - multi_degree(current_level));
          rightrange = std::min (rightrange, radius); // rightrange describes the largest entry (in the current dimension) that results in legal level (<=jmax)

          next_level[current_dim]+= leftrange;

          /* the 2 ball around (5,3) with jmax = 8 leads to a call with current_level = (6,3) and radius 1.
           * The only legal choice for the 2nd dimension would result in (6,2), but with j0 in this dimension being 3
           * => leftrange = 0, rightrange = -1.
           * In 3D: 2 ball around (5,3,4) (j0 = (3,3,3), jmax = 12) leads to 1 Ball around (6,3,4). This leads to the same situation,
           * still this is no error, since in the next dimension the 1 ball around (6,3,4) contains (6,3,3) (legal).
           */
/*
          if (multi_degree(next_level)>maxlevel)
          {
              cout << "lambda " << lambda.number() << " " << lambda << "cur_lvl " << current_level << " next_lvl " << next_level << " lr " << leftrange << " rr " << rightrange << endl;
          }
 * */
          //assert (multi_degree(next_level)<=maxlevel);
          //assert (leftrange<=rightrange);

          for (int k=leftrange; k <=rightrange; k++)
          {
              add_level_recurse(lambda,w,radius-abs(k),factor,next_level,current_dim+1,maxlevel,true,precond);
              ++next_level[current_dim];
              // lambda_it[rec_level]=lambda[rec_level]+k; // fuer den naechsten iteranten versteht sich
          }
          // discs araound sublevels outside jmax can still contain levels below jmax! 
          // If one entry grows another has to get smaller, thus we never need to increase the current index by more than radius/2
          int rightrange2;

          //if (rightrange > 0 )
          {
              rightrange2 = (radius-rightrange)/2;
          }
          if (current_dim < space_dimension -1)
          {
              // if there would be no minimal level the following line would be satisfactory:
              //for (int k(rightrange+1); k<=std::max(rightrange,0)+(radius - abs(rightrange))/2;k++)
              for (int k(std::max(rightrange+1,leftrange)); k<=std::max(rightrange,0)+(radius - abs(rightrange))/2;k++)
              {
                  add_level_recurse(lambda,w,radius-abs(k),factor,next_level,current_dim+1,maxlevel,false,precond);
                  ++next_level[current_dim];
              }
          }
      }
      else
      { // we have iterated over all dimensions and can now add the current level (if it is legal)!
          assert (legal == true);
          Index mu;
          if (current_level == this->basis().j0())
          {
              mu = this->basis().first_generator();
          }
          else
          {
              mu = this->basis().first_wavelet(current_level);
          }
          // the result of the following call is that the cache holds the whole level block corresponing to all wavelets with level = |mu|
          // ideally one should replace this mehtod call with a slightly modified version of the code from a(,)
          a(mu,lambda);


          index_lt level_key(this->basis().j0());
          for (int k=0;k<space_dimension;k++)
          {
              level_key[k] = current_level[k]-level_key[k];
          }
          //PERFORMANCE: store the numbers of all levels up to jmax, do not compute anything here:
          int blocknumber(level_key.number());

          // return the block corresponding to level = |mu|

          // same as:
          // w.add(factor,entries_cache[lambda.number()][blocknumber]);
          Block& block (entries_cache[lambda.number()][blocknumber]);
#if 1
          for (typename Block::const_iterator it(block.begin()), itend(block.end()); it != itend; ++it)
          {
              // high caching strategy
              w[it->first]=w[it->first]+factor*(it->second)/precond/D(this->basis().get_wavelet(it->first));
              // low caching strategy
              //w[it->first]=w[it->first]+factor*(it->second)/precond/(problem->D(this->basis().get_wavelet(it->first)));
          }
#else
          // mid caching strategy
          typename Block::const_iterator it(block.begin()), itend(block.end());
          w[it->first]=w[it->first]+factor*(it->second)/precond/D_diag(this->basis().get_wavelet(it->first));
          it++;
          for (; it != itend; ++it)
          {
              w[it->first]=w[it->first]+factor*(it->second)/precond/D_already(this->basis().get_wavelet(it->first));
          }
#endif
      }

  }
/*
  template <class PROBLEM>
  void
  CachedTProblem<PROBLEM>::add_ball(const Index& lambda,
				    //InfiniteVector<double, Index>& w,
				    Vector<double>& w,
				    const int radius,
				    const double factor,
				    const int maxlevel,
                                    const set<unsigned int> levelwindow) const
  {

      // For dim=1,2 it is easy to describe all levels in the (1-norm) ball of radius around lambda.
      // For higher dimensions this is done recursivly.

      double d1 = D(lambda);
      if (space_dimension == 1)
      {
          int j0 = this->basis().j0()[0];
          Index mu;
          set<unsigned int>:: iterator lwit(levelwindow.begin()),lwitend(levelwindow.end());
          for (int level = (max (j0, lambda.j()[0]-radius) > *lwit)?max (j0, lambda.j()[0]-radius):*lwit; level  < ((min(lambda.j()[0]+radius,maxlevel) < *lwitend)?min(lambda.j()[0]+radius,maxlevel):*lwitend) +1; level++)
          {
              if (levelwindow.find(level) == lwitend)
              {
                  continue;
              }
              mu = this->basis().first_wavelet(level);
              // the result of the following call is that the cache holds the whole level block corresponing to all wavelets with level = |mu|
              // ideally one should replace this mehtod call with a slightly modified version of the code from a(,)
              // using the first generator on level j0 would meand no difference, since block 0 contains wavelets and generators
              a(mu,lambda);
              // add the level
              Block& block (entries_cache[lambda.number()][level-j0]);
# if 1
              for (typename Block::const_iterator it(block.begin()), itend(block.end()); it != itend; ++it)
              {
                  // high caching strategy
                  w[it->first]=w[it->first]+factor*(it->second)/d1/D(this->basis().get_wavelet(it->first));
                  // low caching strategy
                  //w[it->first]=w[it->first]+factor*(it->second)/d1/(problem->D(this->basis().get_wavelet(it->first)));
              }
#else
              // mid caching strategy
              typename Block::const_iterator it(block.begin()), itend(block.end());
              w[it->first]=w[it->first]+factor*(it->second)/d1/D(this->basis().get_wavelet(it->first));
              ++it;
              for (; it != itend; ++it)
              {
                  w[it->first]=w[it->first]+factor*(it->second)/d1/D_already(this->basis().get_wavelet(it->first));
              }
#endif
          }
      }
      else if (space_dimension == 2)
      {

          // The ball can be described of levellines consisting of levels with the same multidegree
          // The first is determined with the distance of lambda.j and j0. The last with maxlevel
          // The first level in a levelline is determined with minx = min(j0[0], lambda.j[0]-radius)
          // The last level in a levelline is determined with miny = min(j0[1], lambda.j[1]-radius)
          Index mu;
          index_lt j0(this->basis().j0());
          int lambdaline = lambda.j()[0]+lambda.j()[1];
          int lowestline = j0[0]+j0[1];
          int dist2j0=lambdaline-lowestline;
          int dist2maxlevel=maxlevel-lambdaline;
          MultiIndex<int,space_dimension> currentlevel;
          int xstart,xend,ystart,blocknumber;
          // iterate the levellines. offset relative to lambdas levelline
          for (int offset = -std::min(dist2j0,radius); offset < std::min(dist2maxlevel,radius)+1; offset++)
          {
              // iterate over the levels on the levelline
              xstart = lambda.j()[0]-radius+ceil((radius+offset)/2.0); //x coordinate of the first level, ignoring restrictions by j0[0]
              xend = lambda.j()[0]+floor((radius+offset)/2.0); // same for the last
              ystart = lambda.j()[1]+floor((radius+offset)/2.0); // and for the second dimension
              // the first level in the ball is denoted with steps=0.
              // The first level on the current levelline in the ball may have steps >0.
              for (int steps = max(0,j0[0]-xstart); steps <= min(xend-xstart,ystart-j0[1]) ; steps++ )
              {
                  currentlevel[0]= xstart+steps;
                  currentlevel[1]= ystart-steps;
                  mu = this->basis().first_wavelet(currentlevel);

                  a(mu,lambda); // computes & stores the block in the cache
                  blocknumber = currentlevel[0]-j0[0] + ((lambdaline-lowestline+offset)*(lambdaline-lowestline+offset+1))/2;
                  // add the level if it is in levelwindow
                  if (levelwindow.find(blocknumber) == levelwindow.end())
                  {
                      continue;
                  }
                  Block& block (entries_cache[lambda.number()][blocknumber]);
#if 1
                  for (typename Block::const_iterator it(block.begin()), itend(block.end()); it != itend; ++it)
                  {
                      // high caching strategy
                      w[it->first]=w[it->first]+factor*(it->second)/d1/D(this->basis().get_wavelet(it->first));
                      // low caching strategy
                      //w[it->first]=w[it->first]+factor*(it->second)/d1/(problem->D(this->basis().get_wavelet(it->first)));
                  }
#else
                  // mid caching strategy
                  typename Block::const_iterator it(block.begin()), itend(block.end());
                  w[it->first]=w[it->first]+factor*(it->second)/d1/D(this->basis().get_wavelet(it->first));
                  ++it;
                  for (; it != itend; ++it)
                  {
                      w[it->first]=w[it->first]+factor*(it->second)/d1/D_already(this->basis().get_wavelet(it->first));
                  }
#endif
              }
          }
      }
      else // dim > 2. iteration over all levels in 'range' is done recursivly for arbitrary dimensions
      {
          add_level_recurse(lambda,w,radius,factor,lambda.j(),0,maxlevel,true,d1,levelwindow);
      }

  }

  template <class PROBLEM>
  void
  CachedTProblem<PROBLEM>::add_level_recurse(const Index& lambda,
                                            Vector<double>& w,
                                            const int radius,
                                            const double factor,
                                            const index_lt & current_level,
                                            const int current_dim,
                                            const int maxlevel,
                                            const bool legal,
                                            const double precond,
                                            const set<unsigned int> levelwindow) const
  {

      if (current_dim < space_dimension)
      {
          // we call a (multiindex) level legal, if its norm is <=maxlevel.
          // let lmin and lmax be the lowest and largest (integer) values
          // that result in a level with norm <= jmax if written at [curdimension]
          // in the current_level. then we have j0[curdim]<= lmin[curdim] (left boundary)
          // For the right boundary we observe that curlevel does not need to be a legal
          // level in order for the disc with radius around it to contain legal levels.
          // In such a case let r2 be the minimal radius such that the r2-ball around
          // curlevel contains a legal level. Denote the level resulting from curlevel, where
          // the entry [curdim] is increased by k by lvl_k. We observe that the minimal radius
          // around lvl_k that gives a ball that contains legal levels has a radius of r2+2*k.
          // r2 can be computed with lmax.
          
          int current_j0 (this->basis().j0()[current_dim]);
          index_lt next_level(current_level);
          int leftrange = - min (next_level[current_dim]-current_j0,radius);
          // compiler complains about the line:
          // int rightrange = std::min(radius, (maxlevel - multi_degree(current_level)));

          int rightrange = (maxlevel - multi_degree(current_level));
          rightrange = std::min (rightrange, radius); // rightrange describes the largest entry (in the current dimension) that results in legal level (<=jmax)

          next_level[current_dim]+= leftrange;

          // the 2 ball around (5,3) with jmax = 8 leads to a call with current_level = (6,3) and radius 1.
          // The only legal choice for the 2nd dimension would result in (6,2), but with j0 in this dimension being 3
          // => leftrange = 0, rightrange = -1.
          // In 3D: 2 ball around (5,3,4) (j0 = (3,3,3), jmax = 12) leads to 1 Ball around (6,3,4). This leads to the same situation,
          // still this is no error, since in the next dimension the 1 ball around (6,3,4) contains (6,3,3) (legal).
          //

          //assert (multi_degree(next_level)<=maxlevel);
          //assert (leftrange<=rightrange);

          for (int k=leftrange; k <=rightrange; k++)
          {
              add_level_recurse(lambda,w,radius-abs(k),factor,next_level,current_dim+1,maxlevel,true,precond);
              ++next_level[current_dim];
              // lambda_it[rec_level]=lambda[rec_level]+k; // fuer den naechsten iteranten versteht sich
          }
          // discs araound sublevels outside jmax can still contain levels below jmax!
          // If one entry grows another has to get smaller, thus we never need to increase the current index by more than radius/2
          int rightrange2;

          //if (rightrange > 0 )
          {
              rightrange2 = (radius-rightrange)/2;
          }
          if (current_dim < space_dimension -1)
          {
              // if there would be no minimal level the following line would be satisfactory:
              //for (int k(rightrange+1); k<=std::max(rightrange,0)+(radius - abs(rightrange))/2;k++)
              for (int k(std::max(rightrange+1,leftrange)); k<=std::max(rightrange,0)+(radius - abs(rightrange))/2;k++)
              {
                  add_level_recurse(lambda,w,radius-abs(k),factor,next_level,current_dim+1,maxlevel,false,precond);
                  ++next_level[current_dim];
              }
          }
      }
      else
      {   // we have iterated over all dimensions and can now add the current level if ...
          // ... it is legal
          // ... and in levelwindow
          assert (legal == true);
          Index mu;

          index_lt level_key(this->basis().j0());
          for (int k=0;k<space_dimension;k++)
          {
              level_key[k] = current_level[k]-level_key[k];
          }
          //PERFORMANCE: store the numbers of all levels up to jmax, do not compute anything here:
          int blocknumber(level_key.number());

          if (levelwindow.find(blocknumber) != levelwindow.end())
          {
              if (current_level == this->basis().j0())
              {
                  mu = this->basis().first_generator();
              }
              else
              {
                  mu = this->basis().first_wavelet(current_level);
              }
              // the result of the following call is that the cache holds the whole level block corresponing to all wavelets with level = |mu|
              // ideally one should replace this mehtod call with a slightly modified version of the code from a(,)
              a(mu,lambda);

              // return the block corresponding to level = |mu|

              // same as:
              // w.add(factor,entries_cache[lambda.number()][blocknumber]);
              Block& block (entries_cache[lambda.number()][blocknumber]);
              typename Block::const_iterator it(block.begin()), itend(block.end());
              w[it->first]=w[it->first]+factor*(it->second)/precond/D(this->basis().get_wavelet(it->first));
              ++it;
              for (; it != itend; ++it)
              {
                  w[it->first]=w[it->first]+factor*(it->second)/precond/D_already(this->basis().get_wavelet(it->first));
              }
          }
      }
  }
*/ // comment for the levelwindow code
  

  template <class PROBLEM>
  double
  CachedTProblem<PROBLEM>::norm_A() const
  {
      if (normA == 0.0)
      {
          //const int offset = 2;
          int offset;
          switch (space_dimension) {
              case 1:
                  offset = 2;
                  break;
              case 2:
                  offset = 1;
                  break;
              default:
                  offset = 0;
          }
        
#if _WAVELETTL_CACHEDPROBLEM_VERBOSITY >= 1
          cout << "CachedTProblem()::norm_A() called..." << endl;
#endif

          set<Index> Lambda;
          //cout << "cached_tproblem.norm_A :: last wavelet = " << (problem->basis().last_wavelet(multi_degree(problem->basis().j0())+offset)) << endl;
          for (Index lambda ( problem->basis().first_generator() ), itend(problem->basis().last_wavelet(multi_degree(problem->basis().j0())+offset));; ++lambda) {
            Lambda.insert(lambda);
            if (lambda == itend) break;
          }
          SparseMatrix<double> A_Lambda;
          setup_stiffness_matrix(*this, Lambda, A_Lambda);
//#if 1
          double help;
          unsigned int iterations;
          LanczosIteration(A_Lambda, 1e-6, help, normA, 200, iterations);
          normAinv = 1./help;
//#else
//        Vector<double> xk(Lambda.size(), false);
//        xk = 1;
//        unsigned int iterations;
//        normA = PowerIteration(A_Lambda, xk, 1e-6, 100, iterations);
//#endif

#if _WAVELETTL_CACHEDPROBLEM_VERBOSITY >= 1
          cout << "... done!" << endl;
#endif
      }
      return normA;
  }

  template <class PROBLEM>
  double
  CachedTProblem<PROBLEM>::norm_Ainv() const
  {
      //cout << " this is norm_Ainv(). Verbosity is " << _WAVELETTL_CACHEDPROBLEM_VERBOSITY << endl;

      if (normAinv == 0.0)
      {
          //const int offset = 2;
          int offset;
          switch (space_dimension) {
              case 1:
                  offset = 2;
                  break;
              case 2:
                  offset = 1;
                  break;
              default:
                  offset = 0;
          }
#if _WAVELETTL_CACHEDPROBLEM_VERBOSITY >= 1
          cout << "CachedTProblem()::norm_Ainv() called..." << endl;
          cout << "verbosity is "<< _WAVELETTL_CACHEDPROBLEM_VERBOSITY << endl;
#endif
#if _WAVELETTL_CACHEDPROBLEM_VERBOSITY >=2
          Index tempindex(problem->basis().last_wavelet(multi_degree(problem->basis().j0())+offset));
          cout << "number of last wavelet is "<<tempindex.number()<<endl;
#endif
          set<Index> Lambda;
          for (Index lambda ( problem->basis().first_generator() ), itend(problem->basis().last_wavelet(multi_degree(problem->basis().j0())+offset));; ++lambda)
          {
              Lambda.insert(lambda);
              if (lambda == itend) break;
          }
          SparseMatrix<double> A_Lambda;
          setup_stiffness_matrix(*this, Lambda, A_Lambda);
#if 1
          double help;
          unsigned int iterations;
          LanczosIteration(A_Lambda, 1e-6, help, normA, 200, iterations);
          normAinv = 1./help;
#else
          Vector<double> xk(Lambda.size(), false);
          xk = 1;
          unsigned int iterations;
          normAinv = InversePowerIteration(A_Lambda, xk, 1e-6, 200, iterations);
#endif

#if _WAVELETTL_CACHEDPROBLEM_VERBOSITY >= 1
          cout << "... done!" << endl;
#endif
      }
      return normAinv;
  }

  template <class PROBLEM>
  void
  CachedTProblem<PROBLEM>::apply(const std::set<int>& window, const Vector<double>& x,
				Vector<double>& res) const
  {
//TODO PERFORMANCE :: reduce get_wavelet calls

    res.resize(x.size());
    typedef typename Index::type_type generator_type;
    // cout << " size = " << entries_cache.size() << endl;
    if (entries_cache.size() == 0)
    {
        unsigned int l = 0;
        for (typename std::set<int>::const_iterator win_it_col = window.begin(); win_it_col != window.end(); win_it_col++, l++)
        {
            const double d1 = this->D(*(problem->basis().get_wavelet(*win_it_col)));
            unsigned int k = 0;
            for (typename std::set<int>::const_iterator win_it_row = window.begin(); win_it_row != window.end(); win_it_row++, k++)
            {
                res[k] += x[l] * this->a(*(problem->basis().get_wavelet(*win_it_row)),*(problem->basis().get_wavelet(*win_it_col)))
                        / (d1*this->D(*(problem->basis().get_wavelet(*win_it_row))));
            }
        }
        return;
    }

    typename std::set<int>::const_iterator win_it_col = window.begin();
    unsigned int l = 0;

    for (typename ColumnCache::iterator col_it(entries_cache.begin()),col_end(entries_cache.end()); col_it != col_end; col_it++)
    {
        if (col_it->first < *win_it_col)
        {
            continue;
        }
        else if (col_it->first > *win_it_col)
        {
            const double d1 = this->D(*(problem->basis().get_wavelet(*win_it_col)));
            unsigned int k = 0;
            for (typename std::set<int>::const_iterator win_it_row(window.begin()), win_it_row_end(window.end()); win_it_row != win_it_row_end; win_it_row++, k++)
            {
                res[k] += x[l] * this->a(*(problem->basis().get_wavelet(*win_it_row)), *(problem->basis().get_wavelet(*win_it_col)))
                        / (d1*this->D(*(problem->basis().get_wavelet(*win_it_row))));
            }
            col_it--;
// TODO PERFORMANCE: col_it zeigt jetzt auf die gerade eingefügte Column, durch das col_it++ am Ende der for Schleife also wieder auf die selbe Column ... ineffizient.
            // 	cout << "column should be " << (*win_it_col) << " and it is " << col_it->first << endl;
            // 	cout << "done inserting column " << (*win_it_col) << endl;
            win_it_col++;
            l++;
            if(l == res.size())
                return;
        }
        else if (col_it->first == *win_it_col)
        {
            typename std::set<int>::const_iterator win_it_row = window.begin();
            unsigned int k = 0;
            index_lt row_key, first_level(basis().j0());
            int rowind_blocknumber;
            for (typename Column::const_iterator block_it = (col_it->second).begin(), block_it_end((col_it->second).end()); block_it != block_it_end; block_it++)
            {
// TODO Performance: replace the following lines with a function call to a nice mapping: Wavelet_number -> Level_number
                const Index* rowind = problem->basis().get_wavelet(*win_it_row);
                row_key = rowind->j(); //,first_level(basis().j0());
                for (int kk=0;kk<space_dimension;kk++)
                {
                    row_key[kk] = row_key[kk]-first_level[kk];
                }
// TODO PERFORMANCE: store the numbers of all levels up to jmax, do not compute anything here:
                rowind_blocknumber = row_key.number();
                if (block_it->first < rowind_blocknumber)
                {
                    continue;
                }
                else if (block_it->first > rowind_blocknumber)
                { // insert this block!
                    //	    cout << "Rowblock number " << rowind_blocknumber << endl;
                    Column& col(col_it->second);
                    typename Column::iterator lb(col.lower_bound(rowind_blocknumber));
                    typename Column::iterator it(lb);
                    typedef typename Column::value_type value_type;
                    it = col.insert(lb, value_type(rowind_blocknumber, Block()));
                    std::list<Index> nus;

                    // WORKS ONLY FOR LOCAL OPERATORS
                    // compute basis functions on the level of the current row that
                    // have intersecting support with the wavelet given by the current column.
                    // Again: Generators and Wavelets are thrown together in one index set.
                    // The structure of "intersect_wavelets" makes things a bit messy (see a() )
                    // so we use intersecting_elements instead
                    // caution: ensure that jmax_ in setup_full_collection is high enough! 
                    // otherwise get_wavelet will fail!

// TODO produce nicer looking code (like in a() )
                    intersecting_elements(basis(),
                                          *(problem->basis().get_wavelet(col_it->first)),
                                          rowind->j(),
                                          nus);                    
                    {
                        Block& block(it->second);
                        const double d1 = this->D(*(problem->basis().get_wavelet(col_it->first)));
                        double entry;
                        typedef typename Block::value_type value_type_block;
                        
                        for (typename std::list<Index>::const_iterator nus_it(nus.begin()), nus_itend(nus.end()); nus_it != nus_itend; ++nus_it)
                        {
                            //	      cout << "nonzero entry at " << (*it).number() << endl;
                            while (*win_it_row < (*nus_it).number() && win_it_row != window.end())
                            {
                                win_it_row++;
                                k++;
                            }
                            entry = problem->a(*nus_it, problem->basis().get_wavelet(col_it->first));
                            
                            if (fabs(entry) > 1e-16 ) //(entry != 0.)
                            {
                                block.insert(block.end(), value_type_block((*nus_it).number(), entry));
                            }
                            // at this point holds either (*win_it_row).number() < (*nus_it).number() or (*win_it_row).number() == (*nus_it).number()

                            if ( *win_it_row == (*nus_it).number() )
                            {
                                if (fabs(entry) > 1e-16 ) //(entry != 0.)
                                {
#if 1
                                    // high caching strategy
                                    res[k] += x[l] * (entry / (d1*D(*nus_it)));
                                    // low caching strategy
                                    //res[k] += x[l] * (entry / (d1*problem->D(*nus_it)));
#else
                                    // mid caching strategy
                                    if (!flag)
                                    {
                                        res[k] += x[l] * (entry / (d1*D_diag(*nus_it)));
                                        flag = true;
                                    }
                                    else
                                    {
                                        res[k] += x[l] * (entry / (d1*D_already(*nus_it)));
                                    }
#endif
                                }
                            }
                        }
                    }
                    if (win_it_row == window.end())
                    {
                        break;
                    }
                    block_it--;

                    // leaving this level
                    // iterate win_it_row till the end of the current sublevel
// TODO PERFORMANCE: as before a mapping IndexNumber->BlockNumber could speed up things
                    // this code is also executed in another branch of the if then else statement. apply any code changes also there
                    {
                        index_lt row_key_temp;
                        int rowind_blocknumber_temp;
                        while (true)
                        {
                            if (win_it_row == window.end())
                            {
                                break;
                            }
                            //compute current blocknumber for win_it_row:
                            row_key_temp = problem->basis().get_wavelet(*win_it_row)->j(); //,first_level(basis().j0());
                            for (int kk=0;kk<space_dimension;kk++)
                            {
                                row_key_temp[kk] = row_key_temp[kk]-first_level[kk];
                            }
                            rowind_blocknumber_temp = row_key_temp.number();
                            if (rowind_blocknumber_temp == rowind_blocknumber)
                            {
                                win_it_row++;
                                k++;
                            }
                            else
                            {
                                break;
                            }
                        }
                    }
                    if (win_it_row == window.end())
                    {
                        break;
                    }
                }
                else if (block_it->first == rowind_blocknumber)
                { // extract existing block!
                    //	    cout << "extracting level " << j << " in column " << col_it->first << endl;
                    const double d1 = this->D(*(problem->basis().get_wavelet(col_it->first)));
                    
                    for (typename Block::const_iterator it(block_it->second.begin()), itend(block_it->second.end()); it != itend; ++it)
                    {
                        //	      cout << "nonzero entry in row " << it->first << endl;
                        while (*win_it_row < it->first && k < res.size()) {
                            win_it_row++;
                            k++;
                        }
                        //	      cout << "k = " << k << endl;
                        if (k == res.size())
                            break;
                        // at this point holds either (*win_it_row).number() < (*it).number() or (*win_it_row).number() == (*it).number()
                        if ( *win_it_row == it->first )
                        {
                            //		cout << "putting row " << *win_it_row << endl;
#if 1
                            // high caching strategy
                            res[k] += x[l] * (it->second / (d1*D(problem->basis().get_wavelet(it->first))));
                            //low caching strategy
                            //res[k] += x[l] * (it->second / (d1*problem->D(problem->basis().get_wavelet(it->first))));
#else
                            // mid caching strategy
                            if (!flag)
                            {
                                res[k] += x[l] * (it->second / (d1*D_diag(problem->basis().get_wavelet(it->first))));
                                flag = true;
                            }
                            else
                            {
                                res[k] += x[l] * (it->second / (d1*D_already(problem->basis().get_wavelet(it->first))));
                            }
#endif
                        }
                    }// end loop over level block

// TODO PERFORMANCE: kick the following line!
                    if (k == res.size())
                        break;

                    //cout << "leaving the level " << endl;

                    // leaving this level
                    // iterate win_it_row till the end of the current sublevel
// TODO PERFORMANCE: as before a mapping IndexNumber->BlockNumber could speed up things
                    // this code is also executed in anther branch of the if then else statement. apply any code changes also there
                    {
                        index_lt row_key_temp;
                        int rowind_blocknumber_temp;
                        while (true)
                        {
                            if (win_it_row == window.end())
                            {
                                break;
                            }
                            //compute current blocknumber for win_it_row:
                            row_key_temp = problem->basis().get_wavelet(*win_it_row)->j(); //,first_level(basis().j0());
                            for (int kk=0;kk<space_dimension;kk++)
                            {
                                row_key_temp[kk] = row_key_temp[kk]-first_level[kk];
                            }
                            rowind_blocknumber_temp = row_key_temp.number();
                            if (rowind_blocknumber_temp == rowind_blocknumber)
                            {
                                win_it_row++;
                                k++;
                            }
                            else
                            {
                                break;
                            }
                        }
                    }

                    if (win_it_row == window.end())
                    {
                        break;
                    }
                }// end else if (block_it->first == j)
            }// end loop over level blocks
            win_it_col++;

            l++;

            if (l == res.size())
                return;
        }// end if existing column found
    }// end loop over all columns in cache

    // insert the remaining columns
    while (l < res.size())
    {
        //      cout << "inserting column " << *win_it_col << endl;
        const double d1 = problem->D(*(problem->basis().get_wavelet(*win_it_col)));
        unsigned int k = 0;
        for (typename std::set<int>::const_iterator win_it_row = window.begin();
	     win_it_row != window.end(); win_it_row++, k++)
        {
            res[k] += x[l] * this->a(*(problem->basis().get_wavelet(*win_it_row)),
				     *(problem->basis().get_wavelet(*win_it_col)))
	            / (d1*problem->D(*(problem->basis().get_wavelet(*win_it_row))));
        }
        win_it_col++;
        l++;
    }
  }

// TODO PERFORMANCE: "Shift this case upwards"
  template <class PROBLEM>
  bool
  CachedTProblem<PROBLEM>::CG(const std::set<int>& window, const Vector<double> &b, Vector<double> &xk,
	  const double tol, const unsigned int maxiter, unsigned int& iterations)
  {
    // see: "Templates for the Solution of Linear Systems: Building Blocks for Iterative Methods"

    Vector<double> rk(window.size(), false),
      zk(window.size(), false),
      pk(window.size(), false),
      Apk(window.size(), false);

    // first (negative) residual
    //cout << "entering windowed matrix vector product..." << endl;
    apply(window, xk, rk);
    //cout << "done windowed matrix vector product..." << endl;
    //cout << rk << endl;
    rk.subtract(b);
    const double normr0 = l2_norm_sqr(rk);
    double normrk = normr0, rhok = 0, oldrhok = 0;
    for (iterations = 1; normrk/normr0 > tol*tol && iterations <= maxiter; iterations++)
      {
	//P.apply_preconditioner(rk, zk);
	zk = rk;
	rhok = rk * zk;

	if (iterations == 1) // TODO: shift this case upwards!
	  pk = zk;
	else
	  pk.sadd(rhok/oldrhok, zk);
	//cout << "entering windowed matrix vector product..." << endl;
	apply(window, pk, Apk);

	//cout << "done windowed matrix vector product..." << endl;
	//cout << Apk << endl;
	const double alpha = rhok/(pk*Apk);
	xk.add(-alpha,  pk);
	rk.add(-alpha, Apk);
	normrk = l2_norm_sqr(rk);

	oldrhok = rhok;
      }

    return (iterations <= maxiter);
  }

  /* Some code for the mid caching strategy:
  template <class PROBLEM>
  double
  CachedTProblem<PROBLEM>::D_already(const Index& lambda) const
  {
#if 1
      index_lt level_key(this->basis().j0()), lambda_level(lambda.j());
      for (int k=0;k<space_dimension;k++)
      {
          level_key[k] = lambda_level[k]-level_key[k];
      }
      //PERFORMANCE: store the numbers of all levels up to jmax, do not compute anything here:
      double is = entries_cache[lambda.number()][level_key.number()][lambda.number()];
      if (isnan(is))
      {
          cout << "error in D_already!" << endl;
          pause;
      }
      double shouldbe = a(lambda,lambda);
      if (is != shouldbe)
      {
          cout << "problem in D_already" << endl;
          pause;
      }
      double returnvalue = sqrt (is);
      int lanu(lambda.number()), lenu(level_key.number());
      return sqrt (entries_cache[lambda.number()][level_key.number()][lambda.number()]);
      //return sqrt(a(lambda,lambda));
#else
      return 1;
#endif
  }
  */
  
}
