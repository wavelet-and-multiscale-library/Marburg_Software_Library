// implementation for a_basis.h

#include <interval/a_basis.h>
#include <cmath>
#include <algebra/polynomial.h>
#include <algebra/piecewise.h>
#include <numerics/gram_schmidt.h>
#include <numerics/ortho_poly.h>

#define ALPERT_CODE // use the port of the code I got from Bradley Alpert

namespace WaveletTL
{
  /* constructors **********************************************************/
  template <int n>
  ABasis<n>::ABasis()
  {
    setup();
  }

  /* setup routine shared by different constructors */
  #ifndef ALPERT_CODE
  /* construction as described in [A] */
  template <int n>
  void
  ABasis<n>::setup()
  {
    int i, j, k;
    Polynomial<double> p;
    Array1D<Piecewise<double> > monom(2*n-1);
    Array1D<double> a(n);

    /* construct generators: normed Legendre polynomials */
    generators = Array1D<Piecewise<double> >(n);
    // start with monomial basis
    // create monomial basis on [-1,1] in monom, create generators on the fly
    for (i = 0; i < 2*n-1; i++) {
      p.set_coefficient(i,1);
      monom[i].set_local_expansion(-1,p);
      monom[i].set_local_expansion(0,p);
      if (i < n)
        generators[i].set_local_expansion(0,p);
      p.set_coefficient(i,0);
    }
    // orthogonalize generators with Gram-Schmidt procedure
    gramSchmidtProcess(generators);

    /* construct wavelets */
    // we first construct them on [-1,1] and then scale it to [0,1]
    // initialize wavelets
    wavelets = Array1D<Piecewise<double> >(n);
    for (i = 0; i < n; i++) {
      wavelets[i].set_local_expansion(-1,monom[i].get_local_expansion(0));
      wavelets[i].set_local_expansion(0,(-1.0)*monom[i].get_local_expansion(0));
    }

    // step 1: Orthogonalize wavelets to the generators via Gram-Schmidt
    for (i = 0; i < n; i++) // loop through all wavelets
      for (j = 0; j < n; j++) // loop through the monoms
        wavelets[i] -= monom[j].inner_product(wavelets[i]) * monom[j];
    cout << "Step 1 finished." << endl;
    cout << wavelets << endl;

    // step 2: Orthogonalize them to x^n,..., x^(2n-2)
    for (i = 0; i < n-1; i++) {
      for (j = i, k = -1; j < n; j++) {
        if ( (a[j] = wavelets[j].inner_product(monom[n+i])) != 0 )
          k = j; // memorise which inner product is != 0
      }
      if (k != -1) { // wavelets are not yet orthogonal
        // swap k-th and i-th entry
        wavelets.swap(i,k);
        a.swap(i,k);
        // orthogonalize
        for (j = i+1; j < n; j++)
          wavelets[j] -= a[j]/a[i] * wavelets[i];
      }
    }
    cout << "Step 2 finished." << endl;
    cout << wavelets << endl;

    /* step 3: Orthogonalize wavelets among themselves, preserving proper
       orthogonality to x^n,...,x^(2n-2), and norm them */
    for (i = n-1; i >= 0; i--) {
      for (j = i; j < n; j++)
        wavelets[i] -= wavelets[i].inner_product(wavelets[j]) * wavelets[j];
      wavelets[i] *= 1.0/sqrt(wavelets[i].inner_product(wavelets[i]));
    }
    cout << "Step 3 finished." << endl;
    cout << wavelets << endl;

    #if 0     
    // scale wavelets to [0,1]
    for (i = 0; i < n; i++) {
      wavelets[i].dilate_me(1);
      wavelets[i].shift_me(1);
    }
    #endif // 0
  }

  #else // ALPERT_CODE
  /* port of the Maple code I got from Bradley Alpert */
  #define ip(g, k) (g.inner_product(legendre[k],0.0,1.0)*(double)(2*(k)+1))
  template <int n>
  void
  ABasis<n>::setup()
  {
    int i, j;
    MathTL::LegendrePolynomial legendre_pol;
    Array1D<Polynomial<double> > legendre(2*n-1);
    Array1D<Polynomial<double> > f(n);
    double a;

    // we first construct them on [-1,1] and then scale them to [0,1]

    /* construct Legendre polynomials (needed later) */
    for (i = 0, a = 1.0; i < 2*n-1; i++) {
      legendre[i] = legendre_pol.assemble(i);
      // renorm with (2i)!/(2^i (i!)^2)
      if (i > 1) {
        a *= ((double)(2*i-1))/(double)i;
        legendre[i] *= a;
      }
    }
    #if 0
    // start with monomial basis
    for (i = 0; i < 2*n-1; i++)
      legendre[i].set_coefficient(i,1);
    // create Legendre Polynomials by Gram-Schmidt process
    // orthogonalization with respect to [-1,1]
    for (i = 0; i < 2*n-1; i++) {
      for (j = 1; j < i; j++) {
        legendre[i] -= legendre[i].inner_product(legendre[j],-1.0,1.0) * legendre[j];
      }
      // norm 
      legendre[i] *= 1.0/sqrt(legendre[i].inner_product(legendre[i],-1.0,1.0));
    }
    cout << "Legendre polynomials constructed." << endl;
    #endif // 0
    
    /* construct generators: normed Legendre polynomials */
    generators = Array1D<Piecewise<double> >(n);
    // copy first n Legendre Polynomials into enerators
    for (i = 0; i < n; i++) {
      Polynomial<double> p;
      p.set_coefficient(1,2);
      p.set_coefficient(0,-1);
      generators[i].set_local_expansion(0,legendre[i].substitute_into(p)); // L_i(2x-1), i.e. scale to [0,1]
      generators[i] *= 1.0/sqrt(generators[i].inner_product(generators[i])); // norm
    }

//    cout << "Constructing functions f ..." << endl;
    /* construct wavelets */
    // step 1: create wavelets and orthogonalize then to the generators via Gram-Schmidt
    for (i = 0; i < n; i++) {
      f[i].set_coefficient(n-i-1,1);
      for (j = (n-i)%2; j <= n-1; j+=2) {
        f[i] -= ip(f[i],j) * legendre[j];
      }
    }
//    cout << "Step 1 finished." << endl;
//    cout << f << endl;

    // step 2: Orthogonalize them to x^n,..., x^(2n-2)
    for (i = 0; i < n-1; i++) {
      a = ip(f[i],i+n);
      for (j = i+2; j < n; j+=2)
        f[j] -= (ip(f[j],i+n) / a) * f[i];
    }
//    cout << "Step 2 finished." << endl;
//    cout << f << endl;

    /* step 3: Orthogonalize wavelets among themselves, preserving proper
       orthogonality to x^n,...,x^(2n-2), and norm them */
    for (i = n-1; i >= 0; i--) {
      for (j = i+2; j < n; j+=2)
        f[i] -= f[i].inner_product(f[j],0.0,1.0)*2 * f[j];
      f[i] *= 1.0/sqrt(f[i].inner_product(f[i],0.0,1.0)*2);
    }
//    cout << "Step 3 finished." << endl;
//    cout << f << endl;

    wavelets = Array1D<Piecewise<double> >(n);
    #if 0 // for tests
    for (i = 0; i < n; i++)
      wavelets[i].set_local_expansion(0, ((i+n)%2 ? (-1.0) : 1.0) * f[i] );
    #else
    /* set local expansions and continue wavelets to [-1,0] */
    for (i = 0; i < n; i++) {
      wavelets[i].set_local_expansion(0, ((i+n)%2 ? (-1.0) : 1.0) * f[i] );
      f[i].scale(-1.0); // f -> f(-x)
      wavelets[i].set_local_expansion(-1, f[i]);
      // scale to [0,1]
      wavelets[i].dilate_me(1);
      wavelets[i].shift_me(1);
    }
    #endif // 0
//    cout << "Construction finished." << endl;
  }
  #endif // ALPERT_CODE

  /* member functions ******************************************************/
  template <int n>
  inline
  typename ABasis<n>::Index
  ABasis<n>::first_generator(const int j) const
  {
    assert(j >= j0());
    return IntervalMultiIndex<ABasis<n> >(j, E_GENERATOR, 0, 0, this);
  }
  
  template <int n>
  inline
  typename ABasis<n>::Index
  ABasis<n>::last_generator(const int j) const
  {
    assert(j >= j0());
    return IntervalMultiIndex<ABasis<n> >(j, E_GENERATOR, (1<<j)-1, n-1, this);
  }
  
  template <int n>
  inline
  typename ABasis<n>::Index
  ABasis<n>::first_wavelet(const int j) const
  {
    assert(j >= j0());
    return IntervalMultiIndex<ABasis<n> >(j, E_WAVELET, 0, 0, this);
  }
  
  template <int n>
  inline
  typename ABasis<n>::Index
  ABasis<n>::last_wavelet(const int j) const
  {
    assert(j >= j0());
    return IntervalMultiIndex<ABasis<n> >(j, E_WAVELET, (1<<j)-1, n-1, this);
  }

  template <int n>
  Piecewise<double>
  ABasis<n>::get_function(const ABasis<n>::Index index) const
  {
    assert(index.c() < n); // component is within range

    Piecewise<double> fkt;

    if (index.e() == 0) // generator
      fkt = generators[index.c()];
    else { // wavelet
      fkt = wavelets[index.c()];
      fkt.dilate_me(index.j());
    }
    fkt.shift_me(index.k());

    return fkt;
  }
  
  #if 0
  template <int n>
  void
  ABasis<n>::reconstruct(const InfiniteVector<double, Index>& c,
                       const int j,
                       InfiniteVector<double, Index>& v) const
  {
    ;
  }
  
  template <int n>
  void
  ABasis<n>::reconstruct_1(const Index& lambda,
                         const int j,
                         InfiniteVector<double, Index>& c) const
  {
    ;
  }
  #endif
  
}
