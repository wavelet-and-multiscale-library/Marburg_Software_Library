// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_TINY_TOOLS_H
#define _MATHTL_TINY_TOOLS_H

#include <cassert>
#include <cmath>
#include <functional>
#include <map>

// some utility functions

//! absolute value of a real number
template <class R>
inline R abs(const R a)
{
  return  (a > 0 ? a : -a);
}

//! signum function
template <class R>
inline R sign(const R a)
{
  return (a > 0 ? 1. : (a < 0 ? -1. : 0.));
}

//! stable computation of hypotenuse
/*!
  stable computation of the hypotenuse length for scalars a and b;
  we do _not_ compute
    sqrt(a*a+b*b)
  but the more stable
    a*sqrt(1+b/a*b/a)
 */
template <class R>
inline R hypot(const R a, const R b)
{
  R r(0);
  if (a == 0)
    r = abs(b);
  else
    {
      R c(b/a);
      r = abs(a) * sqrt(1+c*c);
    }
  return r;
}

//! factorial of a (signed or unsigned) integer
template <class I>
I factorial(const I n)
{
  I r(1);
  for (I i(2); i <= n; r *= i, i++);
  return r;
}

/*!
  binomial coefficient n over k
  By convention, we have

   (0) = 1
   (0)

   (n) = 0 for n < k or k < 0
   (k)
 */

inline int binomial(const int n, const int k)
{
  int r(1);
  
  if (k > n || k < 0)
    return 0;
  
  if (k == 0)
    return 1;
  
  for (int i(k+1); i <= n; i++)
    r *= i;
  
  for (int i(2); i+k <= n; i++)
    r /= i; // always possible without remainder
  
  return r;
}

//! (-1)^k
inline int minus1power(const int k)
{
  return (k%2==0 ? 1 : -1);
}

//! integer power n^k, 0^0=1, k >= 0
template <class I, class J>
int intpower(const I n, const J k)
{
  int r(1);
  for (J j(1); j <= k; j++)
    r *= n;

  return r;
}

inline
double mypow(double base, int power){
    double r=1.0;
    while(power--){
        r*=base;
    }
    return r;
}


inline
long long fast_power(long long base, long long power) {
    #define MOD 1000000007
    long long result = 1;
    while(power > 0) {

        if(power % 2 == 1) { // Can also use (power & 1) to make code even faster
            result = (result*base) % MOD;
        }
        base = (base * base) % MOD;
        power = power / 2; // Can also use power >>= 1; to make code even faster
    }
    return result;
}

inline
long long fast_power2(long base, long power) {
    #define MOD 1000000007
    long long result = 1;
    while(power > 0) {

        if(power & 1) { // Can also use (power & 1) to make code even faster
            result = (result*base) % MOD;
        }
        base = (base * base) % MOD;
        power>>=1; // Can also use power >>= 1; to make code even faster
    }
    return result;
}

inline double fastPow(double a, int b) {
  union {
    double d;
    int x[2];
  } u = { a };
  u.x[1] = (int)(b * (u.x[1] - 1072632447) + 1072632447);
  u.x[0] = 0;
  return u.d;
}

inline double fastPrecisePow(double a, int b) {
  // calculate approximation with fraction of the exponent
  int e = (int) b;
  union {
    double d;
    int x[2];
  } u = { a };
  u.x[1] = (int)((b - e) * (u.x[1] - 1072632447) + 1072632447);
  u.x[0] = 0;

  // exponentiation by squaring with the exponent's integer part
  // double r = u.d makes everything much slower, not sure why
  double r = 1.0;
  while (e) {
    if (e & 1) {
      r *= a;
    }
    a *= a;
    e >>= 1;
  }

  return r * u.d;
}

/*!
  helper function object for thresholding within a std::map<I,C>:
  returns true if argument is strictly less than eta in modulus
*/
template <class I, class C>
class threshold_criterion
  : public std::unary_function<I, C>
{
public:
  threshold_criterion(const double eta) : eta_(eta) {}
  bool operator() (std::pair<const I, C>& p) { return fabs(p.second) < eta_; }

private:
  const double eta_;
};

/*!
  fast routine to compute 2^(j/2) without a sqrt() call
  (we assume that j>=0)
*/
inline
double twotothejhalf(const int j)
{
  return j%2 ? M_SQRT2 * (1<<(j>>1)) : 1<<(j>>1);
}

/*!
  dyadic modulo x |-> x mod 2^j >= 0
  (we assume that j>=0)
*/
inline
int dyadic_modulo (const int x, const int j)
{
  return (x >= 0
	  ? x-((x>>j)<<j)
	  : (x+(((-x)>>j)<<j)) ? x+((((-x)>>j)+1)<<j) : 0);
}


/*
 * computes floor(log2(n)), where n is a positive integer
 * Tested for int n, unsigned int n.
 * 2.5 times faster than log2(double)
 * 
 * Be careful with negative arguments and automatic conversion to unsigned int.
 * 
 * There might be more clever alternatives for this function!
 * 
 * Assumption: 
 * sizeof(int) = sizeof(unsigned int) = 4
 * 
 * Code borrowed from:
 * Bit Twiddling Hacks - http://graphics.stanford.edu/~seander/bithacks.html#IntegerLogLookup
 * By Sean Eron Anderson
 * seander@cs.stanford.edu 
 * (with additions from Andrew Shapira) 
 */
const unsigned int log2Bits[] = {0x2, 0xC, 0xF0, 0xFF00, 0xFFFF0000}; //, 0xFFFFFFFF00000000}; // uncomment for 64 bit integers
const unsigned int log2Sizes[] = {1, 2, 4, 8, 16}; //, 32}; // uncomment for 64 bit integers

inline unsigned int log2(const unsigned int n)
{
    unsigned int v;  // 32-bit value to find the log2 of 
    v = n;
    register unsigned int r = 0; // result of log2(v) will go here
    int i;
    for (i = 4; i >= 0; i--) // unroll for speed...
    {
        if (v & log2Bits[i])
        {
            v >>= log2Sizes[i];
            r |= log2Sizes[i];
        } 
    }
    return r;
}



/* For sorting arrays, e.g., multiindices, by their values from position 
 * 
 * i to DIM-1
 * 
 * T = const Array1D < MultiIndex<int,DIM> >  == j0
 * Assumption: f < DIM, otherwise nothing will be compared!
 *
template<class T> 
struct index_cmp
{
    index_cmp(const T arr, unsigned int f) : arr(arr), from(f) {}
    bool operator()(const size_t a, const size_t b) const
    {
        //return arr[a] > arr[b];
        // compare the MultiIndex stored in arr from f to the end lexicographically
        for (unsigned int i=from; i < arr[a].size(); ++i)
        {
            if (arr[a][i] > arr[b][i])
            {
                return false;
            }
            else if (arr[a][i] < arr[b][i])
            {
                return true;
            }
        }
        return (a<b);
    }
    const T arr;
    const unsigned int from;
};


 * Similar, but the last entry of arr is ignored. This is relevant for the first level j_ with a certain norm \|j_\|. 
 * All but the last entry of such a level are equal to j0()[patch][i]. 
 * However, the last entry is of some value independent of j0()[patch][DIM-1]
 
template<class T> 
struct index_cmp_ignoreLastEntry
{
    index_cmp_ignoreLastEntry(const T arr, unsigned int f) : arr(arr), from(f) {}
    bool operator()(const size_t a, const size_t b) const
    {
        //return arr[a] > arr[b];
        // compare the MultiIndex stored in arr from f to the end lexicographically
        for (unsigned int i=from; i < (arr[a].size()-1); ++i)
        {
            if (arr[a][i] > arr[b][i])
            {
                return false;
            }
            else if (arr[a][i] < arr[b][i])
            {
                return true;
            }
        }
        return (a<b);
    }
    const T arr;
    const unsigned int from;
};
 For sorting arrays, e.g., multiindices, by their values from position 
 * 
 * DIM-1 to 0 (reversed order in the dimensions)
 * 
 * (ordering w.r.t. entry in the array, i.e., a and b, is not reversed)
 * 
 * T = const Array1D < MultiIndex<int,DIM>>  == j0
 *
template<class T> 
struct index_cmp_reversed
{
    index_cmp_reversed(const T arr) : arr(arr) {}
    bool operator()(const size_t a, const size_t b) const
    {
        //return arr[a] > arr[b];
        // compare the MultiIndex stored in arr from f to the end lexicographically
        for (int i=arr[a].size()-1; i >= 0; i--)
        {
            if (arr[a][i] > arr[b][i])
            {
                return false;
            }
            else if (arr[a][i] < arr[b][i])
            {
                return true;
            }
        }
        return (a<b);
    }
    const T arr;
};
*/
#endif
