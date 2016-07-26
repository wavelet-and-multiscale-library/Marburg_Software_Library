/* 
 * File:   Periodic_TestProblem.h
 * Author: keding
 *
 * Created on July 7, 2016, 10:01 AM
 */

#include <galerkin/periodic_gramian.h>
#include <galerkin/periodic_laplacian.h>


#ifndef PERIODIC_TESTPROBLEM_H
#define	PERIODIC_TESTPROBLEM_H

template <unsigned int N>
class PeriodicTestProblem
  : public PeriodicGramianProblem
{
public:
  double g(const double t) const {
    switch(N) {
    case 1:
      //u(t) = 100*(t*t*t-1.5*t*t+0.5*t)
        //return 100*(t*t*t-1.5*t*t+0.5*t);
        return 1;
      //return 100*(-6*t+3);
      break;
    case 2:
      return t*(1-t);
    case 3:
      return cos(2*M_PI*t);      
    case 4:
      return sin(3*M_PI*t);
    case 5:
        if(t<0.5)
           return t;
        else
           return 1-t;
    default:
      return 0;
      break;
    }
  }
};

template <unsigned int N>
class PeriodicTestProblem2
  : public PeriodicLaplacianProblem
{
public:
  double g(const double t) const {
    switch(N) {
    case 1:
      //u(t) = 100*(t*t*t-1.5*t*t+0.5*t)
        //return 100*(t*t*t-1.5*t*t+0.5*t);
        return 2;
      //return 100*(-6*t+3);
      break;
    case 2:
      return 4*M_PI*M_PI*cos(2*M_PI*t);      
    case 3:
      return 9*M_PI*M_PI*sin(3*M_PI*t);
    default:
      return 0;
      break;
    }
  }
};


#endif	/* PERIODIC_TESTPROBLEM_H */

