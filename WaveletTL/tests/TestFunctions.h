#ifndef TESTFUNCTIONS_H
#define	TESTFUNCTIONS_H

#include <utils/function.h>

// some 1-periodic test functions

// f(x)=1
class Function1 : public Function<1> {
public:
  inline double value(const Point<1>& p, const unsigned int component = 0) const {
    return 1.0;
  }
  
  void vector_value(const Point<1> &p, Vector<double>& values) const {
    values.resize(1, false);
    values[0] = value(p);
  }
};

// f(x)=x*(1-x)
class Function2 : public Function<1> {
public:
  inline double value(const Point<1>& p, const unsigned int component = 0) const {
    return p[0]*(1-p[0]);
  }
  
  void vector_value(const Point<1> &p, Vector<double>& values) const {
    values.resize(1, false);
    values[0] = value(p);
  }
};

// f(x)=cos(2*pi*x)
class Function3 : public Function<1> {
public:
  inline double value(const Point<1>& p, const unsigned int component = 0) const {
    return cos(2*M_PI*p[0]);
  }
  
  void vector_value(const Point<1> &p, Vector<double>& values) const {
    values.resize(1, false);
    values[0] = value(p);
  }
};

// f(x)=sin(3*pi*x)
class Function4 : public Function<1> {
public:
  inline double value(const Point<1>& p, const unsigned int component = 0) const {
    return sin(3*M_PI*p[0]);
  }
  
  void vector_value(const Point<1> &p, Vector<double>& values) const {
    values.resize(1, false);
    values[0] = value(p);
  }
//solution to sin(3*pi*x)-projection
//class Function5 : public Function<1> {
//public:
//  inline double value(const Point<1>& p, const unsigned int component = 0) const {
//    return sin(3*M_PI*p[0])+3*M_PI*p[0]*p[0]-3*M_PI*p[0]+0.5*M_PI-2./(3*M_PI);
//  }
//  
//  void vector_value(const Point<1> &p, Vector<double>& values) const {
//    values.resize(1, false);
//    values[0] = value(p);
//  }
//};


//hat function
class Hat : public Function<1>
{
    public:
        inline double value(const Point<1>& p,
                            const unsigned int component = 0) const
        {
            return (std::max(0.0,0.5-abs(p[0]-0.5)));
        }

        void vector_value(const Point<1> &p,
                          Vector<double>& values) const
        {
            values.resize(1, false);
            values[0] = value(p);
        }
};

//hat function
class scaledHat : public Function<1>
{
    public:
        inline double value(const Point<1>& p,
                            const unsigned int component = 0) const
        {
            return (std::max(0.0,2-8*abs(p[0]-0.5)));
        }

        void vector_value(const Point<1> &p,
                          Vector<double>& values) const
        {
            values.resize(1, false);
            values[0] = value(p);
        }
};

#endif	/* PERIODIC_TESTFUNCTIONS_H */

