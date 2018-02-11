#ifndef TESTFUNCTIONS_H
#define	TESTFUNCTIONS_H

#include <utils/function.h>

/* some 1-periodic test functions
1: 
2: 
3: 
 */

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

// f(x)=x*x*(1-x)
class Function2a : public Function<1> {
public:
  inline double value(const Point<1>& p, const unsigned int component = 0) const {
    return p[0]*p[0]*(1-p[0]);
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

// f(x)=-sin(3*pi*x)-4
class Function4 : public Function<1> {
public:
  inline double value(const Point<1>& p, const unsigned int component = 0) const {
    return -sin(3*M_PI*p[0]);
  }


  void vector_value(const Point<1> &p, Vector<double>& values) const {
    values.resize(1, false);
    values[0] = value(p);
  }
};  


// f(x)=-sin(3*pi*x)+2x^2 or 2(1-x)^2)
class SolutionToDeltadis : public Function<1> {
public:
  inline double value(const Point<1>& p, const unsigned int component = 0) const {
    return -sin(3*M_PI*p[0])+(p[0]<0.5 ? 2*p[0]*p[0] : 2*(1-p[0])*(1-p[0]));
  }


  void vector_value(const Point<1> &p, Vector<double>& values) const {
    values.resize(1, false);
    values[0] = value(p);
  }
};  

// f(x)=x
class Function5 : public Function<1> {
public:
  inline double value(const Point<1>& p, const unsigned int component = 0) const {
    return p[0];
  }
  
  void vector_value(const Point<1> &p, Vector<double>& values) const {
    values.resize(1, false);
    values[0] = value(p);
  }  
};  

//solution to sin(3*pi*x)-projection
class Function7 : public Function<1> {
public:
  inline double value(const Point<1>& p, const unsigned int component = 0) const {
    return sin(3*M_PI*p[0])+3*M_PI*p[0]*p[0]-3*M_PI*p[0]+0.5*M_PI-2./(3*M_PI);
  }
  
  void vector_value(const Point<1> &p, Vector<double>& values) const {
    values.resize(1, false);
    values[0] = value(p);
  }
};

class Function8 : public Function<1> {
public:
  inline double value(const Point<1>& p, const unsigned int component = 0) const {
    return sin(2*M_PI*p[0]);
  }
  
  void vector_value(const Point<1> &p, Vector<double>& values) const {
    values.resize(1, false);
    values[0] = value(p);
  }
};

//solution to -u''=x-0.5 with periodic b.c.
class Function6 : public Function<1> {
public:
  inline double value(const Point<1>& p, const unsigned int component = 0) const {
    return -1./6*p[0]*p[0]*p[0]+1./4*p[0]*p[0]-1./12*p[0];
  }
  
  void vector_value(const Point<1> &p, Vector<double>& values) const {
    values.resize(1, false);
    values[0] = value(p);
  }  
};  



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

//scaled hat function
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

//scaled hat function
class scaledHat2a : public Function<1>
{
    public:
        inline double value(const Point<1>& p,
                            const unsigned int component = 0) const
        {
            return std::max(0.0, 4-32*abs(p[0]-0.125))/sqrt(2);
        }

        void vector_value(const Point<1> &p,
                          Vector<double>& values) const
        {
            values.resize(1, false);
            values[0] = value(p);
        }
};
//scaled hat function (d=2, j=3, k=1)
class scaledHat2b : public Function<1>
{
    public:
        inline double value(const Point<1>& p,
                            const unsigned int component = 0) const
        {
            return std::max(0.0, 4-16*abs(p[0]-0.25))/sqrt(2);
        }

        void vector_value(const Point<1> &p,
                          Vector<double>& values) const
        {
            values.resize(1, false);
            values[0] = value(p);
        }
};

//scaled Quark function (d=2, p=1, j=3, k=1)
class scaledQuark : public Function<1>
{
    public:
        inline double value(const Point<1>& p,
                            const unsigned int component = 0) const
        {
            return (8*p[0]-1)*std::max(0.0, 4-32*abs(p[0]-0.125))/sqrt(2);
        }

        void vector_value(const Point<1> &p,
                          Vector<double>& values) const
        {
            values.resize(1, false);
            values[0] = value(p);
        }
};



//f(x)=exp(-100*(x-0.5)^2)
class Function2b : public Function<1>
{
    public:
        inline double value(const Point<1>& p,
                            const unsigned int component = 0) const
        {
            return exp(-100*(p[0]-0.5)*(p[0]-0.5));
        }

        void vector_value(const Point<1> &p,
                          Vector<double>& values) const
        {
            values.resize(1, false);
            values[0] = value(p);
        }
};

class FunctionBBCCDDU : public Function<1>
{
    public:
        inline double value(const Point<1>& p,
                            const unsigned int component = 0) const
        {
            return 4*((exp(5.0*p[0])-1)*(1-(exp(5.0*p[0])-1)/(exp(5.0)-1)))/(exp(5.0)-1);
        }

        void vector_value(const Point<1> &p,
                          Vector<double>& values) const
        {
            values.resize(1, false);
            values[0] = value(p);
        }
};
#endif	/* PERIODIC_TESTFUNCTIONS_H */

