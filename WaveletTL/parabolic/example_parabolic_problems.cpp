/* just a couple of test problems that are used frequently
TestProblem - a SimpleSturmBVP
 * TestRHS - 2D examples
 * Exact_Sol1D - 1D solutions and ...
 * TestRHS1D - 1D rihgthand sides corresponding to them
 * PertubedBVP - pertubed problem class
0 Operator
*/

/*
  different test problems with homogeneous Dirichlet b.c.'s
  1: y(t)=x*(1-x), -y''(t)=2
 */
template <unsigned int N>
class TestProblem
  : public SimpleSturmBVP
{
public:
  double p(const double t) const {
    switch(N) {
    case 1:
      return 1;
      break;
    default:
      return 0;
      break;
    }
  }
  double p_prime(const double t) const {
    switch(N) {
    case 1:
      return 0;
      break;
    default:
      return 0;
      break;
    }
  }
  double q(const double t) const {
    switch(N) {
    case 1:
      return 0;
      break;
    default:
      return 0;
      break;
    }
  }
  double g(const double t) const {
    switch(N) {
    case 1:
      return 2;
      break;
    default:
      return 0;
      break;
    }
  }
  bool bc_left() const { return true; }
  bool bc_right() const { return true; }
};



  /* Some 1D solutions for different problems
   * 1: -delta u = f; hom Dirichlet BC; u=x*(1-x); f = 2;
   * 2: (I-delta) u = f; hom Dirichlet BC; u=x*(1-x); f=x*(1-x)+2;
   * 3: (I-delta) u = f; hom Neumann BC; u = x^3-1.5x^2+d; f=x^3-1.5x^2-6x+3+d;
   * 4: u'=delta u + f; hom Dirichlet BC, f time independent; u=(1-exp(-pi^2t))sin(pix); f = pi^2sin(pix)
   * 5: u'=delta u + f; u = tcos(pix); f = (1+pi^2)cos(pix)
   * 6: u'=delta u + f; u = ??, u_0=hat; f=0;
   * 7: u_0 = \Xi_[0,0.5]; f=0;
   * 8: u'=delta u + f; u= exp(alpha*(1-t)) * sin^2(2pix); f= u_t-u_xx=exp(alpha(1-t))*(-alpha*sin^2+8pi^2sin^2-8pi^2cos^2);
   * 9: u'=delta u + f; u=p(t); f= p'(t); p(t) = exp(2t)
   * 10: u'= delta u + f; hom Neumann BC; u_x = 0; u= t; f= 1
   * 11: u'= delta u + f; hom Neumann BC; u_x = 0; f= 0, u= const
   * 12: u'' = f; hom Dirichlet; u = 4x(1-x)exp(-50(x-0.5)^2)
   * 13: u'= delta u + 0.5: u = 0.5t+u_0; u_0=1 // simple version of projects' problem (time forward case)
   * 14: u'= delta u - u + f; f = 0.5: u = (u_0-0.5)exp(-t)+0.5; u_0=1 // 2nd simple version of projects' problem (time forward case)
   * 15: u'= delta u + f    ; f = g(1-t) = -0.5t+0.5+g_0 with g being the u_13. u = -0.25t^2+t(0.5+g_0); u(0) = 0 // time backward case corresponding to #13
   * 16: u'= delta u - u + f; f = g(1-t) = (g_0-0.5)exp(t-1)+0.5; g=u_14.  u = 0.25*(exp(-t)-exp(t))-0.5*(exp(-t)-1)   ????u = exp(-t)*(-0.5-0.5*exp(-1)*(g_0-0.5))+exp(t)*(0.5*exp(-1)(g_0-0.5))+0.5
   * 17: w = max(1-1/c * abs(0.5-x),0) // Hat function. 1 at 0.5, 0 at 0.5 +- c
   * 18: w = max (0, -1/(c*c) * (0.5+c -x)*(0.5-c -x)) // piecewise quadratic, centered at 0.5, area = [0.5-c,0.5+c]
   * 19: w = 0 on [0,0.5-c]; #17 on [0.5-c,0.5]; #18 on [0.5,0.5+c]; 0 on [0.5+c;1] // piecewise combination of #17 and #18
   */
  template <unsigned int N>
  class Exact_Sol1D
  : public Function<1,double>
  {
  public:
      virtual ~Exact_Sol1D() {};
      double value(const Point<1>& p, const unsigned int component = 0) const {
          switch(N)
          {
              case 1:
                  return p[0]*(1-p[0]);
                  break;
              case 2:
                  return p[0]*(1-p[0]);
                  break;
              case 3:
                  return p[0]*p[0]*(-1.5+p[0]);
                  break;
              case 4:
                  return (1-exp(-M_PI*M_PI*get_time()))*sin(M_PI*p[0]);
                  break;
              case 5:
                  return get_time()*cos(M_PI*p[0]);
                  break;
              case 6:
                  return 0.5-abs(p[0]-0.5);
                  break;
              case 7:
                  return (p[0] <= 0.5)?1.0:0;
                  break;
              case 8:
                  return exp(1*(1-get_time()))*sin(2*M_PI*p[0])*sin(2*M_PI*p[0]);
                  break;
              case 9:
                  return exp(2*get_time());
                  break;
              case 10:
                  return get_time();
                  break;
              case 11:
                  return 1;
                  break;
              case 12:
                  return 4*p[0]*(1-p[0])*exp(-50*(p[0]-0.5)*(p[0]-0.5));
                  break;
              case 13:
                  return 0.5*get_time()+1;
                  break;
              case 14:
                  return 0.5*exp(-get_time())+0.5;
                  break;
              case 15:
                  return (-0.25*get_time()+1.5)*get_time();
                  break;
              case 16:
                  //return exp(-get_time())*(-0.5-0.25*exp(-1))+exp(get_time()-1)*0.25+0.5;
                  return -0.75*exp(-get_time())+0.25*exp(get_time())+0.5;
                  break;
              case 17:
                  return std::max(0.0, 1-8.0* std::abs(p[0]-0.5)); // hat on [0.375,0.625]
                  break;
              case 18:
                  return std::max(0.0,64.0*(0.625-p[0])*(0.375-p[0])); // wrong quadratic araound 0.5
                  // return std::max(0.0,64.0*(0.625-p[0])*(p[0]-0.375)); // quadratic araound 0.5
                  break;
              case 19:
                  return ((p[0] <= 0.5)?(std::max(0.0, 1-8.0* abs(p[0]-0.5))):(std::max(0.0,-64.0*(0.625-p[0])*(p[0]-0.375))));
                  break;
              case 20:
                  return ((p[0] <= 0.75)?1:(-1)); // should reproduce example from the paper
                  break;
              default:
                  abort();
                  break;
          }
      }
      void vector_value(const Point<1>& p, Vector<double>& values) const
      {
          values[0] = value(p);
      }
  };

  /* Some 1D righthandsides for different problems
   *  1: -delta u = f; hom Dirichlet BC; u=x*(1-x); f = 2;
     //1: u(x,y) = 2x^3-3x^2, -u''+u = f = 2x^3+3x^2-12x+g // now shifted to position 3
   *  2: (I-delta) u = f; hom Dirichlet BC; u=x*(1-x); f=x*(1-x)+2;
   *  3: (I-delta) u = f; hom Neumann BC; u = x^3-1.5x^2+d; f=x^3-1.5x^2-6x+3+d;
   *  4: u'= delta u + f; hom Dirichlet BC, f time independent; u=(1-exp(-pi^2t))sin(pix); f = pi^2sin(pix)
   *  5: u'= delta u + f; u = tcos(pix); f = (1+pi^2)cos(pix)
   *  6: u'= delta u + f; u = ??, u_0=hat; f=0;
   *  7: u_0 = \Xi_[0,0.5]; f=0;
   *  8: u'= delta u + f; hom Neumann BC; u= exp(alpha*(1-t)) * sin^2(2pix); f= u_t-u_xx=exp(alpha(1-t))*(-alpha*sin^2+8pi^2sin^2-8pi^2cos^2)
   *  9: u'= delta u + f; hom Neumann BC; u_x = 0; u=p(t); f= p'(t); p(t) = exp(2t)
   * 10: u'= delta u + f; hom Neumann BC; u_x = 0; u= t; f= 1
   * 11: u'= delta u + f; hom Neumann BC; u_x = 0; f= 0, u= const
   * 12: u = 8*x*(1-x)*exp(-50(x-0.5)^2), u'' = f; hom Dirichlet
   * 13: u'= delta u + f    ; f = 0.5: u = 0.5t+u_0; u_0=1 // simple version of projects' problem (time forward case)
   * 14: u'= delta u - u + f; f = 0.5: u = (u_0-0.5)exp(-t)+0.5; u_0=1 // 2nd simple version of projects' problem (time forward case)
   * 15: u'= delta u + f    ; f = g(1-t) = -0.5t+0.5+g_0 with g being the u_13. u = -0.25t^2+t(0.5+g_0); u(0) = 0 // time backward case corresponding to #13
   * 16: u'= delta u - u + f; f = g(1-t) = (g_0-0.5)exp(t-1)+0.5; g=u_14. u = exp(-t)*(-0.5-0.5*exp(-1)*(g_0-0.5))+exp(t)*(0.5*exp(-1)(g_0-0.5))+0.5
   */
  template <unsigned int N>
  class TestRHS1D
  : public Function<1,double>

  {
  public:
      virtual ~TestRHS1D() {};
      double value(const Point<1>& p, const unsigned int component = 0) const
      {
          switch(N) {
              case 1:
                  return 2.0;
                  break;
                  //return p[0]*( p[0]*( p[0]*2+3)-12)-6;
                  //break;
              case 2:
                  return p[0]*(1-p[0])+2.0;
                  break;
              case 3:
                  return 3.0 + p[0]*(-6.0+p[0]*(-1.5+p[0]));
                  break;
              case 4:
                  return M_PI*M_PI*sin(M_PI*p[0]);
                  break;
              case 5:
                  return (1+M_PI*M_PI*get_time())*cos(M_PI*p[0]);
                  break;
              case 6:
              case 7:
                  return 0;
                  break;
              case 8:
                  return exp(1*(1-get_time()))*( (8*M_PI*M_PI-1)*sin(2*M_PI*p[0])*sin(2*M_PI*p[0]) -8*M_PI*M_PI*cos(2*M_PI*p[0])*cos(2*M_PI*p[0]) ) ;
                  break;
              case 9:
                  return 2*exp(2*get_time());
                  break;
              case 10:
                  return 1;
                  break;
              case 11:
                  return 0;
                  break;
              case 12:
                  return 8*(-49-1000*p[0]*(1+p[0]*(-6+p[0]*(10-5*p[0]))))*exp(-50*(p[0]-0.5)*(p[0]-0.5));
                  break;
              case 13:
//                  return 0.5;
//                  break;
              case 14:
                  return 0.5;
                  break;
              case 15:
                  return -0.5*get_time()+1.5;
                  break;
              case 16:
                  return 0.5*exp(get_time()-1)+0.5;
                  break;
              default:
                  abort();
                  break;
          }
      }
      void vector_value(const Point<1>& p, Vector<double>& values) const
      {
          values[0] = value(p);
      }
  };

/*
  Some test problems on the interval
  1: u(x,y) = 2x^3-3x^2, -u''+u = f = 2x^3+3x^2-12x+g
  2: u = 4*x*(1-x)*exp(-50(x-0.5)^2), u'' = f
*/
  /*
template <unsigned int N>
class TestRHS1D
  : public Function<1,double>
{
public:
  virtual ~TestRHS1D() {};
  double value(const Point<1>& p, const unsigned int component = 0) const {
    switch(N) {
    case 1:
        return p[0]*( p[0]*( p[0]*2+3)-12)-6;
      break;
    case 2:
        return -4*(-49-1000*p[0]*(1+p[0]*(-6+p[0]*(10-5*p[0]))))*exp(-50*(p[0]-0.5)*(p[0]-0.5));
        break;
    default:
      abort();
      break;
    }
  }


  void vector_value(const Point<1>& p, Vector<double>& values) const {
    values[0] = value(p);
  }
};
*/


/*
  Some driving terms for the parabolic problem on the cube
incomplete. each entry corresponds to TestRHS

2: u(x,y) = exp(-50*((x-0.5)^2+(y-0.5)^2)), -Delta u(x,y)= (200-(100x-50)^2-(100y-50)^2)*u(x,y); nearly homogeneous Dirichlet BC
  3: u(x,y) = x(1-x)^2y^2(1-y), -Delta u(x,y)= 4*(1-x)*y^2*(1-y)-2*x*y^2*(1-y)-2*x*(1-x)^2*(1-y)+4*x*(1-x)^2*y

  7: u(x,y) = (1-exp(-pi^2 t))*sin(pi*x)*sin(pi*y), f = 2*pi*sin(pi*x)*sin(pi*y); u_t = u'' + f; homogeneous Dirichlet BC

 8: -delta u + u = f : u = (x^2*(x-1)^2*y^2*(y-1)^2, -delta u + u = f = - (12x^2-12*x+2)*y^2*(y-1)^2 - (12y^2-12*y+2)*x^2*(x-1)^2 + x^2*(x-1)^2*y^2*(y-1)^2 homogeneous Neumann BC
 9: u = t*x^2(x-1)^2y^2(y-1)^2+10; u_t = delta u+f; f = - (12x^2-12*x+2)*y^2*(y-1)^2*t - (12y^2-12*y+2)*x^2*(x-1)^2*t + x^2*(x-1)^2*y^2*(y-1)^2
   * 13: u'= delta u + f    ; f = 0.5: u = 0.5t+u_0; u_0=1 // simple version of projects' problem (time forward case)
   * 14: u'= delta u - u + f; f = 0.5: u = (u_0-0.5)exp(-t)+0.5; u_0=1 // 2nd simple version of projects' problem (time forward case)
   * 15: u'= delta u + f    ; f = g(1-t) = -0.5t+0.5+g_0 with g being the u_13. u = -0.25t^2+t(0.5+g_0); u(0) = 0 // time backward case corresponding to #13
   * 16: u'= delta u - u + f; f = g(1-t) = (g_0-0.5)exp(t-1)+0.5; g=u_14. u = exp(-t)*(-0.5-0.5*exp(-1)*(g_0-0.5))+exp(t)*(0.5*exp(-1)(g_0-0.5))+0.5
  // u= sum of 2 gauss kernels
*/

template <unsigned int N>
class Exact_Sol2D
: public Function<2,double>
{
public:
    virtual ~Exact_Sol2D() {};
    double value(const Point<2>& p, const unsigned int component = 0) const {
        switch(N) {
            case 1:
                return 1;
                break;
            case 2:
                return exp(-50*((p[0]-0.5)*(p[0]-0.5)+(p[1]-0.5)*(p[1]-0.5)));
                break;
            case 3:
                return p[0]*(1-p[0])*p[1]*(1-p[1]);
                break;
            case 7:
                return (1-exp(-2*M_PI*M_PI*get_time()))*sin(M_PI*p[0])*sin(M_PI*p[1]);
                break;
            case 8:
                //return (2*p[0]*p[0]*p[0]-3*p[0]*p[0]+2)*(2*p[1]*p[1]*p[1]-3*p[1]*p[1]+2);
                return p[0]*p[0]*(p[0]-1)*(p[0]-1)*p[1]*p[1]*(p[1]-1)*(p[1]-1);
                break;
            case 9:
                return get_time()*p[0]*p[0]*(p[0]-1)*(p[0]-1)*p[1]*p[1]*(p[1]-1)*(p[1]-1)+10.0;
                break;
            case 13:
                return 0.5*get_time()+1;
                break;
            case 14:
                return 0.5*exp(-get_time())+0.5;
                break;
            case 15:
                return (-0.25*get_time()+1.5)*get_time();
                break;
            case 16:
                return exp(-get_time())*(-0.5-0.25*exp(-1))+exp(get_time()-1)*0.25+0.5;
                break;
            case 17:
                return (((p[0] >= 0.375)&&(p[0] <= 0.625))&&((p[1] >= 0.375)&&(p[1] <= 0.625)))
                        ?((1-8.0* std::abs(p[0]-0.5)) * (1-8.0* std::abs(p[1]-0.5))) : (0.0); // hat around 0.5
                break;
            case 18:
                return ((std::max(0.0,64.0*(0.75-p[0])*(p[0]-0.5)))
                       *(std::max(0.0,64.0*(0.75-p[1])*(p[1]-0.5)))); // quadratic on (0.5, 0.75)^2
                break;
            case 19:
                return ((std::max(0.0,64.0*(0.75-p[0])*(p[0]-0.5)))
                       *(std::max(0.0,64.0*(0.5-p[1])*(p[1]-0.25)))); // quadratic on (0.5, 0.75)x(0.25,0.5)
                break;
            case 20:
                return ((std::max(0.0,64.0*(0.5-p[0])*(p[0]-0.25)))
                       *(std::max(0.0,64.0*(0.5-p[1])*(p[1]-0.25)))); // quadratic on (0.25, 0.5)^2
                break;
            case 21:
                return ((std::max(0.0,64.0*(0.5-p[0])*(p[0]-0.25)))
                       *(std::max(0.0,64.0*(0.75-p[1])*(p[1]-0.5)))); // quadratic on (0.25, 0.5)x(0.5,0.75)
                break;
            case 22: // symmetric case
                return ((std::max(0.0,64.0*(0.625-p[0])*(p[0]-0.375)))
                       *(std::max(0.0,64.0*(0.625-p[1])*(p[1]-0.375)))); // quadratic on (0.375,0.625)^2
                break;
                /*
            case 19:
                return ((p[0] >= 0.375)
                         ? ( (p[0] <= 0.5)
                              ? (1-8.0* std::abs(p[0]-0.5))
                              : ( (p[0] <= 0.625)
                                   ? (64.0*(0.625-p[0])*(p[0]-0.375))
                                   : (0.0)
                                )
                           )
                         : (0.0)
                       )
                        *
                       ((p[1] >= 0.375)
                         ? ( (p[1] <= 0.5)
                              ? (1-8.0* std::abs(p[1]-0.5))
                              : ( (p[1] <= 0.625)
                                   ? (64.0*(0.625-p[1])*(p[1]-0.375))
                                   : (0.0)
                                )
                           )
                         : (0.0)
                       );
              break;
                 * */
              /*
            case 20: // a single generator
                return (((p[0] >= 0.5)&&(p[0] <= 0.75))&&((p[1] >= 0.5)&&(p[1] <= 0.75)))?(1.0):(0.0);
                break;
                */
                /*
                 * case X:
                 * return 2*exp(-50*((p[0]-0.2)*(p[0]-0.2)+(p[1]-0.3)*(p[1]-0.3))) + 3*exp(-40*((p[0]-0.8)*(p[0]-0.8)+(p[1]-0.6)*(p[1]-0.6)));
                 * break;
                 */
                /*
            case 21:
                return std::max(0.0,  ((p[0]-0.25)*(0.75-p[0])) )
                       *std::max(0.0,  ((p[1]-0.25)*(0.75-p[1])) )
                       *256.0;
                break;
                 * */
            case 23:
                return p[0];
                break;
            case 24:
                return p[1];
                break;
            case 25:
                return 1-p[0];
                break;
            case 26:
                return 1-p[1];
                break;
            case 27: // bumps: centered in one direction, assymetric in the other
                return (((p[0] >= 0.375)&&(p[0] <= 0.625))&&((p[1] >= 0.625)&&(p[1] <= 0.875)))?(10.0):(0.0);
                break;
            case 28:
                return (((p[0] >= 0.375)&&(p[0] <= 0.625))&&((p[1] >= 0.125)&&(p[1] <= 0.375)))?(10.0):(0.0);
                break;
            case 29:
                return (((p[0] >= 0.125)&&(p[0] <= 0.375))&&((p[1] >= 0.375)&&(p[1] <= 0.625)))?(10.0):(0.0);
                break;
            case 30:
                return (((p[0] >= 0.625)&&(p[0] <= 0.875))&&((p[1] >= 0.375)&&(p[1] <= 0.625)))?(10.0):(0.0);
                break;
            case 31: // haar_wavelet # 57
                return ((p[0] >= 0.5)?
                                     ((p[0]<0.625)? (1.0) :( (p[0]<0.75)?(-1.0):(0.0) ))
                                    :(0.0)
                       )*
                       ((p[1] >= 0.25)?
                                     ((p[1]<0.375)? (1.0) :( (p[1]<0.5)?(-1.0):(0.0) ))
                                    :(0.0)
                       );
                break;
            case 32:
                return std::max(0.0, std::min(4.0*(p[0]-0.25),1.0) )
                        * std::max(0.0, std::min(1.5-4.0*abs(p[1]-0.375),1.0) );
                break;
            default:
                cout << "test_problems::Exact_Sol2D unspecified case. aborting." << endl;
                abort();
                break;
        }
    }
    void vector_value(const Point<2>& p, Vector<double>& values) const {
        values[0] = value(p);
    }
};


/*
  Some test problems for the Poisson equation on the cube with homogeneous Dirichlet b.c.'s:
  1: u(x,y) = x(1-x)y(1-y), -Delta u(x,y) = 2(x(1-x)+y(1-y)); homogeneous Dirichlet BC
  2: u(x,y) = exp(-50*((x-0.5)^2+(y-0.5)^2)), -Delta u(x,y)= (200-(100x-50)^2-(100y-50)^2)*u(x,y); nearly homogeneous Dirichlet BC
  3: u(x,y) = x(1-x)^2y^2(1-y), -Delta u(x,y)= 4*(1-x)*y^2*(1-y)-2*x*y^2*(1-y)-2*x*(1-x)^2*(1-y)+4*x*(1-x)^2*y
  4: f= (1-x)*(1-y)  driving term for Robs example
  5: -delta u + u = f : u = (2x^3-3x^2+1)*(2y^3-3y^2+2), -delta u + u = f =-(12x-6)(2y^3-3y^2+2) -(12y-6)(2x^3-3x^2+2) +  (2x^3-3x^2+2)(2y^3-3y^2+2)
  6: u = x*(1-x)*y*(1-y)*exp(-40*((x-0.8)^2+(y-0.6)^2)), f=-delta u; homogeneous Dirichlet BC
  7: u(x,y) = (1-exp(-pi^2 t))*sin(pi*x)*sin(pi*y), f = 2*pi*sin(pi*x)*sin(pi*y); u_t = u'' + f; homogeneous Dirichlet BC
  8: -delta u + u = f : u = (x^2*(x-1)^2*y^2*(y-1)^2, -delta u + u = f = - (12x^2-12*x+2)*y^2*(y-1)^2 - (12y^2-12*y+2)*x^2*(x-1)^2 + x^2*(x-1)^2*y^2*(y-1)^2 homogeneous Neumann BC
  9: u = t*x^2(x-1)^2y^2(y-1)^2+10; u_t = delta u+f; f = - (12x^2-12*x+2)*y^2*(y-1)^2*t - (12y^2-12*y+2)*x^2*(x-1)^2*t + x^2*(x-1)^2*y^2*(y-1)^2
 13: u'= delta u + f    ; f = 0.5: u = 0.5t+u_0; u_0=1 // simple version of projects' problem (time forward case)
 14: u'= delta u - u + f; f = 0.5: u = (u_0-0.5)exp(-t)+0.5; u_0=1 // 2nd simple version of projects' problem (time forward case)
 15: u'= delta u + f    ; f = g(1-t) = -0.5t+0.5+g_0 with g being the u_13. u = -0.25t^2+t(0.5+g_0); u(0) = 0 // time backward case corresponding to #13
 16: u'= delta u - u + f; f = g(1-t) = (g_0-0.5)exp(t-1)+0.5; g=u_14. u = exp(-t)*(-0.5-0.5*exp(-1)*(g_0-0.5))+exp(t)*(0.5*exp(-1)(g_0-0.5))+0.5
*/
template <unsigned int N>
class TestRHS2D
: public Function<2,double>
{

public:
    virtual ~TestRHS2D() {};
    double value(const Point<2>& p, const unsigned int component = 0) const {
        switch(N) {
            case 2:
                return (200.-(100.*p[0]-50.)*(100.*p[0]-50.)-(100.*p[1]-50.)*(100.*p[1]-50.))
                        * exp(-50.*((p[0]-0.5)*(p[0]-0.5)+(p[1]-0.5)*(p[1]-0.5)));
                break;
            case 3:
                return 4*(1-p[0])*p[1]*p[1]*(1-p[1])
                        - 2*p[0]*p[1]*p[1]*(1-p[1])
                        - 2*p[0]*(1-p[0])*(1-p[0])*(1-p[1])
                        + 4*p[0]*(1-p[0])*(1-p[0])*p[1];
                break;
            case 4:
                return (1-p[0])*(1-p[1]);
                break;
            case 5:
                return -(12*p[0]-6)*(2*p[0]*p[0]*p[0]-3*p[0]*p[0]+2)
                        -(12*p[1]-6)*(2*p[1]*p[1]*p[1]-3*p[1]*p[1]+2)
                        +(2*p[0]*p[0]*p[0]-3*p[0]*p[0]+2)*(2*p[1]*p[1]*p[1]-3*p[1]*p[1]+2);
                break;
            case 6:
                return exp(-40*((p[0]-0.8)*(p[0]-0.8)+(p[1]-0.6)*(p[1]-0.6)))*( 126*p[1]*(-1+p[1]) + p[0]*( -94+p[1]*(-5472+p[1]*(13184+p[1]*(-14080+6400*p[1]))) + p[0]*( 94+p[1]*(15808+p[1]*(-23520+p[1]*(14080-6400*p[1]))) +p[0]*(-16640*p[1]*(1-p[1]) +p[0]*6400*p[1]*(1-p[1])  )))     );
                break;
            case 7:
                return 2*M_PI*M_PI*sin(M_PI*p[0])*sin(M_PI*p[1]);
                break;
            case 8:
                return -(12*p[0]*p[0]-12*p[0]+2)*p[1]*p[1]*(p[1]-1)*(p[1]-1)
                        -(12*p[1]*p[1]-12*p[1]+2)*p[0]*p[0]*(p[0]-1)*(p[0]-1)
                        + p[0]*p[0]*(p[0]-1)*(p[0]-1)*p[1]*p[1]*(p[1]-1)*(p[1]-1); //
                break;
            case 9:
                return -(12*p[0]*p[0]-12*p[0]+2)*p[1]*p[1]*(p[1]-1)*(p[1]-1)*get_time()
                        -(12*p[1]*p[1]-12*p[1]+2)*p[0]*p[0]*(p[0]-1)*(p[0]-1)*get_time()
                        + p[0]*p[0]*(p[0]-1)*(p[0]-1)*p[1]*p[1]*(p[1]-1)*(p[1]-1);
            case 13:
            case 14:
                return 0.5;
                break;
            case 15:
                return -0.5*get_time()+1.5;
                break;
            case 16:
                return 0.5*exp(get_time()-1)+0.5;
                break;
            case 1:
            default:
                return 2*(p[0]*(1-p[0])+p[1]*(1-p[1]));
                break;
        }
    }
    void vector_value(const Point<2>& p, Vector<double>& values) const {
        values[0] = value(p);
    }
};



template <unsigned int DIM>
  class PertubedBVP
    : public EllipticBVP<DIM>
  {
  public:
    /*!
      constructor with given right-hand side
    */
    PertubedBVP(const Function<DIM>* f, const double eta)
    : eta_(eta), EllipticBVP<DIM> (f, f, f)
    {
    }

    /*!
      diffusion coefficient a
     */
    const double a(const Point<DIM>& x) const { return eta_; }

    /*!
      reaction coefficient q
    */
    const double q(const Point<DIM>& x) const { return 1.0; }

    /*!
      flag for constant coefficients
    */
    const bool constant_coefficients() const { return true; }
  protected:
      const double eta_;
  };

/* 0 Operator*/
template <unsigned int DIM>
class ZeroBVP
: public EllipticBVP<DIM>
{
public:
    /*!
      constructor with given right-hand side
    */
    ZeroBVP(const Function<DIM>* f)
    : EllipticBVP<DIM> (f, f, f)
    {
    }

    /*!
      diffusion coefficient a
     */
    const double a(const Point<DIM>& x) const { return 0; }

    /*!
      reaction coefficient q
    */
    const double q(const Point<DIM>& x) const { return 0; }

    /*!
      flag for constant coefficients
    */
    const bool constant_coefficients() const { return true; }
};

/* a problem that aborts if it is used. used to check that indeed all attributes are unused or overwritten*/
template <unsigned int DIM>
class AbortBVP
: public EllipticBVP<DIM>
{
public:
    /*!
      constructor with given right-hand side
    */
    AbortBVP(const Function<DIM>* f)
    : EllipticBVP<DIM> (f, f, f)
    {
    }

    /*!
      diffusion coefficient a
     */
    const double a(const Point<DIM>& x) const { abort(); return 0; }

    /*!
      reaction coefficient q
    */
    const double q(const Point<DIM>& x) const { abort(); return 0; }

    /*!
      flag for constant coefficients
    */
    const bool constant_coefficients() const { abort(); return true; }
};
