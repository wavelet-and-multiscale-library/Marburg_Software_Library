#include <string>

/**
 * TestProblem
 * Different test problems with homogeneous Dirichlet b.c.'s on [0,1]
 *
 *  N =
 *      1: y(t) = x*(1-x),  -y''(t) = 2
 *      2:   y(t) = exp(-100*(x-0.5)^2),
 *        -y''(t) = (200-(200x-100)^2)*exp(-100*(x-0.5)^2)
 *      3: 1D example from [BBCCDDU]
 */
template <unsigned int N>
class TestProblem : public SimpleSturmBVP
{
    public:
        double p(const double t) const
        {
            switch(N)
            {
                case 1:
                    return 1;
                    break;
                case 2:
                    return 1;
                    break;
                case 3:
                    return 1;
                    break;
                case 4:
                    return 1;
                    break;
                case 5:
                    return 0;
                    break;
               case 6:
                    return 0;
                    break;
                case 7:
                    return 0;
                    break;
                case 8:
                    return 0;
                    break;
                case 9:
                    return 1;
                    break;
                case 10:
                    return 0;
                    break;
                case 11:
                    return 1;
                    break;
                case 12:
                    return 1;
                    break;
                case 13:
                    return 0;
                    break;
                default:
                    return 0;
                    break;
            }
        }

        double p_prime(const double t) const
        {
            switch(N)
            {
                case 1:
                    return 0;
                    break;
                case 2:
                    return 0;
                    break;
                case 3:
                    return 0;
                    break;
                case 4:
                    return 0;
                    break;
                case 5:
                    return 0;
                    break;
                case 6:
                    return 0;
                    break;
                case 7:
                    return 0;
                    break;
                case 8:
                    return 0;
                    break;
                case 9:
                    return 0;
                    break;
                case 10:
                    return 0;
                    break; 
                case 11:
                    return 0;
                    break;    
                case 12:
                    return 0;
                    break;  
                case 13:
                    return 0;
                    break;
                default:
                    return 0;
                    break;
            }
        }

        double q(const double t) const
        {
            switch(N)
            {
                case 1:
                    return 0;
                    break;
                case 2:
                    return 0;
                    break;
                case 3:
                    return 0;
                    break;
                case 4:
                    return 0;
                    break;
                case 5:
                    return 1;
                    break;
                case 6:
                    return 1;
                    break;
                case 7:
                    return 1;
                    break;
                case 8:
                    return 1;
                    break;
                case 9:
                    return 0;
                    break;
                case 10:
                    return 1;
                    break; 
                case 11:
                    return 0;
                    break;    
                case 12:
                    return 0;
                    break;
                case 13:
                    return 1;
                    break;
                default:
                    return 0;
                    break;
            }
        }

        double g(const double t) const
        {
            switch(N)
            {
                case 1:
                    return 2;
                    break;
                case 2:
                    return ( (200.0-(200.0*t-100.0)*(200.0*t-100.0)) *
                                            exp(-100.0*(t-0.5)*(t-0.5)) );
                    break;
                case 3:
                    return ( -100*exp(5*t)*(1-(exp(5*t)-1)/(exp(5.)-1)) /
                             (exp(5.)-1)+200*exp(10*t)/((exp(5.)-1) *
                             (exp(5.)-1))+100*(exp(5*t)-1)*exp(5*t) /
                             ((exp(5.)-1)*(exp(5.)-1)) );
                    break;
                case 4:
                    return -9*M_PI*M_PI*sin(3*M_PI*t);  
                    break;
                case 5:
                    return t*t*(1-t);
                case 6:
                    return std::max(0.0, 4-32*abs(t-0.125))/sqrt(2);
                    break;
                case 7:
                    return std::max(0.0, 4-16*abs(t-0.25))/sqrt(2);
                    break;
                case 8:
                    return (8*t-1)*std::max(0.0, 4-32*abs(t-0.125))/sqrt(2);
                    break;
                case 9:
                    return  -sin(3.*M_PI*t)*9.*M_PI*M_PI - 4.;
                    break;                    
                case 10:
                    return -sin(3*M_PI*t)+(t<0.5 ? 2*t*t : 2*(1-t)*(1-t));  
                    break;
                case 11:
                    return -2*(t-1)*(t-1)-8*t*(t-1)-2*t*t;  
                    break;    
                case 12:
                    return -cos(2.*M_PI*t)*4.*M_PI*M_PI - 4.;  
                    break; 
                case 13:
                    return pow(t,0.75)-t;
                default:
                    return 0;
                    break;
            }
        }

        string toString() const
        {
            switch(N)
            {
                case 1:
                    return "y(t) = x*(1-x), -y''(t) = 2, y(0) = y(1) = 0";
                    break;
                case 2:
                    return "y(t) = exp(-100*(x-0.5)^2), -y''(t) = (200-(200x-100)^2)*exp(-100*(x-0.5)^2), y(0) = y(1) = 0";
                    break;
                case 3:
                    return "1D example from [BBCCDDU]";
                    break;
                case 4:
                     return "y(t) = sin(3*M_PI*t), -y''(t) = 9*M_PI*M_PI*sin(3*M_PI*t), y(0) = y(1) = 0";  
                     break;
                case 6:
                     return "Hat function as gramian problem (for test purposes)";  
                     break;
                case 11:
                     return "Neumann BC test problem";  
                     break;     
                case 12:
                     return "Neumann BC test problem 2";  
                     break;  
                case 13:
                    return "identity bvp of x^alpha singularity function";
                    break;
                default:
                    return "TestProblem: N not defined.";
                    break;
            }
        }

  bool bc_left() const { return true; }   /* boundary conditions left: true */
  bool bc_right() const { return true; }  /* boundary conditions right: true */
};