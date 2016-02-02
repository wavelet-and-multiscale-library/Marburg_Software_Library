/**
 * Hat function
 */
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



/**
 * Cut sqrt function
 */
class CutSqrt : public Function<1>
{
     public:
         inline double value(const Point<1>& p,
                             const unsigned int component = 0) const
         {
             return (p[0] < 0.1 || p[0] > 0.75
                     ? 0
                     : sqrt(p[0]-0.1));
         }

         void vector_value(const Point<1> &p,
                           Vector<double>& values) const
         {
             values.resize(1, false);
             values[0] = value(p);
         }
};



/**
 * Cut pwr function
 */
class CutPwr : public Function<1>
{
     public:
         inline double value(const Point<1>& p,
                             const unsigned int component = 0) const
         {
             return (p[0] < 0.1 || p[0] > 0.75
                     ? 0
                     : pow(p[0]-0.1,0.75));
         }

         void vector_value(const Point<1> &p,
                           Vector<double>& values) const
         {
             values.resize(1, false);
             values[0] = value(p);
         }
};



/**
 * Cut lin function
 */
class CutLin : public Function<1>
{
     public:
         inline double value(const Point<1>& p,
                             const unsigned int component = 0) const
         {
             return (p[0] < 0.1 || p[0] > 0.75
                     ? 0
                     : 2*(p[0]-0.1));
         }

         void vector_value(const Point<1> &p,
                           Vector<double>& values) const
         {
             values.resize(1, false);
             values[0] = value(p);
         }
};



/**
 *  Kink at 0 < a < 1
 */
class Kink : public Function<1>
{
     public:
         Kink(const double a = 0.5) : a_(a) {}

         inline double value(const Point<1>& p,
                             const unsigned int component = 0) const
         {
             if (0. <= p[0] && p[0] < a_)
                 return 1/(2*a_*a_)*p[0]*p[0];

             if (a_ <= p[0] && p[0] <= 1.0)
                 return 0.5*(1-(p[0]-a_)/(1-a_))*(1-(p[0]-a_)/(1-a_));

             return 0.;
         }

         void vector_value(const Point<1> &p,
                           Vector<double>& values) const
         {
             values.resize(1, false);
             values[0] = value(p);
         }

     protected:
         double a_;
};