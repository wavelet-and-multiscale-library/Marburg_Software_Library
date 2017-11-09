#include <iostream>
#include <fstream>

#include <geometry/sampled_mapping.h>
#include <utils/fixed_array1d.h>

#include <interval/p_basis.h>
#include <interval/ds_basis.h>

#include <utils/function.h>
#include <cube/tbasis.h>
#include <cube/tbasis_indexplot.h>
#include <cube/tbasis_evaluate.h>

using namespace std;
using namespace WaveletTL;

class Plateau
  : public Function<2>
{
  public:
  inline double value(const Point<2>& p,
		      const unsigned int component = 0) const
  {
    return (p[0] < 0.25 || p[0] > 0.75 ? 0 : 1)
      * (p[1] < 0.25 || p[1] > 0.75 ? 0 : 1);
  }

  void vector_value(const Point<2> &p,
		    Vector<double>& values) const
  {
    values.resize(1, false);
    values[0] = value(p);
  }
};

class Jump_and_Singularity
  : public Function<2>
{
  public:
  inline double value(const Point<2>& p,
		      const unsigned int component = 0) const
  {
    Point<2> x;
    x[0] = 0.5;
    x[1] = 0.5;

    return atan2((p-x)[0],(p-x)[1]);
  }

  void vector_value(const Point<2> &p,
		    Vector<double>& values) const
  {
    values.resize(1, false);
    values[0] = value(p);
  }
};

class My_Constant
:public Function<2>
{
public:
    inline double value(const Point<2>& p,
                        const unsigned int component = 0) const
    {
        return 1;
    }

    void vector_value(const Point<2> &p,
                      Vector<double>& values) const
    {values.resize(1,false);
    values[0] = value(p);
    }
};

int main()
{
    cout << "* expand some interesting functions on the square into their wavelet series..." << endl;

    Function<2>* f  = new Plateau();
    Function<2>* f2 = new Jump_and_Singularity();
    Function<2>* f3 = new My_Constant();

    int res(8);
#if 1
    Grid<2> grid(Point<2>(0.0, 0.0), Point<2>(1.0, 1.0), 1<<5, 1<<5);
    SampledMapping<2> mapping(grid, *f);
    ofstream file("plateau_function.m");
    mapping.matlab_output(file);
    file.close();

    Grid<2> grid2(Point<2>(0.0, 0.0), Point<2>(1.0, 1.0), 1<<5, 1<<5);
    SampledMapping<2> mapping2(grid2, *f2);
    ofstream file2("Jump_and_Singularity_function.m");
    mapping2.matlab_output(file2);
    file2.close();

    Grid<2> grid3(Point<2>(0.0, 0.0), Point<2>(1.0, 1.0), 1<<res, 1<<res);
    SampledMapping<2> abbildung(grid3, *f3);
    ofstream datei("Constant_Function.m");
    abbildung.matlab_output(datei);
    datei.close();
#endif

#if 1
    typedef PBasis<3,5> Basis1D; // GANZ WICHTIG !!! zur rekonstruktion mittels dualer wavelets muessen diese regulaer genug sein!! also bitte NICHT 2,2 hier eintragen!!
    //typedef DSBasis<2,4> Basis1D;
    typedef TensorBasis<Basis1D> Basis;
    typedef Basis::Index Index;

    //Basis1D basis1d(1,1); // PBasis, complementary b.c.'s
    Basis1D basis1d(false,false);
    FixedArray1D<Basis1D*,2> bases1d;
    bases1d[0] = bases1d[1] = &basis1d;


#if 0

    cout << "Testing Tensor Basis 1D functionality" << endl;

    Vector<double> values (1);
    values[0]=1;
    Function<1>* f4 = new ConstantFunction<1,double>(values);

    FixedArray1D<Basis1D*,1> eineBasis;
    eineBasis[0] = &basis1d;
    TensorBasis<Basis1D,1> bas(eineBasis);

    const int maxrange2 = 1;
    MultiIndex<int,1> jmax2;
    jmax2[0]=bas.bases()[0]->j0()+maxrange2;

    InfiniteVector<double,TensorBasis<Basis1D,1>::Index> coeffs4;
    bas.expand(f4, false, jmax2, coeffs4);
    cout << "- f4: (approx.) expansion coefficients in the primal basis:" << endl
            << coeffs4;


    //cout << bas.first_generator()<< endl;
    //TensorBasis<Basis1D,1>::Index temp_ind(bas.first_generator());
    
    // PBasis<2,4>
    //coeffs4[temp_ind]=0.353553;
    // PBasis<3,5>
    //coeffs4[temp_ind]=0.25;
    //++temp_ind;
    //coeffs4[temp_ind]=0.25;
    //cout << coeffs4 <<endl;

    
    SampledMapping<1> mapping6(evaluate(bas, coeffs4,false,8));
    ofstream file6("Constant_Function1d.m");
    cout << "writing Constant_Function1d.m" << endl;
    mapping6.matlab_output(file6);
    file6.close();
    delete f4;

    return 0;
    
#endif
cout << "HERE" << endl;
    Basis basis(bases1d);


    const int maxrange = 2;
    MultiIndex<int,2> jmax;
    jmax[0]=basis.bases()[0]->j0()+maxrange;
    jmax[1]=basis.bases()[1]->j0();

    InfiniteVector<double,Index> coeffs;
    basis.expand(f, false, jmax, coeffs);
    cout << "- (approx.) expansion coefficients in the primal basis:" << endl
            << coeffs;

     // /*
    InfiniteVector<double,Index> coeffs2;
    basis.expand(f2, false, jmax, coeffs2);
    cout << "- (approx.) expansion coefficients in the primal basis:" << endl
            << coeffs2;

    InfiniteVector<double,Index> coeffs3;
    basis.expand(f3, false, jmax, coeffs3);
    cout << "- (approx.) expansion coefficients in the primal basis:" << endl
            << coeffs3;

    std::ofstream plotstream;
    plotstream.open("coefficient_plot.m");
cout << "Writing coefficient_plot.m" << endl;
    plot_indices_tbasis(&basis,coeffs, maxrange, plotstream, "jet", true, true);
    plotstream.close();
    // */
    if (f) delete f;
    // /*
    std::ofstream plotstream2;
    plotstream2.open("coefficient_plot2.m");
cout << "Writing coefficient_plot2.m" << endl;
    plot_indices_tbasis(&basis,coeffs2, maxrange, plotstream2, "jet", false, true);
    plotstream2.close();

    // */
    if (f2) delete f2;

    std::ofstream plotstream3;
    plotstream3.open("coefficient_plot3.m");
cout << "Writing coefficient_plot3.m" << endl;
    plot_indices_tbasis(&basis,coeffs3, maxrange, plotstream3, "jet", false, true);
    plotstream3.close();

    if (f3) delete f3;
#endif
#if 1
    // folgende Zeile funktioniert nicht.
    /*
    evaluate(*(basis.bases()[1]),
             Basis1D::Index(3,1,0,basis.bases()[1]),
             true,
             3).values();
    */
    //Grid<2> grid3(Point<2>(0.0, 0.0), Point<2>(1.0, 1.0), 1<<5, 1<<5);
    //evaluate(basis, coeffs,true,3);
    cout << "HERE 2" << endl;
    SampledMapping<2> mapping3(evaluate(basis, coeffs,false,res));
    //SampledMapping<2> s       (evaluate(homBC(), domainList[dom]->get_uD(), true, 6));

    ofstream file3("plateau_function2.m");
    cout << "Writing plateau_function2.m" << endl;
    mapping3.matlab_output(file3);
    file3.close();

    //Grid<2> grid4(Point<2>(0.0, 0.0), Point<2>(1.0, 1.0), 1<<5, 1<<5);
    SampledMapping<2> mapping4(evaluate(basis, coeffs2,false,res));
    ofstream file4("Jump_and_Singularity_function2.m");
    cout << "Writing Jump_and_Singularity_function2.m" << endl;
    mapping4.matlab_output(file4);
    file4.close();
    
    SampledMapping<2> mapping5(evaluate(basis, coeffs3,false,res));
    ofstream file5("Constant_Function2.m");
    cout << "Writing Constant_Function.m" << endl;
    mapping5.matlab_output(file5);
    file5.close();
    
    return 0;
#endif
}

