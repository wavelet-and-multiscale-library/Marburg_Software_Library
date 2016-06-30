/**
 * Example_PlotCoefficients.cpp (version 1.0)
 *
 * Intended to expand a given function into a 
 * wavelet series and plot wavelet coefficients.
 *
 * This software is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
 * either expressed or implied.
 *
 * Contact: AG Numerik, Philipps University Marburg
 * http://www.mathematik.uni-marburg.de/~numerik/
 */

#include <iostream>
#include <utils/function.h>
#include <interval/p_basis.h>
#include <interval/p_expansion.h>
#include <interval/i_indexplot.h>
#include <./TestFunctions.h>

using namespace std;
using namespace MathTL;
using namespace WaveletTL;

/**
 * Example PlotCoefficients: main()
 *
 *
 */
int main()
{
    cout << "Example [PlotCoefficients]" << endl
         << "\nExpanding some functions into their wavelet series ..." << endl;

    /*
     * Setup wavelet basis
     */
    const int d = 2;
    const int dt = 2;

    typedef PBasis<d,dt> Basis;
    typedef Basis::Index Index;
    Basis basis(0,0);  /* basis with no b.c.'s */

    //const int j0 = basis.j0();
    const int jmax = 9;  /* maximal level of wavelet basis */
    basis.set_jmax(jmax);


    /*
     * Hat function
     */
    Function<1>* f1 = new Hat();
    const char* filenameCoefficients1 = "Example_PlotCoefficients_Hat.m";

    cout << "\nExpanding 'Hat' function ..." << endl;

    InfiniteVector<double,Index> coeffs1;
    expand(f1, basis, false, jmax, coeffs1);

    cout << "\nWavelet coefficients in the primal basis:" << endl
         << coeffs1 << endl
         << "\nWrite coefficients of 'Hat' function to file ..." << endl;

    std::ofstream coeff_stream1;
    coeff_stream1.open(filenameCoefficients1);
    coeff_stream1 << "figure;" << endl;
    plot_indices(&basis, coeffs1, jmax, coeff_stream1, "jet", false, true, -8);
    coeff_stream1 << "title('Example [PlotCoefficients]: Hat function (Primbs basis)');" << endl;
    coeff_stream1.close();

    delete f1;

    cout << "Coefficient plot written to:" << endl
         << filenameCoefficients1 << endl;



    /*
     * Cut lin function
     */
    Function<1>* f2 = new CutLin();
    const char* filenameCoefficients2 = "Example_PlotCoefficients_CutLin.m";

    cout << "\nExpanding 'cut lin' function ..." << endl;

    InfiniteVector<double,Index> coeffs2;
    expand(f2, basis, false, jmax, coeffs2);

    cout << "\nWavelet coefficients in the primal basis:" << endl
         << coeffs2 << endl
         << "\nWrite coefficients of 'cut lin' function to file ..." << endl;

    std::ofstream coeff_stream2;
    coeff_stream2.open(filenameCoefficients2);
    coeff_stream2 << "figure;" << endl;
    plot_indices(&basis, coeffs2, jmax, coeff_stream2, "jet", false, true, -8);
    coeff_stream2 << "title('Example [PlotCoefficients]: cut lin function (Primbs basis)');" << endl;
    coeff_stream2.close();

    delete f2;

    cout << "Coefficient plot written to:" << endl
         << filenameCoefficients2 << endl;



    /*
     * Cut sqrt function
     */
    Function<1>* f3 = new CutSqrt();
    const char* filenameCoefficients3 = "Example_PlotCoefficients_CutSqrt.m";

    cout << "\nExpanding 'cut sqrt' function ..." << endl;

    InfiniteVector<double,Index> coeffs3;
    expand(f3, basis, false, jmax, coeffs3);

    cout << "\nWavelet coefficients in the primal basis:" << endl
         << coeffs3 << endl
         << "\nWrite coefficients of 'cut sqrt' function to file ..." << endl;

    std::ofstream coeff_stream3;
    coeff_stream3.open(filenameCoefficients3);
    coeff_stream3 << "figure;" << endl;
    plot_indices(&basis, coeffs3, jmax, coeff_stream3, "jet", false, true, -8);
    coeff_stream3 << "title('Example [PlotCoefficients]: cut sqrt function (Primbs basis)');" << endl;
    coeff_stream3.close();

    delete f3;

    cout << "Coefficient plot written to:" << endl
         << filenameCoefficients3 << endl;



    /*
     * Cut pwr function
     */
    Function<1>* f4 = new CutPwr();
    const char* filenameCoefficients4 = "Example_PlotCoefficients_CutPwr.m";

    cout << "\nExpanding 'cut pwr' function ..." << endl;

    InfiniteVector<double,Index> coeffs4;
    expand(f4, basis, false, jmax, coeffs4);

    cout << "\nWavelet coefficients in the primal basis:" << endl
         << coeffs4 << endl
         << "\nWrite coefficients of 'cut pwr' function to file ..." << endl;

    std::ofstream coeff_stream4;
    coeff_stream4.open(filenameCoefficients4);
    coeff_stream4 << "figure;" << endl;
    plot_indices(&basis, coeffs4, jmax, coeff_stream4, "jet", false, true, -8);
    coeff_stream4 << "title('Example [PlotCoefficients]: cut pwr function (Primbs basis)');" << endl;
    coeff_stream4.close();

    delete f4;

    cout << "Coefficient plot written to:" << endl
         << filenameCoefficients4 << endl;



    /*
     * Kink function
     */
    Function<1>* f5 = new Kink(0.5);
    const char* filenameCoefficients5 = "Example_PlotCoefficients_Kink.m";

    cout << "\nExpanding 'Kink' function ..." << endl;

    InfiniteVector<double,Index> coeffs5;
    expand(f5, basis, false, jmax, coeffs5);

    cout << "\nWavelet coefficients in the primal basis:" << endl
         << coeffs5 << endl
         << "\nWrite coefficients of 'Kink' function to file ..." << endl;

    std::ofstream coeff_stream5;
    coeff_stream5.open(filenameCoefficients5);
    coeff_stream5 << "figure;" << endl;
    plot_indices(&basis, coeffs5, jmax, coeff_stream5, "jet", false, true, -8);
    coeff_stream5 << "title('Example [PlotCoefficients]: Kink function (Primbs basis)');" << endl;
    coeff_stream5.close();

    delete f5;

    cout << "Coefficient plot written to:" << endl
         << filenameCoefficients5 << endl;



    /**
     * Contact: AG Numerik, Philipps University Marburg
     * http://www.mathematik.uni-marburg.de/~numerik/
     */
    cout << "\nContact: AG Numerik, Philipps University Marburg" << endl
         << "http://www.mathematik.uni-marburg.de/~numerik/" << endl
         << "\nEnd of Example [PlotCoefficients]." << endl;

    return EXIT_SUCCESS;
}
