#include <iostream>
#include <fstream>
#include <sstream>

#include <algebra/infinite_vector.h>
#include <utils/array1d.h>
#include <geometry/sampled_mapping.h>

#include <interval/a_basis.h>
#include <interval/a_evaluate.h>
#include <interval/i_indexplot.h>

using namespace std;
using namespace WaveletTL;

int main()
{
  cout << "Testing the A bases ..." << endl;

  const int n = 4;
  typedef ABasis<n> Basis;
  const unsigned int plot_resolution = 10;
  int i;

  ABasis<n> basis;
  ABasis<n>::Index lambda(&basis);
  SampledMapping<1> map;
  std::ofstream fs;

  cout << "n = " << n << endl;

  cout << "Point evaluation of index " << lambda << " at 0.5: ";
  cout << evaluate(basis, 0, lambda, 0.5) << endl;

  cout << "Generators:" << endl;
  for (i = 0; i < n; i++)
    cout << basis.get_function(Basis::Index(0, E_GENERATOR, 0, i, &basis)) << endl;

  cout << "Wavelets:" << endl;
  for (i = 0; i < n; i++)
    cout << basis.get_function(Basis::Index(0, E_WAVELET, 0, i, &basis)) << endl;

  cout << "- evaluating primal A functions: writing file 'agenerator-i.dat'" << endl;
  for (i = 0; i < n; i++) {
    ostringstream filename;
    filename << "agenerator-" << i << ".dat";
    fs.open(filename.str().c_str());
    lambda = Basis::Index(0, E_GENERATOR, 0, i, &basis);
    map = WaveletTL::evaluate(basis, lambda, plot_resolution);
    map.gnuplot_output(fs);
    fs.close();
  }

  cout << "- evaluating primal A functions: writing file 'awavelet-i.dat'" << endl;
  for (i = 0; i < n; i++) {
    ostringstream filename;
    filename << "awavelet-" << i << ".dat";
    fs.open(filename.str().c_str());
    lambda = Basis::Index(0, E_WAVELET, 0, i, &basis);
    map = evaluate(basis, lambda, plot_resolution);
    map.gnuplot_output(fs);
    fs.close();
  }

  cout << endl << "Testing the A basis for n = 1 (should be the Haar basis)" << endl;
  ABasis<1> haar;
  cout << haar.get_function(ABasis<1>::Index(0, E_GENERATOR, 0, 0, &haar)) << endl; // Haar generator
  cout << haar.get_function(ABasis<1>::Index(0, E_WAVELET, 0, 0, &haar)) << endl; // Haar wavelet
  // plot Haar generator
  cout << "Writing file 'haar-generator.dat'" << endl;
  fs.open("haar-generator.dat");
  map = WaveletTL::evaluate(haar, ABasis<1>::Index(0, E_GENERATOR, 0, 0, &haar), plot_resolution);
  map.gnuplot_output(fs);
  fs.close();
  // plot Haar wavelet
  cout << "Writing file 'haar-wavelet.dat'" << endl;
  fs.open("haar-wavelet.dat");
  map = WaveletTL::evaluate(haar, ABasis<1>::Index(0, E_WAVELET, 0, 0, &haar), plot_resolution);
  map.gnuplot_output(fs);
  fs.close();

  return 0;
}
