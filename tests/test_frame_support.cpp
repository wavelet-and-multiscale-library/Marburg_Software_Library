
#include <iostream>
#include <frame_index.h>
#include <time.h> 

#include <aggregated_frame.h>
#include <interval/ds_basis.h>
#include <cube/cube_basis.h>
#include <frame_support.h>

using std::cout;
using std::endl;


using FrameTL::FrameIndex;
using FrameTL::AggregatedFrame;
using WaveletTL::CubeIndex;

using namespace FrameTL;
using namespace std;
using namespace WaveletTL;

int main()
{

  const int DIM = 2;
  cout << "Testing routines for supports..." << endl;
  typedef DSBasis<2,2> Basis1D;
  typedef CubeBasis<Basis1D> Basis;
  typedef AggregatedFrame<Basis1D,2,2> Frame2D;
  typedef Basis::Index Index;
  Basis basis;
  
  FrameIndex<Basis1D,2,2> ind;


  return 0;
}
