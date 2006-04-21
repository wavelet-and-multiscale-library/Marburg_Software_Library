
#include <iostream>
#include <frame_index.h>
#include <time.h> 

#include <interval/ds_basis.h>
#include <interval/p_basis.h>
#include <cube/cube_basis.h>

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
  cout << "Testing class FrameIndex..." << endl;
  //typedef DSBasis<2,2> Basis1D;
  typedef PBasis<2,2> Basis1D;
  typedef CubeBasis<Basis1D> Basis;
  typedef AggregatedFrame<Basis1D,2,2> Frame2D;
  typedef Basis::Index Index;
  Basis basis;
  
  FrameIndex<Basis1D,2,2> ind;

  //##############################  
  Matrix<double> A(DIM,DIM);
  A(0,0) = 2.0;
  A(1,1) = 1.0;
  Point<2> b;
  b[0] = -1.0;
  b[1] = -1.0;
  AffineLinearMapping<2> affineP(A,b);
  //##############################

  //##############################
  LinearBezierMapping bezierP(Point<2>(-1.,-1),Point<2>(-1.,1),
			      Point<2>(1.,-1.), Point<2>(1,0));
  //##############################
  Array1D<Chart<DIM,DIM>* > charts(2);
  charts[0] = &affineP;
  charts[1] = &bezierP;

  SymmetricMatrix<bool> adj(2);
  adj(0,0) = 1;
  adj(1,1) = 1;
  adj(1,0) = 1;

  //to specify primal boundary the conditions
  Array1D<FixedArray1D<int,2*DIM> > bc(2);

  //primal boundary conditions for first patch: all Dirichlet
  FixedArray1D<int,2*DIM> bound_1;
  bound_1[0] = 1;
  bound_1[1] = 1;
  bound_1[2] = 1;
  bound_1[3] = 1;

  bc[0] = bound_1;

  //primal boundary conditions for second patch: all Dirichlet
  FixedArray1D<int,2*DIM> bound_2;
  bound_2[0] = 1;
  bound_2[1] = 1;
  bound_2[2] = 1;
  bound_2[3] = 1;

  bc[1] = bound_2;

//to specify primal boundary the conditions
  Array1D<FixedArray1D<int,2*DIM> > bcT(2);

  //dual boundary conditions for first patch
  FixedArray1D<int,2*DIM> bound_3;
  bound_3[0] = 0;
  bound_3[1] = 0;
  bound_3[2] = 0;
  bound_3[3] = 0;

  bcT[0] = bound_3;

  //dual boundary conditions for second patch
  FixedArray1D<int,2*DIM> bound_4;
  bound_4[0] = 0;
  bound_4[1] = 0;
  bound_4[2] = 0;
  bound_4[3] = 0;
 
  bcT[1] = bound_4;

  Atlas<DIM,DIM> Lshaped(charts,adj);  
  cout << Lshaped << endl;

  MultiIndex<int, 2> e2;
  e2[0] = 0;
  e2[1] = 0;
  MultiIndex<int, 2> k2;
  k2[0] = 1;
  k2[1] = 1;

  unsigned int p2 = 0;
  int j2 = 3;



  //a frame has finally been constructed
  //AggregatedFrame<Basis1D, DIM, DIM> frame(&Lshaped, bc, bcT);
  AggregatedFrame<Basis1D, DIM, DIM> frame(&Lshaped, bc);

  FrameIndex<Basis1D,2,2> ind1(&frame,j2,e2,p2,k2);
  cout << "one index:" << endl
       << ind1 << endl;

//   FrameIndex<Basis1D,2,2> ind2(&frame,
// 			       CubeIndex<Basis1D,2>(&basis),
// 			       0);
//   cout << "another index:" << endl
//        << ind2 << endl;

//   cout << "testing < operator "
//        << (ind1 < ind2) << endl;

//   cout << "testing == operator "
//        << (ind1 == ind2) << endl;

//   cout << "testing != operator "
//        << (ind1 != ind2) << endl;

//   cout << "testing <= operator "
//        << (ind1 <= ind2) << endl;

#if 0
  for (int i = 0; i < 1000; i++) 
    {
      cout << ind1 << endl;
      ++ind1;
    }
#endif

  cout << "##################################################" << endl;
  int i = 0;
  for (FrameIndex<Basis1D,2,2> lambda = FrameTL::first_generator<Basis1D,2,2,Frame2D>(&frame, frame.j0());
       lambda <= FrameTL::last_wavelet<Basis1D,2,2,Frame2D>(&frame, frame.j0()+4); ++lambda) {
//     if (!((i == 97) || (i == 98))) {
//       ++i;
//       continue;
//     }
    
    FrameIndex<Basis1D,2,2> mu(lambda);
    cout << lambda << " number  = " << mu.number() << endl;

    if (! (FrameIndex<Basis1D,2,2>(&frame, i) == lambda))
      abort();
    cout << "-----------------------------------------" << endl;
    i++;
  }
  
  for (int j = frame.j0(); j < frame.j0()+5; j++ ) {
    if (j == frame.j0()) {
      cout << "number of first generator on coarstest level = " << first_generator_num<Basis1D,2,2,Frame2D>(&frame) << endl;
      cout << "number of last generator on coarstest level = " << last_generator_num<Basis1D,2,2,Frame2D>(&frame) << endl;
    }
    cout << "number of first wavelet on level " << j << " = " << first_wavelet_num<Basis1D,2,2,Frame2D>(&frame, j) << endl;
    cout << "number of last wavelet on level " << j << " = " << last_wavelet_num<Basis1D,2,2,Frame2D>(&frame, j) << endl;    
  }

  //++++++++++++++++++++++++++++++++++ 1D example ++++++++++++++++++++++++++++++++++
  //##############################  
  Matrix<double> A1D(1,1);
  A1D(0,0) = 1.5;
  Point<1> b1D;
  b1D[0] = -1.0;
  AffineLinearMapping<1> affineP1D(A1D,b1D);
  //##############################
  Array1D<Chart<1,1>* > charts1D(2);
  charts1D[0] = &affineP1D;
  charts1D[1] = &affineP1D;

  SymmetricMatrix<bool> adj1D(2);
  adj1D(0,0) = 1;
  adj1D(1,1) = 1;
  adj1D(1,0) = 1;
  
  //to specify primal boundary the conditions
  Array1D<FixedArray1D<int,2*1> > bc_1D(2);

  //primal boundary conditions for first patch: all Dirichlet
  FixedArray1D<int,2*1> bound_1_1D;
  bound_1_1D[0] = 1;
  bound_1_1D[1] = 1;

  bc_1D[0] = bound_1_1D;

  //primal boundary conditions for second patch: all Dirichlet
  FixedArray1D<int,2*1> bound_2_1D;
  bound_2_1D[0] = 1;
  bound_2_1D[1] = 1;

  bc_1D[1] = bound_2_1D;

  //to specify primal boundary the conditions
  Array1D<FixedArray1D<int,2*1> > bcT_1D(2);

  //dual boundary conditions for first patch
  FixedArray1D<int,2*1> bound_1T_1D;
  bound_1T_1D[0] = 0;
  bound_1T_1D[1] = 0;

  bcT_1D[0] = bound_1T_1D;

  //dual boundary conditions for second patch
  FixedArray1D<int,2*1> bound_2T_1D;
  bound_2T_1D[0] = 0;
  bound_2T_1D[1] = 0;

  bcT_1D[1] = bound_2T_1D;

  Atlas<1,1> interv(charts1D,adj1D);  
  cout << interv << endl;


  typedef CubeBasis<Basis1D,1> RefBasis1D;
  RefBasis1D basis1;


  //a frame has finally been constructed
  //AggregatedFrame<Basis1D, 1, 1> frame1D(&interv, bc_1D, bcT_1D);
  AggregatedFrame<Basis1D, 1, 1> frame1D(&interv, bc_1D);


  MultiIndex<unsigned int, 1> e;
  e[0] = 0;
  MultiIndex<int, 1> k;
  k[0] = 1;

  unsigned int p = 0;
  int j = 3;

//   FrameIndex<Basis1D,1,1> ind1D(&frame1D,j,e,p,k);

//   for (int i = 0; i < 1; i++) 
//     {
//       cout << ind1D << endl;
//       ++ind1D;
//     }

 //  cout << "################################" << endl;
//   cout <<  FrameTL::last_wavelet<Basis1D,2,2,Frame2D>(&frame, 5) << endl;
  
//   for (FrameIndex<Basis1D,2,2> ind = FrameTL::first_generator<Basis1D,2,2,frame>(&frame, 5);
//        ind <= FrameTL::last_wavelet<Basis1D,2,2,frame>(&frame, 5); ++ind) 
//     {
//       cout << ind << endl;
//     }

//   CubeIndex<Basis1D,2> cube_index1(&basis);
//   cout << cube_index1 << endl;
//   ++cube_index1;
//   cout << cube_index1 << endl;

//   FrameIndex<Basis1D,2,2> ind2(cube_index1,1,4);
//   cout << ind2 << endl;

//   ++cube_index1;
//   FrameIndex<Basis1D,2> ind3(cube_index1,1,4);

//   cout << "test for equality returns "
//        << (ind2 == ind3) << endl;



  return 0;
}
