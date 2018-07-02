# Marburg_Software_Library
The Marburg Software Library is a collection of several wavelet based C++ algorithms
 to solve partial differential equations on one-dimensional intervals [a,b] as well as on two-dimensional cubes and L-shaped domains. While the Poisson equation in 1D and 2D serves as the prime example, also other PDE's, like Sturm- and Helmholtz boundary value problems (BVP), can be solved. All implemented algorithms realize adaptive wavelet approaches based on aposteriori error estimators and associated refinement strategies. Therefore they allow to solve the BVP's wihout any a priori information on the solution. The library contains algorithms based on classical wavelet bases as well as generalizations to wavelet and quarklet frames.
The core ingredient for all our wavelet bases/frames is the boundary adapted wavelet bases on the unit interval as constructed in [Pr] and [DS]. 
To use this software library the user has to select his favourite problem class, to choose the underlying wavelet bases
together with the order and the amount of vanishing moments, to choose a solution tolerance and to enter a right hand side.
As a good starting point we recommend to take a look at the Examples folder which contains implementations of the algorithms proposed in [CDD1] and [CDD2] in one dimension. There is also available a graphical user interface for the [CDD1] algorithm (Examples/SturmBVP_Solver_v1.0) where the above mentioned input selection is quite convenient.
For the other implemented algorithms, the task for the user stays the same, however the management requires some basic C++ knowledge. During the execution of every algorithm, matlab plots of the solution and the solution indices will be produced and stored in the respective tests/Matlab_outputs folder.  
The in practice very efficient frame schemes based on domain decompositions, e.g., the multiplicative Schwarz algorithm, can be found in the folder FrameTL/tests (see [We], [DFRSW], [DFPRW], [DFR] for the theoretical background).  
In WaveletTL/tests there are implementations of the quarklet frame schemes as outlined in [DKR].

We give a short overview of the available schemes and where to find them:


Wavelet schemes in 1D:  
Examples/Example_CDD1/Example_CDD1.cpp  
Examples/Example_CDD2/Example_CDD2.cpp  
Examples/SturmBVP_Solver_v1.0 (GUI)

Quarklet schemes in 1D:
WaveletTL/tests/(test_sturm_bvp.cpp, test_periodic_laplacian.cpp, test_periodic_gramian.cpp)

Quarklet schemes in 2D:
WaveletTL/tests/(test_ldomain_frame.cpp, test_ringdomain_frame.cpp, test_slitdomain_frame.cpp)


Wavelet frame schemes in 1 and 2D:
FrameTL/tests/(test_additive_Schwarz.cpp, test_multiplicative_Schwarz.cpp, test_richardson.cpp, test_steepest_descent.cpp)  



Sources:

[CDD1] A. Cohen, W. Dahmen and R. DeVore, *Adaptive wavelet methods for elliptic operator eqautions - Convergence rates*, Math. Comput. **70** (2001), no. 233, 27-75

[CDD2] A. Cohen, W. Dahmen and R. DeVore, *Adaptive wavelet methods II: Beyond the elliptic case*, Found. Comput. Math. **2** (2002), no. 3, 203-245

[DFRSW] S. Dahlke, M. Fornasier, T. Raasch, R. Stevenson and M. Werner, *Adaptive frame methods for elliptic operator equations: The steepest descent approach*, IMA J. Numer. Anal. **27** (2007), no. 4, 717-740

[DFPRW] S. Dahlke, M. Fornasier, M. Primbs, T. Raasch and M. Werner, *Nonlinear and adaptive frame approximation schemes for elliptic PDEs: Theory and numerical experiments*, Numer. Methods Partial Differ. Equations **6** (2009), no. 6, 1366-1401

[DFR] S. Dahlke, M. Fornasier and T. Raasch, *Adaptive frame methods for elliptic operator equations*, Adv. Comput. Math. **27** (2007), no.1, 27-63 

[DKR] S. Dahlke, P. Keding and T. Raasch, *Quarkonial frames with compression properties*, Calcolo **27** (2017), no.3, 823-855

[DS]  W. Dahmen and R. Schneider, *Composite wavelet bases for operator equations*, Math. Comput. **68** (1999), 1533-1567 

[Pr]  M. Primbs, *New stable biorthogonal spline-wavelets on the interval*, **57** (2010), no.1, 121-162 

[We]  M. Werner, *Adaptive Wavelet Frame Domain Decomposition Methods for Elliptic Operators*, Ph.D. thesis, Philipps-Universität Marburg, 2009 
 
