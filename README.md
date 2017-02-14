# Marburg_Software_Library
The Marburg Software Library is a collection of several wavelet based c++ algorithms
 to solve partial differential equations on one-dimensional intervals [a,b] as well as on two-dimensional cubes and L-shaped domains. While the Poisson equation in 1D and 2D serves as the prime example, also other PDE's, like Sturm- and Helmholtz boundary value problems (BVP), can be solved. All implemented algorithms have two things in common:
first, the adaptive approach and therefore the ability to solve BVP's wihout any a priori information on the solution
and second, the boundary adapted wavelet bases on the unit interval as constructed in [Pr] and [DS] as the core of more involved constructions of wavelet/quarklet bases and frames on one and two-dimensional domains. 
As a good starting point we recommend to take a look at the Examples folder which consists of implementations of the algorithms proposed in [CDD1] and [CDD2] in one dimension. There is also available a graphical user interface for the [CDD1] algorithm (Examples/SturmBVP_Solver_v1.0) where you can specify a Sturm BVP , the underlying wavelet bases together with the order and the amount of vanishing moments and a solution tolerance.
While for the other implemented algorithms it takes some basic c++ knowledge to specify the latter, these are essentially all the configurations the user needs to handle for the algorithms to work properly.  
The in parctice very efficient frame schemes based on domain decompositions, e.g., the multiplicative Schwarz algorithm, can be found in the folder FrameTL/tests (see [We], [DFRSW], [DFPRW], [DFR] for the theoretical background).  
In WaveletTL/tests there are implementations of the quarklet frame schemes as outlined in [DKR].
We give a short overview of the available schemes and where to find them:


Wavelet schemes in 1D:  
Examples/Example_CDD1/Example_CDD1.cpp  
Examples/Example_CDD2/Example_CDD2.cpp  
Examples/SturmBVP_Solver_v1.0 (GUI)

Quarklet schemes in 1D:
WaveletTL/tests/(test_sturm_bvp.cpp, test_periodic_laplacian.cpp, test_periodic_gramian.cpp)

Wavelet frame schemes in 1 and 2D:
FrameTL/tests/(test_additive_Schwarz.cpp, test_multiplicative_Schwarz.cpp, test_richardson.cpp, test_steepest_descent.cpp)  

...



Sources:

[CDD1] ...  
[CDD2] ...  
[DFRSW] ...  
[DFPRW] ...  
[DFR] ... 
[DKR] ...  
[DS] ...  
[Pr] ...  
[We] ...  
 
