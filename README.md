# Marburg_Software_Library
The Marburg Software Library is a collection of several wavelet based c++ algorithms
 to solve partial differential equations on one-dimensional intervals [a,b] as well as on two-dimensional cubes and L-shaped domains. The Poisson equation in 1D and 2D serves as the prime example but also other PDE's, like Sturm- and Helmholtz boundary value problems (BVP), can be solved. All implemented algorithms have in commen the adaptive approach and therefore are able to solve BVP's wihout any a priori information on the solution.  
As a good starting point we recommend to take a look at the Examples folder which consists of implementations of the algorithms proposed in [CDD1] and [CDD2] in one dimension. For the [CDD1] algorithm there is also available a graphical user interface (Examples/SturmBVP_Solver_v1.0) where you can specify a Sturm BVP , the underlying wavelet bases and a solution tolerance.
In the folder WaveletTL/tests the algorithms  
