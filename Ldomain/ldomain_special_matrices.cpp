// implementation for ldomain_basis.h, template specialization to SplineBasis<d,dT,DS_construction>,
// application of diverse matrices

#include <cmath>
#include <time.h>
#include <iostream>
#include <algorithm>
#include <utils/map_tools.h>

using std::cout;
using std::endl;

namespace WaveletTL
{
  template <int d, int dT>
  void
  LDomainBasis<SplineBasis<d,dT,DS_construction> >::apply_Mj0
  (const int j,
   const std::map<size_type,double>& x, 
   std::map<size_type,double>& y) const
  {
    typedef std::map<size_type,double> V;
    
    const unsigned int Deltaj   = basis1d().Deltasize(j);
    const unsigned int Deltajp1 = basis1d().Deltasize(j+1);
    
    basis1d().Mj0_.set_level(j);

    for (V::const_iterator itx(x.begin()); itx != x.end(); ++itx) {
      // determine patch number
      const unsigned int patch =
	itx->first < 3*(Deltaj-2)*(Deltaj-2)
	? itx->first / ((Deltaj-2)*(Deltaj-2))
	: 3+(itx->first-3*(Deltaj-2)*(Deltaj-2))/(Deltaj-2);
      
      switch(patch) {
      case 0: {
	// psi_lambda decomposes into generators on patch 0 alone
	// apply kron(M#,M#), cf. KroneckerMatrix::apply()
 	// 1. get column of second factor M#
 	V z1, z2, z3;
	const size_type block_index = itx->first % (Deltaj-2);
 	z1.insert(std::pair<size_type,double>(block_index, 1.0));
 	basis1d().Mj0_.apply_central_block(z1, z2);
	
 	// 2. get column of first factor M#
 	z1.clear();
	const size_type block_nr    = itx->first / (Deltaj-2);
 	z1.insert(std::pair<size_type,double>(block_nr, 1.0));
 	basis1d().Mj0_.apply_central_block(z1, z3);
	
 	// 3. combine results
 	for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
 	  for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
	    y[it3->first*(Deltajp1-2)+it2->first]
	      += itx->second * it2->second * it3->second;
      }
	break;
      case 1: {
	// psi_lambda decomposes into generators on patch 1 alone
	// apply kron(M#,M#)
	const size_type block_nr    = (itx->first-(Deltaj-2)*(Deltaj-2)) / (Deltaj-2);
	const size_type block_index = (itx->first-(Deltaj-2)*(Deltaj-2)) % (Deltaj-2);
	
 	V z1, z2, z3;
 	z1.insert(std::pair<size_type,double>(block_index, 1.0));
 	basis1d().Mj0_.apply_central_block(z1, z2);

 	z1.clear();
 	z1.insert(std::pair<size_type,double>(block_nr, 1.0));
 	basis1d().Mj0_.apply_central_block(z1, z3);

 	for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
 	  for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
	    y[(Deltajp1-2)*(Deltajp1-2)+it3->first*(Deltajp1-2)+it2->first]
	      += itx->second * it2->second * it3->second;
      }
      case 2: {
	// psi_lambda decomposes into generators on patch 2 alone
	// apply kron(M#,M#)
	const size_type block_nr    = (itx->first-2*(Deltaj-2)*(Deltaj-2)) / (Deltaj-2);
	const size_type block_index = (itx->first-2*(Deltaj-2)*(Deltaj-2)) % (Deltaj-2);

 	V z1, z2, z3;
 	z1.insert(std::pair<size_type,double>(block_index, 1.0));
 	basis1d().Mj0_.apply_central_block(z1, z2);

 	z1.clear();
 	z1.insert(std::pair<size_type,double>(block_nr, 1.0));
 	basis1d().Mj0_.apply_central_block(z1, z3);

 	for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
 	  for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
	    y[2*(Deltajp1-2)*(Deltajp1-2)+it3->first*(Deltajp1-2)+it2->first]
	      += itx->second * it2->second * it3->second;
      }
	break;
      case 3: {
	// psi_lambda decomposes into generators on patches 0,1 and 3

	// contribution from patch 0
	//
	// apply kron(M#,ML):
 	// 1. get column of second factor ML
	V z1, z2, z3;
	z1.insert(std::pair<size_type,double>(0, 1.0));
	basis1d().Mj0_.apply(z1, z2); // we have to neglect the first entry of z2 later

 	// 2. get column of first factor M#
 	z1.clear();
 	const size_type block_nr = itx->first-3*(Deltaj-2)*(Deltaj-2);
 	z1.insert(std::pair<size_type,double>(block_nr, 1.0));
 	basis1d().Mj0_.apply_central_block(z1, z3);

 	// 3. combine results
 	for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
 	  for (V::const_iterator it2(++(z2.begin())); it2 != z2.end(); ++it2)
	    y[it3->first*(Deltajp1-2)+it2->first-1] // note the "-1"
	      += itx->second * it2->second * it3->second * M_SQRT1_2;

	// contribution from patch 1
	//
	// apply kron(M#,MR):
 	// 1. get column of second factor MR
	z1.clear();
	z1.insert(std::pair<size_type,double>(Deltaj-1, 1.0));
	z2.clear();
	basis1d().Mj0_.apply(z1, z2); // we have to neglect the last entry of z2 later

	// 2. get column of first factor M#: this has already been done above
	
	// 3. combine results
 	for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
 	  for (V::reverse_iterator it2(++(z2.rbegin())); it2 != z2.rend(); ++it2)
	    y[(Deltajp1-2)*(Deltajp1-2)+it3->first*(Deltajp1-2)+it2->first-1] // note the "-1"
	      += itx->second * it2->second * it3->second * M_SQRT1_2;

	// contribution from patch 3
	//
	// apply kron(M#,Mtopleft):
	// 1. get column of first factor M#: this has already been done above

	// 2. get top left entry of Mj0
	const double Mtopleft = basis1d().Mj0_.get_entry(0,0);
	
	// 3. combine results
	for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
	  y[3*(Deltajp1-2)*(Deltajp1-2)+it3->first]
	    += itx->second * it3->second * Mtopleft;
      }
	break;
      case 4: {
	// psi_lambda decomposes into generators on patches 1,2 and 4

	// contribution from patch 1
	//
	// apply kron(MR,M#):
 	// 1. get column of second factor M#
	V z1, z2, z3;
	const size_type block_index = itx->first-(3*(Deltaj-2)+1)*(Deltaj-2);
	z1.insert(std::pair<size_type,double>(block_index, 1.0));
	basis1d().Mj0_.apply_central_block(z1, z2);

 	// 2. get column of first factor MR
	z1.clear();
	z1.insert(std::pair<size_type,double>(Deltaj-1, 1.0));
	basis1d().Mj0_.apply(z1, z3); // we have to neglect the last entry of z3 later
	
	// 3. combine results
 	for (V::reverse_iterator it3(++(z3.rbegin())); it3 != z3.rend(); ++it3)
 	  for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
	    y[(Deltajp1-2)*(Deltajp1-2)+(it3->first-1)*(Deltajp1-2)+it2->first] // note the "-1"
	      += itx->second * it2->second * it3->second * M_SQRT1_2;

	// contribution from patch 2
	//
	// apply kron(ML,M#):
	// 1. get column of second factor M#: this has already been done above

	// 2. get column of first factor ML
	z1.clear();
	z1.insert(std::pair<size_type,double>(0, 1.0));
	z3.clear();
	basis1d().Mj0_.apply(z1, z3); // we have to neglect the first entry of z3 later

	// 3. combine results
 	for (V::const_iterator it3(++(z3.begin())); it3 != z3.end(); ++it3)
 	  for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
	    y[2*(Deltajp1-2)*(Deltajp1-2)+(it3->first-1)*(Deltajp1-2)+it2->first] // note the "-1"
	      += itx->second * it2->second * it3->second * M_SQRT1_2;

	// contribution from patch 4
	//
	// apply kron(Mtopleft,M#):
	// 1. get column of second factor M#: this has already been done above

	// 2. get top left entry of Mj0
	const double Mtopleft = basis1d().Mj0_.get_entry(0,0);
	
	// 3. combine results
	for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
	  y[(3*(Deltajp1-2)+1)*(Deltajp1-2)+it2->first]
	    += itx->second * it2->second * Mtopleft;
      }
	break;
      default:
	break;
      }
    }
  }
  
  template <int d, int dT>
  void
  LDomainBasis<SplineBasis<d,dT,DS_construction> >::apply_Mj0T_transposed
  (const int j,
   const std::map<size_type,double>& x, 
   std::map<size_type,double>& y) const
  {
    typedef std::map<size_type,double> V;
    
    const unsigned int Deltaj   = basis1d().Deltasize(j);
    const unsigned int Deltajp1 = basis1d().Deltasize(j+1);
    
    basis1d().Mj0T_.set_level(j);

    for (V::const_iterator itx(x.begin()); itx != x.end(); ++itx) {
      // determine patch number
      const unsigned int patch =
	itx->first < 3*(Deltajp1-2)*(Deltajp1-2)
	? itx->first / ((Deltajp1-2)*(Deltajp1-2))
	: 3+(itx->first-3*(Deltajp1-2)*(Deltajp1-2))/(Deltajp1-2);
      
      switch(patch) {
      case 0: {
	// apply kron((M#)^T,(M#)^T):
 	// 1. get column of second factor (M#)^T
 	V z1, z2, z3;
	const size_type block_index = itx->first % (Deltajp1-2);
 	z1.insert(std::pair<size_type,double>(block_index, 1.0));
 	basis1d().Mj0T_.apply_central_block_transposed(z1, z2);

 	// 2. get column of first factor (M#)^T
 	z1.clear();
	const size_type block_nr    = itx->first / (Deltajp1-2);
 	z1.insert(std::pair<size_type,double>(block_nr, 1.0));
 	basis1d().Mj0T_.apply_central_block_transposed(z1, z3);
	
	// 3. combine results
 	for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
 	  for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
	    y[it3->first*(Deltaj-2)+it2->first]
	      += itx->second * it2->second * it3->second;

 	// apply kron((M#)^T,(ML)^T)
	// 1. get column of second factor (ML)^T (which is a number)
	const double ML_entry = basis1d().Mj0T_.get_entry(block_index+1, 0);

	// 2. get column of first factor (M#)^T: this has already been done above

	// 3. combine results
	for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
	  y[3*(Deltaj-2)*(Deltaj-2)+it3->first]
	    += itx->second * ML_entry * it3->second * M_SQRT1_2;
      }
	break;
      case 1: {
	// apply kron((M#)^T,(M#)^T):
 	// 1. get column of second factor (M#)^T
 	V z1, z2, z3;
	const size_type block_index = (itx->first-(Deltajp1-2)*(Deltajp1-2)) % (Deltajp1-2);
 	z1.insert(std::pair<size_type,double>(block_index, 1.0));
 	basis1d().Mj0T_.apply_central_block_transposed(z1, z2);

 	// 2. get column of first factor (M#)^T
 	z1.clear();
	const size_type block_nr    = (itx->first-(Deltajp1-2)*(Deltajp1-2)) / (Deltajp1-2);
 	z1.insert(std::pair<size_type,double>(block_nr, 1.0));
 	basis1d().Mj0T_.apply_central_block_transposed(z1, z3);
	
	// 3. combine results
 	for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
 	  for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
	    y[(Deltaj-2)*(Deltaj-2)+it3->first*(Deltaj-2)+it2->first]
	      += itx->second * it2->second * it3->second;

	// apply kron((M#)^T,(MR)^T)
	// 1. get column of second factor (MR)^T (which is a number)
	const double MR_entry = basis1d().Mj0T_.get_entry(block_index+1, Deltaj-1);

	// 2. get column of first factor (M#)^T: this has already been done above
	
	// 3. combine results
	for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
	  y[3*(Deltaj-2)*(Deltaj-2)+it3->first]
	    += itx->second * MR_entry * it3->second * M_SQRT1_2;

	// apply kron((MR)^T,(M#)^T)
	// 1. get column of second factor (M#)^T: this has already been done above

	// 2. get column of first factor (MR)^T (which is a number)
	const double MR_entry2 = basis1d().Mj0T_.get_entry(block_nr+1, Deltaj-1);
	
	// 3. combine results
	for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
	  y[(3*(Deltaj-2)+1)*(Deltaj-2)+it2->first]
	    += itx->second * MR_entry2 * it2->second * M_SQRT1_2;
      }
	break;
      case 2: {
	// apply kron((M#)^T,(M#)^T):
 	// 1. get column of second factor (M#)^T
 	V z1, z2, z3;
	const size_type block_index = (itx->first-2*(Deltajp1-2)*(Deltajp1-2)) % (Deltajp1-2);
 	z1.insert(std::pair<size_type,double>(block_index, 1.0));
 	basis1d().Mj0T_.apply_central_block_transposed(z1, z2);

 	// 2. get column of first factor (M#)
 	z1.clear();
	const size_type block_nr    = (itx->first-2*(Deltajp1-2)*(Deltajp1-2)) / (Deltajp1-2);
 	z1.insert(std::pair<size_type,double>(block_nr, 1.0));
 	basis1d().Mj0T_.apply_central_block_transposed(z1, z3);
	
	// 3. combine results
 	for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
 	  for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
	    y[2*(Deltaj-2)*(Deltaj-2)+it3->first*(Deltaj-2)+it2->first]
	      += itx->second * it2->second * it3->second;

	// apply kron((ML)^T,(M#)^T)
	// 1. get column of second factor (M#)^T: this has already been done above
	
	// 2. get column of first factor (ML)^T (which is a number)
	const double ML_entry = basis1d().Mj0T_.get_entry(block_nr+1, 0);

	// 3. combine results
	for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
	  y[(3*(Deltaj-2)+1)*(Deltaj-2)+it2->first]
	    += itx->second * ML_entry * it2->second * M_SQRT1_2;	
      }
	break;
      case 3: {
	// apply kron(Mtopleft,(M#)^T):
	// 1. get column of second factor (M#)^T
	V z1, z2;
 	const size_type block_nr = itx->first-3*(Deltajp1-2)*(Deltajp1-2);
	z1.insert(std::pair<size_type,double>(block_nr, 1.0));
	basis1d().Mj0T_.apply_central_block_transposed(z1, z2);

	// 2. get top left entry of Mj0T
	const double Mtopleft = basis1d().Mj0T_.get_entry(0,0);

	// 3. combine the results
 	for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
 	  y[3*(Deltaj-2)*(Deltaj-2)+it2->first]
 	    += itx->second * it2->second * Mtopleft;
      }
	break;
      case 4: {
	// apply kron(Mtopleft,(M#)^T):
	V z1, z2;
 	const size_type block_nr = itx->first-(3*(Deltajp1-2)+1)*(Deltajp1-2);
	z1.insert(std::pair<size_type,double>(block_nr, 1.0));
	basis1d().Mj0T_.apply_central_block_transposed(z1, z2);

	// 2. get top left entry of Mj0T
	const double Mtopleft = basis1d().Mj0T_.get_entry(0,0);

	// 3. combine the results
 	for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
 	  y[(3*(Deltaj-2)+1)*(Deltaj-2)+it2->first]
 	    += itx->second * it2->second * Mtopleft;	
      }
	break;
      default:
	break;
      }
    }
  }

  template <int d, int dT>
  void
  LDomainBasis<SplineBasis<d,dT,DS_construction> >::apply_Mj1c_01
  (const int j,
   const std::map<size_type,double>& x, 
   std::map<size_type,double>& y) const
  {
    typedef std::map<size_type,double> V;
    
    const unsigned int Deltaj   = basis1d().Deltasize(j);
    const unsigned int Nablaj   = basis1d().Nablasize(j);
    const unsigned int Deltajp1 = basis1d().Deltasize(j+1);
    
    basis1d().Mj0_.set_level(j);
    basis1d().Mj1c_.set_level(j);

    for (V::const_iterator itx(x.begin()); itx != x.end(); ++itx) {
      // determine patch number
      const unsigned int patch =
 	itx->first < 3*(Deltaj-2)*Nablaj
 	? itx->first / ((Deltaj-2)*Nablaj)
 	: 4;
      
      switch(patch) {
      case 0: {
	// apply kron(M#,N#)
 	// 1. get column of second factor N#
 	V z1, z2, z3;
	const size_type block_index = itx->first % Nablaj;
 	z1.insert(std::pair<size_type,double>(block_index, 1.0));
 	basis1d().Mj1c_.apply(z1, z2);
	
 	// 2. get column of first factor M#
 	z1.clear();
	const size_type block_nr    = itx->first / Nablaj;
 	z1.insert(std::pair<size_type,double>(block_nr, 1.0));
 	basis1d().Mj0_.apply_central_block(z1, z3);
	
 	// 3. combine results
 	for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
 	  for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
	    y[it3->first*(Deltajp1-2)+it2->first]
	      += itx->second * it2->second * it3->second;
      }
	break;
      case 1: {
	// apply kron(M#,N#)
 	// 1. get column of second factor N#
 	V z1, z2, z3;
	const size_type block_index = (itx->first-(Deltaj-2)*Nablaj) % Nablaj;
 	z1.insert(std::pair<size_type,double>(block_index, 1.0));
 	basis1d().Mj1c_.apply(z1, z2);
	
 	// 2. get column of first factor M#
 	z1.clear();
	const size_type block_nr    = (itx->first-(Deltaj-2)*Nablaj) / Nablaj;
 	z1.insert(std::pair<size_type,double>(block_nr, 1.0));
 	basis1d().Mj0_.apply_central_block(z1, z3);
	
 	// 3. combine results
 	for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
 	  for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
	    y[(Deltajp1-2)*(Deltajp1-2)+it3->first*(Deltajp1-2)+it2->first]
	      += itx->second * it2->second * it3->second;
      }
	break;
      case 2: {
	// apply kron(M#,N#)
 	// 1. get column of second factor N#
 	V z1, z2, z3;
	const size_type block_index = (itx->first-2*(Deltaj-2)*Nablaj) % Nablaj;
 	z1.insert(std::pair<size_type,double>(block_index, 1.0));
 	basis1d().Mj1c_.apply(z1, z2);
	
 	// 2. get column of first factor M#
 	z1.clear();
	const size_type block_nr    = (itx->first-2*(Deltaj-2)*Nablaj) / Nablaj;
 	z1.insert(std::pair<size_type,double>(block_nr, 1.0));
 	basis1d().Mj0_.apply_central_block(z1, z3);
	
 	// 3. combine results
 	for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
 	  for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
	    y[2*(Deltajp1-2)*(Deltajp1-2)+it3->first*(Deltajp1-2)+it2->first]
	      += itx->second * it2->second * it3->second;
      }
	break;
      case 4: {
	// contribution from patch 1
	//
	// apply kron(MR,N#):
 	// 1. get column of second factor N#
	V z1, z2, z3;
	const size_type block_index = itx->first-3*(Deltaj-2)*Nablaj;
	z1.insert(std::pair<size_type,double>(block_index, 1.0));
	basis1d().Mj1c_.apply(z1, z2);

 	// 2. get column of first factor MR
	z1.clear();
	z1.insert(std::pair<size_type,double>(Deltaj-1, 1.0));
	basis1d().Mj0_.apply(z1, z3); // we have to neglect the last entry of z3 later
	
	// 3. combine results
 	for (V::reverse_iterator it3(++(z3.rbegin())); it3 != z3.rend(); ++it3)
 	  for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
	    y[(Deltajp1-2)*(Deltajp1-2)+(it3->first-1)*(Deltajp1-2)+it2->first] // note the "-1"
	      += itx->second * it2->second * it3->second * M_SQRT1_2;

 	// contribution from patch 2
 	//
 	// apply kron(ML,N#):
 	// 1. get column of second factor N#: this has already been done above

 	// 2. get column of first factor ML
 	z1.clear();
 	z1.insert(std::pair<size_type,double>(0, 1.0));
 	z3.clear();
 	basis1d().Mj0_.apply(z1, z3); // we have to neglect the first entry of z3 later

 	// 3. combine results
  	for (V::const_iterator it3(++(z3.begin())); it3 != z3.end(); ++it3)
  	  for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
 	    y[2*(Deltajp1-2)*(Deltajp1-2)+(it3->first-1)*(Deltajp1-2)+it2->first] // note the "-1"
 	      += itx->second * it2->second * it3->second * M_SQRT1_2;

 	// contribution from patch 4
 	//
 	// apply kron(Mtopleft,N#):
 	// 1. get column of second factor M#: this has already been done above

 	// 2. get top left entry of Mj0
 	const double Mtopleft = basis1d().Mj0_.get_entry(0,0);
	
 	// 3. combine results
 	for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
 	  y[(3*(Deltajp1-2)+1)*(Deltajp1-2)+it2->first]
 	    += itx->second * it2->second * Mtopleft;
      }
	break;
      default:
	break;
      }
    }
  }
  
  template <int d, int dT>
  void
  LDomainBasis<SplineBasis<d,dT,DS_construction> >::apply_Mj1c_10
  (const int j,
   const std::map<size_type,double>& x, 
   std::map<size_type,double>& y) const
  {
    typedef std::map<size_type,double> V;
    
    const unsigned int Deltaj   = basis1d().Deltasize(j);
    const unsigned int Nablaj   = basis1d().Nablasize(j);
    const unsigned int Deltajp1 = basis1d().Deltasize(j+1);
    
    basis1d().Mj0_.set_level(j);
    basis1d().Mj1c_.set_level(j);

    for (V::const_iterator itx(x.begin()); itx != x.end(); ++itx) {
      // determine patch number
      const unsigned int patch =
 	itx->first < 3*(Deltaj-2)*Nablaj
 	? itx->first / ((Deltaj-2)*Nablaj)
 	: 3;

      switch(patch) {
      case 0: {
	// apply kron(N#,M#)
 	// 1. get column of second factor M#
 	V z1, z2, z3;
	const size_type block_index = itx->first % (Deltaj-2);
 	z1.insert(std::pair<size_type,double>(block_index, 1.0));
	basis1d().Mj0_.apply_central_block(z1, z2);
	
 	// 2. get column of first factor N#
 	z1.clear();
	const size_type block_nr    = itx->first / (Deltaj-2);
 	z1.insert(std::pair<size_type,double>(block_nr, 1.0));
 	basis1d().Mj1c_.apply(z1, z3);
	
 	// 3. combine results
 	for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
 	  for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
	    y[it3->first*(Deltajp1-2)+it2->first]
	      += itx->second * it2->second * it3->second;
      }
	break;
      case 1: {
	// apply kron(N#,M#)
 	// 1. get column of second factor M#
 	V z1, z2, z3;
	const size_type block_index = (itx->first-(Deltaj-2)*Nablaj) % (Deltaj-2);
 	z1.insert(std::pair<size_type,double>(block_index, 1.0));
 	basis1d().Mj0_.apply_central_block(z1, z2);
	
 	// 2. get column of first factor N#
 	z1.clear();
	const size_type block_nr    = (itx->first-(Deltaj-2)*Nablaj) / (Deltaj-2);
 	z1.insert(std::pair<size_type,double>(block_nr, 1.0));
 	basis1d().Mj1c_.apply(z1, z3);
	
 	// 3. combine results
 	for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
 	  for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
	    y[(Deltajp1-2)*(Deltajp1-2)+it3->first*(Deltajp1-2)+it2->first]
	      += itx->second * it2->second * it3->second;
      }
	break;
      case 2: {
	// apply kron(N#,M#)
 	// 1. get column of second factor M#
 	V z1, z2, z3;
	const size_type block_index = (itx->first-2*(Deltaj-2)*Nablaj) % (Deltaj-2);
 	z1.insert(std::pair<size_type,double>(block_index, 1.0));
 	basis1d().Mj0_.apply_central_block(z1, z2);
	
 	// 2. get column of first factor N#
 	z1.clear();
	const size_type block_nr    = (itx->first-2*(Deltaj-2)*Nablaj) / (Deltaj-2);
 	z1.insert(std::pair<size_type,double>(block_nr, 1.0));
 	basis1d().Mj1c_.apply(z1, z3);
	
 	// 3. combine results
 	for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
 	  for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
	    y[2*(Deltajp1-2)*(Deltajp1-2)+it3->first*(Deltajp1-2)+it2->first]
	      += itx->second * it2->second * it3->second;
      }
	break;
      case 3: {
	// contribution from patch 0
	//
	// apply kron(N#,ML):
 	// 1. get column of second factor ML
	V z1, z2, z3;
	z1.insert(std::pair<size_type,double>(0, 1.0));
	basis1d().Mj0_.apply(z1, z2); // we have to neglect the first entry of z2 later

 	// 2. get column of first factor N#
 	z1.clear();
 	const size_type block_nr = itx->first-3*(Deltaj-2)*Nablaj;
 	z1.insert(std::pair<size_type,double>(block_nr, 1.0));
 	basis1d().Mj1c_.apply(z1, z3);
	
 	// 3. combine results
 	for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
 	  for (V::const_iterator it2(++(z2.begin())); it2 != z2.end(); ++it2)
	    y[it3->first*(Deltajp1-2)+it2->first-1] // note the "-1"
	      += itx->second * it2->second * it3->second * M_SQRT1_2;

 	// contribution from patch 1
 	//
 	// apply kron(N#,MR):
  	// 1. get column of second factor MR
 	z1.clear();
 	z1.insert(std::pair<size_type,double>(Deltaj-1, 1.0));
 	z2.clear();
 	basis1d().Mj0_.apply(z1, z2); // we have to neglect the last entry of z2 later

 	// 2. get column of first factor N#: this has already been done above
	
 	// 3. combine results
  	for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
  	  for (V::reverse_iterator it2(++(z2.rbegin())); it2 != z2.rend(); ++it2)
 	    y[(Deltajp1-2)*(Deltajp1-2)+it3->first*(Deltajp1-2)+it2->first-1] // note the "-1"
 	      += itx->second * it2->second * it3->second * M_SQRT1_2;

 	// contribution from patch 3
 	//
 	// apply kron(N#,Mtopleft):
 	// 1. get column of first factor N#: this has already been done above

 	// 2. get top left entry of Mj0
 	const double Mtopleft = basis1d().Mj0_.get_entry(0,0);
	
 	// 3. combine results
 	for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
 	  y[3*(Deltajp1-2)*(Deltajp1-2)+it3->first]
 	    += itx->second * it3->second * Mtopleft;
      }
	break;
      default:
	break;
      }
    }      
  }
  
  template <int d, int dT>
  void
  LDomainBasis<SplineBasis<d,dT,DS_construction> >::apply_Mj1c_11
  (const int j,
   const std::map<size_type,double>& x, 
   std::map<size_type,double>& y) const
  {
    typedef std::map<size_type,double> V;
    
    const unsigned int Nablaj   = basis1d().Nablasize(j);
    const unsigned int Deltajp1 = basis1d().Deltasize(j+1);
    
    basis1d().Mj1c_.set_level(j);

    for (V::const_iterator itx(x.begin()); itx != x.end(); ++itx) {
      // determine patch number
      const unsigned int patch = itx->first / (Nablaj*Nablaj);

      switch(patch) {
      case 0: {
	// apply kron(N#,N#)
 	// 1. get column of second factor N#
 	V z1, z2, z3;
	const size_type block_index = itx->first % Nablaj;
 	z1.insert(std::pair<size_type,double>(block_index, 1.0));
 	basis1d().Mj1c_.apply(z1, z2);
	
 	// 2. get column of first factor N#
 	z1.clear();
	const size_type block_nr    = itx->first / Nablaj;
 	z1.insert(std::pair<size_type,double>(block_nr, 1.0));
 	basis1d().Mj1c_.apply(z1, z3);
	
 	// 3. combine results
 	for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
 	  for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
	    y[it3->first*(Deltajp1-2)+it2->first]
	      += itx->second * it2->second * it3->second;
      }
	break;
      case 1: {
	// apply kron(N#,N#)
	const size_type block_nr    = (itx->first-Nablaj*Nablaj) / Nablaj;
	const size_type block_index = (itx->first-Nablaj*Nablaj) % Nablaj;
	
 	V z1, z2, z3;
 	z1.insert(std::pair<size_type,double>(block_index, 1.0));
 	basis1d().Mj1c_.apply(z1, z2);

 	z1.clear();
 	z1.insert(std::pair<size_type,double>(block_nr, 1.0));
 	basis1d().Mj1c_.apply(z1, z3);

 	for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
 	  for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
	    y[(Deltajp1-2)*(Deltajp1-2)+it3->first*(Deltajp1-2)+it2->first]
	      += itx->second * it2->second * it3->second;
      }
	break;
      case 2: {
	// apply kron(N#,N#)
	const size_type block_nr    = (itx->first-2*Nablaj*Nablaj) / Nablaj;
	const size_type block_index = (itx->first-2*Nablaj*Nablaj) % Nablaj;
	
 	V z1, z2, z3;
 	z1.insert(std::pair<size_type,double>(block_index, 1.0));
 	basis1d().Mj1c_.apply(z1, z2);

 	z1.clear();
 	z1.insert(std::pair<size_type,double>(block_nr, 1.0));
 	basis1d().Mj1c_.apply(z1, z3);

 	for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
 	  for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
	    y[2*(Deltajp1-2)*(Deltajp1-2)+it3->first*(Deltajp1-2)+it2->first]
	      += itx->second * it2->second * it3->second;
      }
	break;
      default:
	break;
      }
    }
  }  
}
