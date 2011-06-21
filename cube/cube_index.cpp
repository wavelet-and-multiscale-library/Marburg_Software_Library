// implementation for cube_index.h
#include <cmath>

namespace WaveletTL
{
  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  CubeIndex<IBASIS,DIM,CUBEBASIS>::CubeIndex(const CUBEBASIS* basis)
    : basis_(basis)
  {
    if (basis_ == 0) {
      j_ = 0; // invalid (e and k are initialized by zero automatically)
    } else {
      j_ = basis_->j0(); // coarsest level;
      // e_ is zero by default: generator
      for (unsigned int i = 0; i < DIM; i++)
	k_[i] = basis_->bases()[i]->DeltaLmin();
      
      num_ = -1;
    }
  }

  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  CubeIndex<IBASIS,DIM,CUBEBASIS>::CubeIndex(const int j,
					     const type_type& e,
					     const translation_type& k,
					     const CUBEBASIS* basis)
    : basis_(basis), j_(j), e_(e), k_(k), num_(-1)
  {
    set_number();
  }

  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  CubeIndex<IBASIS,DIM,CUBEBASIS>::CubeIndex(const CubeIndex& lambda)
    : basis_(lambda.basis_), j_(lambda.j_), e_(lambda.e_), k_(lambda.k_), num_(lambda.num_)
  {
    set_number();
  }

  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  CubeIndex<IBASIS,DIM,CUBEBASIS>::CubeIndex(const CubeIndex* lambda)
    : basis_(lambda->basis_), j_(lambda->j_), e_(lambda->e_), k_(lambda->k_), num_(lambda->num_)
  { 
    set_number();
  }

  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  CubeIndex<IBASIS,DIM,CUBEBASIS>::CubeIndex(const int num,
					     const CUBEBASIS* basis)
    : basis_(basis)
  {
    num_ = num; 

  if(DIM == 2 && false)
  {
    int j0 = basis->j0();
    int Nabla1 = basis->bases()[0]->Nablasize(j0);
    int Nabla2 = basis->bases()[1]->Nablasize(j0);
    int Delta1 = basis->bases()[0]->Deltasize(j0);
    int Delta2 = basis->bases()[1]->Deltasize(j0);
    double tmp2 = (double)(Delta1*Nabla2 + Delta2*Nabla1)/(2.*Nabla1*Nabla2);
    double tmp;
    int act_num = num_;

    if (act_num < Delta1*Delta2) 
    {
      j_ = j0;
      e_[0]=0;
      e_[1]=0;
      k_[0] = basis->bases()[0]->DeltaLmin() + act_num/Delta2;
      k_[1] = basis->bases()[1]->DeltaLmin() + (act_num % Delta2);
    }
    else
    {
      tmp = std::sqrt(tmp2*tmp2 + (double)(num-Delta1*Delta2)/(Nabla1*Nabla2) ) +1. -tmp2;
      j_ = floor(log2(tmp)) + j0;
      act_num -= (Delta1 - Nabla1 + Nabla1*(1<<(j_-j0)))*(Delta2 - Nabla2 + Nabla2*(1<<(j_-j0)));
      if(act_num < (Delta1 - Nabla1 +(1<<(j_-j0))*Nabla1)*(1<<(j_-j0))*Nabla2)
      {
        e_[0]=0;
        e_[1]=1;
        k_[0] = basis->bases()[0]->DeltaLmin() + act_num/((1<<(j_-j0))*Nabla2);
        k_[1] = basis->bases()[1]->Nablamin() + (act_num % ((1<<(j_-j0))*Nabla2));
      }
      else if (act_num < (Delta1 - Nabla1 + (1<<(j_-j0))*Nabla1)*(1<<(j_-j0))*Nabla2 + (1<<(j_-j0))*Nabla1*(Delta2 - Nabla2 +(1<<(j_-j0))*Nabla2))
      {
        e_[0]=1;
        e_[1]=0;
        act_num -= (Delta1 - Nabla1 +(1<<(j_-j0))*Nabla1)*(1<<(j_-j0))*Nabla2;
        k_[0] = basis->bases()[0]->Nablamin() + act_num/(Delta2 - Nabla2 +(1<<(j_-j0))*Nabla2);
        k_[1] = basis->bases()[1]->DeltaLmin() + (act_num % (Delta2 - Nabla2 +(1<<(j_-j0))*Nabla2));
      }
      else
      {
        e_[0]=1;
        e_[1]=1;
        act_num -= (Delta1 - Nabla1 +(1<<(j_-j0))*Nabla1)*(1<<(j_-j0))*Nabla2 + (1<<(j_-j0))*Nabla1*(Delta2 - Nabla2 +(1<<(j_-j0))*Nabla2);
        k_[0] = basis->bases()[0]->Nablamin() + act_num/((1<<(j_-j0))*Nabla2);
        k_[1] = basis->bases()[1]->Nablamin() + (act_num % ((1<<(j_-j0))*Nabla2));
      }
    }
  } //end if (DIM == 2)
  else if (DIM == 3)
  {
    int j0 = basis->j0();
    int Nabla1 = basis->bases()[0]->Nablasize(j0);
    int Nabla2 = basis->bases()[1]->Nablasize(j0);
    int Nabla3 = basis->bases()[2]->Nablasize(j0);
    int Delta1 = basis->bases()[0]->Deltasize(j0);
    int Delta2 = basis->bases()[1]->Deltasize(j0);
    int Delta3 = basis->bases()[2]->Deltasize(j0);
    int act_num = num_;

    //Benötigt fürs Lösen der qubischen Gleichung schon normiert
    double a = (double)(Delta1*Nabla2*Nabla3 + Delta2*Nabla1*Nabla3 + Delta3*Nabla2*Nabla1)/(Nabla1*Nabla2*Nabla3);
    double b = (double)(Delta1*Delta2*Nabla3 + Delta3*Delta1*Nabla2 + Delta3*Delta2*Nabla1)/(Nabla1*Nabla2*Nabla3);
    double c = (double) (Delta1*Delta2*Delta3 - num)/(Nabla1*Nabla2*Nabla3);
    //Hilfsvariablen der Cardanischen Formeln
    double p = b - (a*a)/3.;
    double q = (2.*a*a*a)/27. - (a*b)/3. + c;
    //Diskriminante
    double D = (q*q)/4. + (p*p*p)/27.;
    //std::cout << "D = " << D << std::endl;
    //Hilfsvariablen der Lösung
    double u = pow(sqrt(D) -q/2.,1./3.);
    double v = pow(-sqrt(D) -q/2.,1./3.);
    double z = u + v;
    double x = z - a/3. +1.00001 ;


    if (num < Delta1*Delta2*Delta3) 
    {
      j_ = j0;
      e_[0]=0;
      e_[1]=0;
      e_[2]=0;
      k_[0] = basis->bases()[0]->DeltaLmin() + act_num/(Delta3*Delta2);
      k_[1] = basis->bases()[1]->DeltaLmin() + (act_num/Delta3 % Delta2);
      k_[2] = basis->bases()[2]->DeltaLmin() + (act_num % (Delta3));
    }
    else
    {
      j_ = floor(log2(x)) + j0;
      act_num -= (Delta1 - Nabla1 +(1<<(j_-j0))*Nabla1)*(Delta2 - Nabla2 +(1<<(j_-j0))*Nabla2)*(Delta3 - Nabla3 +(1<<(j_-j0))*Nabla3);
      if (act_num < (Delta1 - Nabla1 +(1<<(j_-j0))*Nabla1)*(Delta2 - Nabla2 +(1<<(j_-j0))*Nabla2)*(1<<(j_-j0))*Nabla3)
      {
        e_[0]=0;
        e_[1]=0;
        e_[2]=1;
        k_[0] = basis->bases()[0]->DeltaLmin() + act_num/((1<<(j_-j0))*Nabla3*(Delta2 - Nabla2 +(1<<(j_-j0))*Nabla2));
        k_[1] = basis->bases()[1]->DeltaLmin() + (act_num/((1<<(j_-j0))*Nabla3) % (Delta2 - Nabla2 +(1<<(j_-j0))*Nabla2));
        k_[2] = basis->bases()[2]->Nablamin() +  (act_num % ((1<<(j_-j0))*Nabla3));
      }
      else if (act_num < (Delta1 - Nabla1 +(1<<(j_-j0))*Nabla1)*(Delta2 - Nabla2 +(1<<(j_-j0))*Nabla2)*(1<<(j_-j0))*Nabla3
                          + (Delta1 - Nabla1 +(1<<(j_-j0))*Nabla1)*(1<<(j_-j0))*Nabla2*(Delta3 - Nabla3 +(1<<(j_-j0))*Nabla3))
      {
        e_[0]=0;
        e_[1]=1;
        e_[2]=0;
        act_num -= (Delta1 - Nabla1 +(1<<(j_-j0))*Nabla1)*(Delta2 - Nabla2 +(1<<(j_-j0))*Nabla2)*(1<<(j_-j0))*Nabla3;
        k_[0] = basis->bases()[0]->DeltaLmin() + act_num/((Delta3 - Nabla3 +(1<<(j_-j0))*Nabla3)*(1<<(j_-j0))*Nabla2);
        k_[1] = basis->bases()[1]->Nablamin() + act_num/((Delta3 - Nabla3 +(1<<(j_-j0))*Nabla3)) % ((1<<(j_-j0))*Nabla2);
        k_[2] = basis->bases()[2]->DeltaLmin() +  (act_num % ((Delta3 - Nabla3 +(1<<(j_-j0))*Nabla3)));
      }
      else if (act_num < (Delta1 - Nabla1 +(1<<(j_-j0))*Nabla1)*(Delta2 - Nabla2 +(1<<(j_-j0))*Nabla2)*(1<<(j_-j0))*Nabla3
                          + (Delta1 - Nabla1 +(1<<(j_-j0))*Nabla1)*(1<<(j_-j0))*Nabla2*(Delta3 - Nabla3 +(1<<(j_-j0))*Nabla3)
                          + (Delta1 - Nabla1 +(1<<(j_-j0))*Nabla1)*(1<<(j_-j0))*Nabla2*(1<<(j_-j0))*Nabla3)
      {
        e_[0]=0;
        e_[1]=1;
        e_[2]=1;
        act_num -= (Delta1 - Nabla1 +(1<<(j_-j0))*Nabla1)*(Delta2 - Nabla2 +(1<<(j_-j0))*Nabla2)*(1<<(j_-j0))*Nabla3
                    + (Delta1 - Nabla1 +(1<<(j_-j0))*Nabla1)*(1<<(j_-j0))*Nabla2*(Delta3 - Nabla3 +(1<<(j_-j0))*Nabla3);
        k_[0] = basis->bases()[0]->DeltaLmin() + act_num/((1<<(j_-j0))*Nabla3*(1<<(j_-j0))*Nabla2);
        k_[1] = basis->bases()[1]->Nablamin() + act_num/((1<<(j_-j0))*Nabla3) % ((1<<(j_-j0))*Nabla2);
        k_[2] = basis->bases()[2]->Nablamin() +  (act_num % ((1<<(j_-j0))*Nabla3));
      }
      else if (act_num < (Delta1 - Nabla1 +(1<<(j_-j0))*Nabla1)*(Delta2 - Nabla2 +(1<<(j_-j0))*Nabla2)*(1<<(j_-j0))*Nabla3
                          + (Delta1 - Nabla1 +(1<<(j_-j0))*Nabla1)*(1<<(j_-j0))*Nabla2*(Delta3 - Nabla3 +(1<<(j_-j0))*Nabla3)
                          + (Delta1 - Nabla1 +(1<<(j_-j0))*Nabla1)*(1<<(j_-j0))*Nabla2*(1<<(j_-j0))*Nabla3
                          + (1<<(j_-j0))*Nabla1*(Delta2 - Nabla2 +(1<<(j_-j0))*Nabla2)*(Delta3 - Nabla3 +(1<<(j_-j0))*Nabla3))
      {
        e_[0]=1;
        e_[1]=0;
        e_[2]=0;
        act_num -= (Delta1 - Nabla1 +(1<<(j_-j0))*Nabla1)*(Delta2 - Nabla2 +(1<<(j_-j0))*Nabla2)*(1<<(j_-j0))*Nabla3
                    + (Delta1 - Nabla1 +(1<<(j_-j0))*Nabla1)*(1<<(j_-j0))*Nabla2*(Delta3 - Nabla3 +(1<<(j_-j0))*Nabla3)
                    + (Delta1 - Nabla1 +(1<<(j_-j0))*Nabla1)*(1<<(j_-j0))*Nabla2*(1<<(j_-j0))*Nabla3;
        k_[0] = basis->bases()[0]->Nablamin() + act_num/((Delta3 - Nabla3 +(1<<(j_-j0))*Nabla3)*(Delta2 - Nabla2 +(1<<(j_-j0))*Nabla2));
        k_[1] = basis->bases()[1]->DeltaLmin() + act_num/((Delta3 - Nabla3 +(1<<(j_-j0))*Nabla3)) % ((Delta2 - Nabla2 +(1<<(j_-j0))*Nabla2));
        k_[2] = basis->bases()[2]->DeltaLmin() +  (act_num % ((Delta3 - Nabla3 +(1<<(j_-j0))*Nabla3)));
      }
      else if (act_num < (Delta1 - Nabla1 +(1<<(j_-j0))*Nabla1)*(Delta2 - Nabla2 +(1<<(j_-j0))*Nabla2)*(1<<(j_-j0))*Nabla3
                          + (Delta1 - Nabla1 +(1<<(j_-j0))*Nabla1)*(1<<(j_-j0))*Nabla2*(Delta3 - Nabla3 +(1<<(j_-j0))*Nabla3)
                          + (Delta1 - Nabla1 +(1<<(j_-j0))*Nabla1)*(1<<(j_-j0))*Nabla2*(1<<(j_-j0))*Nabla3
                          + (1<<(j_-j0))*Nabla1*(Delta2 - Nabla2 +(1<<(j_-j0))*Nabla2)*(Delta3 - Nabla3 +(1<<(j_-j0))*Nabla3)
                          + (1<<(j_-j0))*Nabla1*(Delta2 - Nabla2 +(1<<(j_-j0))*Nabla2)*(1<<(j_-j0))*Nabla3)
      {
        e_[0]=1;
        e_[1]=0;
        e_[2]=1;
        act_num -= (Delta1 - Nabla1 +(1<<(j_-j0))*Nabla1)*(Delta2 - Nabla2 +(1<<(j_-j0))*Nabla2)*(1<<(j_-j0))*Nabla3
                    + (Delta1 - Nabla1 +(1<<(j_-j0))*Nabla1)*(1<<(j_-j0))*Nabla2*(Delta3 - Nabla3 +(1<<(j_-j0))*Nabla3)
                    + (Delta1 - Nabla1 +(1<<(j_-j0))*Nabla1)*(1<<(j_-j0))*Nabla2*(1<<(j_-j0))*Nabla3
                    + (1<<(j_-j0))*Nabla1*(Delta2 - Nabla2 +(1<<(j_-j0))*Nabla2)*(Delta3 - Nabla3 +(1<<(j_-j0))*Nabla3);
        k_[0] = basis->bases()[0]->Nablamin() + act_num/((1<<(j_-j0))*Nabla3*(Delta2 - Nabla2 +(1<<(j_-j0))*Nabla2));
        k_[1] = basis->bases()[1]->DeltaLmin() + act_num/((1<<(j_-j0))*Nabla3) % ((Delta2 - Nabla2 +(1<<(j_-j0))*Nabla2));
        k_[2] = basis->bases()[2]->Nablamin() +  (act_num % ((1<<(j_-j0))*Nabla3));
      }
      else if (act_num < (Delta1 - Nabla1 +(1<<(j_-j0))*Nabla1)*(Delta2 - Nabla2 +(1<<(j_-j0))*Nabla2)*(1<<(j_-j0))*Nabla3
                          + (Delta1 - Nabla1 +(1<<(j_-j0))*Nabla1)*(1<<(j_-j0))*Nabla2*(Delta3 - Nabla3 +(1<<(j_-j0))*Nabla3)
                          + (Delta1 - Nabla1 +(1<<(j_-j0))*Nabla1)*(1<<(j_-j0))*Nabla2*(1<<(j_-j0))*Nabla3
                          + (1<<(j_-j0))*Nabla1*(Delta2 - Nabla2 +(1<<(j_-j0))*Nabla2)*(Delta3 - Nabla3 +(1<<(j_-j0))*Nabla3)
                          + (1<<(j_-j0))*Nabla1*(Delta2 - Nabla2 +(1<<(j_-j0))*Nabla2)*(1<<(j_-j0))*Nabla3
                          + (1<<(j_-j0))*Nabla1*(1<<(j_-j0))*Nabla2*(Delta3 - Nabla3 +(1<<(j_-j0))*Nabla3))
      {
        e_[0]=1;
        e_[1]=1;
        e_[2]=0;
        act_num -= (Delta1 - Nabla1 +(1<<(j_-j0))*Nabla1)*(Delta2 - Nabla2 +(1<<(j_-j0))*Nabla2)*(1<<(j_-j0))*Nabla3
                    + (Delta1 - Nabla1 +(1<<(j_-j0))*Nabla1)*(1<<(j_-j0))*Nabla2*(Delta3 - Nabla3 +(1<<(j_-j0))*Nabla3)
                    + (Delta1 - Nabla1 +(1<<(j_-j0))*Nabla1)*(1<<(j_-j0))*Nabla2*(1<<(j_-j0))*Nabla3
                    + (1<<(j_-j0))*Nabla1*(Delta2 - Nabla2 +(1<<(j_-j0))*Nabla2)*(Delta3 - Nabla3 +(1<<(j_-j0))*Nabla3)
                    + (1<<(j_-j0))*Nabla1*(Delta2 - Nabla2 +(1<<(j_-j0))*Nabla2)*(1<<(j_-j0))*Nabla3;
        k_[0] = basis->bases()[0]->Nablamin() + act_num/((Delta3 - Nabla3 +(1<<(j_-j0))*Nabla3)*(1<<(j_-j0))*Nabla2);
        k_[1] = basis->bases()[1]->Nablamin() + act_num/((Delta3 - Nabla3 +(1<<(j_-j0))*Nabla3)) % ((1<<(j_-j0))*Nabla2);
        k_[2] = basis->bases()[2]->DeltaLmin() +  (act_num % ((Delta3 - Nabla3 +(1<<(j_-j0))*Nabla3)));
      }
      else 
      {
        e_[0]=1;
        e_[1]=1;
        e_[2]=1;
        act_num -= (Delta1 - Nabla1 +(1<<(j_-j0))*Nabla1)*(Delta2 - Nabla2 +(1<<(j_-j0))*Nabla2)*(1<<(j_-j0))*Nabla3
                    + (Delta1 - Nabla1 +(1<<(j_-j0))*Nabla1)*(1<<(j_-j0))*Nabla2*(Delta3 - Nabla3 +(1<<(j_-j0))*Nabla3)
                    + (Delta1 - Nabla1 +(1<<(j_-j0))*Nabla1)*(1<<(j_-j0))*Nabla2*(1<<(j_-j0))*Nabla3
                    + (1<<(j_-j0))*Nabla1*(Delta2 - Nabla2 +(1<<(j_-j0))*Nabla2)*(Delta3 - Nabla3 +(1<<(j_-j0))*Nabla3)
                    + (1<<(j_-j0))*Nabla1*(Delta2 - Nabla2 +(1<<(j_-j0))*Nabla2)*(1<<(j_-j0))*Nabla3
                    + (1<<(j_-j0))*Nabla1*(1<<(j_-j0))*Nabla2*(Delta3 - Nabla3 +(1<<(j_-j0))*Nabla3);
        k_[0] = basis->bases()[0]->Nablamin() + act_num/((1<<(j_-j0))*Nabla3*(1<<(j_-j0))*Nabla2);
        k_[1] = basis->bases()[1]->Nablamin() + act_num/((1<<(j_-j0))*Nabla3) % ((1<<(j_-j0))*Nabla2);
        k_[2] = basis->bases()[2]->Nablamin() +  (act_num % ((1<<(j_-j0))*Nabla3));
      }
    }

    //act_num -= tmp2;

  }//else DIM == 3
  else  
  {

    // to be decreased successively
    int act_num = num_;

    int j = basis_->j0();

    int tmp2 = 0;

    bool gen = 0;

    // determine the level of the index
    while (true) {
      int tmp = 1;
      for (unsigned int i = 0; i < DIM; i++)
	tmp *= (basis_->bases()[i])->Deltasize(j);
      if (tmp > num)
	break;
      tmp2 = tmp;
      j++;
    }
    if (j == basis_->j0()) {
      j_ = j;
      gen = 1;
    }
    else
      j_ = j-1;

    act_num -= tmp2;

    tmp2 = 0;

    // determine the type of the index 
    int tmp2_old = 0;
    if (gen)
      e_ = type_type();
    else {
      MultiIndex<int,DIM> type;
      type[DIM-1] = 1;
      bool exit = 0;
      // loop over types
      while (true) {
	int tmp = 1;
	for (unsigned int i = 0; i < DIM; i++) {
	  if (type[i] == 0)
	    tmp *= (basis_->bases()[i])->Deltasize(j_);
	  else
	    tmp *= (basis_->bases()[i])->Nablasize(j_);
	}
	tmp2_old = tmp2;
	tmp2 += tmp;
	
	// right type found
	if (tmp2 > act_num)
	  break;

	// determine next type
	for (unsigned int i = DIM-1; i >= 0; i--) {
	  if ( type[i] == 1 ) {
	    type[i] = 0;
	    exit = (i == 0);
	  }
	  else {
	    type[i]++;
	    break;
	  }
	}
	if (exit)
	  break;
      }//end while
      e_ = type;
    }// end else

    act_num -= tmp2_old;
    tmp2 = 1;

    // determine the position of the index
    for (unsigned int i = DIM-1; i > 0; i--) {
      if (e_[i] == 0) {
	tmp2 *= (basis_->bases()[i])->Deltasize(j_);
      }
      else {
	tmp2 *= (basis_->bases()[i])->Nablasize(j_);
      }
    }


    act_num += 1;
    for (unsigned int i = 0; i < DIM; i++) {
      int tmp = 0;
      if (act_num <= tmp2) {
 	if (e_[i] == 0)
	  k_[i] = (basis_->bases()[i])->DeltaLmin();
 	else
 	  k_[i] = (basis_->bases()[i])->Nablamin();
      }
      else {// act_num > tmp
	if (tmp2 == 1)
	  tmp = act_num-1;
	else if ((tmp2 != 1) && ((act_num % tmp2) != 0))
	  tmp = act_num / tmp2;
	else if ( (act_num % tmp2) == 0 ) 
	  tmp = act_num / tmp2 - 1;

 	if (e_[i] == 0)
 	  k_[i] = tmp + (basis_->bases()[i])->DeltaLmin();
 	else
 	  k_[i] = tmp + (basis_->bases()[i])->Nablamin();
	
	if ( (act_num % tmp2) != 0 )
	  act_num = act_num % tmp2;
	else
	  act_num = tmp2;
      }
      if ((i+1) < DIM) {
	if (e_[i+1] == 0) {
	  tmp2 /= (basis_->bases()[i+1])->Deltasize(j_);
	}
	else {
	  tmp2 /= (basis_->bases()[i+1])->Nablasize(j_);
	}
      }
    }
   }//end DIM >3
  }
    
  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  CubeIndex<IBASIS,DIM,CUBEBASIS>&
  CubeIndex<IBASIS,DIM,CUBEBASIS>::operator = (const CubeIndex<IBASIS,DIM,CUBEBASIS>& lambda)
  {
    j_ = lambda.j();
    e_ = lambda.e();
    k_ = lambda.k();
    basis_ = lambda.basis();
    num_ = lambda.number();

    return *this;
  }


  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  bool
  CubeIndex<IBASIS,DIM,CUBEBASIS>::operator == (const CubeIndex& lambda) const
  {
    return (j_ == lambda.j() &&
	    e_ == lambda.e() &&
	    k_ == lambda.k());
  }

  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  bool
  CubeIndex<IBASIS,DIM,CUBEBASIS>::operator < (const CubeIndex& lambda) const
  {
    // standard lexicographic order on (j,e,k),
    // we assume that e and k are already lexicographically ordered (cf. MultiIndex)
    return (j_ < lambda.j() ||
	    (j_ == lambda.j() && (e_ < lambda.e() ||
				  (e_ == lambda.e() && k_ < lambda.k()))));
  }
  
  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  CubeIndex<IBASIS,DIM,CUBEBASIS>&
  CubeIndex<IBASIS,DIM,CUBEBASIS>::operator ++ ()
  {
    if (num_ > -1)
      num_++;

    bool eplusplus = false;
    for (int i = DIM-1; i >= 0; i--) {
      const int last_index = (e_[i] == 0
			      ? basis_->bases()[i]->DeltaRmax(j_)
			      : basis_->bases()[i]->Nablamax(j_));
      if (k_[i] == last_index) {
	k_[i] = (e_[i] == 0
		 ? basis_->bases()[i]->DeltaLmin()
		 : basis_->bases()[i]->Nablamin());
	eplusplus = (i == 0);
      } else {
	++k_[i];
	break;
      }
    }

    bool jplusplus = false;
    if (eplusplus) {
      for (int i = DIM-1; i >= 0; i--) {
	if (e_[i] == 1) {
	  e_[i] = 0;
	  jplusplus = (i == 0);
	} else {
	  ++e_[i];
	  break;
	}
      }
     
      if (!jplusplus) // then choose lowest translation index k=k(j,e)
	for (unsigned int i = 0; i < DIM; i++)
	  k_[i] = (e_[i] == 0
		   ? basis_->bases()[i]->DeltaLmin()
		   : basis_->bases()[i]->Nablamin());
    }

    if (jplusplus) {
      ++j_;
      // choose lowest type e=(0,...,0,1) and lowest translation index k=k(j,e)
      for (unsigned int i = 0; i < DIM-1; i++) {
	e_[i] = 0;
	k_[i] = basis_->bases()[i]->DeltaLmin();
      }
      e_[DIM-1] = 1;
      k_[DIM-1] = basis_->bases()[DIM-1]->Nablamin();
    }
    
    return *this;
  }

  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  void
  CubeIndex<IBASIS,DIM,CUBEBASIS>::get_Index_Parameter(const CUBEBASIS* basis, const int num, int& j, type_type& e, translation_type& k) const
  {

  if(DIM == 2)
  {
    int j0 = basis->j0();
    int Nabla1 = basis->bases()[0]->Nablasize(j0);
    int Nabla2 = basis->bases()[1]->Nablasize(j0);
    int Delta1 = basis->bases()[0]->Deltasize(j0);
    int Delta2 = basis->bases()[1]->Deltasize(j0);
    double tmp2 = (double)(Delta1*Nabla2 + Delta2*Nabla1)/(2.*Nabla1*Nabla2);
    double tmp;
    int act_num = num;

    if (num < Delta1*Delta2) 
    {
      j = j0;
      e[0]=0;
      e[1]=0;
      k[0] = basis->bases()[0]->DeltaLmin() + act_num/Delta2;
      k[1] = basis->bases()[1]->DeltaLmin() + (act_num % Delta2);
    }
    else
    {
      tmp = std::sqrt(tmp2*tmp2 + (double)(num-Delta1*Delta2)/(Nabla1*Nabla2) ) +1. -tmp2;
      j = floor(log2(tmp)) + j0;
      act_num -= (Delta1 - Nabla1 +(1<<(j-j0))*Nabla1)*(Delta2 - Nabla2 +(1<<(j-j0))*Nabla2);
      if(act_num < (Delta1 - Nabla1 +(1<<(j-j0))*Nabla1)*(1<<(j-j0))*Nabla2)
      {
        e[0]=0;
        e[1]=1;
        k[0] = basis->bases()[0]->DeltaLmin() + act_num/((1<<(j-j0))*Nabla2);
        k[1] = basis->bases()[1]->Nablamin() + (act_num % ((1<<(j-j0))*Nabla2));
      }
      else if (act_num < (Delta1 - Nabla1 + (1<<(j-j0))*Nabla1)*(1<<(j-j0))*Nabla2 + (1<<(j-j0))*Nabla1*(Delta2 - Nabla2 +(1<<(j-j0))*Nabla2))
      {
        e[0]=1;
        e[1]=0;
        act_num -= (Delta1 - Nabla1 +(1<<(j-j0))*Nabla1)*(1<<(j-j0))*Nabla2;
        k[0] = basis->bases()[0]->Nablamin() + act_num/(Delta2 - Nabla2 +(1<<(j-j0))*Nabla2);
        k[1] = basis->bases()[1]->DeltaLmin() + (act_num % (Delta2 - Nabla2 +(1<<(j-j0))*Nabla2));
      }
      else
      {
        e[0]=1;
        e[1]=1;
        act_num -= (Delta1 - Nabla1 +(1<<(j-j0))*Nabla1)*(1<<(j-j0))*Nabla2 + (1<<(j-j0))*Nabla1*(Delta2 - Nabla2 +(1<<(j-j0))*Nabla2);
        k[0] = basis->bases()[0]->Nablamin() + act_num/((1<<(j-j0))*Nabla2);
        k[1] = basis->bases()[1]->Nablamin() + (act_num % ((1<<(j-j0))*Nabla2));
      }
    }
  } //end if (DIM == 2)
  else if (DIM == 3)
  {
    int j0 = basis->j0();
    int Nabla1 = basis->bases()[0]->Nablasize(j0);
    int Nabla2 = basis->bases()[1]->Nablasize(j0);
    int Nabla3 = basis->bases()[2]->Nablasize(j0);
    int Delta1 = basis->bases()[0]->Deltasize(j0);
    int Delta2 = basis->bases()[1]->Deltasize(j0);
    int Delta3 = basis->bases()[2]->Deltasize(j0);
    int act_num = num;

    //Benötigt fürs Lösen der qubischen Gleichung schon normiert
    double a = (double)(Delta1*Nabla2*Nabla3 + Delta2*Nabla1*Nabla3 + Delta3*Nabla2*Nabla1)/(Nabla1*Nabla2*Nabla3);
    double b = (double)(Delta1*Delta2*Nabla3 + Delta3*Delta1*Nabla2 + Delta3*Delta2*Nabla1)/(Nabla1*Nabla2*Nabla3);
    double c = (double) (Delta1*Delta2*Delta3 - num)/(Nabla1*Nabla2*Nabla3);
    //Hilfsvariablen der Cardanischen Formeln
    double p = b - (a*a)/3.;
    double q = (2.*a*a*a)/27. - (a*b)/3. + c;
    //Diskriminante
    double D = (q*q)/4. + (p*p*p)/27.;
    //std::cout << "D = " << D << std::endl;
    //Hilfsvariablen der Lösung
    double u = pow(sqrt(D) -q/2.,1./3.);
    double v = pow(-sqrt(D) -q/2.,1./3.);
    double z = u + v;
    double x = z - a/3. +1.00001 ;


    if (num < Delta1*Delta2*Delta3) 
    {
      j = j0;
      e[0]=0;
      e[1]=0;
      e[2]=0;
      k[0] = basis->bases()[0]->DeltaLmin() + act_num/(Delta3*Delta2);
      k[1] = basis->bases()[1]->DeltaLmin() + (act_num/Delta3 % Delta2);
      k[2] = basis->bases()[2]->DeltaLmin() + (act_num % (Delta3));
    }
    else
    {
      j = floor(log2(x)) + j0;
      act_num -= (Delta1 - Nabla1 +(1<<(j-j0))*Nabla1)*(Delta2 - Nabla2 +(1<<(j-j0))*Nabla2)*(Delta3 - Nabla3 +(1<<(j-j0))*Nabla3);
      if (act_num < (Delta1 - Nabla1 +(1<<(j-j0))*Nabla1)*(Delta2 - Nabla2 +(1<<(j-j0))*Nabla2)*(1<<(j-j0))*Nabla3)
      {
        e[0]=0;
        e[1]=0;
        e[2]=1;
        k[0] = basis->bases()[0]->DeltaLmin() + act_num/((1<<(j-j0))*Nabla3*(Delta2 - Nabla2 +(1<<(j-j0))*Nabla2));
        k[1] = basis->bases()[1]->DeltaLmin() + (act_num/((1<<(j-j0))*Nabla3) % (Delta2 - Nabla2 +(1<<(j-j0))*Nabla2));
        k[2] = basis->bases()[2]->Nablamin() +  (act_num % ((1<<(j-j0))*Nabla3));
      }
      else if (act_num < (Delta1 - Nabla1 +(1<<(j-j0))*Nabla1)*(Delta2 - Nabla2 +(1<<(j-j0))*Nabla2)*(1<<(j-j0))*Nabla3
                          + (Delta1 - Nabla1 +(1<<(j-j0))*Nabla1)*(1<<(j-j0))*Nabla2*(Delta3 - Nabla3 +(1<<(j-j0))*Nabla3))
      {
        e[0]=0;
        e[1]=1;
        e[2]=0;
        act_num -= (Delta1 - Nabla1 +(1<<(j-j0))*Nabla1)*(Delta2 - Nabla2 +(1<<(j-j0))*Nabla2)*(1<<(j-j0))*Nabla3;
        k[0] = basis->bases()[0]->DeltaLmin() + act_num/((Delta3 - Nabla3 +(1<<(j-j0))*Nabla3)*(1<<(j-j0))*Nabla2);
        k[1] = basis->bases()[1]->Nablamin() + act_num/((Delta3 - Nabla3 +(1<<(j-j0))*Nabla3)) % ((1<<(j-j0))*Nabla2);
        k[2] = basis->bases()[2]->DeltaLmin() +  (act_num % ((Delta3 - Nabla3 +(1<<(j-j0))*Nabla3)));
      }
      else if (act_num < (Delta1 - Nabla1 +(1<<(j-j0))*Nabla1)*(Delta2 - Nabla2 +(1<<(j-j0))*Nabla2)*(1<<(j-j0))*Nabla3
                          + (Delta1 - Nabla1 +(1<<(j-j0))*Nabla1)*(1<<(j-j0))*Nabla2*(Delta3 - Nabla3 +(1<<(j-j0))*Nabla3)
                          + (Delta1 - Nabla1 +(1<<(j-j0))*Nabla1)*(1<<(j-j0))*Nabla2*(1<<(j-j0))*Nabla3)
      {
        e[0]=0;
        e[1]=1;
        e[2]=1;
        act_num -= (Delta1 - Nabla1 +(1<<(j-j0))*Nabla1)*(Delta2 - Nabla2 +(1<<(j-j0))*Nabla2)*(1<<(j-j0))*Nabla3
                    + (Delta1 - Nabla1 +(1<<(j-j0))*Nabla1)*(1<<(j-j0))*Nabla2*(Delta3 - Nabla3 +(1<<(j-j0))*Nabla3);
        k[0] = basis->bases()[0]->DeltaLmin() + act_num/((1<<(j-j0))*Nabla3*(1<<(j-j0))*Nabla2);
        k[1] = basis->bases()[1]->Nablamin() + act_num/((1<<(j-j0))*Nabla3) % ((1<<(j-j0))*Nabla2);
        k[2] = basis->bases()[2]->Nablamin() +  (act_num % ((1<<(j-j0))*Nabla3));
      }
      else if (act_num < (Delta1 - Nabla1 +(1<<(j-j0))*Nabla1)*(Delta2 - Nabla2 +(1<<(j-j0))*Nabla2)*(1<<(j-j0))*Nabla3
                          + (Delta1 - Nabla1 +(1<<(j-j0))*Nabla1)*(1<<(j-j0))*Nabla2*(Delta3 - Nabla3 +(1<<(j-j0))*Nabla3)
                          + (Delta1 - Nabla1 +(1<<(j-j0))*Nabla1)*(1<<(j-j0))*Nabla2*(1<<(j-j0))*Nabla3
                          + (1<<(j-j0))*Nabla1*(Delta2 - Nabla2 +(1<<(j-j0))*Nabla2)*(Delta3 - Nabla3 +(1<<(j-j0))*Nabla3))
      {
        e[0]=1;
        e[1]=0;
        e[2]=0;
        act_num -= (Delta1 - Nabla1 +(1<<(j-j0))*Nabla1)*(Delta2 - Nabla2 +(1<<(j-j0))*Nabla2)*(1<<(j-j0))*Nabla3
                    + (Delta1 - Nabla1 +(1<<(j-j0))*Nabla1)*(1<<(j-j0))*Nabla2*(Delta3 - Nabla3 +(1<<(j-j0))*Nabla3)
                    + (Delta1 - Nabla1 +(1<<(j-j0))*Nabla1)*(1<<(j-j0))*Nabla2*(1<<(j-j0))*Nabla3;
        k[0] = basis->bases()[0]->Nablamin() + act_num/((Delta3 - Nabla3 +(1<<(j-j0))*Nabla3)*(Delta2 - Nabla2 +(1<<(j-j0))*Nabla2));
        k[1] = basis->bases()[1]->DeltaLmin() + act_num/((Delta3 - Nabla3 +(1<<(j-j0))*Nabla3)) % ((Delta2 - Nabla2 +(1<<(j-j0))*Nabla2));
        k[2] = basis->bases()[2]->DeltaLmin() +  (act_num % ((Delta3 - Nabla3 +(1<<(j-j0))*Nabla3)));
      }
      else if (act_num < (Delta1 - Nabla1 +(1<<(j-j0))*Nabla1)*(Delta2 - Nabla2 +(1<<(j-j0))*Nabla2)*(1<<(j-j0))*Nabla3
                          + (Delta1 - Nabla1 +(1<<(j-j0))*Nabla1)*(1<<(j-j0))*Nabla2*(Delta3 - Nabla3 +(1<<(j-j0))*Nabla3)
                          + (Delta1 - Nabla1 +(1<<(j-j0))*Nabla1)*(1<<(j-j0))*Nabla2*(1<<(j-j0))*Nabla3
                          + (1<<(j-j0))*Nabla1*(Delta2 - Nabla2 +(1<<(j-j0))*Nabla2)*(Delta3 - Nabla3 +(1<<(j-j0))*Nabla3)
                          + (1<<(j-j0))*Nabla1*(Delta2 - Nabla2 +(1<<(j-j0))*Nabla2)*(1<<(j-j0))*Nabla3)
      {
        e[0]=1;
        e[1]=0;
        e[2]=1;
        act_num -= (Delta1 - Nabla1 +(1<<(j-j0))*Nabla1)*(Delta2 - Nabla2 +(1<<(j-j0))*Nabla2)*(1<<(j-j0))*Nabla3
                    + (Delta1 - Nabla1 +(1<<(j-j0))*Nabla1)*(1<<(j-j0))*Nabla2*(Delta3 - Nabla3 +(1<<(j-j0))*Nabla3)
                    + (Delta1 - Nabla1 +(1<<(j-j0))*Nabla1)*(1<<(j-j0))*Nabla2*(1<<(j-j0))*Nabla3
                    + (1<<(j-j0))*Nabla1*(Delta2 - Nabla2 +(1<<(j-j0))*Nabla2)*(Delta3 - Nabla3 +(1<<(j-j0))*Nabla3);
        k[0] = basis->bases()[0]->Nablamin() + act_num/((1<<(j-j0))*Nabla3*(Delta2 - Nabla2 +(1<<(j-j0))*Nabla2));
        k[1] = basis->bases()[1]->DeltaLmin() + act_num/((1<<(j-j0))*Nabla3) % ((Delta2 - Nabla2 +(1<<(j-j0))*Nabla2));
        k[2] = basis->bases()[2]->Nablamin() +  (act_num % ((1<<(j-j0))*Nabla3));
      }
      else if (act_num < (Delta1 - Nabla1 +(1<<(j-j0))*Nabla1)*(Delta2 - Nabla2 +(1<<(j-j0))*Nabla2)*(1<<(j-j0))*Nabla3
                          + (Delta1 - Nabla1 +(1<<(j-j0))*Nabla1)*(1<<(j-j0))*Nabla2*(Delta3 - Nabla3 +(1<<(j-j0))*Nabla3)
                          + (Delta1 - Nabla1 +(1<<(j-j0))*Nabla1)*(1<<(j-j0))*Nabla2*(1<<(j-j0))*Nabla3
                          + (1<<(j-j0))*Nabla1*(Delta2 - Nabla2 +(1<<(j-j0))*Nabla2)*(Delta3 - Nabla3 +(1<<(j-j0))*Nabla3)
                          + (1<<(j-j0))*Nabla1*(Delta2 - Nabla2 +(1<<(j-j0))*Nabla2)*(1<<(j-j0))*Nabla3
                          + (1<<(j-j0))*Nabla1*(1<<(j-j0))*Nabla2*(Delta3 - Nabla3 +(1<<(j-j0))*Nabla3))
      {
        e[0]=1;
        e[1]=1;
        e[2]=0;
        act_num -= (Delta1 - Nabla1 +(1<<(j-j0))*Nabla1)*(Delta2 - Nabla2 +(1<<(j-j0))*Nabla2)*(1<<(j-j0))*Nabla3
                    + (Delta1 - Nabla1 +(1<<(j-j0))*Nabla1)*(1<<(j-j0))*Nabla2*(Delta3 - Nabla3 +(1<<(j-j0))*Nabla3)
                    + (Delta1 - Nabla1 +(1<<(j-j0))*Nabla1)*(1<<(j-j0))*Nabla2*(1<<(j-j0))*Nabla3
                    + (1<<(j-j0))*Nabla1*(Delta2 - Nabla2 +(1<<(j-j0))*Nabla2)*(Delta3 - Nabla3 +(1<<(j-j0))*Nabla3)
                    + (1<<(j-j0))*Nabla1*(Delta2 - Nabla2 +(1<<(j-j0))*Nabla2)*(1<<(j-j0))*Nabla3;
        k[0] = basis->bases()[0]->Nablamin() + act_num/((Delta3 - Nabla3 +(1<<(j-j0))*Nabla3)*(1<<(j-j0))*Nabla2);
        k[1] = basis->bases()[1]->Nablamin() + act_num/((Delta3 - Nabla3 +(1<<(j-j0))*Nabla3)) % ((1<<(j-j0))*Nabla2);
        k[2] = basis->bases()[2]->DeltaLmin() +  (act_num % ((Delta3 - Nabla3 +(1<<(j-j0))*Nabla3)));
      }
      else 
      {
        e[0]=1;
        e[1]=1;
        e[2]=1;
        act_num -= (Delta1 - Nabla1 +(1<<(j-j0))*Nabla1)*(Delta2 - Nabla2 +(1<<(j-j0))*Nabla2)*(1<<(j-j0))*Nabla3
                    + (Delta1 - Nabla1 +(1<<(j-j0))*Nabla1)*(1<<(j-j0))*Nabla2*(Delta3 - Nabla3 +(1<<(j-j0))*Nabla3)
                    + (Delta1 - Nabla1 +(1<<(j-j0))*Nabla1)*(1<<(j-j0))*Nabla2*(1<<(j-j0))*Nabla3
                    + (1<<(j-j0))*Nabla1*(Delta2 - Nabla2 +(1<<(j-j0))*Nabla2)*(Delta3 - Nabla3 +(1<<(j-j0))*Nabla3)
                    + (1<<(j-j0))*Nabla1*(Delta2 - Nabla2 +(1<<(j-j0))*Nabla2)*(1<<(j-j0))*Nabla3
                    + (1<<(j-j0))*Nabla1*(1<<(j-j0))*Nabla2*(Delta3 - Nabla3 +(1<<(j-j0))*Nabla3);
        k[0] = basis->bases()[0]->Nablamin() + act_num/((1<<(j-j0))*Nabla3*(1<<(j-j0))*Nabla2);
        k[1] = basis->bases()[1]->Nablamin() + act_num/((1<<(j-j0))*Nabla3) % ((1<<(j-j0))*Nabla2);
        k[2] = basis->bases()[2]->Nablamin() +  (act_num % ((1<<(j-j0))*Nabla3));
      }
    }

    //act_num -= tmp2;

  }//else DIM == 3
  }

  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  void
  CubeIndex<IBASIS,DIM,CUBEBASIS>::set_number()
  {
    if (basis_ == 0) {
      num_ = -1;
      return;
    }
      
    type_type gen_type;
    assert( e_ != gen_type || j_ == basis_->j0() );
    int result = 0;
    int genfstlvl = 0;
    bool gen = 1;

    // check if wavelet is a generator
    for (unsigned int i = 0; i < DIM; i++) {
      if (e_[i] == 1) {
	gen = 0;
	break;
      }
    }

    // determine how many wavelets there are on all the levels
    // below the level of this index
    if (! gen) {
      result = 0;
      genfstlvl =1;
      //generators on level j0
      for (unsigned int i = 0; i< DIM; i++)                     
	genfstlvl *= (basis_->bases()[i])->Deltasize((basis_->bases()[i])->j0());
      //additional wavelets on level j
      //            =(#Gen[1]+#Wav[1])*...*(#Gen[Dim-1]+#Wav[Dim-1])
      //             -#Gen[1]*...*#Gen[Dim-1]
      for (int lvl= 0 ; 
	   lvl < (j_-basis_->j0());
	   lvl++){
	int genCurLvl = 1;
	int addWav = 1;
	for (unsigned int i = 0; i< DIM; i++) {
	  unsigned int curJ = basis_->bases()[i]->j0()+lvl;
	  int genCurDim = (basis_->bases()[i])->Deltasize(curJ);
	  genCurLvl *= genCurDim;
	  addWav *= genCurDim+ (basis_->bases()[i])->Nablasize(curJ);
	}
	result += addWav-genCurLvl;
      }
      result += genfstlvl;
    }

    // now determine how many wavelets there are on the same level
    // that belong to another type less than the type of this wavelet,
    // add the result to res afterwards
    if (! gen) {
      MultiIndex<int,DIM> type;
      type[DIM-1] = 1;
      bool exit = 0;
      // loop over all ''smaller'' types
      while (type < e_) {
	int tmp = 1;
	for (unsigned int i = 0; i < DIM; i++) {
	  if (type[i] == 0)
	    tmp *= (basis_->bases()[i])->Deltasize(j_);
	  else
	    tmp *= (basis_->bases()[i])->Nablasize(j_);
	}
	result += tmp;
	// determine next type
	for (unsigned int i = DIM-1; i >= 0; i--) {
	  if ( type[i] == 1 ) {
	    type[i] = 0;
	    exit = (i == 0);
	  }
	  else {
	    type[i]++;
	    break;
	  }
	}
	if (exit)
	  break;
      }//end while
    }// end if
    // count the wavelets that belong to the same level and to the same type
    // whose indices are smaller than this index,
    // add the result to res
    for (unsigned int i = 0; i < DIM; i++) {
      int tmp = 1;
      if (e_[i] == 0) {
	if (k_[i] == (basis_->bases()[i])->DeltaLmin())
	  continue;
      }
      else
	if (k_[i] == (basis_->bases()[i])->Nablamin())
	  continue;
      
      if (e_[i] == 0)
	tmp *= k_[i]-(basis_->bases()[i])->DeltaLmin();
      else
	tmp *= k_[i]-(basis_->bases()[i])->Nablamin();
      
      for (unsigned int l = i+1; l < DIM; l++) {
	if (e_[l] == 0) {
	  tmp *= (basis_->bases()[l])->Deltasize(j_);
	}
	else {
	  tmp *= (basis_->bases()[l])->Nablasize(j_);
	}
      }
      result += tmp;
    }
    num_ = result;
  }


  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  CubeIndex<IBASIS,DIM,CUBEBASIS>
  first_generator(const CUBEBASIS* basis, const int j)
  {
    assert(j >= basis->j0());

    typename CubeIndex<IBASIS,DIM,CUBEBASIS>::type_type e;
    typename CubeIndex<IBASIS,DIM,CUBEBASIS>::translation_type k;
    for (unsigned int i = 0; i < DIM; i++) {
      typename IBASIS::Index lambda(first_generator<IBASIS>(basis->bases()[i], j));
      k[i] = lambda.k();
    }

    return CubeIndex<IBASIS,DIM,CUBEBASIS>(j, e, k, basis);
  }

  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  CubeIndex<IBASIS,DIM,CUBEBASIS>
  last_generator(const CUBEBASIS* basis, const int j)
  {
    assert(j >= basis->j0());

    typename CubeIndex<IBASIS,DIM,CUBEBASIS>::type_type e;
    typename CubeIndex<IBASIS,DIM,CUBEBASIS>::translation_type k;
    for (unsigned int i = 0; i < DIM; i++) {
      typename IBASIS::Index lambda(last_generator<IBASIS>(basis->bases()[i], j));
      k[i] = lambda.k();
    }

    return CubeIndex<IBASIS,DIM,CUBEBASIS>(j, e, k, basis);
  }

  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  CubeIndex<IBASIS,DIM,CUBEBASIS>
  first_wavelet(const CUBEBASIS* basis, const int j)
  {
    assert(j >= basis->j0());
    
    typename CubeIndex<IBASIS,DIM,CUBEBASIS>::type_type e;
    typename CubeIndex<IBASIS,DIM,CUBEBASIS>::translation_type k;
    for (unsigned int i = 0; i < DIM-1; i++) {
      typename IBASIS::Index lambda(first_generator<IBASIS>(basis->bases()[i], j));
      k[i] = lambda.k();
    }
    typename IBASIS::Index lambda(first_wavelet<IBASIS>(basis->bases()[DIM-1], j));
    k[DIM-1] = lambda.k();
    e[DIM-1] = 1;

    return CubeIndex<IBASIS,DIM,CUBEBASIS>(j, e, k, basis);
  }

  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  CubeIndex<IBASIS,DIM,CUBEBASIS>
  last_wavelet(const CUBEBASIS* basis, const int j)
  {
    assert(j >= basis->j0());
    
    typename CubeIndex<IBASIS,DIM,CUBEBASIS>::type_type e;
    typename CubeIndex<IBASIS,DIM,CUBEBASIS>::translation_type k;
    for (unsigned int i = 0; i < DIM; i++) {
      typename IBASIS::Index lambda(last_wavelet<IBASIS>(basis->bases()[i], j));
      k[i] = lambda.k();
      e[i] = 1;
    }

    return CubeIndex<IBASIS,DIM,CUBEBASIS>(j, e, k, basis);
  }

  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  const int
  first_generator_num(const CUBEBASIS* basis)
  {
    CubeIndex<IBASIS,DIM,CUBEBASIS> ind(first_generator<IBASIS,DIM,CUBEBASIS>(basis, basis->j0()));
    ind.set_number();
    return ind.number();
  }

  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  const int
  last_generator_num(const CUBEBASIS* basis)
  {
    CubeIndex<IBASIS,DIM,CUBEBASIS> ind(last_generator<IBASIS,DIM,CUBEBASIS>(basis, basis->j0()));
    ind.set_number();
    return ind.number();

  }

  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  const int
  first_wavelet_num(const CUBEBASIS* basis, const int j)
  {
    assert(j >= basis->j0());
    CubeIndex<IBASIS,DIM,CUBEBASIS> ind(first_wavelet<IBASIS,DIM,CUBEBASIS>(basis, j));
    ind.set_number();
    return ind.number();

  }

  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  const int 
  last_wavelet_num(const CUBEBASIS* basis, const int j)
  {
    assert(j >= basis->j0());
    CubeIndex<IBASIS,DIM,CUBEBASIS> ind(last_wavelet<IBASIS,DIM,CUBEBASIS>(basis, j));
    ind.set_number();
    return ind.number();

  }


}
