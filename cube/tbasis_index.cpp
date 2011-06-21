// implementation for tbasis_index.h

namespace WaveletTL
{
	template <class IBASIS, unsigned int DIM, class TENSORBASIS>
	TensorIndex<IBASIS,DIM,TENSORBASIS>::TensorIndex(const TENSORBASIS* basis)
	:basis_(basis)
	{
		if (basis_ == 0) {
      		//j_ = 0; // invalid (e and k are initialized by zero automatically)
      		num_ = -1;
    	} else {
    		for (unsigned int i = 0; i < DIM; i++) {
    			j_[i] = basis_->bases()[i]->j0(); // coarsest level;
      			// e_ is zero by default: generator
				k_[i] = basis_->bases()[i]->DeltaLmin();
    		}
      		num_ = 1;
		}
	}

	template <class IBASIS, unsigned int DIM, class TENSORBASIS>
	TensorIndex<IBASIS,DIM,TENSORBASIS>:: TensorIndex(const level_type& j, const type_type& e, const translation_type& k, const TENSORBASIS* basis)
	: basis_(basis), j_(j), e_(e), k_(k)
	{
		MathTL::FixedArray1D<std::map<int, int>,DIM> sizes; // store number of basis elements. Generators on level j0 (0), wavelets on level j0 (1), j0+1 (2), ...
                int uptothislevel(0); // number of basis function below the current level
		int oncurrentlevel(1); // number of base elements on current level j
		int range(0); // determines wether a new entry has to be computed in "sizes"
		level_type currentlevel, j0(basis_->j0());
		type_type currenttype;		
		for (unsigned int i=0; i < DIM; i++) {
                    currentlevel[i]=j0[i];
                    currenttype[i]=0;
                    sizes[i][0] = basis_->bases()[i]->Deltasize(j0[i]); // Number of generators on level j0
                    sizes[i][1] = basis_->bases()[i]->Nablasize(j0[i]); // Number of Wavelets on level j0
		}
		// iterate over all level up to j_ and add up number of basis functions up to the level
                if (currentlevel != j_ || currenttype != e_)
                while (true)
		{
                    // compute number of functions on this level
                    oncurrentlevel = 1;
                    for (unsigned int i = 0; i < DIM; i++)
                    {
                        oncurrentlevel *= sizes[i][currentlevel[i]+currenttype[i]-j0[i]];
                    }
                    uptothislevel += oncurrentlevel;
                    // increase index = (currentlevel,currenttype)
                    // "small loop"
                    // iterate over all combinations of generators/wavelets for all dimensions with currentlevel[i]=j0[i]
                    // this looks like binary addition: (in 3 Dim:) gwg is followed by gww (g=0,w=1)
                    bool done = true;
                    for (int i(DIM-1); i >= 0; i--)
                    {
                        // find first position on level j0
                        if (currentlevel[i] == j0[i])
                        {
                            if (currenttype[i] == 1)
                            {
                                currenttype[i]=0;
                            } else
                            {
                                currenttype[i]=1;
                                done = false;
                                break;
                            }
                        }
                    }
                    // done == true means that all components with currentlevel[i]=j0[i] were wavelets.
                    // the level j has to be increased.
                    // iterate "big loop", meaning: "currentlevel++"

                    if (done == true)
                    {
                    for (int i(DIM-1); i >= 0; i--)
                    {
                        if (i != 0)
                        {
                            if (currentlevel[i] != j0[i])
                            {
                                // increase left neighbor
                                currentlevel[i-1]=currentlevel[i-1]+1;
                                if (currentlevel[i-1]-j0[i] == range) sizes[i-1][range+1]=basis ->bases()[i]->Nablasize(currentlevel[i-1]); // if needed compute and store new size information
                                currenttype[i-1]=1;
                                int temp = currentlevel[i]-j0[i];
                                currentlevel[i]=j0[i];
                                currenttype[i]=0;
                                currentlevel[DIM-1]=j0[DIM-1]+temp-1;
                                currenttype[DIM-1]= (temp == 1?0:1);
  				break;
                            }
  			} else // i == 0. "big loop" arrived at the last index. We have to increase the level!
  			{
                            range = range +1;
                            if (DIM == 1)
                            {
                                currenttype[i] = 1;
                                currentlevel[i]=currentlevel[i]+1;
                                sizes[0][range+1]=basis ->bases()[0]->Nablasize(currentlevel[0]); // if needed compute and store new size information
                            }
                            else
                            {
                                //currentlevel[DIM-1]=j0[DIM-1]+currentlevel[0]-j0[0]+1; currenttype[DIM-1]=1;
                                currentlevel[DIM-1]=j0[DIM-1]+range; currenttype[DIM-1]=1;
				currenttype[0]=0; currentlevel[0]=j0[0];
                                sizes[DIM-1][range+1]=basis ->bases()[DIM-1]->Nablasize(currentlevel[DIM-1]); // if needed compute and store new size information
                            }
  			break; // unnoetig, da i==0 gilt.
  			}
                    } // end of "big loop"
                    }
                    if (currentlevel == j_ && currenttype == e_) break;
		}
                // determine number of the function described by k in (currentlevel,currenttype)
		translation_type ktemp; // we only use that this is a multiindex
		for (unsigned int i = 0; i < DIM; i++)
                    if (e_[i] == 0) ktemp[i]=k_[i]-basis_->bases()[i]->DeltaLmin();
                    else ktemp[i]=k_[i]-basis_->bases()[i]->Nablamin();
                oncurrentlevel = ktemp[0];
		for (unsigned int i=1; i <= DIM-1; i++)
		{
                    oncurrentlevel *= sizes[i][currentlevel[i]+currenttype[i]-j0[i]];
                    oncurrentlevel += ktemp[i];
		}
		num_ =uptothislevel+oncurrentlevel;
	}

	template <class IBASIS, unsigned int DIM, class TENSORBASIS>
	TensorIndex<IBASIS, DIM, TENSORBASIS>::TensorIndex(const TensorIndex& lambda)
	: basis_(lambda.basis_), j_(lambda.j_), e_(lambda.e_), k_(lambda.k_), num_(lambda.num_) {}

	template <class IBASIS, unsigned int DIM, class TENSORBASIS>
	TensorIndex<IBASIS, DIM, TENSORBASIS>::TensorIndex(const TensorIndex* lambda)
	: basis_(lambda->basis_), j_(lambda->j_), e_(lambda->e_), k_(lambda->k_), num_(lambda->num_) {}

	template <class IBASIS, unsigned int DIM, class TENSORBASIS>
	TensorIndex<IBASIS, DIM, TENSORBASIS>::TensorIndex(const int number, const TENSORBASIS* basis)
	: basis_(basis), num_(number)
	{
          /* 
            This implementation for DIM = 2,3 assumes that Nablasize(j) = Deltasize(j0) + sum_{l=0}^{j-j0} 2^l * Nablasize(j0)
            and uses a lot of div and mod
          */
          if(DIM == 2)
          {
          level_type j0 = basis->j0();
          int Nabla1 = basis->bases()[0]->Nablasize(j0[0]);
          int Nabla2 = basis->bases()[1]->Nablasize(j0[1]);
          int Delta1 = basis->bases()[0]->Deltasize(j0[0]);
          int Delta2 = basis->bases()[1]->Deltasize(j0[1]);
          int act_num = num_;
          int i = 0; 
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
            act_num -= Delta1*Delta2;
            // computes the outer level
            while(true)
            {
              int tmp = (Delta1*Nabla2 + Delta2*Nabla1 + (i+1)*Nabla1*Nabla2)<<i;
              if(act_num < tmp)
                break;
              act_num -= tmp;
              i++;
            }
            if(act_num < (Delta1*Nabla2)<<i)
            {
              j_[0] = j0[0];
              j_[1] = j0[1] + i;
              e_[0]=0;
              e_[1]=1;
              k_[0] = basis->bases()[0]->DeltaLmin() + act_num/(Nabla2<<i);
              k_[1] = basis->bases()[1]->Nablamin() + (act_num % (Nabla2<<i));
            }
            else if ((act_num - ((Delta1*Nabla2)<<i)) / ((Nabla1*Nabla2)<<i) < i)
            {
              unsigned int tmp2 = (act_num - ((Delta1*Nabla2)<<i)) / ((Nabla1*Nabla2)<<i);
              j_[0] = j0[0] + tmp2;
              j_[1] = j0[1] + i - tmp2;
              e_[0]=1;
              e_[1]=1;
              act_num -= (Delta1*Nabla2 + tmp2 * (Nabla1*Nabla2))<<i;
              k_[0] = basis->bases()[0]->Nablamin() + act_num/(Nabla2<<(i-tmp2));
              k_[1] = basis->bases()[1]->Nablamin() + (act_num % (Nabla2<<(i-tmp2)));
            }
            else if (act_num < (Delta1*Nabla2 + (i) * Nabla1*Nabla2 + Delta2*Nabla1)<<i ) 
            {
              j_[0] = j0[0] + i;
              j_[1] = j0[1];
              e_[0]=1;
              e_[1]=0;
              act_num -= (Delta1*Nabla2 + (i) * (Nabla1*Nabla2))<<i;
              k_[0] = basis->bases()[0]->Nablamin() + act_num/Delta2;
              k_[1] = basis->bases()[1]->DeltaLmin() + (act_num % Delta2);
            }
            else
            {
              j_[0] = j0[0] + i;
              j_[1] = j0[1];
              e_[0]=1;
              e_[1]=1;
              act_num -= (Delta1*Nabla2 + (i) * (Nabla1*Nabla2) + Delta2*Nabla1)<<i;
              k_[0] = basis->bases()[0]->Nablamin() + act_num/Nabla2;
              k_[1] = basis->bases()[1]->Nablamin() + (act_num % Nabla2);
            }
          }
          }
          else if(DIM == 3)
          {
          level_type j0 = basis->j0();
          int Nabla1 = basis->bases()[0]->Nablasize(j0[0]);
          int Nabla2 = basis->bases()[1]->Nablasize(j0[1]);
          int Nabla3 = basis->bases()[2]->Nablasize(j0[2]);
          int Delta1 = basis->bases()[0]->Deltasize(j0[0]);
          int Delta2 = basis->bases()[1]->Deltasize(j0[1]);
          int Delta3 = basis->bases()[2]->Deltasize(j0[2]);
          int act_num = num_;
          int i = 0; 
          if (act_num < Delta1*Delta2*Delta3) 
          {
            j_ = j0;
            e_[0]=0;
            e_[1]=0;
            e_[2]=0;
            k_[0] = basis->bases()[0]->DeltaLmin() + (act_num/(Delta3*Delta2));
            k_[1] = basis->bases()[1]->DeltaLmin() + (act_num/Delta3) % Delta2;
            k_[2] = basis->bases()[2]->DeltaLmin() + (act_num % Delta3);
          }
          else
          {
            act_num -= Delta1*Delta2*Delta3;
            while(true)
            {
              int tmp = (Delta1*Nabla2*Delta3 + Delta2*Nabla1*Delta3 + Delta2*Nabla3*Delta1 
                        + (i+1) * (Nabla1*Nabla2*Delta3 + Nabla1*Nabla3*Delta2 + Nabla3*Nabla2*Delta1)
                        + ((i+1)*(i+2)/2)*Nabla1*Nabla2*Nabla3) <<i;
              if(act_num < tmp)
                break;
              act_num -= tmp;
              i++;
            }
            if(i == 0)
            {
              if (act_num < Delta1*Delta2*Nabla3)
              {
                j_ = j0;
                e_[0] = 0;
                e_[1] = 0;
                e_[2] = 1;
                k_[0] = basis->bases()[0]->DeltaLmin() + (act_num/(Nabla3*Delta2));
                k_[1] = basis->bases()[1]->DeltaLmin() + (act_num/Nabla3) % Delta2;
                k_[2] = basis->bases()[2]->Nablamin() + (act_num % Nabla3);
              }
              else
              {
              act_num -= Delta1*Delta2*Nabla3;
              if(act_num < Delta1*Nabla2*Delta3)
              {
                j_ = j0;
                e_[0] = 0;
                e_[1] = 1;
                e_[2] = 0;
                k_[0] = basis->bases()[0]->DeltaLmin() + (act_num/(Delta3*Nabla2));
                k_[1] = basis->bases()[1]->Nablamin() + (act_num/Delta3) % Nabla2;
                k_[2] = basis->bases()[2]->DeltaLmin() + (act_num % Delta3);
              }
              else
              {
              act_num -= Delta1*Nabla2*Delta3;
              if(act_num < Delta1*Nabla2*Nabla3)
              {
                j_ = j0;
                e_[0] = 0;
                e_[1] = 1;
                e_[2] = 1;
                k_[0] = basis->bases()[0]->DeltaLmin() + (act_num/(Nabla3*Nabla2));
                k_[1] = basis->bases()[1]->Nablamin() + (act_num/Nabla3) % Nabla2;
                k_[2] = basis->bases()[2]->Nablamin() + (act_num % Nabla3);
              }
              else
              {
              act_num -= Delta1*Nabla2*Nabla3;
              if(act_num < Nabla1*Delta2*Delta3)
              {
                j_ = j0;
                e_[0] = 1;
                e_[1] = 0;
                e_[2] = 0;
                k_[0] = basis->bases()[0]->Nablamin() + (act_num/(Delta3*Delta2));
                k_[1] = basis->bases()[1]->DeltaLmin() + (act_num/Delta3) % Delta2;
                k_[2] = basis->bases()[2]->DeltaLmin() + (act_num % Delta3);
              }
              else
              {
              act_num -= Nabla1*Delta2*Delta3;
              if(act_num < Nabla1*Delta2*Nabla3)
              {
                j_ = j0;
                e_[0] = 1;
                e_[1] = 0;
                e_[2] = 1;
                k_[0] = basis->bases()[0]->Nablamin() + (act_num/(Nabla3*Delta2));
                k_[1] = basis->bases()[1]->DeltaLmin() + (act_num/Nabla3) % Delta2;
                k_[2] = basis->bases()[2]->Nablamin() + (act_num % Nabla3);
              }
              else
              {
              act_num -= Nabla1*Delta2*Nabla3;
              if(act_num < Nabla1*Nabla2*Delta3)
              {
                j_ = j0;
                e_[0] = 1;
                e_[1] = 1;
                e_[2] = 0;
                k_[0] = basis->bases()[0]->Nablamin() + (act_num/(Delta3*Nabla2));
                k_[1] = basis->bases()[1]->Nablamin() + (act_num/Delta3) % Nabla2;
                k_[2] = basis->bases()[2]->DeltaLmin() + (act_num % Delta3);
              }
              else
              {
              act_num -= Nabla1*Nabla2*Delta3;
                j_ = j0;
                e_[0] = 1;
                e_[1] = 1;
                e_[2] = 1;
                k_[0] = basis->bases()[0]->Nablamin() + (act_num/(Nabla3*Nabla2));
                k_[1] = basis->bases()[1]->Nablamin() + (act_num/Nabla3) % Nabla2;
                k_[2] = basis->bases()[2]->Nablamin() + (act_num % Nabla3);
              }}}}}}
            }
            else
            {            
              if(act_num < ((Delta1+Nabla1) * (Delta2 + Nabla2) * Nabla3)<<i)
              {
                if(act_num < (Delta1 * Delta2 * Nabla3)<<i)
                {
                  j_[0] = j0[0];
                  j_[1] = j0[1];
                  j_[2] = j0[2] + i;
                  e_[0]=0;
                  e_[1]=0;
                  e_[2]=1;
                  k_[0] = basis->bases()[0]->DeltaLmin() + act_num/((Delta2*Nabla3)<<i);
                  k_[1] = basis->bases()[1]->DeltaLmin() + (act_num/(Nabla3<<i)) % Delta2;
                  k_[2] = basis->bases()[2]->Nablamin() + (act_num % (Nabla3<<i));
                }
                else
                {
                act_num -= (Delta1 * Delta2 * Nabla3)<<i;
                if(act_num < (Delta1 * Nabla2 * Nabla3)<<i)
                {
                  j_[0] = j0[0];
                  j_[1] = j0[1];
                  j_[2] = j0[2] + i;
                  e_[0]=0;
                  e_[1]=1;
                  e_[2]=1;
                  k_[0] = basis->bases()[0]->DeltaLmin() + act_num/((Nabla2*Nabla3)<<i);
                  k_[1] = basis->bases()[1]->Nablamin() + (act_num/(Nabla3<<i)) % Nabla2;
                  k_[2] = basis->bases()[2]->Nablamin() + (act_num % (Nabla3<<i));
                }
                else
                {
                act_num -= (Delta1 * Nabla2 * Nabla3)<<i;
                if(act_num < (Nabla1 * Delta2 * Nabla3)<<i)
                {
                  j_[0] = j0[0];
                  j_[1] = j0[1];
                  j_[2] = j0[2] + i;
                  e_[0]=1;
                  e_[1]=0;
                  e_[2]=1;
                  k_[0] = basis->bases()[0]->Nablamin() + act_num/((Delta2*Nabla3)<<i);
                  k_[1] = basis->bases()[1]->DeltaLmin() + (act_num/(Nabla3<<i)) % Delta2;
                  k_[2] = basis->bases()[2]->Nablamin() + (act_num % (Nabla3<<i));
                }
                else
                {
                act_num -= (Nabla1 * Delta2 * Nabla3)<<i;
                  j_[0] = j0[0];
                  j_[1] = j0[1];
                  j_[2] = j0[2] + i;
                  e_[0]=1;
                  e_[1]=1;
                  e_[2]=1;
                  k_[0] = basis->bases()[0]->Nablamin() + act_num/((Nabla2*Nabla3)<<i);
                  k_[1] = basis->bases()[1]->Nablamin() + (act_num/(Nabla3<<i)) % Nabla2;
                  k_[2] = basis->bases()[2]->Nablamin() + (act_num % (Nabla3<<i));
                }}}
              }
              else if(((act_num - (((Delta1+Nabla1) * (Delta2 + Nabla2) * Nabla3)<<i))/(((Delta1+Nabla1)*Nabla2*Nabla3)<<i)) < (i-1))
              {
                unsigned int tmp2 = (act_num - (((Delta1+Nabla1) * (Delta2 + Nabla2) * Nabla3)<<i))/(((Delta1+Nabla1)*Nabla2*Nabla3)<<i);
                j_[0] = j0[0];
                j_[1] = j0[1] + tmp2+1;
                j_[2] = j0[2] + i - tmp2 -1;
                act_num -= ((Delta1+Nabla1) * (Delta2 + Nabla2) * Nabla3 + tmp2 * (Delta1+Nabla1)*Nabla2*Nabla3)<<i;
                if(act_num < Delta1*Nabla2*Nabla3 <<i)
                {
                  e_[0] = 0;
                  k_[0] = basis->bases()[0]->DeltaLmin() + act_num/((Nabla2*Nabla3)<<i);
                }
                else
                {
                  e_[0] = 1;
                  act_num -= Delta1*Nabla2*Nabla3 <<i;
                  k_[0] = basis->bases()[0]->Nablamin() + act_num/((Nabla2*Nabla3)<<i);
                }
                e_[1] = 1;
                e_[2] = 1;
                k_[1] = basis->bases()[1]->Nablamin() + (act_num/(Nabla3<<(i-tmp2-1))) % (Nabla2 <<(tmp2+1));
                k_[2] = basis->bases()[2]->Nablamin() + (act_num % (Nabla3<<(i-tmp2-1)));
              } 
              else if((act_num -(((Delta1+Nabla1) * (Delta2 + Nabla2) * Nabla3 + (i-1) * (Delta1+Nabla1)*Nabla2*Nabla3)<<i)) 
                          < ((Nabla2*(Delta1+Nabla1)*(Delta3+Nabla3))<<i ) )
              {
                j_[0] = j0[0];
                j_[1] = j0[1] + i;
                j_[2] = j0[2];
                e_[1] = 1;
                act_num -= ((Delta1+Nabla1) * (Delta2 + Nabla2) * Nabla3 + (i-1)*((Delta1+Nabla1)*Nabla2*Nabla3))<<i;
                if(act_num < (Delta1*Nabla2*Delta3) <<i)
                {
                  e_[0] = 0;
                  e_[2] = 0;
                  k_[0] = basis->bases()[0]->DeltaLmin() + act_num/((Nabla2*Delta3)<<i);
                  k_[1] = basis->bases()[1]->Nablamin() + (act_num/Delta3) % (Nabla2 <<i);
                  k_[2] = basis->bases()[2]->DeltaLmin() + (act_num % Delta3);
                }
                else 
                {
                act_num -= (Delta1*Nabla2*Delta3) <<i;
                if(act_num < (Delta1*Nabla2*Nabla3) <<i)
                {
                  e_[0] = 0;
                  e_[2] = 1;
                  k_[0] = basis->bases()[0]->DeltaLmin() + act_num/((Nabla2*Nabla3)<<i);
                  k_[1] = basis->bases()[1]->Nablamin() + (act_num/Nabla3) % (Nabla2 <<i);
                  k_[2] = basis->bases()[2]->Nablamin() + (act_num % Nabla3);
                }
                else 
                {
                act_num -= (Delta1*Nabla2*Nabla3) <<i;
                if(act_num < (Nabla1*Nabla2*Delta3) <<i)
                {
                  e_[0] = 1;
                  e_[2] = 0;
                  k_[0] = basis->bases()[0]->Nablamin() + act_num/((Nabla2*Delta3)<<i);
                  k_[1] = basis->bases()[1]->Nablamin() + (act_num/Delta3) % (Nabla2 <<i);
                  k_[2] = basis->bases()[2]->DeltaLmin() + (act_num % Delta3);
                }
                else 
                {
                act_num -= (Nabla1*Nabla2*Delta3) <<i;
                  e_[0] = 1;
                  e_[2] = 1;
                  k_[0] = basis->bases()[0]->Nablamin() + act_num/((Nabla2*Nabla3)<<i);
                  k_[1] = basis->bases()[1]->Nablamin() + (act_num/Nabla3) % (Nabla2 <<i);
                  k_[2] = basis->bases()[2]->Nablamin() + (act_num % Nabla3);
                }}}
              }
              else
              {
                act_num -= ((Delta1+Nabla1) * (Delta2 + Nabla2) * Nabla3 
                           + (i-1) * (Delta1+Nabla1)*Nabla2*Nabla3 
                           + Nabla2*(Delta1+Nabla1)*(Delta3+Nabla3))<<i;
                int l = 1;
                while(true)
                {
                  int tmp = (Nabla1*Nabla2*Delta3 + Nabla1*Delta2*Nabla3 + (i+1-l)*Nabla1*Nabla2*Nabla3)<<i;
                  if (l == i)
                    tmp += Nabla1*Delta2*Delta3<<i;
                  if (act_num < tmp)
                    break;
                  act_num -= tmp;
                  l++;
                }
                j_[0] = j0[0] + l;
                e_[0] = 1;
                if(l < i)
                {
                  if(act_num < (Nabla1*Delta2*Nabla3)<<i)
                  {
                    j_[1] = j0[1];
                    j_[2] = j0[2] + i - l;
                    e_[1] = 0;
                    e_[2] = 1;
                    k_[0] = basis->bases()[0]->Nablamin() + act_num/((Delta2*Nabla3)<<(i-l));
                    k_[1] = basis->bases()[1]->DeltaLmin() + (act_num/(Nabla3<<(i-l))) % Delta2;
                    k_[2] = basis->bases()[2]->Nablamin() + (act_num % (Nabla3<<(i-l)));
                  }
                  else if( ((act_num - ((Nabla1*Delta2*Nabla3)<<i))/((Nabla1*Nabla2*Nabla3)<<i)) < (i-l) )
                  {
                    unsigned int tmp2 = ((act_num - ((Nabla1*Delta2*Nabla3)<<i))/((Nabla1*Nabla2*Nabla3)<<i));
                    j_[1] = j0[1] + tmp2;
                    j_[2] = j0[2] + i - l - tmp2;
                    e_[1] = 1;
                    e_[2] = 1;
                    act_num -= (Nabla1*Delta2*Nabla3 + tmp2*Nabla1*Nabla2*Nabla3)<<i;
                    k_[0] = basis->bases()[0]->Nablamin() + act_num/((Nabla2*Nabla3)<<(i-l));
                    k_[1] = basis->bases()[1]->Nablamin() + (act_num/(Nabla3<<(i-l-tmp2))) % (Nabla2<<tmp2);
                    k_[2] = basis->bases()[2]->Nablamin() + (act_num % (Nabla3<<(i-l-tmp2)));
                  } 
                  else
                  {
                    act_num -= (Nabla1*Delta2*Nabla3 + (i-l)*Nabla1*Nabla2*Nabla3)<<i;
                    if(act_num < Nabla1*Nabla2*Delta3<<i)
                    {
                      j_[1] = j0[1] + i - l;
                      j_[2] = j0[2];
                      e_[1] = 1;
                      e_[2] = 0;
                      k_[0] = basis->bases()[0]->Nablamin() + act_num/((Nabla2*Delta3)<<(i-l));
                      k_[1] = basis->bases()[1]->Nablamin() + (act_num/Delta3) % (Nabla2<<(i-l));
                      k_[2] = basis->bases()[2]->DeltaLmin() + (act_num % Delta3);
                    }
                    else
                    {
                      j_[1] = j0[1] + i - l;
                      j_[2] = j0[2];
                      e_[1] = 1;
                      e_[2] = 1;
                      act_num -= Nabla1*Nabla2*Delta3<<i;
                      k_[0] = basis->bases()[0]->Nablamin() + act_num/((Nabla2*Nabla3)<<(i-l));
                      k_[1] = basis->bases()[1]->Nablamin() + (act_num/Nabla3) % (Nabla2<<(i-l));
                      k_[2] = basis->bases()[2]->Nablamin() + (act_num % Nabla3);
                    }
                  }
                }
                else
                {
                  if(act_num < (Nabla1*Delta2*Delta3)<<i)
                  {
                    j_[1] = j0[1];
                    j_[2] = j0[2];
                    e_[1] = 0;
                    e_[2] = 0;
                    k_[0] = basis->bases()[0]->Nablamin() + act_num/(Delta2*Delta3);
                    k_[1] = basis->bases()[1]->DeltaLmin() + (act_num/Delta3) % Delta2;
                    k_[2] = basis->bases()[2]->DeltaLmin() + (act_num % Delta3);
                  }
                  else 
                  {
                  act_num -= Nabla1*Delta2*Delta3<<i;
                  if(act_num < (Nabla1*Delta2*Nabla3)<<i)
                  {
                    j_[1] = j0[1];
                    j_[2] = j0[2];
                    e_[1] = 0;
                    e_[2] = 1;
                    k_[0] = basis->bases()[0]->Nablamin() + act_num/(Delta2*Nabla3);
                    k_[1] = basis->bases()[1]->DeltaLmin() + (act_num/Nabla3) % Delta2;
                    k_[2] = basis->bases()[2]->Nablamin() + (act_num % Nabla3);
                  }
                  else
                  {
                  act_num -= Nabla1*Delta2*Nabla3<<i;
                  if(act_num < (Nabla1*Nabla2*Delta3)<<i)
                  {
                    j_[1] = j0[1];
                    j_[2] = j0[2];
                    e_[1] = 1;
                    e_[2] = 0;
                    k_[0] = basis->bases()[0]->Nablamin() + act_num/(Nabla2*Delta3);
                    k_[1] = basis->bases()[1]->Nablamin() + (act_num/Delta3) % Nabla2;
                    k_[2] = basis->bases()[2]->DeltaLmin() + (act_num % Delta3);
                  }
                  else
                  {
                  act_num -= Nabla1*Nabla2*Delta3<<i;
                    j_[1] = j0[1];
                    j_[2] = j0[2];
                    e_[1] = 1;
                    e_[2] = 1;
                    k_[0] = basis->bases()[0]->Nablamin() + act_num/(Nabla2*Nabla3);
                    k_[1] = basis->bases()[1]->Nablamin() + (act_num/Nabla3) % Nabla2;
                    k_[2] = basis->bases()[2]->Nablamin() + (act_num % Nabla3);
                  }}}
                }
              }
            }
          }
          }
          else //DIM >3
          {

                MathTL::FixedArray1D<std::map<int, int>,DIM> sizes; // Store number of basis elements. Generators on level j0 (0), wavelets on level j0 (1), j0+1 (2), ...
		int remains = number+1; // numbering begins at 0
		int oncurrentlevel(1); // number of base elements on current level j
		int range(0); // number of level steps we have climbed so far. determines wether a new entry has to be computed in "sizes"
		level_type currentlevel, j0(basis_->j0());
		type_type currenttype;
		currentlevel = j0;
		for (unsigned int i=0; i < DIM; i++) {
			currenttype[i]=0;
			sizes[i][0] = basis_->bases()[i]->Deltasize(j0[i]); // Number of generators on level j0
                        sizes[i][1] = basis_->bases()[i]->Nablasize(j0[i]); // N o Wavelets
                        oncurrentlevel *= sizes[i][0];
		}

                // iterate over all levels. Add up number of basis functions till looked for level is reached
		while (remains > oncurrentlevel) // break if we are at the right level
		{
                    // else substract number of basis functions on current_index=(currentlevel,currentindex) and increase current_index
                    remains -= oncurrentlevel;
                    // increase index = (currentlevel,currenttype)
                    // "small loop" "currenttype++" (currentlevel fest)
                    // iterate over all combinations of generators/wavelets for all dimensions with currentlevel[i]=j0[i]
                    // this looks like binary addition: (in 3 Dim:) gwg is followed by gww (g=0,w=1)
                    bool done = true;
                    for (int i(DIM-1); i >= 0; i--)
                    {
                        // find first position on level j0
                        if (currentlevel[i] == j0[i])
                        {
                            if (currenttype[i] == 1)
                            {
                                currenttype[i]=0;
                            } else
                            {
                                currenttype[i]=1;
                                done = false;
                                break;
                            }
                        }
                    }
                    // done == true means that all components with currentlevel[i]=j0[i] were wavelets.
                    // the level j has to be increased.
                    // iterate "big loop", meaning: "currentlevel++"
                    if (done == true)
                    {
                        for (int i(DIM-1); i >= 0; i--)
                        {
                            if (i != 0)
                            {
                                if (currentlevel[i] != j0[i])
                                {
                                    // increase left neighbor
                                    currentlevel[i-1]=currentlevel[i-1]+1;
                                    if (currentlevel[i-1]-j0[i-1] == range)
                                    {
                                      sizes[i-1][range+1]=basis ->bases()[i]->Nablasize(currentlevel[i-1]); // if needed compute and store new size information
                                    }
                                    currenttype[i-1]=1;
                                    int temp = currentlevel[i]-j0[i];
                                    currentlevel[i]=j0[i];
                                    currenttype[i]=0;
                                    currentlevel[DIM-1]=j0[DIM-1]+temp-1;
                                    currenttype[DIM-1]= (temp == 1?0:1);
                                    break;
                                }

                            } else // i == 0. "big loop" arrived at the last index. We have to increase range!
                            {
                                range = range +1;
                                if (DIM == 1)
                                {
                                    currenttype[i] = 1;
                                    currentlevel[i]=currentlevel[i]+1;
                                    sizes[0][range+1]=basis ->bases()[0]->Nablasize(currentlevel[0]); // if needed compute and store new size information
                                    //double temp1(currenttype[i]), temp2(currentlevel[i]), temp3(basis ->bases()[0]->Nablasize(currentlevel[0])), temp4(sizes[0][range+1]);
                                    //cout << temp1 << " " << temp2 << " " << temp3 << " " << temp4 << endl;
                                }
                                else
                                {
                                    //currentlevel[DIM-1]=j0[DIM-1]+currentlevel[0]-j0[0]+1; currenttype[DIM-1]=1;
                                    currentlevel[DIM-1]=j0[DIM-1]+range; currenttype[DIM-1]=1;
                                    for(int i(DIM-2); i>=0; i--)
                                    {
                                      currenttype[i]=0; currentlevel[i]=j0[i];
                                    }
                                    sizes[DIM-1][range+1]=basis ->bases()[DIM-1]->Nablasize(currentlevel[DIM-1]); // if needed compute and store new size information
                                }
                            break; // unnoetig, da i==0 gilt.
                            }
                        } // end of "big loop"
                    }
                    // compute number of functions on this level
                    oncurrentlevel = 1;
                    for (unsigned int i = 0; i < DIM; i++)
                    {
                        oncurrentlevel *= sizes[i][currentlevel[i]+currenttype[i]-j0[i]];
                    }
		} // end of while

                // determine k corresponding to the number of the basis function given by "remains" (on level (currentlevel,currenttype) )
		unsigned int modul;                
		j_ = currentlevel;
		e_ = currenttype;
                remains -= 1; // numbering begins at 0
		for (int i = DIM-1; i > 0; i--)
		{
			modul = sizes[i][currentlevel[i]+currenttype[i]-j0[i]];
			k_[i]= remains % modul+(e_[i] == 0 ? basis->bases()[i]->DeltaLmin():basis->bases()[i]->Nablamin());
			remains = remains / modul;
		}
                k_[0]=remains +(e_[0] == 0 ? basis->bases()[0]->DeltaLmin():basis->bases()[0]->Nablamin());

          } //end DIM >3
        }

	template <class IBASIS, unsigned int DIM, class TENSORBASIS>
  	TensorIndex<IBASIS,DIM,TENSORBASIS>&
  	TensorIndex<IBASIS,DIM,TENSORBASIS>::operator = (const TensorIndex<IBASIS,DIM,TENSORBASIS>& lambda)
  	{
    	j_ = lambda.j();
    	e_ = lambda.e();
    	k_ = lambda.k();
    	basis_ = lambda.basis();
    	num_ = lambda.number();
    	return *this;
  	}

	template <class IBASIS, unsigned int DIM, class TENSORBASIS>
	bool
	TensorIndex<IBASIS, DIM, TENSORBASIS>::operator == (const TensorIndex<IBASIS, DIM, TENSORBASIS>& lambda) const
	{
            return (j_ == lambda.j() &&
                    e_ == lambda.e() &&
                    k_ == lambda.k());
	}

        // PERFORMANCE: eventuell die vielen Aufrufe von bases_->j0() durch Hilfsvariable auf einen reduzieren
	template <class IBASIS, unsigned int DIM, class TENSORBASIS>
  	TensorIndex<IBASIS,DIM,TENSORBASIS>&
  	TensorIndex<IBASIS,DIM,TENSORBASIS>::operator ++ ()
  	{
            if (num_ == -1) return *this;
            level_type j0(basis_->j0());
            num_++;
            // determine next translation index
            bool jplusplus = false;
            for (int i = DIM-1; i >= 0; i--) {
      		const int last_index = (e_[i] == 0 ? basis_->bases()[i]->DeltaRmax(j_[i])
                                                   : basis_->bases()[i]->Nablamax(j_[i]));
      		if (k_[i] == last_index)
                {
                    k_[i] = (e_[i] == 0 ? basis_->bases()[i]->DeltaLmin()
	 				: basis_->bases()[i]->Nablamin());
                    jplusplus = (i == 0);
      		} else
                {
                    ++k_[i];
                    break;
      		}
            }
            if (jplusplus == false) return *this;
            // else: determine next level index
            // "small loop" "e_++" (j_ fest)
            // iterate over all combinations of generators/wavelets for all dimensions with j_[i]=j0[i]
            // this looks like binary addition: (in 3 Dim:) gwg is followed by gww (g=0,w=1)
            bool done = true;
            for (int i(DIM-1); i >= 0; i--)
            {
                // find first position on level j0
                if (j_[i] == j0[i])
                {
                    if (e_[i] == 1)
                    {
                        e_[i]=0;
                        k_[i]=basis_->bases()[i]->DeltaLmin();
                    } else
                    {
                        e_[i]=1;
                        k_[i]=basis_->bases()[i]->Nablamin();
                        done = false;
                        break;
                    }
                }
            }
            // done == true bedeutet, dass alle Komponenten auf level j0() wavelets waren.
            // "big loop" "j_++"
            if (done == true)
            {
                for (int i(DIM-1); i >= 0; i--)
                {
                    if (i != 0)
                    {
                        if (j_[i] != j0[i])
                        {
                            // increase left neighbor
                            j_[i-1]=j_[i-1]+1;
                            e_[i-1]=1;
                            k_[i-1]=basis_->bases()[i-1]->Nablamin();
                            int temp = j_[i]-j0[i];
                            j_[i]=j0[i];
                            e_[i]=0;
                            k_[i]=basis_->bases()[i]->DeltaLmin();
                            j_[DIM-1]=j0[DIM-1]+temp-1;
                            e_[DIM-1]= (temp == 1?0:1);
                            k_[DIM-1]= (temp == 1?basis_->bases()[i]->DeltaLmin():basis_->bases()[i]->Nablamin());
                            break;
                        }
                    } else // i == 0. "big loop" arrived at the last index. We have to increase the level!
                    {
                        if (DIM == 1)
                        {
                            e_[i] = 1; // diese Zeile erfüllt nur in der allerersten Iteration einen Zweck
                            k_[i]=basis_->bases()[i]->Nablamin(); // diese Zeile erfüllt nur in der allerersten Iteration einen Zweck
                            j_[i]=j_[i]+1;
                        }
                        else
                        {
                            j_[DIM-1]=j0[DIM-1]+j_[0]-j0[0]+1;
                            e_[DIM-1]=1;
                            k_[DIM-1]=basis_->bases()[i]->Nablamin();
                            j_[0]=j0[0];
                            e_[0]=0;
                            k_[0]=basis_->bases()[i]->DeltaLmin();
                        }
                        break; // unnoetig, da i==0 gilt.
                    }
                } 
            } // end of "big loop"
            return *this;
         
  	}

	template <class IBASIS, unsigned int DIM, class TENSORBASIS>
	bool
	TensorIndex<IBASIS,DIM,TENSORBASIS>::operator < (const TensorIndex& lambda) const
	{
            // Ordering primary by level j as in MultiIndex
            // (ordering of \N^dim, that is the distance from 0, that is the same as
            // ordering first by 1-norm of j and in the case of equal norms lexicographical in j.)
            // secondly and tertiaryly lexicographical in e and k.
            return ( multi_degree(j_) < multi_degree(lambda.j()) ||
                    (multi_degree(j_) == multi_degree(lambda.j()) &&
                     (j_ < lambda.j() ||
                      (j_ == lambda.j() && 
                       (e_.lex(lambda.e()) ||
                        (e_ == lambda.e() && k_.lex(lambda.k()))
                       )
                      )
                     )
                    )
                   );
	}

	template <class IBASIS, unsigned int DIM, class TENSORBASIS>
	TensorIndex<IBASIS,DIM,TENSORBASIS>
	first_generator(const TENSORBASIS* basis)
	{
    	return TensorIndex<IBASIS,DIM,TENSORBASIS>(0, basis);
	}

	template <class IBASIS, unsigned int DIM, class TENSORBASIS>
	TensorIndex<IBASIS,DIM,TENSORBASIS>
	last_generator(const TENSORBASIS* basis)
	{
		typename TensorIndex<IBASIS,DIM,TENSORBASIS>::level_type j;
		typename TensorIndex<IBASIS,DIM,TENSORBASIS>::type_type e;
		typename TensorIndex<IBASIS,DIM,TENSORBASIS>::translation_type k;
		j = basis->j0();
		for (unsigned int i = 0; i < DIM; i++)
		{
			e[i] = 0;
			k[i] = basis->bases()[i]->DeltaRmax(j[i]);
		}
		return TensorIndex<IBASIS,DIM,TENSORBASIS>(j, e, k, basis);
	}

  	template <class IBASIS, unsigned int DIM, class TENSORBASIS>
  	TensorIndex<IBASIS,DIM,TENSORBASIS>
  	first_wavelet(const TENSORBASIS* basis, const typename TensorIndex<IBASIS,DIM,TENSORBASIS>::level_type j)
  	{
            assert(multi_degree(j) > multi_degree(basis->j0())
                    || ((multi_degree(j) == multi_degree(basis->j0())) && (basis->j0() <= j))
                  );
            typename TensorIndex<IBASIS,DIM,TENSORBASIS>::type_type e;
            typename TensorIndex<IBASIS,DIM,TENSORBASIS>::translation_type k;
            bool first_level = true;
            for (unsigned int i = 0; i < DIM; i++) {
                if (j[i] == basis->bases()[i]->j0())
                {
                    e[i] = 0;
                    k[i] = basis->bases()[i]->DeltaLmin();
                } else
                {
                    e[i] = 1;
                    k[i] = basis->bases()[i]->Nablamin();
                    first_level = false;
                }
            }
            if (first_level == true)
            {
                e[DIM-1] = 1;
                k[DIM-1] = basis->bases()[DIM-1]->Nablamin();
            }
            return TensorIndex<IBASIS,DIM,TENSORBASIS>(j, e, k, basis);
        }


        template <class IBASIS, unsigned int DIM, class TENSORBASIS>
        TensorIndex<IBASIS,DIM,TENSORBASIS>
        first_wavelet(const TENSORBASIS* basis, const int level)
        {
            assert(level >= multi_degree(basis->j0()) );
            typename TensorIndex<IBASIS,DIM,TENSORBASIS>::level_type j;
            typename TensorIndex<IBASIS,DIM,TENSORBASIS>::type_type e;
            typename TensorIndex<IBASIS,DIM,TENSORBASIS>::translation_type k;

            j[DIM-1] = basis->jo()[DIM-1] + level - multi_degree(basis->j0());
            e[DIM-1] = 1;
            k[DIM-1] = basis->basis()[DIM-1]->Nablamin();

            for (unsigned int i = 0; i < DIM-1; i++)
            {
                j[i] = basis->j0()[i];
                e[i] = 0;
                k[i] = basis->bases[i]->DeltaLmin();
            }
            return TensorIndex<IBASIS,DIM,TensorBasis<IBASIS,DIM> >(j, e, k, basis);
        }

  	template <class IBASIS, unsigned int DIM, class TENSORBASIS>
  	TensorIndex<IBASIS,DIM,TENSORBASIS>
  	last_wavelet(const TENSORBASIS* basis, const typename TensorIndex<IBASIS,DIM,TENSORBASIS>::level_type j)
  	{
            assert(multi_degree(j) > multi_degree(basis->j0())
                    || ((multi_degree(j) == multi_degree(basis->j0())) && (basis->j0() <= j))
                  );
            typename TensorIndex<IBASIS,DIM,TENSORBASIS>::type_type e;
            typename TensorIndex<IBASIS,DIM,TENSORBASIS>::translation_type k;
            for (unsigned int i = 0; i < DIM; i++) {
                e[i] = 1;
                k[i] = basis->bases()[i]->Nablamax(j[i]);
            }
            return TensorIndex<IBASIS,DIM,TENSORBASIS>(j, e, k, basis);
	}

        // PERFORMANCE !! viele Aufrufe von j0()
        template <class IBASIS, unsigned int DIM, class TENSORBASIS>
  	TensorIndex<IBASIS,DIM,TENSORBASIS>
  	last_wavelet(const TENSORBASIS* basis, const unsigned int level)
  	{
            assert(level >= multi_degree(basis->j0()) );
            typename TensorIndex<IBASIS,DIM,TensorBasis<IBASIS,DIM> >::level_type j;
            typename TensorIndex<IBASIS,DIM,TensorBasis<IBASIS,DIM> >::type_type e;
            typename TensorIndex<IBASIS,DIM,TensorBasis<IBASIS,DIM> >::translation_type k;
            j[0]= basis->j0()[0]+level-multi_degree(basis->j0());
            e[0]=1;
            k[0]= basis->bases()[0]->Nablamax(j[0]);
            for (unsigned int i = 1; i < DIM; i++) {
                j[i] = basis->j0()[i];
                e[i] = 1;
                k[i] = basis->bases()[i]->Nablamax(basis->bases()[i]->j0());
            }
            return TensorIndex<IBASIS,DIM,TensorBasis<IBASIS,DIM> >(j, e, k, basis);
	}

	template <class IBASIS, unsigned int DIM, class TENSORBASIS>
	const int
	first_generator_num(const TENSORBASIS* basis)
	{
		return 0;
	}

	template <class IBASIS, unsigned int DIM, class TENSORBASIS>
	const int
	last_generator_num(const TENSORBASIS* basis)
	{
		int res=1;
		for (unsigned int i = 0; i < DIM; i++)
			res *= basis->bases()[i]->Deltasize(basis->bases()[i]->j0());
		return res-1;
	}

	template <class IBASIS, unsigned int DIM, class TENSORBASIS>
	const int
	first_wavelet_num(const TENSORBASIS* basis, const typename TensorIndex<IBASIS,DIM,TENSORBASIS>::level_type j)
        {
            // Assertions are checked in first_wavelet
            // for (unsigned int i = 0; i < DIM; i++)
            //     assert(j[i] >= (basis->bases()[i]->j0()));
            TensorIndex<IBASIS,DIM,TENSORBASIS> temp (first_wavelet<IBASIS,DIM,TENSORBASIS>(basis,j));
            return temp.number();
	}

	template <class IBASIS, unsigned int DIM, class TENSORBASIS>
	const int
	last_wavelet_num(const TENSORBASIS* basis, const typename TensorIndex<IBASIS,DIM,TENSORBASIS>::level_type j)
	{
            // Assertions are checked in last_wavelet
            // for (unsigned int i = 0; i < DIM; i++)
            //     assert(j[i] >= (basis->bases()[i]->j0()));
            TensorIndex<IBASIS,DIM,TENSORBASIS> temp (last_wavelet<IBASIS,DIM,TENSORBASIS>(basis,j));
            return temp.number();
	}

        template <class IBASIS, unsigned int DIM, class TENSORBASIS>
	const int
	last_wavelet_num(const TENSORBASIS* basis, const unsigned int level)
	{
            // Assertions are checked in last_wavelet
            // for (unsigned int i = 0; i < DIM; i++)
            //     assert(j[i] >= (basis->bases()[i]->j0()));
            TensorIndex<IBASIS,DIM,TENSORBASIS> temp (last_wavelet<IBASIS,DIM,TENSORBASIS>(basis,level));
            return temp.number();
	}
}
