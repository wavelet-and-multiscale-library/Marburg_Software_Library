// implementation for tframe_index.h

namespace WaveletTL
{
	template <class IFRAME, unsigned int DIM, class TENSORFRAME>
	TensorQIndex<IFRAME,DIM,TENSORFRAME>::TensorQIndex(const TENSORFRAME* frame)
	:frame_(frame)
        {
            if (frame_ != 0)
            {
    		for (unsigned int i = 0; i < DIM; i++) 
                {
                    p_[i] = 0;
                    j_[i] = frame_->frames()[i]->j0(); // coarsest level;
                    // e_ is zero by default: generator
                    k_[i] = frame_->frames()[i]->DeltaLmin();
    		}
      		num_ = 0;
#if _TFRAME_DEBUGLEVEL_ >= 1
                cout << "TensorQIndex: Constructor from frame: num_ = 0 instead of 1" << endl;
                abort();
#endif
            }
            else 
            {
      		//j_ = 0; // invalid (e and k are initialized by zero automatically)
                num_ = -1;
            }
	}

        template <class IFRAME, unsigned int DIM, class TENSORFRAME>
	TensorQIndex<IFRAME,DIM,TENSORFRAME>:: TensorQIndex(const polynomial_type& p, const level_type& j, const type_type& e, const translation_type& k, const TENSORFRAME* frame)
	: frame_(frame), p_(p), j_(j), e_(e), k_(k)
	{
            if (frame != 0)
            {
                MathTL::FixedArray1D<std::map<int, int>,DIM> sizes; // store number of frame elements. Generators on level j0 (0), quarklets on level j0 (1), j0+1 (2), ...
                int uptothislevel(0); // number of frame function below the current level
		int oncurrentlevel(1); // number of base elements on current level j
                level_type j0(frame->j0());
                MultiIndex<int,DIM> j_minus_j0(j);
                MultiIndex<bool,DIM> min_type;
                unsigned int levelrange(0); // multi_degree (j) - multi_degree(j0)
                for (unsigned int i=0; i<DIM; ++i)
                {
                    j_minus_j0 [i] -= j0[i];
                    levelrange += j_minus_j0[i];
                    min_type[i] = (j_minus_j0[i] != 0); // min_type = 0 => generators are allowed in this direction
                }
                
                for (unsigned int i=0; i < DIM; i++) {
                    sizes[i][0] = frame->frames()[i]->Deltasize(j0[i]); // Number of generators on level j0
                    for (unsigned int k = 0; k <= levelrange; ++k)
                    {
                        sizes[i][k+1] = frame->frames()[i]->Nablasize(j0[i]+k); // Number of Wavelets on level j0+k
                    }
		}
                // sizes[0][j0[0]+levelrange] is never used
                
                unsigned int l(0);
                unsigned int maxnum(j_minus_j0.number());
                
                MultiIndex<int,DIM> mult_it;
                for (unsigned int i=0; i<DIM; ++i)
                {
                    mult_it[i] = 0;
                }
                
                while (l < maxnum)
                {
                    oncurrentlevel = 1;
                    for (unsigned int i=0; i<DIM; ++i)
                    {
                        if (mult_it[i] == 0)
                        {
                            oncurrentlevel *= (sizes[i][0] + sizes[i][1]);
                        }
                        else
                        {
                            oncurrentlevel *= sizes[i][mult_it[i]+1];
                        }
                    }
                    uptothislevel += oncurrentlevel;
                    ++l;
                    ++mult_it;
                }
                // count the number of (j,e,k) with j = j_, but with e < e_
                // iterate over all allowed types <= e_

                
                l=0;
                for (unsigned int i=0; i<DIM; ++i)
                {
                    mult_it[i] = min_type[i];
                }
                //cout << "target_e = " << e_ << endl;
                while (mult_it.lex(e_))
                {
                    //cout << "current_e at beginning of loop = " << mult_it << endl;
                    oncurrentlevel = 1;
                    for (unsigned int i=0; i<DIM; ++i)
                    {
                        if (mult_it[i] == 0)
                        {
                            oncurrentlevel *= sizes[i][0];
                        }
                        else
                        {
                            oncurrentlevel *= sizes[i][j_minus_j0[i] +1];
                        }
                    }
                    uptothislevel += oncurrentlevel;
                    
                    // mult_it ++ (type++)
                    bool done = true;
                    for (int i(DIM-1); i >= 0; i--)
                    {
                        // find first position on level j0
                        if (min_type[i] == false)
                        {
                            if (mult_it[i] == 1)
                            {
                                mult_it[i]=0;
                                //k_[i]=frame_->frames()[i]->DeltaLmin();
                            } else
                            {
                                mult_it[i]=1;
                                //k_[i]=frame_->frames()[i]->Nablamin();
                                done = false;
                                break;
                            }
                        }
                    }
                    assert (!done); // while loop assumes mult_it < e_
                    //cout << "current_e at end of loop = " << mult_it << "; " << mult_it << " < " << e_ << " = " << (mult_it < e_) << endl;
                    
                }
                //cout << "current_e after loop = " << mult_it << endl;
                // count quarklets with same j,e but lower k
                oncurrentlevel = 1; // product of ksize(DIM-1) * ...* ksize(current_dim)
                for (int i=DIM-1; i>0; --i)
                {
                    if (e_[i] == 0)
                    {
                        uptothislevel += (k_[i] - frame->frames()[i]->DeltaLmin())*oncurrentlevel;
                        oncurrentlevel *= sizes[i][0];
                    }
                    else
                    {
                        uptothislevel += (k_[i] - frame->frames()[i]->Nablamin())*oncurrentlevel;
                        oncurrentlevel *= sizes[i][j_[i]-j0[i]+1];
                    }
                }
                // dimension 0 : no blocksizes need to be computed
                if (e_[0] == 0)
                {
                    uptothislevel += (k_[0] - frame->frames()[0]->DeltaLmin())*oncurrentlevel;
                }
                else
                {
                    uptothislevel += (k_[0] - frame->frames()[0]->Nablamin())*oncurrentlevel;
                }
                num_ = uptothislevel;
            }
            else
            {
                num_ = -1;
#if _TFRAME_DEBUGLEVEL_ >= 1
                cout << "TensorQIndex(p,j,e,k, frame) with frame = 0" << endl;
                abort();
#endif
            }
        }

        template <class IFRAME, unsigned int DIM, class TENSORFRAME>
	TensorQIndex<IFRAME, DIM, TENSORFRAME>::TensorQIndex(const int& p, const int& j, const int& e, const int& k, const TENSORFRAME* frame)
	: frame_(frame), num_(0)//, j_(j), e_(e), k_(k)
        {
            assert (DIM == 1);
            p_[0]=p;
            j_[0]=j;
            e_[0]=e;
            k_[0]=k;
        }

	template <class IFRAME, unsigned int DIM, class TENSORFRAME>
	TensorQIndex<IFRAME, DIM, TENSORFRAME>::TensorQIndex(const TensorQIndex& lambda)
	: frame_(lambda.frame_), p_(lambda.p_), j_(lambda.j_), e_(lambda.e_), k_(lambda.k_), num_(lambda.num_) {}

	template <class IFRAME, unsigned int DIM, class TENSORFRAME>
	TensorQIndex<IFRAME, DIM, TENSORFRAME>::TensorQIndex(const TensorQIndex* lambda)
	: frame_(lambda->frame_), p_(lambda.p_), j_(lambda->j_), e_(lambda->e_), k_(lambda->k_), num_(lambda->num_) {}

        template <class IFRAME, unsigned int DIM, class TENSORFRAME>
        TensorQIndex<IFRAME,DIM,TENSORFRAME>:: TensorQIndex(const polynomial_type& p, const level_type& j, const type_type& e, const translation_type& k, const int number, const TENSORFRAME* frame)
        : frame_(frame), p_(p), j_(j), e_(e), k_(k), num_(number) {}
        
	template <class IFRAME, unsigned int DIM, class TENSORFRAME>
	TensorQIndex<IFRAME, DIM, TENSORFRAME>::TensorQIndex(const int number, const TENSORFRAME* frame)
	: frame_(frame), num_(number)
	{
#if _TFRAME_DEBUGLEVEL_ >= 1
            cout << "TensorQIndex(number) is called. Use full_collection[number] instead!" << endl;
#endif
          /* 
            This implementation for DIM = 2,3 assumes that Nablasize(j) = Deltasize(j0) + sum_{l=0}^{j-j0} 2^l * Nablasize(j0)
            and uses a lot of div and mod
          */
          if(DIM == 2)
          {
          level_type j0 = frame->j0();
          int Nabla1 = frame->frames()[0]->Nablasize(j0[0]);
          int Nabla2 = frame->frames()[1]->Nablasize(j0[1]);
          int Delta1 = frame->frames()[0]->Deltasize(j0[0]);
          int Delta2 = frame->frames()[1]->Deltasize(j0[1]);
          int act_num = num_;
          int i = 0; 
          if (act_num < Delta1*Delta2) 
          {
            j_ = j0;
            e_[0]=0;
            e_[1]=0;
            k_[0] = frame->frames()[0]->DeltaLmin() + act_num/Delta2;
            k_[1] = frame->frames()[1]->DeltaLmin() + (act_num % Delta2);
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
              k_[0] = frame->frames()[0]->DeltaLmin() + act_num/(Nabla2<<i);
              k_[1] = frame->frames()[1]->Nablamin() + (act_num % (Nabla2<<i));
            }
            else if ((act_num - ((Delta1*Nabla2)<<i)) / ((Nabla1*Nabla2)<<i) < i)
            {
              unsigned int tmp2 = (act_num - ((Delta1*Nabla2)<<i)) / ((Nabla1*Nabla2)<<i);
              j_[0] = j0[0] + tmp2;
              j_[1] = j0[1] + i - tmp2;
              e_[0]=1;
              e_[1]=1;
              act_num -= (Delta1*Nabla2 + tmp2 * (Nabla1*Nabla2))<<i;
              k_[0] = frame->frames()[0]->Nablamin() + act_num/(Nabla2<<(i-tmp2));
              k_[1] = frame->frames()[1]->Nablamin() + (act_num % (Nabla2<<(i-tmp2)));
            }
            else if (act_num < (Delta1*Nabla2 + (i) * Nabla1*Nabla2 + Delta2*Nabla1)<<i ) 
            {
              j_[0] = j0[0] + i;
              j_[1] = j0[1];
              e_[0]=1;
              e_[1]=0;
              act_num -= (Delta1*Nabla2 + (i) * (Nabla1*Nabla2))<<i;
              k_[0] = frame->frames()[0]->Nablamin() + act_num/Delta2;
              k_[1] = frame->frames()[1]->DeltaLmin() + (act_num % Delta2);
            }
            else
            {
              j_[0] = j0[0] + i;
              j_[1] = j0[1];
              e_[0]=1;
              e_[1]=1;
              act_num -= (Delta1*Nabla2 + (i) * (Nabla1*Nabla2) + Delta2*Nabla1)<<i;
              k_[0] = frame->frames()[0]->Nablamin() + act_num/Nabla2;
              k_[1] = frame->frames()[1]->Nablamin() + (act_num % Nabla2);
            }
          }
          }
          else if(DIM == 3)
          {
          level_type j0 = frame->j0();
          int Nabla1 = frame->frames()[0]->Nablasize(j0[0]);
          int Nabla2 = frame->frames()[1]->Nablasize(j0[1]);
          int Nabla3 = frame->frames()[2]->Nablasize(j0[2]);
          int Delta1 = frame->frames()[0]->Deltasize(j0[0]);
          int Delta2 = frame->frames()[1]->Deltasize(j0[1]);
          int Delta3 = frame->frames()[2]->Deltasize(j0[2]);
          int act_num = num_;
          int i = 0; 
          if (act_num < Delta1*Delta2*Delta3) 
          {
            j_ = j0;
            e_[0]=0;
            e_[1]=0;
            e_[2]=0;
            k_[0] = frame->frames()[0]->DeltaLmin() + (act_num/(Delta3*Delta2));
            k_[1] = frame->frames()[1]->DeltaLmin() + (act_num/Delta3) % Delta2;
            k_[2] = frame->frames()[2]->DeltaLmin() + (act_num % Delta3);
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
                k_[0] = frame->frames()[0]->DeltaLmin() + (act_num/(Nabla3*Delta2));
                k_[1] = frame->frames()[1]->DeltaLmin() + (act_num/Nabla3) % Delta2;
                k_[2] = frame->frames()[2]->Nablamin() + (act_num % Nabla3);
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
                k_[0] = frame->frames()[0]->DeltaLmin() + (act_num/(Delta3*Nabla2));
                k_[1] = frame->frames()[1]->Nablamin() + (act_num/Delta3) % Nabla2;
                k_[2] = frame->frames()[2]->DeltaLmin() + (act_num % Delta3);
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
                k_[0] = frame->frames()[0]->DeltaLmin() + (act_num/(Nabla3*Nabla2));
                k_[1] = frame->frames()[1]->Nablamin() + (act_num/Nabla3) % Nabla2;
                k_[2] = frame->frames()[2]->Nablamin() + (act_num % Nabla3);
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
                k_[0] = frame->frames()[0]->Nablamin() + (act_num/(Delta3*Delta2));
                k_[1] = frame->frames()[1]->DeltaLmin() + (act_num/Delta3) % Delta2;
                k_[2] = frame->frames()[2]->DeltaLmin() + (act_num % Delta3);
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
                k_[0] = frame->frames()[0]->Nablamin() + (act_num/(Nabla3*Delta2));
                k_[1] = frame->frames()[1]->DeltaLmin() + (act_num/Nabla3) % Delta2;
                k_[2] = frame->frames()[2]->Nablamin() + (act_num % Nabla3);
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
                k_[0] = frame->frames()[0]->Nablamin() + (act_num/(Delta3*Nabla2));
                k_[1] = frame->frames()[1]->Nablamin() + (act_num/Delta3) % Nabla2;
                k_[2] = frame->frames()[2]->DeltaLmin() + (act_num % Delta3);
              }
              else
              {
              act_num -= Nabla1*Nabla2*Delta3;
                j_ = j0;
                e_[0] = 1;
                e_[1] = 1;
                e_[2] = 1;
                k_[0] = frame->frames()[0]->Nablamin() + (act_num/(Nabla3*Nabla2));
                k_[1] = frame->frames()[1]->Nablamin() + (act_num/Nabla3) % Nabla2;
                k_[2] = frame->frames()[2]->Nablamin() + (act_num % Nabla3);
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
                  k_[0] = frame->frames()[0]->DeltaLmin() + act_num/((Delta2*Nabla3)<<i);
                  k_[1] = frame->frames()[1]->DeltaLmin() + (act_num/(Nabla3<<i)) % Delta2;
                  k_[2] = frame->frames()[2]->Nablamin() + (act_num % (Nabla3<<i));
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
                  k_[0] = frame->frames()[0]->DeltaLmin() + act_num/((Nabla2*Nabla3)<<i);
                  k_[1] = frame->frames()[1]->Nablamin() + (act_num/(Nabla3<<i)) % Nabla2;
                  k_[2] = frame->frames()[2]->Nablamin() + (act_num % (Nabla3<<i));
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
                  k_[0] = frame->frames()[0]->Nablamin() + act_num/((Delta2*Nabla3)<<i);
                  k_[1] = frame->frames()[1]->DeltaLmin() + (act_num/(Nabla3<<i)) % Delta2;
                  k_[2] = frame->frames()[2]->Nablamin() + (act_num % (Nabla3<<i));
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
                  k_[0] = frame->frames()[0]->Nablamin() + act_num/((Nabla2*Nabla3)<<i);
                  k_[1] = frame->frames()[1]->Nablamin() + (act_num/(Nabla3<<i)) % Nabla2;
                  k_[2] = frame->frames()[2]->Nablamin() + (act_num % (Nabla3<<i));
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
                  k_[0] = frame->frames()[0]->DeltaLmin() + act_num/((Nabla2*Nabla3)<<i);
                }
                else
                {
                  e_[0] = 1;
                  act_num -= Delta1*Nabla2*Nabla3 <<i;
                  k_[0] = frame->frames()[0]->Nablamin() + act_num/((Nabla2*Nabla3)<<i);
                }
                e_[1] = 1;
                e_[2] = 1;
                k_[1] = frame->frames()[1]->Nablamin() + (act_num/(Nabla3<<(i-tmp2-1))) % (Nabla2 <<(tmp2+1));
                k_[2] = frame->frames()[2]->Nablamin() + (act_num % (Nabla3<<(i-tmp2-1)));
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
                  k_[0] = frame->frames()[0]->DeltaLmin() + act_num/((Nabla2*Delta3)<<i);
                  k_[1] = frame->frames()[1]->Nablamin() + (act_num/Delta3) % (Nabla2 <<i);
                  k_[2] = frame->frames()[2]->DeltaLmin() + (act_num % Delta3);
                }
                else 
                {
                act_num -= (Delta1*Nabla2*Delta3) <<i;
                if(act_num < (Delta1*Nabla2*Nabla3) <<i)
                {
                  e_[0] = 0;
                  e_[2] = 1;
                  k_[0] = frame->frames()[0]->DeltaLmin() + act_num/((Nabla2*Nabla3)<<i);
                  k_[1] = frame->frames()[1]->Nablamin() + (act_num/Nabla3) % (Nabla2 <<i);
                  k_[2] = frame->frames()[2]->Nablamin() + (act_num % Nabla3);
                }
                else 
                {
                act_num -= (Delta1*Nabla2*Nabla3) <<i;
                if(act_num < (Nabla1*Nabla2*Delta3) <<i)
                {
                  e_[0] = 1;
                  e_[2] = 0;
                  k_[0] = frame->frames()[0]->Nablamin() + act_num/((Nabla2*Delta3)<<i);
                  k_[1] = frame->frames()[1]->Nablamin() + (act_num/Delta3) % (Nabla2 <<i);
                  k_[2] = frame->frames()[2]->DeltaLmin() + (act_num % Delta3);
                }
                else 
                {
                act_num -= (Nabla1*Nabla2*Delta3) <<i;
                  e_[0] = 1;
                  e_[2] = 1;
                  k_[0] = frame->frames()[0]->Nablamin() + act_num/((Nabla2*Nabla3)<<i);
                  k_[1] = frame->frames()[1]->Nablamin() + (act_num/Nabla3) % (Nabla2 <<i);
                  k_[2] = frame->frames()[2]->Nablamin() + (act_num % Nabla3);
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
                    k_[0] = frame->frames()[0]->Nablamin() + act_num/((Delta2*Nabla3)<<(i-l));
                    k_[1] = frame->frames()[1]->DeltaLmin() + (act_num/(Nabla3<<(i-l))) % Delta2;
                    k_[2] = frame->frames()[2]->Nablamin() + (act_num % (Nabla3<<(i-l)));
                  }
                  else if( ((act_num - ((Nabla1*Delta2*Nabla3)<<i))/((Nabla1*Nabla2*Nabla3)<<i)) < (i-l) )
                  {
                    unsigned int tmp2 = ((act_num - ((Nabla1*Delta2*Nabla3)<<i))/((Nabla1*Nabla2*Nabla3)<<i));
                    j_[1] = j0[1] + tmp2;
                    j_[2] = j0[2] + i - l - tmp2;
                    e_[1] = 1;
                    e_[2] = 1;
                    act_num -= (Nabla1*Delta2*Nabla3 + tmp2*Nabla1*Nabla2*Nabla3)<<i;
                    k_[0] = frame->frames()[0]->Nablamin() + act_num/((Nabla2*Nabla3)<<(i-l));
                    k_[1] = frame->frames()[1]->Nablamin() + (act_num/(Nabla3<<(i-l-tmp2))) % (Nabla2<<tmp2);
                    k_[2] = frame->frames()[2]->Nablamin() + (act_num % (Nabla3<<(i-l-tmp2)));
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
                      k_[0] = frame->frames()[0]->Nablamin() + act_num/((Nabla2*Delta3)<<(i-l));
                      k_[1] = frame->frames()[1]->Nablamin() + (act_num/Delta3) % (Nabla2<<(i-l));
                      k_[2] = frame->frames()[2]->DeltaLmin() + (act_num % Delta3);
                    }
                    else
                    {
                      j_[1] = j0[1] + i - l;
                      j_[2] = j0[2];
                      e_[1] = 1;
                      e_[2] = 1;
                      act_num -= Nabla1*Nabla2*Delta3<<i;
                      k_[0] = frame->frames()[0]->Nablamin() + act_num/((Nabla2*Nabla3)<<(i-l));
                      k_[1] = frame->frames()[1]->Nablamin() + (act_num/Nabla3) % (Nabla2<<(i-l));
                      k_[2] = frame->frames()[2]->Nablamin() + (act_num % Nabla3);
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
                    k_[0] = frame->frames()[0]->Nablamin() + act_num/(Delta2*Delta3);
                    k_[1] = frame->frames()[1]->DeltaLmin() + (act_num/Delta3) % Delta2;
                    k_[2] = frame->frames()[2]->DeltaLmin() + (act_num % Delta3);
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
                    k_[0] = frame->frames()[0]->Nablamin() + act_num/(Delta2*Nabla3);
                    k_[1] = frame->frames()[1]->DeltaLmin() + (act_num/Nabla3) % Delta2;
                    k_[2] = frame->frames()[2]->Nablamin() + (act_num % Nabla3);
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
                    k_[0] = frame->frames()[0]->Nablamin() + act_num/(Nabla2*Delta3);
                    k_[1] = frame->frames()[1]->Nablamin() + (act_num/Delta3) % Nabla2;
                    k_[2] = frame->frames()[2]->DeltaLmin() + (act_num % Delta3);
                  }
                  else
                  {
                  act_num -= Nabla1*Nabla2*Delta3<<i;
                    j_[1] = j0[1];
                    j_[2] = j0[2];
                    e_[1] = 1;
                    e_[2] = 1;
                    k_[0] = frame->frames()[0]->Nablamin() + act_num/(Nabla2*Nabla3);
                    k_[1] = frame->frames()[1]->Nablamin() + (act_num/Nabla3) % Nabla2;
                    k_[2] = frame->frames()[2]->Nablamin() + (act_num % Nabla3);
                  }}}
                }
              }
            }
          }
          }
          else //DIM >3
          {
                MathTL::FixedArray1D<std::map<int, int>,DIM> sizes; // Store number of frame elements. Generators on level j0 (0), quarklets on level j0 (1), j0+1 (2), ...
		int remains = number+1; // numbering begins at 0
		int oncurrentlevel(1); // number of base elements on current level j
		int range(0); // number of level steps we have climbed so far. determines wether a new entry has to be computed in "sizes"
		level_type currentlevel, j0(frame_->j0());
		type_type currenttype;
		currentlevel = j0;
		for (unsigned int i=0; i < DIM; i++) {
			currenttype[i]=0;
			sizes[i][0] = frame_->frames()[i]->Deltasize(j0[i]); // Number of generators on level j0
                        sizes[i][1] = frame_->frames()[i]->Nablasize(j0[i]); // N o Wavelets
                        oncurrentlevel *= sizes[i][0];
		}
                
                // iterate over all levels. Add up number of frame functions till looked for level is reached
		while (remains > oncurrentlevel) // break if we are at the right level
		{
                    // else substract number of frame functions on current_index=(currentlevel,currentindex) and increase current_index
                    remains -= oncurrentlevel;
                    
                    // increase index = (currentlevel,currenttype)
                    // "small loop" "currenttype++" (currentlevel fest)
                    // iterate over all combinations of generators/quarklets for all dimensions with currentlevel[i]=j0[i]
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
                    // done == true means that all components with currentlevel[i]=j0[i] were quarklets.
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
                                        sizes[i-1][range+1]=frame ->frames()[i-1]->Nablasize(currentlevel[i-1]); // if needed compute and store new size information
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
                                    currenttype[0] = 1;
                                    currentlevel[0]=currentlevel[0]+1;
                                    sizes[0][range+1]=frame ->frames()[0]->Nablasize(currentlevel[0]); // if needed compute and store new size information
                                }
                                else
                                {
                                    //currentlevel[DIM-1]=j0[DIM-1]+currentlevel[0]-j0[0]+1; currenttype[DIM-1]=1;
                                    currentlevel[DIM-1]=j0[DIM-1]+range; currenttype[DIM-1]=1;
                                    for(int i(DIM-2); i>=0; i--)
                                    {
                                      currenttype[i]=0; currentlevel[i]=j0[i];
                                    }
                                    sizes[DIM-1][range+1]=frame ->frames()[DIM-1]->Nablasize(currentlevel[DIM-1]); // if needed compute and store new size information
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

                // determine k corresponding to the number of the frame function given by "remains" (on level (currentlevel,currenttype) )
		unsigned int modul;                
		j_ = currentlevel;
		e_ = currenttype;
                remains -= 1; // numbering begins at 0
		for (int i = DIM-1; i > 0; i--)
		{
			modul = sizes[i][currentlevel[i]+currenttype[i]-j0[i]];
			k_[i]= remains % modul+(e_[i] == 0 ? frame->frames()[i]->DeltaLmin():frame->frames()[i]->Nablamin());
			remains = remains / modul;
		}
                k_[0]=remains +(e_[0] == 0 ? frame->frames()[0]->DeltaLmin():frame->frames()[0]->Nablamin());

          } //end DIM >3
        }

	template <class IFRAME, unsigned int DIM, class TENSORFRAME>
  	TensorQIndex<IFRAME,DIM,TENSORFRAME>&
  	TensorQIndex<IFRAME,DIM,TENSORFRAME>::operator = (const TensorQIndex<IFRAME,DIM,TENSORFRAME>& lambda)
  	{
            p_ = lambda.p();
            j_ = lambda.j();
            e_ = lambda.e();
            k_ = lambda.k();
            frame_ = lambda.frame();
            num_ = lambda.number();
            return *this;
  	}

	template <class IFRAME, unsigned int DIM, class TENSORFRAME>
	bool
	TensorQIndex<IFRAME, DIM, TENSORFRAME>::operator == (const TensorQIndex<IFRAME, DIM, TENSORFRAME>& lambda) const
	{
            return (p_ == lambda.p() &&
                    j_ == lambda.j() &&
                    e_ == lambda.e() &&
                    k_ == lambda.k());
	}

        // PERFORMANCE: eventuell die vielen Aufrufe von frames_->j0() durch Hilfsvariable auf einen reduzieren
	template <class IFRAME, unsigned int DIM, class TENSORFRAME>
  	TensorQIndex<IFRAME,DIM,TENSORFRAME>&
  	TensorQIndex<IFRAME,DIM,TENSORFRAME>::operator ++ ()
  	{
            if ((int)num_ == -1) return *this;
            level_type j0(frame_->j0());
            num_++;
            // determine next translation index
            bool jplusplus = false;
            for (int i = DIM-1; i >= 0; i--) {
      		const int last_index = (e_[i] == 0 ? frame_->frames()[i]->DeltaRmax(j_[i])
                                                   : frame_->frames()[i]->Nablamax(j_[i]));
      		if (k_[i] == last_index)
                {
                    k_[i] = (e_[i] == 0 ? frame_->frames()[i]->DeltaLmin()
	 				: frame_->frames()[i]->Nablamin());
                    jplusplus = (i == 0);
      		} else
                {
                    ++k_[i];
                    break;
      		}
            }
            if (jplusplus == false) return *this;
            // else: determine next level index
            // "small loop" "e_++" (j_ is not changed)
            // iterate over all combinations of generators/quarklets for all dimensions with j_[i]=j0[i]
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
                        k_[i]=frame_->frames()[i]->DeltaLmin();
                    } else
                    {
                        e_[i]=1;
                        k_[i]=frame_->frames()[i]->Nablamin();
                        done = false;
                        break;
                    }
                }
            }
            // done == true bedeutet, dass alle Komponenten auf level j0() quarklets waren.
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
                            k_[i-1]=frame_->frames()[i-1]->Nablamin();
                            int temp = j_[i]-j0[i];
                            j_[i]=j0[i];
                            e_[i]=0;
                            k_[i]=frame_->frames()[i]->DeltaLmin();
                            j_[DIM-1]=j0[DIM-1]+temp-1;
                            e_[DIM-1]= (temp == 1?0:1);
                            k_[DIM-1]= (temp == 1?frame_->frames()[DIM-1]->DeltaLmin():frame_->frames()[DIM-1]->Nablamin());
                            break;
                        }
                    } else // i == 0. "big loop" arrived at the last index. We have to increase the level!
                    {
                        if (DIM == 1)
                        {
                            e_[i] = 1; // diese Zeile erfüllt nur in der allerersten Iteration einen Zweck
                            k_[i]=frame_->frames()[i]->Nablamin(); // diese Zeile erfüllt nur in der allerersten Iteration einen Zweck
                            j_[i]=j_[i]+1;
                        }
                        else
                        {
                            j_[DIM-1]=j0[DIM-1]+j_[0]-j0[0]+1;
                            e_[DIM-1]=1;
                            k_[DIM-1]=frame_->frames()[DIM-1]->Nablamin();
                            j_[0]=j0[0];
                            e_[0]=0;
                            k_[0]=frame_->frames()[0]->DeltaLmin();
                        }
                        break; // unnoetig, da i==0 gilt.
                    }
                } 
            } // end of "big loop"
            return *this;
         
  	}

	template <class IFRAME, unsigned int DIM, class TENSORFRAME>
	bool
	TensorQIndex<IFRAME,DIM,TENSORFRAME>::operator < (const TensorQIndex& lambda) const
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

	template <class IFRAME, unsigned int DIM, class TENSORFRAME>
	TensorQIndex<IFRAME,DIM,TENSORFRAME>
	first_q_generator(const TENSORFRAME* frame)
	{
            return TensorQIndex<IFRAME,DIM,TENSORFRAME>(0, frame);
	}

	template <class IFRAME, unsigned int DIM, class TENSORFRAME>
	TensorQIndex<IFRAME,DIM,TENSORFRAME>
	last_q_generator(const TENSORFRAME* frame)
	{
		typename TensorQIndex<IFRAME,DIM,TENSORFRAME>::polynomial_type p;
                typename TensorQIndex<IFRAME,DIM,TENSORFRAME>::level_type j;
		typename TensorQIndex<IFRAME,DIM,TENSORFRAME>::type_type e;
		typename TensorQIndex<IFRAME,DIM,TENSORFRAME>::translation_type k;
		j = frame->j0();
		for (unsigned int i = 0; i < DIM; i++)
		{
			e[i] = 0;
			k[i] = frame->frames()[i]->DeltaRmax(j[i]);
		}
                return TensorQIndex<IFRAME,DIM,TENSORFRAME>(p, j, e, k, last_q_generator_num<IFRAME,DIM,TENSORFRAME>(frame), frame);
	}

  	template <class IFRAME, unsigned int DIM, class TENSORFRAME>
  	TensorQIndex<IFRAME,DIM,TENSORFRAME>
  	first_quarklet(const TENSORFRAME* frame, const typename TensorQIndex<IFRAME,DIM,TENSORFRAME>::level_type j, 
                       const typename TensorQIndex<IFRAME,DIM,TENSORFRAME>::polynomial_type p)
  	{
#if _TFRAME_DEBUGLEVEL_ >= 1
            assert(multi_degree(j) > multi_degree(frame->j0())
                    || ((multi_degree(j) == multi_degree(frame->j0())) && (frame->j0() <= j))
                  );
#endif
            typename TensorQIndex<IFRAME,DIM,TENSORFRAME>::type_type e;
            typename TensorQIndex<IFRAME,DIM,TENSORFRAME>::translation_type k;
            
            bool sofar_only_generators = true;
            for (unsigned int i = 0; i < DIM-1; i++) {
                if (j[i] == frame->frames()[i]->j0())
                {
                    e[i] = 0;
                    k[i] = frame->frames()[i]->DeltaLmin();
                } else
                {
                    e[i] = 1;
                    k[i] = frame->frames()[i]->Nablamin();
                    sofar_only_generators = false;
                }
            }
            if ( (sofar_only_generators == true) || (j[DIM-1] != frame->frames()[DIM-1]->j0()) )
            {
                e[DIM-1] = 1;
                k[DIM-1] = frame->frames()[DIM-1]->Nablamin();
                //sofar_only_generators = false;
            } else
            {
                e[DIM-1] = 0;
                k[DIM-1] = frame->frames()[DIM-1]->DeltaLmin();
            }
            return TensorQIndex<IFRAME,DIM,TENSORFRAME>(j, e, k, frame);
        }


        template <class IFRAME, unsigned int DIM, class TENSORFRAME>
        TensorQIndex<IFRAME,DIM,TENSORFRAME>
        first_quarklet(const TENSORFRAME* frame, const int level, typename TensorQIndex<IFRAME,DIM,TENSORFRAME>::polynomial_type p)
        {
#if _TFRAME_DEBUGLEVEL_ >= 1
            assert(level >= (int)multi_degree(frame->j0()) );
#endif
            typename TensorQIndex<IFRAME,DIM,TENSORFRAME>::level_type j;
            typename TensorQIndex<IFRAME,DIM,TENSORFRAME>::type_type e;
            typename TensorQIndex<IFRAME,DIM,TENSORFRAME>::translation_type k;

            int temp_i;
            j[DIM-1] = level;
            for (unsigned int i = 0; i < DIM-1; i++)
            {
                temp_i = frame->j0()[i];
                j[DIM-1] -= temp_i;
                j[i] = temp_i;
                e[i] = 0;
                k[i] = frame->frames()[i]->DeltaLmin();
            }
            e[DIM-1] = 1;
            k[DIM-1] = frame->frames()[DIM-1]->Nablamin();
            return TensorQIndex<IFRAME,DIM,TensorFrame<IFRAME,DIM> >(p, j, e, k, frame);
        }

  	template <class IFRAME, unsigned int DIM, class TENSORFRAME>
  	TensorQIndex<IFRAME,DIM,TENSORFRAME>
  	last_quarklet(const TENSORFRAME* frame, const typename TensorQIndex<IFRAME,DIM,TENSORFRAME>::level_type j,
                      const typename TensorQIndex<IFRAME,DIM,TENSORFRAME>::polynomial_type p)
  	{
#if _TFRAME_DEBUGLEVEL_ >= 1
            assert(multi_degree(j) > multi_degree(frame->j0())
                    || ((multi_degree(j) == multi_degree(frame->j0())) && (frame->j0() <= j))
                  );
#endif
            typename TensorQIndex<IFRAME,DIM,TENSORFRAME>::type_type e;
            typename TensorQIndex<IFRAME,DIM,TENSORFRAME>::translation_type k;
            for (unsigned int i = 0; i < DIM; i++) {
                e[i] = 1;
                k[i] = frame->frames()[i]->Nablamax(j[i]);
            }
            return TensorQIndex<IFRAME,DIM,TENSORFRAME>(p, j, e, k, frame);
	}

        template <class IFRAME, unsigned int DIM, class TENSORFRAME>
  	TensorQIndex<IFRAME,DIM,TENSORFRAME>
  	last_quarklet(const TENSORFRAME* frame, const unsigned int level, typename TensorQIndex<IFRAME,DIM,TensorFrame<IFRAME,DIM> >::polynomial_type p)
  	{
#if _TFRAME_DEBUGLEVEL_ >= 1
            assert(level >= multi_degree(frame->j0()));
#endif
            
            typename TensorQIndex<IFRAME,DIM,TensorFrame<IFRAME,DIM> >::level_type j;
            typename TensorQIndex<IFRAME,DIM,TensorFrame<IFRAME,DIM> >::type_type e;
            typename TensorQIndex<IFRAME,DIM,TensorFrame<IFRAME,DIM> >::translation_type k;
            int temp_i;
            j[0] = level;
            for (unsigned int i = DIM -1 ; i > 0 ; i--)
            {
                temp_i = frame ->j0()[i];
                j[0] -= temp_i;
                j[i] = temp_i;
                e[i] = 1;
                k[i] = frame->frames()[i]->Nablamax(temp_i);
                assert(frame->frames()[i]->j0() == frame -> j0()[i] );
            }
            e[0]=1;
            k[0]= frame->frames()[0]->Nablamax(j[0]);
            return TensorQIndex<IFRAME,DIM,TensorFrame<IFRAME,DIM> >(p, j, e, k, frame);
	}
        
	template <class IFRAME, unsigned int DIM, class TENSORFRAME>
	const int
	first_q_generator_num(const TENSORFRAME* frame)
	{
		return 0;
	}

	template <class IFRAME, unsigned int DIM, class TENSORFRAME>
	const int
	last_q_generator_num(const TENSORFRAME* frame)
	{
		int res=1;
		for (unsigned int i = 0; i < DIM; i++)
			res *= frame->frames()[i]->Deltasize(frame->frames()[i]->j0());
		return res-1;
	}

	template <class IFRAME, unsigned int DIM, class TENSORFRAME>
	const int
	first_quarklet_num(const TENSORFRAME* frame, const typename TensorQIndex<IFRAME,DIM,TENSORFRAME>::level_type j,
                           const typename TensorQIndex<IFRAME,DIM,TENSORFRAME>::polynomial_type p)
        {
            // Assertions are checked in first_quarklet
            // for (unsigned int i = 0; i < DIM; i++)
            //     assert(j[i] >= (frame->frames()[i]->j0()));
            TensorQIndex<IFRAME,DIM,TENSORFRAME> temp (first_quarklet<IFRAME,DIM,TENSORFRAME>(frame,j,p));
            return temp.number();
	}

	template <class IFRAME, unsigned int DIM, class TENSORFRAME>
	const int
	last_quarklet_num(const TENSORFRAME* frame, const typename TensorQIndex<IFRAME,DIM,TENSORFRAME>::level_type j,
                          const typename TensorQIndex<IFRAME,DIM,TENSORFRAME>::level_type p)
	{
            // Assertions are checked in last_quarklet
            // for (unsigned int i = 0; i < DIM; i++)
            //     assert(j[i] >= (frame->frames()[i]->j0()));
            TensorQIndex<IFRAME,DIM,TENSORFRAME> temp (last_quarklet<IFRAME,DIM,TENSORFRAME>(frame,j,p));
            return temp.number();
	}

        template <class IFRAME, unsigned int DIM, class TENSORFRAME>
	const int
	last_quarklet_num(const TENSORFRAME* frame, const unsigned int level)
	{
            // Assertions are checked in last_quarklet
            // for (unsigned int i = 0; i < DIM; i++)
            //     assert(j[i] >= (frame->frames()[i]->j0()));
            TensorQIndex<IFRAME,DIM,TENSORFRAME> temp (last_quarklet<IFRAME,DIM,TENSORFRAME>(frame,level));
            return temp.number();
	}
}
