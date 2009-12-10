// implementation for tbasis_index.h
// #include "tbasis_index.h"

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
		MathTL::FixedArray1D<std::map<int, int>,DIM> sizes; // Store number of basis elements. Generators on level j0 (0), wavelets on level j0 (1), j0+1 (2), ...
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
                    //cout << "sizes["<<i<<"][0]"<<sizes[i][0] <<endl;
                    //cout << "sizes["<<i<<"][1]"<<sizes[i][1] <<endl;
		}
		// iterate over all level up to j_ and add up number of basis functions up to the level
                if (currentlevel != j_ || currenttype != e_)
                while (true)
		{
                    // compute number of functions on this level
                    oncurrentlevel = 1;
                    for (unsigned int i = 0; i < DIM; i++)
                    {
                        //if (currentlevel[i]+currenttype[i]-j0[i] == range) sizes[i][range] = (currenttype[i] == 0 ? basis_->bases()[i]->Deltasize(j0[i]) : basis_->bases()[i]->Nablasize(j0[i]+range-1)); // store new information
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
                            if (currenttype[i] == 1)
                            {
                                currenttype[i]=0;
                                //if (i == 0) done = true;
                            } else
                            {
                                currenttype[i]=1;
                                done = false;
                                break;
                            }
                    }
                    // done == true bedeutet, dass alle Komponenten auf level j0() wavelets waren.
                    // Jetzt "big loop" also einmal currentlevel erhoehen
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
                                sizes[0][range]=basis ->bases()[0]->Nablasize(currentlevel[0]); // if needed compute and store new size information
                            }
                            else
                            {
                                currentlevel[DIM-1]=j0[DIM-1]+currentlevel[0]-j0[0]+1; currenttype[DIM-1]=1;
				currenttype[0]=0; currentlevel[0]=j0[0];
                                sizes[DIM-1][range+1]=basis ->bases()[DIM-1]->Nablasize(currentlevel[DIM-1]); // if needed compute and store new size information
                            }
  			break; // unnoetig, da i==0 gilt.
  			}
                    } // end of "big loop"
                    }
                    //cout << "(currentlevel, currenttype), (j_, e_) = ("<< currentlevel << ", " << currenttype << "), ("<<j_<<", "<<e_<<")"<< endl;
                    //cout << "(currentlevel == j_ == " << (currentlevel == j_) << "currenttype == e_ == " << (currenttype == e_) << endl;
                    if (currentlevel == j_ && currenttype == e_) break;
		}
                // determine number of the function described by k in (currentlevel,currenttype)
                //cout << "currentlevel, currenttype = "<< currentlevel << ", " << currenttype << endl;
		translation_type ktemp; // we only use that this is a multiindex
		for (int i = 0; i < DIM; i++)
                    if (e_[i] == 0) ktemp[i]=k_[i]-basis_->bases()[i]->DeltaLmin();
                    else ktemp[i]=k_[i]-basis_->bases()[i]->Nablamin();
                //cout << "ktemp "<<ktemp<<endl;
                oncurrentlevel = ktemp[0];
		for (int i=1; i <= DIM-1; i++)
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

	// ACHTUNG k_ wird mit 0 initialisiert, nicht mit DeltaLmin, Nablamin!!

	template <class IBASIS, unsigned int DIM, class TENSORBASIS>
	TensorIndex<IBASIS, DIM, TENSORBASIS>::TensorIndex(const int number, const TENSORBASIS* basis)
	: basis_(basis), num_(number)
	{
                MathTL::FixedArray1D<std::map<int, int>,DIM> sizes; // Store number of basis elements. Generators on level j0 (0), wavelets on level j0 (1), j0+1 (2), ...
		int remains = number+1; // numbering begins at 0
		int uptothislevel(0); // number of basis function below the current level
		int oncurrentlevel(1); // number of base elements on current level j
		int range(0); // number of level steps we have climbed so far. determines wether a new entry has to be computed in "sizes"
		level_type currentlevel, j0(basis_->j0());
		type_type currenttype;
		currentlevel = j0;
		for (unsigned int i=0; i < DIM; i++) {
			currenttype[i]=0;
			sizes[i][0] = basis_->bases()[i]->Deltasize(j0[i]); // Number of generators on level j0
                        sizes[i][1] = basis_->bases()[i]->Nablasize(j0[i]); // N o Wavelets
                        // PROBLEM! Hier kommt je nach verwendeter Basis nicht das selbe raus!
                        //cout << "sizes["<<i<<"][0]" << sizes[i][0] << endl;
                        //cout << "sizes["<<i<<"][1]" << sizes[i][1] << endl;
                        oncurrentlevel *= sizes[i][0];
		}

                // iterate over all levels. Add up number of basis functions till looked for level is reached
		while (remains > oncurrentlevel) // break if we are at the right level
		{
                    // else substract number of basis functions on current_index=(currentlevel,currentindex) and increase current_index
                    //cout << "st of while ct= "<<currenttype<<" cl"<<currentlevel<<endl;
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
                            if (currenttype[i] == 1)
                            {
                                currenttype[i]=0;
                                //if (i == 0) done = true;
                            } else
                            {
                                currenttype[i]=1;
                                done = false;
                                break;
                            }
                    }
                    // done == true bedeutet, dass alle Komponenten auf level j0() wavelets waren.
                    // "big loop" "currentlevel++"
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
                                if (currentlevel[i-1]-j0[i] == range) 
                                {
                                    sizes[i-1][range+1]=basis ->bases()[i]->Nablasize(currentlevel[i-1]); // if needed compute and store new size information
                                    //cout << "sizes[" << i-1 << "][" << range+1 << "]" << sizes[i-1][range+1] << endl;
                                }
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
                                sizes[0][range]=basis ->bases()[0]->Nablasize(currentlevel[0]); // if needed compute and store new size information
                            }
                            else
                            {
                                //cout << "ct= "<<currenttype<<" cl"<<currentlevel<<endl;
                                currentlevel[DIM-1]=j0[DIM-1]+currentlevel[0]-j0[0]+1; currenttype[DIM-1]=1;
				currenttype[0]=0; currentlevel[0]=j0[0];
                                sizes[DIM-1][range+1]=basis ->bases()[DIM-1]->Nablasize(currentlevel[DIM-1]); // if needed compute and store new size information
                                //cout << "ct= "<<currenttype<<" cl"<<currentlevel<<endl;
                                //cout << "sizes["<<DIM-1<<"]["<<range+1<<"]="<<sizes[DIM-1][range+1]<<endl;
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
                    //cout << "eowhile ct= "<<currenttype<<" cl"<<currentlevel<<endl;
		} // end of while

                // determine k corresponding to the number of the basis function given by "remains" (on level (currentlevel,currenttype) )
		unsigned int modul;                
		j_ = currentlevel;
		e_ = currenttype;
                //cout << "ct= "<<currenttype<<" cl"<<currentlevel<<endl;
                //cout << "j_= "<<j_<<" e_= "<<e_<<endl;
                remains -= 1; // numbering begins at 0
		for (int i = DIM-1; i > 0; i--)
		{
			modul = sizes[i][currentlevel[i]+currenttype[i]-j0[i]];
			k_[i]= remains % modul+(e_[i] == 0 ? basis->bases()[i]->DeltaLmin():basis->bases()[i]->Nablamin());
			remains = remains / modul;
		}
                k_[0]=remains +(e_[0] == 0 ? basis->bases()[0]->DeltaLmin():basis->bases()[0]->Nablamin());
        }
/*
            MathTL::FixedArray1D<std::map<int, int>,DIM> sizes; // Store number of basis elements. Generators on level j0 (0), wavelets on level j0 (1), j0+1 (2), ...
		int remains = number+1; // numbering begins at 0
		int uptothislevel(0); // number of basis function below the current level
		int oncurrentlevel(1); // number of base elements on current level j
		int range(0); // number of level steps we have climbed so far. determines wether a new entry has to be computed in "sizes"
		level_type currentlevel, j0(basis_->j0());
		type_type currenttype;
		currentlevel = j0;
		for (unsigned int i=0; i < DIM; i++) {
			currenttype[i]=0;
			sizes[i][0] = basis_->bases()[i]->Deltasize(j0[i]); // Number of generators on level j0
                        // PROBLEM! Hier kommt je nach verwendeter Basis nicht das selbe raus!
                        // cout << "sizes["<<i<<"][0]" << sizes[i][0] << endl;
			oncurrentlevel *= sizes[i][0];
		}
		// iterate over all levels. Add up number of basis functions till looked for level is reached
		while (true)
		{
                    if (remains <= oncurrentlevel) break; // we are in the right level
                    else //increase sublevel (j++)
                    {
                        remains -= oncurrentlevel;
                        // determine next sublevel
  			for (int i(DIM-1); i >= 0; i--)
  			{
  				if (i != 0)
  				{
  					if (currenttype[i] != 0)
  					{
                                                // if possible increase index "left" from the current one. if there is no index on the "left" we are at the end of this level
		  				if (currenttype[i-1] == 0) currenttype[i-1] = 1;
  						else currentlevel[i-1]=currentlevel[i-1]+1;
  						int temp = currentlevel[i]-j0[i];
  						currentlevel[i]=0;
  						currenttype[i]=0;
  						if (temp == 0)
  						{
		  					currentlevel[DIM-1]=j0[DIM-1];
  							currenttype[DIM-1]=0;
  						} else
  						{
		  					currentlevel[DIM-1]=j0[DIM-1]+temp-1;
  							currenttype[DIM-1]=1;
  						}
  						break;
  					}
  				} else
  				{
  					range = range +1;
		  			if (DIM == 1)
		  				if (currenttype[i] == 0) currenttype[i] = 1;
		  				else currentlevel[i]=currentlevel[i]+1;
		  			else
  					{
		  				currentlevel[DIM-1]=j0[DIM-1]+currentlevel[0]-j0[0]+currenttype[0]; currenttype[DIM-1]=1;
						currenttype[0]=0; currentlevel[0]=j0[0];
  					}
  					break; // eigentlich unnoetig, da wir schon bei i=0 sind.
  				}
  			}
  			// compute number of basis functions on current level
  			oncurrentlevel = 1;
 			for (unsigned int i=0; i < DIM; i++) {
                            
                            cout << "currentlevel[i]" << currentlevel[i] <<endl
                            << " j0[i]" << j0[i]<<endl
                            <<" currenttype[i]"<< currenttype[i] <<endl
                            << "sizes[i][range]"<<sizes[i][range] <<endl
                            <<"sizes[i][range+1]"<<sizes[i][range+1] <<endl
                            << " basis_->bases()[i]->Nablasize(j0[i]+range)"<<basis_->bases()[i]->Nablasize(j0[i]+range)<<endl;
                            if (currentlevel[i]+currenttype[i]-j0[i] == range) sizes[i][range] = basis_->bases()[i]->Nablasize(j0[i]+range-1); // store new information
                            oncurrentlevel *= sizes[i][currentlevel[i]+currenttype[i]-j0[i]];
			}			
                    }
		}
                // now currentlevel is the looked for level!
		unsigned int modul;
		j_ = currentlevel;
		e_ = currenttype;
                remains -= 1; // this is the number in the current sublevel
		for (int i = DIM-1; i > 0; i--)
		{
			modul = sizes[i][currentlevel[i]+currenttype[i]-j0[i]];
			k_[i]= remains % modul;
			remains = remains / modul;
		}
                k_[0]=remains;
	}
*/
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

        // eventuell die vielen Aufrufe von bases_->j0() durch Hilfsvariable auf einen reduzieren
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
                        //range = range +1;
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
/*
		for (int i(DIM-1); i >= 0; i--)
  		{
                    if (i != 0)
                    {
                        if (e_[i] != 0)
  			{
                            if (e_[i-1] == 0) e_[i-1] = 1;
                            else j_[i-1]=j_[i-1]+1;
                            int temp = j_[i]-basis_->j0()[i];
                            j_[i]=0;
                            e_[i]=0;
                            if (temp == 0)
                            {
                                j_[DIM-1]=basis_->j0()[DIM-1];
  				e_[DIM-1]=0;
                            } else
                            {
                                j_[DIM-1]=basis_->j0()[DIM-1]+temp-1;
  				e_[DIM-1]=1;
                            }
                            break;
  			}
                    } else
                    {
                        if (DIM == 1)
                            if (e_[i] == 0) e_[i] = 1;
                            else j_[i]=j_[i]+1;
		 	else
                        {
                            //int temp = j_[0]-basis_->j0()[i];
                            j_[DIM-1]=basis_->j0()[DIM-1]+j_[0]-basis_->j0()[0]+e_[0]; e_[DIM-1]=1;
                            e_[0]=0; j_[0]=basis_->j0()[0];
  			}
                        break;
                    }
  		}
  		// then choose lowest translation index k=k(j,e)
		for (unsigned int i = 0; i < DIM; i++)
			k_[i] = (e_[i] == 0
				? basis_->bases()[i]->DeltaLmin()
	   			: basis_->bases()[i]->Nablamin());
  		return *this;
                */
  	}

	template <class IBASIS, unsigned int DIM, class TENSORBASIS>
	bool
	TensorIndex<IBASIS,DIM,TENSORBASIS>::operator < (const TensorIndex& lambda) const
	{
            // Ordering by level (1-norm of) j
            // then standard lexicographic order on (j,e,k),
            // we assume that e and k are already lexicographically ordered (cf. MultiIndex)
            return ( multi_degree(j_) < multi_degree(lambda.j()) ||
                    (multi_degree(j_) == multi_degree(lambda.j()) &&
                     (j_ < lambda.j() ||
                      (j_ == lambda.j() &&
                       (e_ < lambda.e() ||
                        (e_ == lambda.e() && k_ < lambda.k())
                       )
                      )
                     )
                    )
                   );
/*
            return (j_ < lambda.j() ||
                    (
                        j_ == lambda.j() &&
                        (e_ < lambda.e() || (e_ == lambda.e() && k_ < lambda.k()))
                    )
                    );
             * */
	}

	// numbering begins at 0?
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
            //cout << "in first_wavelet. j,e,k = "<<j<<", "<<e<<" ,"<<k<<endl;
            return TensorIndex<IBASIS,DIM,TENSORBASIS>(j, e, k, basis);
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
}
