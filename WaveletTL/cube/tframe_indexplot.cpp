// implementation for tframe_indexplot.h

namespace WaveletTL
{
    template <class TENSORFRAME>
    void plot_indices(const TENSORFRAME* frame,
            const InfiniteVector<double, typename TENSORFRAME::Index>& coeffs,
            const int maxrange,
            std::ostream& os,
            const typename TENSORFRAME::Index::polynomial_type p,
            const char* colormap,
            bool boxed,
            bool colorbar,
            const double aa)
    {
        typedef typename TENSORFRAME::Index Index;
//        typedef typename TENSORFRAME::IntervalFrame Frame1D;
//        typedef typename Index::polynomial_type polynomial_type;
        typedef typename Index::level_type level_type;
        typedef typename Index::type_type type_type;
        typedef typename Index::translation_type translation_type;
        
        const int DIM = 2;
        /*
        // we dont want to be limited to the 2 dimensional case!
        const Frame1D* xframe = frame->frames()[0];
        const Frame1D* yframe = frame->frames()[1];
         */

        const level_type j0 = frame->j0();
        const double maxnorm = linfty_norm(coeffs);

        // determine number of sublevels on each level
        // this gives the number of rows in the multiplot
        //FixedArray1D<int, maxrange+1> number_of_rows;
        Array1D<int> number_of_rows;
        number_of_rows.resize(maxrange+1);
        number_of_rows[0]=2;
        for (int i=1;i<=maxrange;i++)
            {
                number_of_rows[i]=1;
            }
        for (int d=1; d<DIM; d++)
        {
            Array1D<int> temp_nor;
            temp_nor.resize(maxrange+1);
            temp_nor[0]= 1<<(d+1);
            for (int i=1;i<=maxrange;i++)
            {
                temp_nor[i]=0;
                for (int j=0; j<i;j++)
                {
                    temp_nor[i]+=number_of_rows[j];
                }
                temp_nor[i]+=2*number_of_rows[i];
            }
            for (int i=0;i<=maxrange;i++)
                number_of_rows[i]=temp_nor[i];
        }

        const double threshold = 1e-15;
        cout << "linfinity norm of coefficients = " << maxnorm << endl;
        cout << "number of coefficients" << coeffs.size() << endl;

        // initialize m-file
        os << "colormap(" << colormap << ");" << endl;
        os << "X=" << colormap << ";" << endl;

        // iterate over all coeffs.        
        // if ++lambda leads to an increased sublevel we need to plot a new row
        // if ++lambda leads to an increased level we need to start a new column
        MathTL::FixedArray1D<std::map<int, int>,DIM> sizes;
        level_type currentlevel;
        type_type currenttype;
        translation_type current_k;
        for (int i=0; i < DIM; i++)
        {
            currentlevel[i]=j0[i];
            currenttype[i]=0;
            current_k[i]=frame->frames()[i]->DeltaLmin();
            sizes[i][0] = frame->frames()[i]->Deltasize(j0[i]); // Number of generators on level j0
            sizes[i][1] = frame->frames()[i]->Nablasize(j0[i]); // Number of Wavelets on level j0
        }

        int row(1),column(1);
        bool atmaxrange = false;
        int range(0); // determines wether a new entry has to be computed in "sizes"

        while (!atmaxrange)
        {
            // start a new box
            //os << "subplot("<< number_of_rows[maxrange]<<"," << (maxrange +1) << ","<< (((boxstep-1)% number_of_rows[maxrange])*(maxrange+1)+ ((boxstep-1)/ number_of_rows[maxrange]) +1)<<");" << endl;
            os << "subplot("<< number_of_rows[maxrange]<<"," << (maxrange +1) << ","<< ((row-1)*(maxrange+1)+column)<<");" << endl;
            os << "box on" << endl;
            os << "set(gca,'Layer','top')" << endl;
            // turn off x ticks
            os << "set(gca,'XTick',[])" << endl;
            // turn off y ticks
            os << "set(gca,'YTick',[])" << endl;
            os << "ylabel('"<<currentlevel<<", "<<currenttype << ":               ','rotation',0)" << endl;
            os << "title 'level " << multi_degree(currentlevel) << "'" << endl;
            os << "axis([0,1,0,1]);" << endl;

            // based on tframe_index :: operator ++
            FixedArray1D<double,DIM> pos;
            bool jplusplus = false;
            while (!jplusplus)
            {
                // determine current position
                for (int i=0; i<DIM;i++)
                    pos[i] = (current_k[i]-((currenttype[i] == 0) ? frame->frames()[i]->DeltaLmin() : frame->frames()[i]->Nablamin()))
                             * (1. / sizes[i][currentlevel[i]-j0[i]+currenttype[i]]);
                const double c = coeffs.get_coefficient(Index(p,currentlevel, currenttype, current_k, frame));
                if (fabs(c) < threshold) {
                } else
                {
                    double col = (std::max(log10(fabs(c)/maxnorm),aa)+(-aa))/(-aa);
                    // this part is for 2 space dimensions
                    os << "Y=ceil(" << col << " * length(colormap));" << endl;
                    os << "if Y==0" << endl << "Y=Y+1;" << endl << "end;" << endl;
                    os << "patch(" << "[" << pos[0] << ";" << pos[0] + (1. / sizes[0][currentlevel[0]-j0[0]+currenttype[0]]) << ";" << pos[0] + (1. / sizes[0][currentlevel[0]-j0[0]+currenttype[0]]) << ";" << pos[0] << "],"
                       << "[" << pos[1] << ";" << pos[1] << ";" << pos[1] + (1. / sizes[1][currentlevel[1]-j0[1]+currenttype[1]]) << ";" << pos[1] + (1. / sizes[1][currentlevel[1]-j0[1]+currenttype[1]]) << "],"
                       << "X(Y,:)";
                    if (!boxed)
                        os << ",'EdgeColor','none'";
                    os
                            // 	    << ",'LineStyle','none'"
                            // 	    << ",'LineWidth',0.125"
                            << ")" << endl;
                }

                // determine next translation index
                for (int i = DIM-1; i >= 0; i--) {
                    const int last_index = (currenttype[i] == 0 ? frame->frames()[i]->DeltaRmax(currentlevel[i])
                                                                : frame->frames()[i]->Nablamax(currentlevel[i]));
                    if (current_k[i] == last_index)
                    {
                        current_k[i] = (currenttype[i] == 0 ? frame->frames()[i]->DeltaLmin()
                                : frame->frames()[i]->Nablamin());
                        jplusplus = (i == 0);
                    } else
                    {
                        ++current_k[i];
                        break;
                    }
                }
            } // end of while (!jplusplus)
            // current sublevel has been drawn. Increase sublevel or level if needed

            // determine next (sub)level index
            // "small loop" "currenttype++" (currentlevel is fixed)
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
                        current_k[i]=frame->frames()[i]->DeltaLmin();
                    } else
                    {
                        currenttype[i]=1;
                        current_k[i]=frame->frames()[i]->Nablamin();
                        row++;
                        done = false;
                        break;
                    }
                }
            }
            // done == true bedeutet, dass alle Komponenten auf level j0() wavelets waren.
            // "big loop" "currentlevel++"
            if (done == true)
            {
                for (int i(DIM-1); i >= 0; i--)
                {
                    if (i != 0)
                    { // try to increase sublevel
                        if (currentlevel[i] != j0[i])
                        {
                            // increase left neighbor
                            currentlevel[i-1]=currentlevel[i-1]+1;
                            if (currentlevel[i-1]-j0[i] == range) sizes[i-1][range+1]=frame ->frames()[i]->Nablasize(currentlevel[i-1]); // if needed compute and store new size information
                            currenttype[i-1]=1;
                            current_k[i-1]=frame->frames()[i-1]->Nablamin();
                            int temp = currentlevel[i]-j0[i];
                            currentlevel[i]=j0[i];
                            currenttype[i]=0;
                            current_k[i]=frame->frames()[i]->DeltaLmin();
                            currentlevel[DIM-1]=j0[DIM-1]+temp-1;
                            currenttype[DIM-1]= (temp == 1?0:1);
                            current_k[DIM-1]= (temp == 1?frame->frames()[i]->DeltaLmin():frame->frames()[i]->Nablamin());
                            row++;
                            break;
                        }
                    } else // i == 0. "big loop" arrived at the last index. We have to increase the level!
                    {
                        range = range +1;
                        row=1;
                        column++;
                        if (DIM == 1)
                        {
                            currenttype[i] = 1; // diese Zeile erfüllt nur in der allerersten Iteration einen Zweck
                            current_k[i]=frame->frames()[i]->Nablamin(); // diese Zeile erfüllt nur in der allerersten Iteration einen Zweck
                            currentlevel[i]=currentlevel[i]+1;
                            sizes[0][range]=frame->frames()[0]->Nablasize(currentlevel[0]); // if needed compute and store new size information
                        }
                        else
                        {
                            currentlevel[DIM-1]=j0[DIM-1]+currentlevel[0]-j0[0]+1;
                            currenttype[DIM-1]=1;
                            current_k[DIM-1]=frame->frames()[i]->Nablamin();
                            currentlevel[0]=j0[0];
                            currenttype[0]=0;
                            current_k[0]=frame->frames()[i]->DeltaLmin();
                            sizes[DIM-1][range+1]=frame ->frames()[DIM-1]->Nablasize(currentlevel[DIM-1]); // if needed compute and store new size information
                        }
                        atmaxrange = (range > maxrange);
                        break; // unnoetig, da i==0 gilt.
                    }
                }
            } // end of "big loop"
        } // end of while(!atmaxrange)
    }
}
