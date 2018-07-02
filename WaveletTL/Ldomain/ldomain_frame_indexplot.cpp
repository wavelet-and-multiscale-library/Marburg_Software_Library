// implementation for tframe_indexplot.h

namespace WaveletTL
{
    
    template <class LDOMAINFRAME>
    void plot_indices_ldomain(const LDOMAINFRAME* frame,
                       const InfiniteVector<double, typename LDOMAINFRAME::Index>& coeffs,
                       std::ostream& os,
                       const typename LDOMAINFRAME::Index::polynomial_type p,
                       const typename LDOMAINFRAME::Index::level_type j,
                       const typename LDOMAINFRAME::Index::type_type e,
                       const char* colormap,
                       bool boxed,
                       bool colorbar,
                       const double lowerclim)
    {
        typedef typename LDOMAINFRAME::Index Index;
        
        const double maxnorm = linfty_norm(coeffs);
        
        //setup figure and properties
        //os<<"figure;"<<endl;
        os<<"clf;"<<endl;
        //os<<"box on;"<<endl;
        os << "colormap( " << colormap<<");" << endl;
//        os << "set(gca,'CLim',[" << lowerclim << " 0])" << endl;
        //convention: we denote generators as wavelets on level j_0 -1
        
        //set size of subplots
        os<<"axlabels='on';"<<endl;
        os<<"psize=0.35;"<<endl;
        os<<"smallpsize=0.05;"<<endl;
        os<<"imargin=0.025;"<<endl;
        os<<"omargin=0.05;"<<endl;
        
        //empty subplot for centered title
//        os<<"subplot(3,3,2)"<<endl;
        os<<"axis('off')"<<endl;
        os<<"title(sprintf('coefficients on degree (%i,%i) and level (%i,%i)',"<<p[0]<<","<<p[1]<<","<<j[0]-1+e[0]<<","<<j[1]-1+e[1]<<"))"<<endl;
               
        int startx, starty, columns, rows;
        Vector<int> patches(5,"3 4 0 1 2"); //order of subplots is important to octave plot
        for(Vector<int>::const_iterator it(patches.begin()); it!=patches.end(); ++it ){
            int patch=*it;
            //set subplot properties for patches
            switch(patch){
                case 0:
                    os<<"#patch 0:"<<endl;
                    os<<"axes('position',[omargin omargin+2*imargin+psize+smallpsize psize psize])"<<endl; 
                    startx= e[0] ? frame->frame1d_11()->first_wavelet(j[0],p[0]).k() : frame->frame1d_11()->first_generator(j[0],p[0]).k();
                    starty= e[1] ? frame->frame1d_01()->first_wavelet(j[1],p[1]).k()+1 : frame->frame1d_01()->first_generator(j[1],p[1]).k()+1; //shift because first generator/wavelet lies on patch 3
                    columns= e[0] ? frame->frame1d_11()->Nablasize(j[0],p[0]) : frame->frame1d_11()->Deltasize(j[0],p[0]);
                    rows=e[1] ? frame->frame1d_01()->Nablasize(j[1],p[1])-1 : frame->frame1d_01()->Deltasize(j[1],p[1])-1;//shift because first generator/wavelet lies on patch 3
                    //os<<"axis(axlabels);"<<endl;
                    os<<"axis(["<<(double)startx-0.5<<" "<<(double)startx-0.5+columns<<" "<<(double)starty-0.5<<" "<<(double)starty-0.5+rows<<"],axlabels);"<<endl;
                    os<<"set(gca,'xaxislocation','top')"<<endl;
                    break;
                case 1:
                    os<<"#patch 1:"<<endl;
                    os<<"axes('position',[omargin omargin psize psize])"<<endl;
                    startx= e[0] ? frame->frame1d_11()->first_wavelet(j[0],p[0]).k() : frame->frame1d_11()->first_generator(j[0],p[0]).k();
                    starty= e[1] ? frame->frame1d_11()->first_wavelet(j[1],p[1]).k() : frame->frame1d_11()->first_generator(j[1],p[1]).k();
                    columns= e[0] ? frame->frame1d_11()->Nablasize(j[0],p[0]) : frame->frame1d_11()->Deltasize(j[0],p[0]);
                    rows=e[1] ? frame->frame1d_11()->Nablasize(j[1],p[1]) : frame->frame1d_11()->Deltasize(j[1],p[1]);
                    os<<"axis(["<<(double)startx-0.5<<" "<<(double)startx-0.5+columns<<" "<<(double)starty-0.5<<" "<<(double)starty-0.5+rows<<"],axlabels);"<<endl;
                    break;
                case 2:
                    os<<"#patch 2:"<<endl;
                    os<<"axes('position',[omargin+2*imargin+psize+smallpsize omargin psize psize])"<<endl;
                    
                    startx= e[0] ? frame->frame1d_01()->first_wavelet(j[0],p[0]).k()+1 : frame->frame1d_01()->first_generator(j[0],p[0]).k()+1;   //shift because first generator/wavelet lies on patch 4
                    starty= e[1] ? frame->frame1d_11()->first_wavelet(j[1],p[1]).k() : frame->frame1d_11()->first_generator(j[1],p[1]).k(); 
                    columns= e[0] ? frame->frame1d_01()->Nablasize(j[0],p[0])-1 : frame->frame1d_01()->Deltasize(j[0],p[0])-1;
                    rows=e[1] ? frame->frame1d_11()->Nablasize(j[1],p[1]) : frame->frame1d_11()->Deltasize(j[1],p[1]);
                    os<<"axis(["<<(double)startx-0.5<<" "<<(double)startx-0.5+columns<<" "<<(double)starty-0.5<<" "<<(double)starty-0.5+rows<<"],axlabels);"<<endl;
                    os<<"set(gca,'yaxislocation','right')"<<endl;
                    break;
                case 3:
                    os<<"#patch 3:"<<endl;
                    os<<"axes('position',[omargin omargin+psize+imargin psize smallpsize])"<<endl;                    
                    startx= e[0] ? frame->frame1d_11()->first_wavelet(j[0],p[0]).k() : frame->frame1d_11()->first_generator(j[0],p[0]).k();
                    starty= e[1] ? frame->frame1d_01()->first_wavelet(j[1],p[1]).k() : frame->frame1d_01()->first_generator(j[1],p[1]).k();
                    columns= e[0] ? frame->frame1d_11()->Nablasize(j[0],p[0]) : frame->frame1d_11()->Deltasize(j[0],p[0]);
                    rows=1;
                    os<<"axis(["<<(double)startx-0.5<<" "<<(double)startx-0.5+columns<<" "<<(double)starty-0.5<<" "<<(double)starty-0.5+rows<<"],'nolabel','ticx');"<<endl;
                    break;
                case 4:
                    os<<"#patch 4:"<<endl;
                    os<<"axes('position',[omargin+psize+imargin omargin smallpsize psize])"<<endl;       
                    startx= e[0] ? frame->frame1d_01()->first_wavelet(j[0],p[0]).k() : frame->frame1d_01()->first_generator(j[0],p[0]).k();
                    starty= e[1] ? frame->frame1d_11()->first_wavelet(j[1],p[1]).k() : frame->frame1d_11()->first_generator(j[1],p[1]).k();
                    columns=1;
                    rows= e[1] ? frame->frame1d_11()->Nablasize(j[1],p[1]) : frame->frame1d_11()->Deltasize(j[1],p[1]);
                    os<<"axis(["<<(double)startx-0.5<<" "<<(double)startx-0.5+columns<<" "<<(double)starty-0.5<<" "<<(double)starty-0.5+rows<<"],'nolabel','ticy');"<<endl;
                    break;
            }
            os<<"box on;"<<endl;
            //plot coefficients
            for (typename InfiniteVector<double, Index>::const_iterator it(coeffs.begin()); it != coeffs.end(); ++it){
                Index ind=it.index();
                if(ind.j()==j && ind.e()==e && ind.p()==p && ind.patch()==patch){              
                    const double ax=(double)ind.k()[0]-0.5;
                    const double bx=(double)ind.k()[0]+0.5;
                    const double ay=(double)ind.k()[1]-0.5;
                    const double by=(double)ind.k()[1]+0.5;
                    const double val=std::max(log10(fabs(*it)/maxnorm),lowerclim);
                    //const double val=log10(fabs(*it)/maxnorm);
                    os << "patch(["<<ax<<","<<bx<<","<<bx<<","<<ax<<"],["<<ay<<","<<ay<<","<<by<<","<<by<<"],"<<val;
                    if(!boxed){
                        os<<",'edgecolor','none'";
                    }
                    os<<");"<<endl;
                }           
            }
        }
        
        //correct colorbar labels
        if (colorbar) {
            //empty subplot for colorbar
            os<<"axes('position',[omargin+2*imargin+psize+smallpsize omargin+2*imargin+psize+smallpsize psize psize])"<<endl;
            os<<"axis('off');"<<endl;
            os << "set(gca,'CLim',[" << lowerclim << " 0])" << endl;
            os<<"h=colorbar('position',[0.9 0.05 0.025 0.85]);"<<endl;
            //os << "h=colorbar;" << endl;
            os << "ytick=get(h,'ytick');"<<endl;
            os << "labels = {};"<<endl;
            os << "for v=ytick(1:end), labels{end+1} = sprintf('10^%i',v); end"<<endl;
            os << "set(h, 'yticklabel',labels);"<<endl;
        }
    }
}
