// implementation for recring_frame_indexplot.h

namespace WaveletTL
{
    
    template <class RECRINGFRAME>
    void plot_indices_recring(const RECRINGFRAME* frame,
                       const InfiniteVector<double, typename RECRINGFRAME::Index>& coeffs,
                       std::ostream& os,
                       const typename RECRINGFRAME::Index::polynomial_type p,
                       const typename RECRINGFRAME::Index::level_type j,
                       const typename RECRINGFRAME::Index::type_type e,
                       const char* colormap,
                       bool boxed,
                       bool colorbar,
                       const double lowerclim)
    {
        typedef typename RECRINGFRAME::Index Index;
        
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
        os<<"psize=0.2;"<<endl;
        os<<"smallpsize=0.05;"<<endl;
        os<<"imargin=0.025;"<<endl;
        os<<"omargin=0.05;"<<endl;
        
        //empty subplot for centered title
//        os<<"subplot(3,3,2)"<<endl;
        os<<"axis('off')"<<endl;
        os<<"title(sprintf('coefficients on degree (%i,%i) and level (%i,%i)',"<<p[0]<<","<<p[1]<<","<<j[0]-1+e[0]<<","<<j[1]-1+e[1]<<"))"<<endl;
          
        
        
        int startx, starty, columns, rows;
//        Vector<int> patches(7,"4 5 6 0 1 2 3"); //order of subplots is important to octave plot
//        for(Vector<int>::const_iterator it(patches.begin()); it!=patches.end(); ++it ){
        for(int patch=0; patch<=15; ++patch ){

//            int patch=*it;
            //set subplot properties for patches
            switch(patch){
                case 1:
                case 5:  
                    //new: axes seems to be better than subplot
//                    os<<"subplot(3,3,1,'position',[omargin omargin+2*imargin+psize+smallpsize psize psize])"<<endl; 
                    switch(patch){
                       case 1:
                           os<<"#patch 1:"<<endl;
                           os<<"axes('position',[omargin omargin+2*imargin+psize+smallpsize psize psize])"<<endl;
                           break;
                       case 5:
                           os<<"#patch 5:"<<endl;
                           os<<"axes('position',[omargin+4*imargin+2*psize+2*smallpsize omargin+2*imargin+psize+smallpsize psize psize])"<<endl;
                           break;
                    } 
                    startx= e[0] ? frame->frame1d_11()->first_wavelet(j[0],p[0]).k() : frame->frame1d_11()->first_generator(j[0],p[0]).k();
                    starty= e[1] ? frame->frame1d()->first_wavelet(j[1],p[1]).k()+1 : frame->frame1d()->first_generator(j[1],p[1]).k()+1; //shift because first generator/wavelet lies on patch 3
                    columns= e[0] ? frame->frame1d_11()->Nablasize(j[0],p[0]) : frame->frame1d_11()->Deltasize(j[0],p[0]);
                    rows=e[1] ? frame->frame1d()->Nablasize(j[1],p[1])-2 : frame->frame1d()->Deltasize(j[1],p[1])-2;//shift because first generator/wavelet lies on patch 3
                    //os<<"axis(axlabels);"<<endl;
                    if(patch==5){
                        os<<"axis(["<<(double)startx-0.5<<" "<<(double)startx-0.5+columns<<" "<<(double)starty-0.5<<" "<<(double)starty-0.5+rows<<"],'nolabel');"<<endl;
                    }
                    else{
                        os<<"axis(["<<(double)startx-0.5<<" "<<(double)startx-0.5+columns<<" "<<(double)starty-0.5<<" "<<(double)starty-0.5+rows<<"],axlabels);"<<endl;
                    }
                    os<<"set(gca,'XTickLabel',{})"<<endl;
                    break;
                case 0:
                case 2:
                case 4:
                case 6:    
//                    os<<"subplot(3,3,7,'position',[omargin omargin psize psize])"<<endl;
                    switch(patch){
                       case 0:
                           os<<"#patch 0:"<<endl;
                           os<<"axes('position',[omargin omargin+4*imargin+2*psize+2*smallpsize psize psize])"<<endl;
                           break;
                       case 2:
                           os<<"#patch 2:"<<endl;
                           os<<"axes('position',[omargin omargin psize psize])"<<endl;
                           break;
                       case 4:
                           os<<"#patch 4:"<<endl;
                           os<<"axes('position',[omargin+4*imargin+2*psize+2*smallpsize omargin psize psize])"<<endl;
                           break;
                       case 6:
                           os<<"#patch 6:"<<endl;
                           os<<"axes('position',[omargin+4*imargin+2*psize+2*smallpsize omargin+4*imargin+2*psize+2*smallpsize psize psize])"<<endl;
                           break;    
                    }
//                    os<<"axes('position',[omargin omargin psize psize])"<<endl;
                    startx= e[0] ? frame->frame1d_11()->first_wavelet(j[0],p[0]).k() : frame->frame1d_11()->first_generator(j[0],p[0]).k();
                    starty= e[1] ? frame->frame1d_11()->first_wavelet(j[1],p[1]).k() : frame->frame1d_11()->first_generator(j[1],p[1]).k();
                    columns= e[0] ? frame->frame1d_11()->Nablasize(j[0],p[0]) : frame->frame1d_11()->Deltasize(j[0],p[0]);
                    rows=e[1] ? frame->frame1d_11()->Nablasize(j[1],p[1]) : frame->frame1d_11()->Deltasize(j[1],p[1]);
                    if(patch==6){
                        os<<"axis(["<<(double)startx-0.5<<" "<<(double)startx-0.5+columns<<" "<<(double)starty-0.5<<" "<<(double)starty-0.5+rows<<"],'nolabel');"<<endl;
                    }
//                    else if(patch==0){
//                        os<<"axis(["<<(double)startx-0.5<<" "<<(double)startx-0.5+columns<<" "<<(double)starty-0.5<<" "<<(double)starty-0.5+rows<<"],ylabel);"<<endl;
//                    }
                    else{
                       os<<"axis(["<<(double)startx-0.5<<" "<<(double)startx-0.5+columns<<" "<<(double)starty-0.5<<" "<<(double)starty-0.5+rows<<"],axlabels);"<<endl; 
                    }
                    if(patch==0){
                        os<<"set(gca,'XTickLabel',{})"<<endl;
                    } 
                    if(patch==4){
                        os<<"set(gca,'YTickLabel',{})"<<endl;
                    } 
                    break;
                case 3:
                case 7:    
//                    os<<"subplot(3,3,9,'position',[omargin+2*imargin+psize+smallpsize omargin psize psize])"<<endl;
                    switch(patch){
                       case 3:
                           os<<"#patch 3:"<<endl;
                           os<<"axes('position',[omargin+2*imargin+psize+smallpsize omargin psize psize])"<<endl;
                           break;
                       case 7:
                           os<<"#patch 7:"<<endl;
                           os<<"axes('position',[omargin+2*imargin+psize+smallpsize omargin+4*imargin+2*psize+2*smallpsize psize psize])"<<endl;
                           break;
                    }
                    startx= e[0] ? frame->frame1d()->first_wavelet(j[0],p[0]).k()+1 : frame->frame1d()->first_generator(j[0],p[0]).k()+1;   //shift because first generator/wavelet lies on patch 4
                    starty= e[1] ? frame->frame1d_11()->first_wavelet(j[1],p[1]).k() : frame->frame1d_11()->first_generator(j[1],p[1]).k(); 
                    columns= e[0] ? frame->frame1d()->Nablasize(j[0],p[0])-2 : frame->frame1d()->Deltasize(j[0],p[0])-2;
                    rows=e[1] ? frame->frame1d_11()->Nablasize(j[1],p[1]) : frame->frame1d_11()->Deltasize(j[1],p[1]);
                    if(patch==7){
                        os<<"axis(["<<(double)startx-0.5<<" "<<(double)startx-0.5+columns<<" "<<(double)starty-0.5<<" "<<(double)starty-0.5+rows<<"],'nolabel');"<<endl;
                    }
                    else{
                        os<<"axis(["<<(double)startx-0.5<<" "<<(double)startx-0.5+columns<<" "<<(double)starty-0.5<<" "<<(double)starty-0.5+rows<<"],axlabels);"<<endl;
                        os<<"set(gca,'YTickLabel',{})"<<endl;
                    }                    
                    break;
                case 8:
                case 13:    
//                    os<<"subplot(3,3,3,'position',[omargin+2*imargin+psize+smallpsize omargin+2*imargin+psize+smallpsize  psize psize])"<<endl; 
                    switch(patch){
                       case 8:
                           os<<"#patch 8:"<<endl;
                           os<<"axes('position',[omargin omargin+3*imargin+2*psize+smallpsize psize smallpsize])"<<endl;
                           break;
                       case 13:
                           os<<"#patch 13:"<<endl;
                           os<<"axes('position',[omargin+4*imargin+2*psize+2*smallpsize omargin+3*imargin+2*psize+smallpsize psize smallpsize])"<<endl;
                           break;
                    }
                    startx= e[0] ? frame->frame1d_11()->first_wavelet(j[0],p[0]).k() : frame->frame1d_11()->first_generator(j[0],p[0]).k();
                    starty= e[1] ? frame->frame1d()->last_wavelet(j[1],p[1]).k() : frame->frame1d()->last_generator(j[1],p[1]).k(); //shift because first generator/wavelet lies on patch 3
                    columns= e[0] ? frame->frame1d_11()->Nablasize(j[0],p[0]) : frame->frame1d_11()->Deltasize(j[0],p[0]);
                    rows=1;//shift because first generator/wavelet lies on patch 3
                    //os<<"axis(axlabels);"<<endl;
                    os<<"axis(["<<(double)startx-0.5<<" "<<(double)startx-0.5+columns<<" "<<(double)starty-0.5<<" "<<(double)starty-0.5+rows<<"],'nolabel','ticx');"<<endl;
                    break;
                case 9:
                case 12:
                    
//                    os<<"subplot(3,3,4,'position',[omargin omargin+psize+imargin psize smallpsize])"<<endl;   
                    switch(patch){
                       case 9:
                           os<<"#patch 9:"<<endl;
                           os<<"axes('position',[omargin omargin+psize+imargin psize smallpsize])"<<endl;
                           break;
                       case 12:
                           os<<"#patch 12:"<<endl;
                           os<<"axes('position',[omargin+4*imargin+2*psize+2*smallpsize omargin+imargin+psize psize smallpsize])"<<endl;
                           break;
                    } 
                    startx= e[0] ? frame->frame1d_11()->first_wavelet(j[0],p[0]).k() : frame->frame1d_11()->first_generator(j[0],p[0]).k();
                    starty= e[1] ? frame->frame1d()->first_wavelet(j[1],p[1]).k() : frame->frame1d()->first_generator(j[1],p[1]).k();
                    columns= e[0] ? frame->frame1d_11()->Nablasize(j[0],p[0]) : frame->frame1d_11()->Deltasize(j[0],p[0]);
                    rows=1;
                    os<<"axis(["<<(double)startx-0.5<<" "<<(double)startx-0.5+columns<<" "<<(double)starty-0.5<<" "<<(double)starty-0.5+rows<<"],'nolabel','ticx');"<<endl;
                    break;
                case 10:
                case 15:                    
//                    os<<"subplot(3,3,8,'position',[omargin+psize+imargin omargin smallpsize psize])"<<endl;  
                    switch(patch){
                       case 10:
                           os<<"#patch 10:"<<endl;
                           os<<"axes('position',[omargin+imargin+psize omargin smallpsize psize])"<<endl;
                           break;
                       case 15:
                           os<<"#patch 15:"<<endl;
                           os<<"axes('position',[omargin+imargin+psize omargin+4*imargin+2*psize+2*smallpsize smallpsize psize])"<<endl;
                           break;
                    }  
                    startx= e[0] ? frame->frame1d()->first_wavelet(j[0],p[0]).k() : frame->frame1d()->first_generator(j[0],p[0]).k();
                    starty= e[1] ? frame->frame1d_11()->first_wavelet(j[1],p[1]).k() : frame->frame1d_11()->first_generator(j[1],p[1]).k();
                    columns=1;
                    rows= e[1] ? frame->frame1d_11()->Nablasize(j[1],p[1]) : frame->frame1d_11()->Deltasize(j[1],p[1]);
                    os<<"axis(["<<(double)startx-0.5<<" "<<(double)startx-0.5+columns<<" "<<(double)starty-0.5<<" "<<(double)starty-0.5+rows<<"],'nolabel','ticy');"<<endl;
                    break;
                case 11:
                case 14:    
//                    os<<"subplot(3,3,6,'position',[omargin+2*imargin+psize+smallpsize omargin+psize+imargin psize smallpsize])"<<endl;   
                    switch(patch){
                       case 11:
                           os<<"#patch 11:"<<endl;
                           os<<"axes('position',[omargin+3*imargin+2*psize+smallpsize omargin smallpsize psize])"<<endl;
                           break;
                       case 14:
                           os<<"#patch 14:"<<endl;
                           os<<"axes('position',[omargin+3*imargin+2*psize+smallpsize omargin+4*imargin+2*psize+2*smallpsize smallpsize psize])"<<endl;
                           break;
                    }
                    startx= e[0] ? frame->frame1d()->last_wavelet(j[0],p[0]).k() : frame->frame1d()->last_generator(j[0],p[0]).k();
                    starty= e[1] ? frame->frame1d_11()->first_wavelet(j[1],p[1]).k() : frame->frame1d_11()->first_generator(j[1],p[1]).k();
                    columns= 1;
                    rows=e[1] ? frame->frame1d_11()->Nablasize(j[1],p[1]) : frame->frame1d_11()->Deltasize(j[1],p[1]);
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
                    double val=std::max(log10(fabs(*it)/maxnorm),lowerclim);
//                    val=(ind.k()[0]==0 && ind.k()[1]==0) ? val+1e-3 : val;
//                    val==lowerclim ? val=val+1e-2*rand()/RAND_MAX : val=val;
                    //const double val=log10(fabs(*it)/maxnorm);
                    os << "patch(["<<ax<<","<<bx<<","<<bx<<","<<ax<<"],["<<ay<<","<<ay<<","<<by<<","<<by<<"],"<<val;
                    if(!boxed){
                        os<<",'edgecolor','none'";
                    }
                    os<<");"<<endl;
                }           
            }
            os << "set(gca,'CLim',[" << lowerclim << " 0])" << endl;
        }
        
        //correct colorbar labels
        if (colorbar) {
            //empty subplot for colorbar
//            os<<"subplot(3,3,5,'position',[omargin+2*imargin+psize+smallpsize omargin+2*imargin+psize+smallpsize psize psize])"<<endl;
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
