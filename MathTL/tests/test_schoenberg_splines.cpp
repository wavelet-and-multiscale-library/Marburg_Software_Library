/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <iostream>
#include <geometry/point.h>
#include <geometry/grid.h>
#include <geometry/sampled_mapping.h>
#include <numerics/schoenberg_splines.h>
#include <utils/array1d.h>

using namespace std;
using namespace MathTL;

int main()
{
    clock_t tic, toc;
    double time;
    
    const int d=3;
    const int k=1;
    
    cout << "testing schoenberg splines" << endl;
    
    Array1D<double> points(400);
    Array1D<double> values(400);
    tic = clock();
    for(int h = 0; h < 400; h++){
        points[h]=(double) h/100;
        values[h]=MathTL::EvaluateSchoenbergBSpline<d>(k, (double) h/100);
    }
    toc = clock();
    time = (double)(toc-tic);
    cout << "Time taken: " << (time/CLOCKS_PER_SEC) << " s\n"<<endl;
    
    Grid<1> gitter(points);
    SampledMapping<1> sm(gitter, values);
    std::ofstream stream("schoenbergspline.m");
    sm.matlab_output(stream);
    stream << "plot(x,y)"<<endl;
    stream << "title('d="<<d<<", k="<<k<<"')"<<endl;
    stream.close();  
    cout << "matlab output plotted" << endl;
    
    cout << "measuring time for recursive call" << endl;
    tic = clock();
    for(int h = 0; h < 400; h++){
        points[h]=(double) h/100;
        values[h]=MathTL::EvaluateCardinalBSpline<d>(1, (double) h/100);
    }
    toc = clock();
    time = (double)(toc-tic);
    cout << "Time taken: " << (time/CLOCKS_PER_SEC) << " s\n"<<endl;
    
    cout << "measuring time for non-recursive call" << endl;
    tic = clock();
    for(int h = 0; h < 400; h++){
        points[h]=(double) h/100;
        values[h]=MathTL::EvaluateCardinalBSpline<d>((double) h/100);
    }
    toc = clock();
    time = (double)(toc-tic);
    cout << "Time taken: " << (time/CLOCKS_PER_SEC) << " s\n"<<endl;
    
    return 0;
}