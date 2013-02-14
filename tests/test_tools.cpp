#include <cmath>
#include <iostream>
#include <utils/tiny_tools.h>

using namespace std;

int main()
{
  cout << "Testing some of the tiny tools..." << endl;

  for (int j = 0; j <= 10; j++) {
    double x1 = sqrt((double) (1<<j));
    double x2 = twotothejhalf(j);
    cout << "* 2^{" << j << "/2}= " << x2 << ", error: " << fabs(x2-x1) << endl;
  }

  cout << endl
       << "- checking positive dyadic modulo routine:" << endl;
  
  const int expo = 3;
  for (int j = -10; j <= 10; j++) {
    cout << j << " modulo 2^" << expo << "="
	 << dyadic_modulo(j,expo)
	 << endl;
  }

#if 1
    
    cout << "improving the log2 Function" << endl;
    
    unsigned int v;  // 32-bit value to find the log2 of 
    const unsigned int b[] = {0x2, 0xC, 0xF0, 0xFF00, 0xFFFF0000}; //, 0xFFFFFFFF00000000};
    const unsigned int S[] = {1, 2, 4, 8, 16}; //, 32};
    int i;
    
#if 1
    // first check what the new method does
    for (i = 4; i >= 0; i--)
    {
        cout << "b[" << i <<"] = " << b[i] << endl;
    }
        
    for (unsigned int num = 1; num < -1; ++num)
    {
        v = num;
        register unsigned int r = 0; // result of log2(v) will go here
        for (i = 4; i >= 0; i--) // unroll for speed...
        {
          if (v & b[i])
          {
            v >>= S[i];
            r |= S[i];
          } 
        }
        assert (r == log2(num));
        //assert (r == floor(log2((double)num))); // classical log2 for double
        //cout << "num = " << num << "; r = " << r << "; v = " << v << "; log2(num) = " << log2(num) << "; floor(log2(num)) = " << floor(log2(num)) << endl;
    }
    
    //for (unsigned int num = 1; num < -1; ++num)
    for (int num = 1; num > 0; ++num)
    {
        v = num;
        register unsigned int r = 0; // result of log2(v) will go here
        for (i = 4; i >= 0; i--) // unroll for speed...
        {
          if (v & b[i])
          {
            v >>= S[i];
            r |= S[i];
          } 
        }
        assert (r == log2((unsigned int)num));
        //assert (r == floor(log2((double)num))); // classical log2 for double
        //cout << "num = " << num << "; r = " << r << "; v = " << v << "; log2(num) = " << log2(num) << "; floor(log2(num)) = " << floor(log2(num)) << endl;
    }
    cout << "sizeof(int) = " << sizeof(int) << endl;
    cout << "sizeof(unsigned int) = " << sizeof(unsigned int) << endl;
    cout << "done" << endl;
    abort();
#endif
    
    cout << "Compare speed of floor(log2(double)) and log2(int)" << endl;
    int repetitions(1); // about 6 minutes per repetition
    clock_t tstart, tend;
    double time1(0), time2(0);
    //
    cout << "performing first speed test" << endl;
    tstart = clock();
    for (unsigned int rep = 0; rep < repetitions; ++rep)
    {
        for (unsigned int num = 1; num< -1; ++num)
        {
            v = num;
            register unsigned int r = 0; // result of log2(v) will go here
            for (i = 4; i >= 0; i--) // unroll for speed...
            {
                if (v & b[i])
                {
                    v >>= S[i];
                    r |= S[i];
                } 
            }
            //v = log2(num);
            //cout << "num = " << num << "; r = " << r << "; v = " << v << "; log2(num) = " << log2(num) << "; floor(log2(num)) = " << floor(log2(num)) << endl;
        }
    }
    tend = clock();
    time1 = (double)(tend-tstart)/CLOCKS_PER_SEC;
    //
    cout << "performing second speed test" << endl;
    tstart = clock();
    for (unsigned int rep = 0; rep < repetitions; ++rep)
    {
        for (unsigned int num = 1; num< -1; ++num)
        {
            //v = floor(log2(num));
            v = log2(num);
        }
    }
    tend = clock();
    time2 = (double)(tend-tstart)/CLOCKS_PER_SEC;
    cout << "Repetitions = " << repetitions << endl
            << "test1 " << time1 << "sec; test2 " << time2 << "sec" << endl;
#endif  
    return 0;
}
