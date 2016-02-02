// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2015                                            |
// | Thorsten Raasch, Manuel Werner, Ulrich Friedrich                   |
// +--------------------------------------------------------------------+


#ifndef _MATHTL_FIXED_VECTOR_IO_H
#define _MATHTL_FIXED_VECTOR_IO_H

#include <algebra/fixed_vector.h>
#include <iostream>
#include <fstream>


using std::cout;
using std::endl;
namespace MathTL
{
    template<unsigned int SIZE>
    void
    fixed_vector_writeToFile(const char *filename, const FixedVector<double, SIZE> & v)
    {
        std::ofstream bin_file(filename);
        //bin_file.is_open();
        if (bin_file.is_open())
        {
            try
            {
                int size = SIZE;
                double temp_d;
                //for (typename InfiniteVector<C,I>::const_iterator it((*this).begin()); it != (*this).end(); ++it)
                bin_file.write((char*)(&size), sizeof(int));
                for (unsigned int i=0; i< SIZE; ++i)
                {
                    //cout << "writing: it = " << it.index() << " value = " << (*it) << endl;
                    temp_d = v[i];
                    //cout << "temp_d="<<temp_d<<endl;
                    bin_file.write(reinterpret_cast<char*>(&temp_d), sizeof(double));
                }
            }
            catch (...)
            {
                cout << "fixed_vector_writeToFile: Error while writing to " << filename << endl;
            }
            bin_file.close();
        }
        else
        {
            cout << "fixed_vector_writeToFile: Could not write to " << filename << endl;
        }
    }

    template <unsigned int SIZE>
    void fixed_vector_readFromFile(const char *filename, FixedVector<double, SIZE>& v)
    {
        std::ifstream bin_file(filename);
        if (bin_file.is_open())
        {
            //while (!bin_file.eof()) // doesn't work, see comment after .eof()
            int temp_i;
            double temp_d;
            try
            {
                bin_file.read((char*)(&temp_i), sizeof(int));
                assert(temp_i == SIZE);
                temp_i=0;
                while(true)
                {
                    bin_file.read(reinterpret_cast<char*>(&temp_d), sizeof(double));
                    if (bin_file.eof()) // the last meaningful .read operation sets the last bit to the last bit of the last variable read. that means that we are NOT at the end of file. So the whole try-block is executed again, resulting in one wrong entry in *this.
                    {
                        break;
                    }
                    //cout << "d="<<temp_d<<endl;
                    v[temp_i] = temp_d;
                    ++temp_i;
                }
                assert(temp_i == SIZE);
            }
            catch (...)
            {
                cout << "fixed_vector_readFromFile: Read error in " << filename << endl;
            }
            bin_file.close();
        }
        else
        {
            cout << "fixed_vector_readFromFile: Could not read from " << filename << endl;
        }
    }
    
}

#endif	/* _MATHTL_FIXED_VECTOR_IO_H */