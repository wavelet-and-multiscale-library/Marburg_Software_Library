// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2014                                            |
// | Ulrich Friedrich
// +--------------------------------------------------------------------+

#ifndef _MATHTL_INFINITE_VECTOR_IO_H
#define _MATHTL_INFINITE_VECTOR_IO_H

#include <algebra/infinite_vector.h>

namespace MathTL
{
    //! write InfiniteVector<double,int> to stream
    void writeIVToFile(const InfiniteVector<double,int >& v, std::ofstream& ofs)
    {
        if (ofs.is_open())
        {
            for ( InfiniteVector<double,int >::const_iterator it(v.begin()), itend(v.end()); it !=itend; ++it)
            {
                try
                {
                    double temp_d(*it);
                    ofs.write(reinterpret_cast<char*>(&temp_d), sizeof(double));
                    int temp_num (it.index());
                    ofs.write(reinterpret_cast<char*>(&temp_num), sizeof(int));
                }
                catch (...)
                {
                    cout << "writeIVToFile: Error while writing" << endl;
                }
            }
            ofs.close();
        }
        else
        {
            cout << "writeIVToFile: Could not write file" << endl;
        }
    }
    
    //! write InfiniteVector<double,int> to stream
    void writeIVToFile(const InfiniteVector<double,int >& v, const char* filename)
    {
        std::ofstream ofs(filename,std::ofstream::binary);
        if (ofs.is_open())
        {
            for ( InfiniteVector<double,int >::const_iterator it(v.begin()), itend(v.end()); it !=itend; ++it)
            {
                try
                {
                    double temp_d(*it);
                    ofs.write(reinterpret_cast<char*>(&temp_d), sizeof(double));
                    int temp_num (it.index());
                    ofs.write(reinterpret_cast<char*>(&temp_num), sizeof(int));
                }
                catch (...)
                {
                    cout << "writeIVToFile: Error while writing" << endl;
                }
            }
            ofs.close();
        }
        else
        {
            cout << "writeIVToFile: Could not write file" << endl;
        }
    }

    void readIVFromFile(InfiniteVector<double, int>& v, std::ifstream& ifs)
    {
        v.clear();
        if (ifs.is_open())
        {
            ifs.seekg (0, ifs.end);
            int length = ifs.tellg()/(sizeof(double)+sizeof(int));
            ifs.seekg (0, ifs.beg);
            double temp_d;
            int temp_num;
            for (unsigned int i=0; i<length;++i)
            {
                try
                {
                    ifs.read(reinterpret_cast<char*>(&temp_d), sizeof(double));
                    ifs.read(reinterpret_cast<char*>(&temp_num), sizeof(int));
                    v.set_coefficient(temp_num, temp_d);
                }
                catch (...)
                {
                    cout << "readIVFromFile: Read error" << endl;
                }
            }
            ifs.close();
        }
        else
        {
            cout << "readIVFromFile: Could not read at all" << endl;
        }
    }
    
    void readIVFromFile(InfiniteVector<double, int>& v, const char* filename)
    {
        std::ifstream ifs(filename, std::ifstream::binary);
        v.clear();
        if (ifs.is_open())
        {
            ifs.seekg (0, ifs.end);
            int length = ifs.tellg()/(sizeof(double)+sizeof(int));
            ifs.seekg (0, ifs.beg);
            double temp_d;
            int temp_num;
            for (unsigned int i=0; i<length;++i)
            {
                try
                {
                    ifs.read(reinterpret_cast<char*>(&temp_d), sizeof(double));
                    ifs.read(reinterpret_cast<char*>(&temp_num), sizeof(int));
                    v.set_coefficient(temp_num, temp_d);
                }
                catch (...)
                {
                    cout << "readIVFromFile: Read error" << endl;
                }
            }
            ifs.close();
        }
        else
        {
            cout << "readIVFromFile: Could not read at all" << endl;
        }
    }
   
}

#endif
