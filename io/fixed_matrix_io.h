// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2015                                            |
// | Thorsten Raasch, Manuel Werner, Ulrich Friedrich                   |
// +--------------------------------------------------------------------+

#ifndef FIXED_MATRIX_IO_H
#define	FIXED_MATRIX_IO_H

#include <iostream>
#include <iomanip>

namespace MathTL
{
    
    template <class C, unsigned int ROW_DIM, unsigned int COL_DIM>
    void
    writeFixedMatrixToFile(const FixedMatrix< C,ROW_DIM,COL_DIM> & A, std::ofstream& ofs)
    {
        typedef typename FixedMatrix< C,ROW_DIM,COL_DIM>::size_type size_type;
        if (ofs.is_open())
        {
            try
            {
                size_type temp_i1(A.row_dimension());
                size_type temp_i2(A.column_dimension());
                ofs.write(reinterpret_cast<char*>(&temp_i1), sizeof(size_type));
                ofs.write(reinterpret_cast<char*>(&temp_i2), sizeof(size_type));
                for (int row=0; row<A.row_dimension(); ++row)
                {
                    for (int column=0; column<A.column_dimension(); ++column)
                    {
                        C temp_d(A(row,column));
                        ofs.write(reinterpret_cast<char*>(&temp_d), sizeof(C));
                    }
                }
            }
            catch (...)
            {
                cout << "writeFixedMatrixToFile: Error while writing" << endl;
            }
            ofs.close();
        }
        else
        {
            cout << "writeFixedMatrixToFile: Could not write file" << endl;
        }
    }
    
    template <class C, unsigned int ROW_DIM, unsigned int COL_DIM>
    void
    writeFixedMatrixToFile(const FixedMatrix< C,ROW_DIM,COL_DIM> & A, const char* filename)
    {
        std::ofstream ofs(filename,std::ofstream::binary);
        typedef typename FixedMatrix< C,ROW_DIM,COL_DIM>::size_type size_type;
        if (ofs.is_open())
        {
            try
            {
                size_type temp_i1(A.row_dimension());
                size_type temp_i2(A.column_dimension());
                ofs.write(reinterpret_cast<char*>(&temp_i1), sizeof(size_type));
                ofs.write(reinterpret_cast<char*>(&temp_i2), sizeof(size_type));
                for (int row=0; row<A.row_dimension(); ++row)
                {
                    for (int column=0; column<A.column_dimension(); ++column)
                    {
                        C temp_d(A(row,column));
                        ofs.write(reinterpret_cast<char*>(&temp_d), sizeof(C));
                    }
                }
            }
            catch (...)
            {
                cout << "writeFixedMatrixToFile: Error while writing" << endl;
            }
            ofs.close();
        }
        else
        {
            cout << "writeFixedMatrixToFile: Could not write file" << endl;
        }
    }
    
    template <class C, unsigned int ROW_DIM, unsigned int COL_DIM>
    void
    readFixedMatrixFromFile(FixedMatrix< C,ROW_DIM,COL_DIM> & A, std::ifstream& ifs)
    {
        typedef typename FixedMatrix< C,ROW_DIM,COL_DIM>::size_type size_type;
        if (ifs.is_open())
        {
            ifs.seekg (0, ifs.end);
            int length = ifs.tellg();
            assert (length > 2*sizeof(size_type));
            ifs.seekg (0, ifs.beg);
            size_type row_dim, col_dim;
            C entry;
            try
            {
                ifs.read(reinterpret_cast<char*>(&row_dim), sizeof(size_type));
                ifs.read(reinterpret_cast<char*>(&col_dim), sizeof(size_type));
                assert ( ((length - 2*sizeof(size_type))/sizeof(C)) == (row_dim * col_dim));
                for (size_type row(0); row < row_dim; ++row)
                {
                    for (size_type col(0); col < col_dim; ++col)
                    {
                        ifs.read(reinterpret_cast<char*>(&entry), sizeof(C));
                        A.set_entry(row,col,entry);
                    }
                }
            }
            catch (...)
            {
                cout << "readFixedMatrixFromFile: Read error" << endl;
            }
            ifs.close();
        }
        else
        {
            cout << "readFixedMatrixFromFile: Could not read at all" << endl;
        }
    }
    
    template <class C, unsigned int ROW_DIM, unsigned int COL_DIM>
    void
    readFixedMatrixFromFile(FixedMatrix< C,ROW_DIM,COL_DIM> & A, const char* filename)
    {
        std::ifstream ifs(filename, std::ifstream::binary);
        typedef typename FixedMatrix< C,ROW_DIM,COL_DIM>::size_type size_type;
        if (ifs.is_open())
        {
            ifs.seekg (0, ifs.end);
            int length = ifs.tellg();
            assert (length > 2*sizeof(size_type));
            ifs.seekg (0, ifs.beg);
            size_type row_dim, col_dim;
            C entry;
            try
            {
                ifs.read(reinterpret_cast<char*>(&row_dim), sizeof(size_type));
                ifs.read(reinterpret_cast<char*>(&col_dim), sizeof(size_type));
                assert ( ((length - 2*sizeof(size_type))/sizeof(C)) == (row_dim * col_dim));
                for (size_type row(0); row < row_dim; ++row)
                {
                    for (size_type col(0); col < col_dim; ++col)
                    {
                        ifs.read(reinterpret_cast<char*>(&entry), sizeof(C));
                        A.set_entry(row,col,entry);
                    }
                }
            }
            catch (...)
            {
                cout << "readFixedMatrixFromFile: Read error" << endl;
            }
            ifs.close();
        }
        else
        {
            cout << "readFixedMatrixFromFile: Could not read at all" << endl;
        }
    }

    
}




#endif	/* FIXED_MATRIX_IO_H */

