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
    fixed_matrix_writeToFile(const FixedMatrix< C,ROW_DIM,COL_DIM> & A, std::ofstream& ofs)
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
                cout << "fixed_matrix_writeToFile: Error while writing" << endl;
            }
            ofs.close();
        }
        else
        {
            cout << "fixed_matrix_writeToFile: Could not write file" << endl;
        }
    }
    
    template <class C, unsigned int ROW_DIM, unsigned int COL_DIM>
    void
    fixed_matrix_writeToFile(const FixedMatrix< C,ROW_DIM,COL_DIM> & A, const char* filename)
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
                cout << "fixed_matrix_writeToFile: Error while writing" << endl;
            }
            ofs.close();
        }
        else
        {
            cout << "fixed_matrix_writeToFile: Could not write file" << endl;
        }
    }
    
    template <class C, unsigned int ROW_DIM, unsigned int COL_DIM>
    void
    fixed_matrix_readFromFile(FixedMatrix< C,ROW_DIM,COL_DIM> & A, std::ifstream& ifs)
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
                cout << "fixed_matrix_readFromFile: Read error" << endl;
            }
            ifs.close();
        }
        else
        {
            cout << "fixed_matrix_readFromFile: Could not read at all" << endl;
        }
    }
    
    template <class C, unsigned int ROW_DIM, unsigned int COL_DIM>
    void
    fixed_matrix_readFromFile(FixedMatrix< C,ROW_DIM,COL_DIM> & A, const char* filename)
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
                cout << "fixed_matrix_readFromFile: Read error" << endl;
            }
            ifs.close();
        }
        else
        {
            cout << "fixed_matrix_readFromFile: Could not read at all" << endl;
        }
    }


/* specifications of the above routines with C = double
 */        
    /*
    template <unsigned int ROWSIZE, unsigned int COLUMNSIZE>
    void fixed_matrix_writeToFile(const char *filename, const FixedMatrix<double, ROWSIZE, COLUMNSIZE> & A)
    {
        //char filename[200];
        //filename[0] = '\x0';
        //strcat(filename, file);
        //strcat(filename, ".bin");
        std::ofstream bin_file(filename);
        if (bin_file.is_open())
        {
            try
            {
                int temp_i;
                temp_i = ROWSIZE;
                bin_file.write((char*)(&temp_i), sizeof(int));
                temp_i = COLUMNSIZE;
                bin_file.write((char*)(&temp_i), sizeof(int));
                double temp_d;
                for (unsigned int i=0; i < ROWSIZE; ++i)
                {
                    for (unsigned int j=0; j < COLUMNSIZE; ++j)
                    {
                        temp_d = A.get_entry(i,j);
                        bin_file.write(reinterpret_cast<char*>(&temp_d), sizeof(double));
                    }
                }
            }
            catch (...)
            {
                cout << "fixed_matrix_writeToFile: Error while writing to " << filename << endl;
            }
            bin_file.close();
        }
        else
        {
            cout << "fixed_matrix_writeToFile: Could not write to " << filename << endl;
        }
    }


    template <unsigned int ROWSIZE, unsigned int COLUMNSIZE>
    void fixed_matrix_readFromFile(const char *filename, FixedMatrix<double, ROWSIZE, COLUMNSIZE>& A)
    {
        //char filename[200];
        //filename[0] = '\x0';
        //strcat(filename, file);
        //strcat(filename, ".bin");
        std::ifstream bin_file(filename);

        if (bin_file.is_open())
        {
            //while (!bin_file.eof()) // doesn't work, see comment after .eof()
            int temp_i, temp_j;
            double temp_d;
            try
            {
                bin_file.read((char*)(&temp_i), sizeof(int));
                assert (temp_i == ROWSIZE);
                bin_file.read((char*)(&temp_j), sizeof(int));
                assert (temp_j == COLUMNSIZE);
                temp_i=0; temp_j=0;
                while(true)
                {
                    bin_file.read(reinterpret_cast<char*>(&temp_d), sizeof(double));
                    if (bin_file.eof()) // the last meaningful .read operation sets the last bit to the last bit of the last variable read. that means that we are NOT at the end of file. So the whole try-block is executed again, resulting in one wrong entry in *this.
                    {
                        break;
                    }
                    //cout << "d="<<temp_d<<endl;
                    A.set_entry(temp_i,temp_j, temp_d);
                    if (temp_j < COLUMNSIZE-1)
                    {
                        ++temp_j;
                    }
                    else
                    {
                        ++temp_i;
                        temp_j = 0;
                    }
                }
                assert ((temp_i == ROWSIZE) && (temp_j == 0));
            }
            catch (...)
            {
                cout << "fixed_matrix_readFromFile: Read error in " << filename << endl;
            }
            bin_file.close();
        }
        else
        {
            cout << "fixed_matrix_readFromFile: Could not read from " << filename << endl;
        }
    }
*/

}




#endif	/* FIXED_MATRIX_IO_H */

