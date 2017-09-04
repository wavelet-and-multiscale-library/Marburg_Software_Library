// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Ulrich Friedrich                   |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_TBASIS_INDEX_H_
#define _WAVELETTL_TBASIS_INDEX_H_

#include <iostream>
using std::cout;
using std::endl;

#include <utils/multiindex.h>
#include <utils/fixed_array1d.h>
#include <algebra/infinite_vector.h>
using MathTL::InfiniteVector;
using MathTL::MultiIndex;

namespace WaveletTL
{
    template <class IBASIS, unsigned int DIM> class TensorBasis;
    /*
     * An index class for tensor product wavelet bases over the d-dimensional unit cube [0,1]^d (or mapped versions thereof) of the type
     * (V_0 \oplus W_0 \oplus W_1 \oplus W_2 \oplus ...)\times (V_0 \oplus W_0 \oplus W_1 \oplus W_2 \oplus ...)
     * as modeled by TensorBasis
    */
    template <class IBASIS, unsigned int DIM, class TENSORBASIS = TensorBasis<IBASIS,DIM> >
    class TensorIndex
    {
    public:
        // level index type
        typedef MultiIndex<int,DIM> level_type;
        // type index type
        typedef MultiIndex<int,DIM> type_type;
        // translation index type
        typedef MultiIndex<int,DIM> translation_type;

        /*
         * Constructor with a given tensor basis
         * (also serves as a default constructor, but yields an invalid index
         * in this case, because the underlying bases must be specified to work correctly)
         */
        TensorIndex(const TENSORBASIS* basis = 0);

        /*
         * Constructor with given j,e,k
         * Code is similar to Tensorindex(int,Basis).
         * 
         * PERFORMANCE: This constructor can be improved by using explicit formulas as in
         * TensorIndex(const int number, const TENSORBASIS* basis)
        */
        TensorIndex(const level_type& j, const type_type& e, const translation_type& k, const TENSORBASIS* basis);

        /*
         * 1D dummy constructor. does not yield a complete TensorIndex
         * number is always set to 0
         */
        TensorIndex(const int& j, const int& e, const int& k, const TENSORBASIS* empty);

        // Copy constructor
        TensorIndex(const TensorIndex& lambda);

        // Copy index from const pointer
        TensorIndex(const TensorIndex* lambda);

        // Constructor with all parameters given
        TensorIndex(const level_type& j, const type_type& e, const translation_type& k, const int number, const TENSORBASIS* basis);
        
        /*
         * Constructor for given number of wavelet index.
         * The number is determined by the ordering "<", i.e.,
         * it is given first by \|level\|_1 
         * and on the same levelnorm lexicographical in the level and then in 
         * type and then in the translation index.
         * 
         * This constructor cannot be used for generators on levels j != j0 (multi index)
         * 
         * Example for level ordering for dim=3: (j,e) =
         * (333,000)->(333,001)->(333,010)->(333,100)->(333,101)->(333,110)->(333,111)-> (range =0)
         * (334,001)->(334,011)->(334,101)->(334,111)->(343,010)->...->(433,111)->  (range=1)
         * (335,001)-> ... (range=2)
         *
         * Code is speed up by using explicit formulas for DIM=2,3
         * DIM=1 is not optimized as it shouldn't be used anyways
         * 
         * full_collection[i] is an alternative! less CPU, more memory usage
        */
        TensorIndex(const int number, const TENSORBASIS* basis);

        // Assignment
        TensorIndex& operator = (const TensorIndex& lambda);

        /*
         * Check equality.
         * Only (j,e,k) are testet. NOT the basis or number
         */
        bool operator == (const TensorIndex& lambda) const;

        // Check non-equality
        inline bool operator != (const TensorIndex& lambda) const
        { return !(*this == lambda); }

        // Preincrement
        TensorIndex& operator ++ ();

        /* Ordering <
         * First by level, i.e. the 1-norm of j,
         * then lexicographically w.r.t. j,e,k
         */
        bool operator < (const TensorIndex& lambda) const;

        // Ordering <=
        bool operator <= (const TensorIndex& lambda) const
        { return (*this < lambda || *this == lambda); }

        // Scale j
        const level_type& j() const { return j_; }

        // Type e
        const type_type& e() const { return e_; }

        // Translation index k
        const translation_type& k() const { return k_; }

        // Underlying basis
        const TENSORBASIS* basis() const { return basis_; }

        const unsigned long int number() const { return num_; }

    protected:

        // Pointer to the underlying basis
        const TENSORBASIS* basis_;

        // Scale
        MultiIndex<int,DIM> j_;

        // Type
        MultiIndex<int,DIM> e_;

        // Translation
        MultiIndex<int,DIM> k_;

        // Number of the index, only for the elements of a wavelet bases
        unsigned int num_;

    };

    //! stream output
    template <class IBASIS, unsigned int DIM, class TENSORBASIS>
    inline std::ostream& operator << (std::ostream& os, const TensorIndex<IBASIS,DIM,TENSORBASIS>& lambda)
    {
        using namespace std;
        os << "("
        << lambda.j()
        << ","
        << lambda.e()
        << ","
        << lambda.k()
        << ")" << " number: " << lambda.number();
        return os;
    }

    /*
     * index of first generator
     */
    template <class IBASIS, unsigned int DIM, class TENSORBASIS>
    TensorIndex<IBASIS,DIM,TENSORBASIS>
    first_generator(const TENSORBASIS* basis);

    /*
     * index of last generator
     */
    template <class IBASIS, unsigned int DIM, class TENSORBASIS>
    TensorIndex<IBASIS,DIM,TENSORBASIS>
    last_generator(const TENSORBASIS* basis);

    /*
     * Index of first wavelet on level j.
     * It is not checket whether j is a valid level, i.e., j-j0 is nonnegative
     */
    template <class IBASIS, unsigned int DIM, class TENSORBASIS>
    TensorIndex<IBASIS,DIM,TENSORBASIS>
    first_wavelet(const TENSORBASIS* basis, const typename TensorIndex<IBASIS,DIM,TENSORBASIS>::level_type j);

    /*
     * Index of first wavelet on level j .
     * The test j >= ||j0||_1 is only made if _TBASIS_DEBUGLEVEL_ >= 1
     */
    template <class IBASIS, unsigned int DIM, class TENSORBASIS>
    TensorIndex<IBASIS,DIM,TENSORBASIS>
    first_wavelet(const TENSORBASIS* basis, const int level);

    /*
     * Index of last wavelet on level j.
     * It is not checket whether j is a valid level, i.e., j-j0 is nonnegative
     */
    template <class IBASIS, unsigned int DIM, class TENSORBASIS>
    TensorIndex<IBASIS,DIM,TENSORBASIS>
    last_wavelet(const TENSORBASIS* basis, const typename TensorIndex<IBASIS,DIM,TENSORBASIS>::level_type j);

    /*
     * Index of last wavelet on level j >= ||j0||_1
     * The test j >= ||j0||_1 is only made if _TBASIS_DEBUGLEVEL_ >= 1
     */
    template <class IBASIS, unsigned int DIM, class TENSORBASIS>
    TensorIndex<IBASIS,DIM,TENSORBASIS>
    last_wavelet(const TENSORBASIS* basis, const unsigned int level);
    
    /*
     * number of first generator on level j0
     */
    template <class IBASIS, unsigned int DIM, class TENSORBASIS>
    const int
    first_generator_num(const TENSORBASIS* basis);

    /*
     * number of last generator on level j0
     */
    template <class IBASIS, unsigned int DIM, class TENSORBASIS>
    const int
    last_generator_num(const TENSORBASIS* basis);

    /*
     * Number of first wavelet on level j >= j0.
     * Note: This methods calls the constructor of TENSORBASIS. It it thus 
     * inefficient to use the so computed number to create an Tensor Index. 
     * This is the same for all ..._wavelet_num methods.
     */
    template <class IBASIS, unsigned int DIM, class TENSORBASIS>
    const int
    first_wavelet_num(const TENSORBASIS* basis, const typename TensorIndex<IBASIS,DIM,TENSORBASIS>::level_type j);

    /*
     * number of last wavelet on level j >= j0
     */
    template <class IBASIS, unsigned int DIM, class TENSORBASIS>
    const int
    last_wavelet_num(const TENSORBASIS* basis, const typename TensorIndex<IBASIS,DIM,TENSORBASIS>::level_type j);

    template <class IBASIS, unsigned int DIM, class TENSORBASIS>
    const int
    last_wavelet_num(const TENSORBASIS* basis, const unsigned int level);
    
    


    /*
     *  write InfiniteVector<double,Index> to stream
     * The only differences to the IO routines from qtbasis.h are that
     *  - this was not tested
     *  - there is no patchnumer p
     */
    
    template <class IBASIS, unsigned int DIM, class TENSORBASIS>
    void writeIVToFile(const InfiniteVector<double, TensorIndex<IBASIS,DIM,TENSORBASIS> >& v, std::ofstream& ofs)
    {
        ofs.is_open();
        if (ofs.is_open())
        {
            for (typename InfiniteVector<double,TensorIndex<IBASIS,DIM,TENSORBASIS> >::const_iterator it(v.begin()), itend(v.end()); it !=itend; ++it)
            {
                try
                {
                    double temp_d(*it);
                    ofs.write(reinterpret_cast<char*>(&temp_d), sizeof(double));
                    for (int i=0; i<DIM; ++i)
                    {
                        int temp_j (it.index().j()[i]);
                        int temp_e (it.index().e()[i]);
                        int temp_k (it.index().k()[i]);
                        ofs.write(reinterpret_cast<char*>(&temp_j), sizeof(int));
                        ofs.write(reinterpret_cast<char*>(&temp_e), sizeof(int));
                        ofs.write(reinterpret_cast<char*>(&temp_k), sizeof(int));
                    }
                    int temp_num (it.index().number());
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
    
    //! write InfiniteVector<double,Index> to stream
    template <class IBASIS, unsigned int DIM, class TENSORBASIS>
    void writeIVToFile(const InfiniteVector<double,TensorIndex<IBASIS,DIM,TENSORBASIS> >& v, const char* filename)
    {
        std::ofstream ofs(filename,std::ofstream::binary);
        if (ofs.is_open())
        {
            for (typename InfiniteVector<double,TensorIndex<IBASIS,DIM,TENSORBASIS> >::const_iterator it(v.begin()), itend(v.end()); it !=itend; ++it)
            {
                try
                {
                    double temp_d(*it);
                    ofs.write(reinterpret_cast<char*>(&temp_d), sizeof(double));
                    for (int i=0; i<DIM; ++i)
                    {
                        int temp_j (it.index().j()[i]);
                        int temp_e (it.index().e()[i]);
                        int temp_k (it.index().k()[i]);
                        ofs.write(reinterpret_cast<char*>(&temp_j), sizeof(int));
                        ofs.write(reinterpret_cast<char*>(&temp_e), sizeof(int));
                        ofs.write(reinterpret_cast<char*>(&temp_k), sizeof(int));
                    }
                    int temp_num (it.index().number());
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
    
    //! open a file written by writeVectorToFile
    template <class IBASIS, unsigned int DIM, class TENSORBASIS>
    void readIVFromFile(const TENSORBASIS* basis, InfiniteVector<double, TensorIndex<IBASIS,DIM,TENSORBASIS> >& v, std::ifstream& ifs)
    {
        v.clear();
        MultiIndex<int,DIM> temp_j, temp_e, temp_k;
        if (ifs.is_open())
        {
            ifs.seekg (0, ifs.end);
            int length = ifs.tellg()/(sizeof(double)+(3*DIM+1)*sizeof(int));
            ifs.seekg (0, ifs.beg);
            double temp_d;
            int temp_i;
            for (unsigned int i=0; i<length;++i)
            {
                try
                {
                    ifs.read(reinterpret_cast<char*>(&temp_d), sizeof(double));
                    for (int i=0; i<DIM; ++i)
                    {
                        ifs.read(reinterpret_cast<char*>(&temp_i), sizeof(int));
                        temp_j[i] = temp_i;
                        ifs.read(reinterpret_cast<char*>(&temp_i), sizeof(int));
                        temp_e[i] = temp_i;
                        ifs.read(reinterpret_cast<char*>(&temp_i), sizeof(int));
                        temp_k[i] = temp_i;
                    }
                    ifs.read(reinterpret_cast<char*>(&temp_i), sizeof(int));
                    v.set_coefficient(TensorIndex<IBASIS,DIM,TENSORBASIS>(temp_j,temp_e, temp_k, temp_i, basis), temp_d);
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
    
    template <class IBASIS, unsigned int DIM, class TENSORBASIS>
    void readIVFromFile(const TENSORBASIS* basis, InfiniteVector<double,TensorIndex<IBASIS,DIM,TENSORBASIS> >& v, const char* filename)
    {
        v.clear();
        std::ifstream ifs(filename, std::ifstream::binary);
        MultiIndex<int,DIM> temp_j, temp_e, temp_k;
        if (ifs.is_open())
        {
            ifs.seekg (0, ifs.end);
            int length = ifs.tellg()/(sizeof(double)+(3*DIM+1)*sizeof(int));
            ifs.seekg (0, ifs.beg);
            double temp_d;
            int temp_i, temp_p;
            for (unsigned int i=0; i<length;++i)
            {
                try
                {
                    ifs.read(reinterpret_cast<char*>(&temp_d), sizeof(double));
                    for (int i=0; i<DIM; ++i)
                    {
                        ifs.read(reinterpret_cast<char*>(&temp_i), sizeof(int));
                        temp_j[i] = temp_i;
                        ifs.read(reinterpret_cast<char*>(&temp_i), sizeof(int));
                        temp_e[i] = temp_i;
                        ifs.read(reinterpret_cast<char*>(&temp_i), sizeof(int));
                        temp_k[i] = temp_i;
                    }
                    ifs.read(reinterpret_cast<char*>(&temp_i), sizeof(int));
                    v.set_coefficient(TensorIndex<IBASIS,DIM,TENSORBASIS>(temp_j,temp_e, temp_k, temp_i, basis), temp_d);
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

#include <cube/tbasis_index.cpp>

#endif /*_WAVELETTL_TBASIS_INDEX_H_*/

