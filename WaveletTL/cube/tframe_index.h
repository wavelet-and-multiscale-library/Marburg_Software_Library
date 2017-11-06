// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Ulrich Friedrich                   |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_TFRAME_INDEX_H_
#define _WAVELETTL_TFRAME_INDEX_H_

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
    template <class IFRAME, unsigned int DIM> class TensorFrame;
    /*
     * An index class for tensor product quarklet frames over the d-dimensional unit cube [0,1]^d (or mapped versions thereof) of the type
     * (V_0 \oplus W_0 \oplus W_1 \oplus W_2 \oplus ...)\times (V_0 \oplus W_0 \oplus W_1 \oplus W_2 \oplus ...)
     * as modeled by Tensorframe
    */
    template <class IFRAME, unsigned int DIM, class TENSORFRAME = TensorFrame<IFRAME,DIM> >
    class TensorQIndex
    {
    public:
        // polynomial index type
        typedef MultiIndex<int,DIM> polynomial_type;
        // level index type
        typedef MultiIndex<int,DIM> level_type;
        // type index type
        typedef MultiIndex<int,DIM> type_type;
        // translation index type
        typedef MultiIndex<int,DIM> translation_type;

        /*
         * Constructor with a given tensor frame
         * (also serves as a default constructor, but yields an invalid index
         * in this case, because the underlying bases must be specified to work correctly)
         */
        TensorQIndex(const TENSORFRAME* frame = 0);

        /*
         * Constructor with given p,j,e,k
         * Code is similar to Tensorindex(int,Frame).
         * 
         * PERFORMANCE: This constructor can be improved by using explicit formulas as in
         * TensorQIndex(const int number, const TENSORFRAME* frame)
        */
        TensorQIndex(const polynomial_type& p, const level_type& j, const type_type& e, const translation_type& k, const TENSORFRAME* frame);

        /*
         * 1D dummy constructor. does not yield a complete TensorQIndex
         * number is always set to 0
         */
        TensorQIndex(const int& p, const int& j, const int& e, const int& k, const TENSORFRAME* empty);

        // Copy constructor
        TensorQIndex(const TensorQIndex& lambda);

        // Copy index from const pointer
        TensorQIndex(const TensorQIndex* lambda);

        // Constructor with all parameters given
        TensorQIndex(const polynomial_type& p, const level_type& j, const type_type& e, const translation_type& k, const int number, const TENSORFRAME* frame);
        
        /*
         * Constructor for given number of quarklet index.
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
        TensorQIndex(const int number, const TENSORFRAME* frame);

        // Assignment
        TensorQIndex& operator = (const TensorQIndex& lambda);

        /*
         * Check equality.
         * Only (p,j,e,k) are tested. NOT the frame or number
         */
        bool operator == (const TensorQIndex& lambda) const;

        // Check non-equality
        inline bool operator != (const TensorQIndex& lambda) const
        { return !(*this == lambda); }

        // Preincrement
        TensorQIndex& operator ++ ();

        /* Ordering <
         * First by polynomial, i.e. the 1-norm of p, 
         * second by level, i.e. the 1-norm of j,
         * then lexicographically w.r.t. p,j,e,k
         */
        bool operator < (const TensorQIndex& lambda) const;

        // Ordering <=
        bool operator <= (const TensorQIndex& lambda) const
        { return (*this < lambda || *this == lambda); }

        // Scale j
        const polynomial_type& p() const { return p_; }
        
        // Scale j
        const level_type& j() const { return j_; }

        // Type e
        const type_type& e() const { return e_; }

        // Translation index k
        const translation_type& k() const { return k_; }

        // Underlying frame
        const TENSORFRAME* frame() const { return frame_; }

        const unsigned long int number() const { return num_; }

    protected:

        // Pointer to the underlying frame
        const TENSORFRAME* frame_;

        // Polynomial
        MultiIndex<int,DIM> p_;
        
        // Scale
        MultiIndex<int,DIM> j_;

        // Type
        MultiIndex<int,DIM> e_;

        // Translation
        MultiIndex<int,DIM> k_;

        // Number of the index, only for the elements of a quarklet frame
        unsigned int num_;

    };

    //! stream output
    template <class IFRAME, unsigned int DIM, class TENSORFRAME>
    inline std::ostream& operator << (std::ostream& os, const TensorQIndex<IFRAME,DIM,TENSORFRAME>& lambda)
    {
        using namespace std;
        os << "("
        << lambda.p()
        << ","
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
    template <class IFRAME, unsigned int DIM, class TENSORFRAME>
    TensorQIndex<IFRAME,DIM,TENSORFRAME>
    first_q_generator(const TENSORFRAME* frame);

    /*
     * index of last generator
     */
    template <class IFRAME, unsigned int DIM, class TENSORFRAME>
    TensorQIndex<IFRAME,DIM,TENSORFRAME>
    last_q_generator(const TENSORFRAME* frame);

    /*
     * Index of first quarklet on level j.
     * It is not checked whether j is a valid level, i.e., j-j0 is nonnegative
     */
    template <class IFRAME, unsigned int DIM, class TENSORFRAME>
    TensorQIndex<IFRAME,DIM,TENSORFRAME>
    first_quarklet(const TENSORFRAME* frame, const typename TensorQIndex<IFRAME,DIM,TENSORFRAME>::level_type j, 
                   const typename TensorQIndex<IFRAME,DIM,TENSORFRAME>::polynomial_type p = typename TensorQIndex<IFRAME,DIM,TensorFrame<IFRAME,DIM> >::polynomial_type());

    /*
     * Index of first quarklet on level j .
     * The test j >= ||j0||_1 is only made if _TFRAME_DEBUGLEVEL_ >= 1
     */
    template <class IFRAME, unsigned int DIM, class TENSORFRAME>
    TensorQIndex<IFRAME,DIM,TENSORFRAME>
    first_quarklet(const TENSORFRAME* frame, const int level, 
                   const typename TensorQIndex<IFRAME,DIM, TensorFrame<IFRAME,DIM> >::polynomial_type p = typename TensorQIndex<IFRAME,DIM,TensorFrame<IFRAME,DIM> >::polynomial_type());

    /*
     * Index of last quarklet on level j.
     * It is not checket whether j is a valid level, i.e., j-j0 is nonnegative
     */
    template <class IFRAME, unsigned int DIM, class TENSORFRAME>
    TensorQIndex<IFRAME,DIM,TENSORFRAME>
    last_quarklet(const TENSORFRAME* frame, const typename TensorQIndex<IFRAME,DIM,TENSORFRAME>::level_type j, 
                  const typename TensorQIndex<IFRAME,DIM,TensorFrame<IFRAME,DIM> >::polynomial_type p = typename TensorQIndex<IFRAME,DIM,TensorFrame<IFRAME,DIM> >::polynomial_type());

    /*
     * Index of last quarklet on level j >= ||j0||_1
     * The test j >= ||j0||_1 is only made if _TFRAME_DEBUGLEVEL_ >= 1
     */
    template <class IFRAME, unsigned int DIM, class TENSORFRAME>
    TensorQIndex<IFRAME,DIM,TENSORFRAME>
    last_quarklet(const TENSORFRAME* frame, const unsigned int level, const typename TensorQIndex<IFRAME,DIM,
                  TensorFrame<IFRAME,DIM> >::polynomial_type p = typename TensorQIndex<IFRAME,DIM,TensorFrame<IFRAME,DIM> >::polynomial_type());
    
    /*
     * number of first generator on level j0
     */
    template <class IFRAME, unsigned int DIM, class TENSORFRAME>
    const int
    first_q_generator_num(const TENSORFRAME* frame);

    /*
     * number of last generator on level j0
     */
    template <class IFRAME, unsigned int DIM, class TENSORFRAME>
    const int
    last_q_generator_num(const TENSORFRAME* frame);

    /*
     * Number of first quarklet on level j >= j0.
     * Note: This methods calls the constructor of TENSORFRAME. It it thus 
     * inefficient to use the so computed number to create an Tensor Index. 
     * This is the same for all ..._quarklet_num methods.
     */
    template <class IFRAME, unsigned int DIM, class TENSORFRAME>
    const int
    first_quarklet_num(const TENSORFRAME* frame, const typename TensorQIndex<IFRAME,DIM,TENSORFRAME>::level_type j,
                       const typename TensorQIndex<IFRAME,DIM,TENSORFRAME>::polynomial_type p);

    /*
     * number of last quarklet on level j >= j0
     */
    template <class IFRAME, unsigned int DIM, class TENSORFRAME>
    const int
    last_quarklet_num(const TENSORFRAME* frame, const typename TensorQIndex<IFRAME,DIM,TENSORFRAME>::level_type j,
                      const typename TensorQIndex<IFRAME,DIM,TENSORFRAME>::polynomial_type p);

    template <class IFRAME, unsigned int DIM, class TENSORFRAME>
    const int
    last_quarklet_num(const TENSORFRAME* frame, const unsigned int level);
    
    


    /*
     *  write InfiniteVector<double,Index> to stream
     * The only differences to the IO routines from qtframe.h are that
     *  - this was not tested
     *  - there is no patchnumer p
     */
    
    template <class IFRAME, unsigned int DIM, class TENSORFRAME>
    void writeIVToFile(const InfiniteVector<double, TensorQIndex<IFRAME,DIM,TENSORFRAME> >& v, std::ofstream& ofs)
    {
        ofs.is_open();
        if (ofs.is_open())
        {
            for (typename InfiniteVector<double,TensorQIndex<IFRAME,DIM,TENSORFRAME> >::const_iterator it(v.begin()), itend(v.end()); it !=itend; ++it)
            {
                try
                {
                    double temp_d(*it);
                    ofs.write(reinterpret_cast<char*>(&temp_d), sizeof(double));
                    for (int i=0; i<DIM; ++i)
                    {
                        int temp_p (it.index().p()[i]);
                        int temp_j (it.index().j()[i]);
                        int temp_e (it.index().e()[i]);
                        int temp_k (it.index().k()[i]);
                        ofs.write(reinterpret_cast<char*>(&temp_p), sizeof(int));
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
    template <class IFRAME, unsigned int DIM, class TENSORFRAME>
    void writeIVToFile(const InfiniteVector<double,TensorQIndex<IFRAME,DIM,TENSORFRAME> >& v, const char* filename)
    {
        std::ofstream ofs(filename,std::ofstream::binary);
        if (ofs.is_open())
        {
            for (typename InfiniteVector<double,TensorQIndex<IFRAME,DIM,TENSORFRAME> >::const_iterator it(v.begin()), itend(v.end()); it !=itend; ++it)
            {
                try
                {
                    double temp_d(*it);
                    ofs.write(reinterpret_cast<char*>(&temp_d), sizeof(double));
                    for (int i=0; i<DIM; ++i)
                    {
                        int temp_p (it.index().p()[i]);
                        int temp_j (it.index().j()[i]);
                        int temp_e (it.index().e()[i]);
                        int temp_k (it.index().k()[i]);
                        ofs.write(reinterpret_cast<char*>(&temp_p), sizeof(int));
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
    template <class IFRAME, unsigned int DIM, class TENSORFRAME>
    void readIVFromFile(const TENSORFRAME* frame, InfiniteVector<double, TensorQIndex<IFRAME,DIM,TENSORFRAME> >& v, std::ifstream& ifs)
    {
        v.clear();
        MultiIndex<int,DIM> temp_p, temp_j, temp_e, temp_k;
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
                        temp_p[i] = temp_i;                       
                        ifs.read(reinterpret_cast<char*>(&temp_i), sizeof(int));
                        temp_j[i] = temp_i;
                        ifs.read(reinterpret_cast<char*>(&temp_i), sizeof(int));
                        temp_e[i] = temp_i;
                        ifs.read(reinterpret_cast<char*>(&temp_i), sizeof(int));
                        temp_k[i] = temp_i;
                    }
                    ifs.read(reinterpret_cast<char*>(&temp_i), sizeof(int));
                    v.set_coefficient(TensorQIndex<IFRAME,DIM,TENSORFRAME>(temp_p, temp_j, temp_e, temp_k, temp_i, frame), temp_d);
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
    
    template <class IFRAME, unsigned int DIM, class TENSORFRAME>
    void readIVFromFile(const TENSORFRAME* frame, InfiniteVector<double,TensorQIndex<IFRAME,DIM,TENSORFRAME> >& v, const char* filename)
    {
        v.clear();
        std::ifstream ifs(filename, std::ifstream::binary);
        MultiIndex<int,DIM> temp_p, temp_j, temp_e, temp_k;
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
                        temp_p[i] = temp_i;
                        ifs.read(reinterpret_cast<char*>(&temp_i), sizeof(int));
                        temp_j[i] = temp_i;
                        ifs.read(reinterpret_cast<char*>(&temp_i), sizeof(int));
                        temp_e[i] = temp_i;
                        ifs.read(reinterpret_cast<char*>(&temp_i), sizeof(int));
                        temp_k[i] = temp_i;
                    }
                    ifs.read(reinterpret_cast<char*>(&temp_i), sizeof(int));
                    v.set_coefficient(TensorQIndex<IFRAME,DIM,TENSORFRAME>(temp_p, temp_j, temp_e, temp_k, temp_i, frame), temp_d);
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

#include <cube/tframe_index.cpp>

#endif /*_WAVELETTL_TFRAME_INDEX_H_*/

