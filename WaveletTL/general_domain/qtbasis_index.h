// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Ulrich Friedrich                                                   |
// +--------------------------------------------------------------------+


#ifndef _WAVELETTL_QTIndex_H
#define	_WAVELETTL_QTIndex_H

#include <utils/multiindex.h>
#include <list>
using MathTL::MultiIndex;
using std::list;
using std::map;
using std::set;


namespace WaveletTL
{
    template <class IBASIS, unsigned int DIM> class QTBasis;
    /*! This class models indices for qtbasis. This basis consists of spatial extended
     * and translated anisotropic tensor wavelet bases on the unit cube. Consequently
     * an index consists of the wavelet type e, level j, translation parameter k and patchnumber p.
     * Ordering is first by the 1-norm of j and then lexicographically in j,p,e,k
     * In each direction generators only appear on the lowest level. They are put in
     * the same index set as the wavelets on the lowest level. Thats the reason for the lexicographical ordering.
     * Numbering begins at 0.
     * Caution: 
     * The first index is related to a GENERATOR on some patch. 
     * After all generators on this level and patch come all WAVELETS with this level on this patch. 
     * Observe that GENERATORS on different patches come after those wavelets on the starting patch.
     * As a result it makes no sense to iterate from the "first wavelet" to some later index. 
     * Consequently the output of the method "first_wavelet(minimal level)" equals "first_generator()"
     */
    template <class IBASIS, unsigned int DIM, class QTBASIS = QTBasis<IBASIS,DIM> >
    class QTIndex
    {
    public:
        // level index type
        typedef MultiIndex<int,DIM> level_type;
        // type index type
        typedef MultiIndex<int,DIM> type_type;
        // translation index type
        typedef MultiIndex<int,DIM> translation_type;

        /*
         * Constructor with a given qtbasis
         * (also serves as a default constructor, but yields an invalid index
         * in this case, because the underlying bases must be specified to work correctly)
         */
        QTIndex(const QTBASIS* basis = 0);

        /*
         * Constructor with given j,e,k,p
         */
        //QTIndex(const level_type& j, const type_type& e, const translation_type& k, const int p, const QTBASIS* basis);

        /*
         * Constructor with given j,e,k,p,num
         * May yield an invalid Index because correctness of num is not checked
         */
        QTIndex(const level_type& j, const type_type& e, const translation_type& k, const int p, const int num, const QTBASIS* basis);
        
        // Copy constructor
        QTIndex(const QTIndex& lambda);

        // Copy index from const pointer
        QTIndex(const QTIndex* lambda);

        /*!
         * Constructor. Returns the index (according to the ordering of the basis 
         * functions) that corresponds to number
        */
        // use full_collection instead!
        //QTIndex(const int number, const QTBASIS* basis);

        // Assignment
        QTIndex& operator = (const QTIndex& lambda);

        // Check equality
        bool operator == (const QTIndex& lambda) const;

        // Check non-equality
        inline bool operator != (const QTIndex& lambda) const
        { return !(*this == lambda); }

        /* Preincrement
         * code works only if minimal levels differ at most by one (see comments)
         */
        QTIndex& operator ++ ();

        /* Ordering <
         * First by level, i.e. the 1-norm of j,
         * then lexicographically w.r.t. j,p,e,k
         */
        bool operator < (const QTIndex& lambda) const;

        bool operator <= (const QTIndex& lambda) const
        { return (*this < lambda || *this == lambda); }

        // Scale j
        const level_type& j() const { return j_; }

        // Type e
        const type_type& e() const { return e_; }

        // Translation index k
        const translation_type& k() const { return k_; }

        // Patchnumber p
        const int& p() const { return p_; }

        // Underlying basis
        const QTBASIS* basis() const { return basis_; }

        // number of index
        const unsigned long int number() const { return num_; }

    protected:

        // Pointer to the underlying basis
        const QTBASIS* basis_;

        // Scale
        MultiIndex<int,DIM> j_;

        // Type
        MultiIndex<int,DIM> e_;

        // Translation
        MultiIndex<int,DIM> k_;

        // Patch
        int p_;

        // Number of the index, only for the elements of a wavelet bases
        unsigned long int num_;

    };
    
    template <class IBASIS, unsigned int DIM, class QTBASIS>
    inline std::ostream& operator << (std::ostream& os, const QTIndex<IBASIS,DIM,QTBASIS>& lambda)
    {
            using namespace std;
            os << "("
            << lambda.j()
            << ","
            << lambda.e()
            << ","
            << lambda.k()
            << ","
            << lambda.p()
            << ")";
    return os;
    }

    //! write InfiniteVector<double,Index> to stream
    template <class IBASIS, unsigned int DIM, class QTBASIS>
    void writeIVToFile(const InfiniteVector<double, QTIndex<IBASIS,DIM,QTBASIS> >& v, std::ofstream& ofs)
    {
        ofs.is_open();
        if (ofs.is_open())
        {
            for (typename InfiniteVector<double,QTIndex<IBASIS,DIM,QTBASIS> >::const_iterator it(v.begin()), itend(v.end()); it !=itend; ++it)
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
                    int temp_p (it.index().p());
                    int temp_num (it.index().number());
                    ofs.write(reinterpret_cast<char*>(&temp_p), sizeof(int));
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
    template <class IBASIS, unsigned int DIM, class QTBASIS>
    void writeIVToFile(const InfiniteVector<double,QTIndex<IBASIS,DIM,QTBASIS> >& v, const char* filename)
    {
        std::ofstream ofs(filename,std::ofstream::binary);
        if (ofs.is_open())
        {
            for (typename InfiniteVector<double,QTIndex<IBASIS,DIM,QTBASIS> >::const_iterator it(v.begin()), itend(v.end()); it !=itend; ++it)
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
                    int temp_p (it.index().p());
                    int temp_num (it.index().number());
                    ofs.write(reinterpret_cast<char*>(&temp_p), sizeof(int));
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
    template <class IBASIS, unsigned int DIM, class QTBASIS>
    void readIVFromFile(const QTBASIS* basis, InfiniteVector<double, QTIndex<IBASIS,DIM,QTBASIS> >& v, std::ifstream& ifs)
    {
        v.clear();
        MultiIndex<int,DIM> temp_j, temp_e, temp_k;
        if (ifs.is_open())
        {
            ifs.seekg (0, ifs.end);
            int length = ifs.tellg()/(sizeof(double)+(3*DIM+2)*sizeof(int));
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
                    temp_p = temp_i;
                    ifs.read(reinterpret_cast<char*>(&temp_i), sizeof(int));
                    v.set_coefficient(QTIndex<IBASIS,DIM,QTBASIS>(temp_j,temp_e, temp_k, temp_p, temp_i, basis), temp_d);
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
    
    template <class IBASIS, unsigned int DIM, class QTBASIS>
    void readIVFromFile(const QTBASIS* basis, InfiniteVector<double,QTIndex<IBASIS,DIM,QTBASIS> >& v, const char* filename)
    {
        v.clear();
        std::ifstream ifs(filename, std::ifstream::binary);
        MultiIndex<int,DIM> temp_j, temp_e, temp_k;
        if (ifs.is_open())
        {
            ifs.seekg (0, ifs.end);
            int length = ifs.tellg()/(sizeof(double)+(3*DIM+2)*sizeof(int));
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
                    temp_p = temp_i;
                    ifs.read(reinterpret_cast<char*>(&temp_i), sizeof(int));
                    v.set_coefficient(QTIndex<IBASIS,DIM,QTBASIS>(temp_j,temp_e, temp_k, temp_p, temp_i, basis), temp_d);
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

#include "qtbasis_index.cpp"

#endif	/* _WAVELETTL_QTIndex_H */

