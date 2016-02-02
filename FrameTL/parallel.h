// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of FrameTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _FRAME_TL_PARALLEL_H
#define _FRAME_TL_PARALLEL_H

#define NBLOCKS 2
#define MASTER 0

namespace FrameTL
{

  /*!
    \file parallel.h
    Routines for the parallelization of adaptive frame
    domain decomposition methods.
  */


  //MPI::Datatype coefficient_datatype;
  MPI_Datatype coefficient_datatype;

  /*!
    Create mpi datatype, consisting of an int and a double.
    Don't forget to call this routine at the very beginning of the mpi program.
   */
  void setup_coefficient_datatype ();

  /*!
    This routine sends an InfiniteVector<double, typename PROBLEM::Index>
    from the current processor to the master (the one with pid 0).
   */
  template <class PROBLEM>
  void send_to_Master (const InfiniteVector<double, typename PROBLEM::Index>&);

  /*!
    All processors, apart from the master, receive an
    InfiniteVector<double, typename PROBLEM::Index>. This routine is intended to work
    hand in hand with broadcast_vec_from_Master. In a parallel adaptive Schwarz
    frame algorithm, the master broadcasts the new global discrete iterate
    in an InfiniteVector<double, typename PROBLEM::Index> while the others
    are receiving this message with receive_all_parts.
   */
  template <class PROBLEM>
  void receive_all_parts (const PROBLEM&,
			  InfiniteVector<double, typename PROBLEM::Index>&);

  /*!
    The master processor sends an InfiniteVector<double, typename PROBLEM::Index>
    to all other (slave) processors.
   */
  template <class PROBLEM>
  void broadcast_vec_from_Master (const PROBLEM&,
				  InfiniteVector<double, typename PROBLEM::Index>&);

  /*!
    The master processor sends a double to all other (slave) processors.
   */
  void broadcast_double_from_Master (double& d);
  

}

#include <parallel.cpp>

#endif
