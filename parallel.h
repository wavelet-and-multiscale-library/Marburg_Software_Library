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


  //MPI::Datatype coefficient_datatype;
  MPI_Datatype coefficient_datatype;

  /*
    create mpi datatype, consisting of 
    an int and a double
   */
  void setup_coefficient_datatype ();


  /*
   */
  template <class PROBLEM>
  void send_to_Master (const InfiniteVector<double, typename PROBLEM::Index>&);

  /*
   */
  template <class PROBLEM>
  void receive_all_parts (const PROBLEM&,
			  InfiniteVector<double, typename PROBLEM::Index>&);

  /*
   */
  template <class PROBLEM>
  void broadcast_vec_from_Master (const PROBLEM&,
				  InfiniteVector<double, typename PROBLEM::Index>&);

  /*
   */
  void broadcast_double_from_Master (double& d);
  

}

#include <parallel.cpp>

#endif
