// implementation for parallel.h

#include <parallel.h>


namespace FrameTL
{


  /*
    create mpi datatype, consisting of 
    an int and a double
   */
  void setup_coefficient_datatype () {
    int array_of_block_lenghts[NBLOCKS] = {1, 1};
    
    MPI::Aint array_of_displacements[NBLOCKS];
    //MPI::Datatype array_of_types[NBLOCKS] = {MPI::INT, MPI::DOUBLE};
    MPI_Datatype array_of_types[NBLOCKS] = {MPI::INT, MPI::DOUBLE};
    
    cout << "Groesse eines int = " << sizeof(int) << endl;
    
//     Coefficient dummy[10];
//     array_of_displacements[0] = 0;
//     MPI_Address (&dummy[0].val, &array_of_displacements[1]);
//     array_of_displacements[1]-=array_of_displacements[0];
//     cout << "disp = " << array_of_displacements[0] << " " << array_of_displacements[1] << endl;
    
    // determine the lower bound and extent of an int
    MPI_Aint extent = 0;
    MPI::Aint lb = 0;
    MPI::INT.Get_extent(lb, extent);
    array_of_displacements[0] = lb;
    //MPI::DOUBLE.Get_extent(lb, extent);
    //array_of_displacements[1] = extent;
    array_of_displacements[1] = 8;
    
    cout << "extent = " << extent << endl;
    cout << "NBLOCKS = "  << NBLOCKS <<  endl;


//     coefficient_datatype = MPI::Datatype::Create_struct(2,
// 							array_of_block_lenghts,
// 							array_of_displacements,
// 							array_of_types
// 							);
//     coefficient_datatype.Commit();

    MPI_Type_struct (2, array_of_block_lenghts, array_of_displacements, array_of_types, &coefficient_datatype);
    MPI_Type_commit(&coefficient_datatype);

  }

  template <class PROBLEM>
  inline
  void send_to_Master (const InfiniteVector<double, typename PROBLEM::Index>& v) {
    int size = v.size();
    MPI::COMM_WORLD.Send( &size, 1, MPI::INT, MASTER, 0);
    
    Coefficient out_vec[size];
    to_array(v, out_vec);
    MPI::COMM_WORLD.Send( &out_vec, size, coefficient_datatype, MASTER, 0);
  }

  template <class PROBLEM>
  inline
  void receive_all_parts (const PROBLEM& P,
			  InfiniteVector<double, typename PROBLEM::Index>& v) {
    const int number_patches = P.basis().n_p();
    int buffer_sizes[number_patches-1];
    for (int i = 1; i < number_patches; i++)
      MPI::COMM_WORLD.Recv( &buffer_sizes[i-1], 1, MPI::INT, i, 0);
    
    for (int i = 1; i < number_patches; i++) {
      Coefficient in_vec[buffer_sizes[i-1]];
      MPI::COMM_WORLD.Recv( &in_vec, buffer_sizes[i-1], coefficient_datatype, i, 0);
      InfiniteVector<double, typename PROBLEM::Index> help;	  
      array_to_map (in_vec, &P.basis(), help, buffer_sizes[i-1]);
      v += help;
      //v.merge(help);
    }
  }

  /*
   */
  template <class PROBLEM>
  inline
  void broadcast_vec_from_Master (const PROBLEM& P,
			      InfiniteVector<double, typename PROBLEM::Index>& v) {
    int p = MPI::COMM_WORLD.Get_rank();
    int current_support_size  = 0;
    if (p == MASTER) {
      current_support_size = v.size();
      //cout << v << endl;
    }
    
    // tell the slaves about the support size of the current approximation
    MPI::COMM_WORLD.Bcast( &current_support_size, 1, MPI::INT, MASTER );

    // send the current approximation to all the slaves
    Coefficient v_vec[current_support_size];
    if (p == MASTER)
      to_array(v, v_vec);
    MPI::COMM_WORLD.Bcast( &v_vec, current_support_size, coefficient_datatype, MASTER);
    if (p != MASTER) {
      array_to_map (v_vec, &P.basis(), v, current_support_size);
    }
  }

  /*
   */
  void broadcast_double_from_Master (double& d)
  {
    //MPI::DOUBLE dd(d);
    MPI::COMM_WORLD.Bcast( &d, 1, MPI::DOUBLE, MASTER );
  }


}
