//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sergiu Ruta 2022. All rights reserved.
//
//   Email: sergiu.ruta@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "spinwaves.hpp"
#include "vmpi.hpp"
// sw module headers
#include "internal.hpp"
#include "iostream"
#include <iomanip>
#include <sstream>

//sergiu for SW
#include "unitcell.hpp"
#include "vector"
#include <cmath>
#include "fstream"
#include "atoms.hpp"


namespace spinwaves {

   // Local variables
   double kx;
   double ky;
   double kz;
   double rx;
   double ry;
   double rz;
   double arg;
   double cosK; 
   double sinK;
   double sx;
   int atom;

   //----------------------------------------------------------------------------
   // Function to initialize sw module
   //----------------------------------------------------------------------------
   void fft_in_space( const std::vector<double>& atom_coords_x,
                   const std::vector<double>& atom_coords_y,
                   const std::vector<double>& atom_coords_z,
                   const int time ){

      for(int k=0;k<internal::nk;k++){

         // kx=spinwaves::internal::kx_list[k];
         // ky=spinwaves::internal::ky_list[k];
         // kz=spinwaves::internal::kz_list[k];
         spinwaves::skx_r[k] = 0.0;
         spinwaves::skx_i[k] = 0.0;

         #ifdef MPICF

            for(int j=0;j<internal::mask.size();j++){
               
               atom=internal::mask[j];
               // rx=atom_coords_x[atom];
               // ry=atom_coords_y[atom];
               // rz=atom_coords_z[atom];
               // arg=-kx*rx - ky*ry - kz*rz ;  
               // cosK= cos(arg);
               // sinK= sin(arg);
               sx = (*internal::sw_array)[atom];

               // JRH cos_k and sin_k I dont think works for mpi
               skx_r[k] += sx*spinwaves::internal::cos_k[k*(internal::mask.size())+atom];
               skx_i[k] += sx*spinwaves::internal::sin_k[k*(internal::mask.size())+atom];
            }



            if (internal::reduc_ver == "direct_scatter"){
               // trying different approach where the transpose isnt needed
               MPI_Reduce(&skx_r[k], &skx_r_scatter[(k % nk_per_rank)*internal::nt+time], 1, MPI_DOUBLE, MPI_SUM, k / nk_per_rank, MPI_COMM_WORLD);
               MPI_Reduce(&skx_i[k], &skx_i_scatter[(k % nk_per_rank)*internal::nt+time], 1, MPI_DOUBLE, MPI_SUM, k / nk_per_rank, MPI_COMM_WORLD);
            }

         #else

            for(int j=0;j<internal::mask.size();j++){
               
               atom=internal::mask[j];
               // rx=atom_coords_x[atom];
               // ry=atom_coords_y[atom];
               // rz=atom_coords_z[atom];
               // arg=-kx*rx - ky*ry - kz*rz ;  
               // cosK= cos(arg);
               // sinK= sin(arg);
               sx = (*internal::sw_array)[atom];

               // testing whether predefing the cos(k) makes much of a difference to speed
               // skx_r[uca] += sx*cosK;
               // skx_i[uca] += sx*sinK;
               skx_r_node[time*internal::nk + k] += sx*spinwaves::internal::cos_k[k*internal::mask.size()+atom];
               skx_i_node[time*internal::nk + k] += sx*spinwaves::internal::sin_k[k*internal::mask.size()+atom];

            }

         #endif

      } 

      #ifdef MPICF
         if (internal::reduc_ver == "rank0"){
            std::cout << "TEST" << std::endl;
               MPI_Reduce(&skx_r[0], &skx_r_node[time*internal::nk], internal::nk, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
               MPI_Reduce(&skx_i[0], &skx_i_node[time*internal::nk], internal::nk, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
         }
      #endif  
          
	   return;

   }

} // end of sw namespace

