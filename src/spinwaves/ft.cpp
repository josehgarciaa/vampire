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
#include <fftw3-mpi.h>


namespace spinwaves {

   // Local variables
   double sx;
   int atom;


   fftw_plan plan;
   fftw_complex *data;
   ptrdiff_t local_no, local_o_start, local_ni, local_i_start;

   //----------------------------------------------------------------------------
   // Function to initialize sw module
   //----------------------------------------------------------------------------
   void fft_in_space( const std::vector<double>& rx,
                   const std::vector<double>& ry,
                   const std::vector<double>& rz,
                   const int time ){

      for(int k=0;k<internal::nk;k++){

         #ifdef MPICF

            for(int j=0;j<internal::mask.size();j++){
               
               atom=internal::mask[j];
               sx = (*internal::sw_array)[atom];

               if (internal::prefactor == true){
                  skx_r[k] += sx*spinwaves::internal::cos_k[k*(internal::mask.size())+atom];
                  skx_i[k] += sx*spinwaves::internal::sin_k[k*(internal::mask.size())+atom];
               }
               else if (internal::prefactor == false){
                  skx_r[k] += sx*cos(rx[atom]*internal::kx[k] + ry[atom]*internal::kx[k] + rz[atom]*internal::kx[k]);
                  skx_i[k] += sx*sin(rx[atom]*internal::kx[k] + ry[atom]*internal::ky[k] + rz[atom]*internal::ky[k]);
               }
            }

            for(int j=0;j<internal::mask.size();j++){
               
               atom=internal::mask2[j];
               sx = (*internal::sw_array)[atom];

               if (internal::prefactor == true){
                  skx_r[k] *= sx*spinwaves::internal::cos_k2[k*(internal::mask2.size())+atom];
                  skx_i[k] *= sx*spinwaves::internal::sin_k2[k*(internal::mask2.size())+atom];
               }
               else if (internal::prefactor == false){
                  skx_r[k] *= sx*cos(rx[atom]*internal::kx[k] + ry[atom]*internal::kx[k] + rz[atom]*internal::kx[k]);
                  skx_i[k] *= sx*sin(rx[atom]*internal::kx[k] + ry[atom]*internal::ky[k] + rz[atom]*internal::ky[k]);
               }
            }




            if (internal::reduc_ver == "direct_scatter"){
               // trying different approach where the transpose isnt needed
               MPI_Reduce(&skx_r[k], &skx_r_scatter[(k % nk_per_rank)*internal::nt+time], 1, MPI_DOUBLE, MPI_SUM, k / nk_per_rank, MPI_COMM_WORLD);
               MPI_Reduce(&skx_i[k], &skx_i_scatter[(k % nk_per_rank)*internal::nt+time], 1, MPI_DOUBLE, MPI_SUM, k / nk_per_rank, MPI_COMM_WORLD);
            }

         #else

            for(int j=0;j<internal::mask.size();j++){
               
               atom=internal::mask[j];
               sx = (*internal::sw_array)[atom];

               if (internal::prefactor == true){
                  skx_r_node[time*internal::nk + k] += sx*spinwaves::internal::cos_k[k*(internal::mask.size())+atom];
                  skx_i_node[time*internal::nk + k] += sx*spinwaves::internal::sin_k[k*(internal::mask.size())+atom];
               }
               else if (internal::prefactor == false){
                  skx_r_node[time*internal::nk + k] += sx*cos(rx[atom]*internal::kx[k] + ry[atom]*internal::kx[k] + rz[atom]*internal::kx[k]);
                  skx_i_node[time*internal::nk + k] += sx*sin(rx[atom]*internal::kx[k] + ry[atom]*internal::ky[k] + rz[atom]*internal::ky[k]);
               }
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

