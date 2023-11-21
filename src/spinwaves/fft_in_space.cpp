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

   //----------------------------------------------------------------------------
   // Function to initialize sw module
   //----------------------------------------------------------------------------
   void fft_in_space( const std::vector<double>& atom_coords_x,
                   const std::vector<double>& atom_coords_y,
                   const std::vector<double>& atom_coords_z,
                   const int time ){

      int nk=spinwaves::internal::kx_FFT_array.size();

      for(int k=0;k<nk;k++){

         kx=spinwaves::internal::kx_FFT_array[k];
         ky=spinwaves::internal::ky_FFT_array[k];
         kz=spinwaves::internal::kz_FFT_array[k];
         spinwaves::Skx_FFT_array_R[k] = 0.0;
         spinwaves::Skx_FFT_array_I[k] = 0.0;

         // double skx_R=0.0;
         // double skx_I=0.0;

         // double Skx_FFT_mat1=0;
         // double Skx_FFT_mat2=0;

         // std::cout << "==============================\n";


         #ifdef MPICF

            for(int atom=0;atom<vmpi::num_core_atoms+vmpi::num_bdry_atoms;atom++){
               
               rx=atom_coords_x[atom];
               ry=atom_coords_y[atom];
               rz=atom_coords_z[atom];
               arg=-kx*rx - ky*ry - kz*rz ;  
               cosK= cos(arg);
               sinK= sin(arg);
               sx=atoms::x_spin_array[atom];


               // JRH cos_k and sin_k I dont think works for mpi
               Skx_FFT_array_R[k] += sx*spinwaves::internal::cos_k[k*(vmpi::num_core_atoms+vmpi::num_bdry_atoms)+atom];
               Skx_FFT_array_I[k] += sx*spinwaves::internal::sin_k[k*(vmpi::num_core_atoms+vmpi::num_bdry_atoms)+atom];
            }

            // // Add value of kpoints for all cpus JRH
            // MPI_Reduce(MPI_IN_PLACE, &Skx_FFT_array_R[0], nk, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            // MPI_Reduce(MPI_IN_PLACE, &Skx_FFT_array_I[0], nk, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            // MPI_Reduce(&Skx_FFT_array_R[0], &Skx_FFT_array_R_node[time*nk], nk, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            // MPI_Reduce(&Skx_FFT_array_I[0], &Skx_FFT_array_I_node[time*nk], nk, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

         #else

            for(int atom=0;atom<atoms::num_atoms;atom++){
               // rx=atom_coords_x[atom];
               // ry=atom_coords_y[atom];
               // rz=atom_coords_z[atom];
               // arg=-kx*rx - ky*ry - kz*rz ;  
               // cosK= cos(arg);
               // sinK= sin(arg);
               sx=atoms::x_spin_array[atom];

               // testing whether predefing the cos(k) makes much of a difference to speed
               // Skx_FFT_array_R[uca] += sx*cosK;
               // Skx_FFT_array_I[uca] += sx*sinK;
               Skx_FFT_array_R_node[time*nk + k] += sx*spinwaves::internal::cos_k[k*atoms::num_atoms+atom];
               Skx_FFT_array_I_node[time*nk + k] += sx*spinwaves::internal::sin_k[k*atoms::num_atoms+atom];

            }

            // file_K_time << time <<" "<<  spinwaves::Skx_FFT_array_R[uca]<<" "<<  spinwaves::Skx_FFT_array_I[uca] << "\n";
            // file_K_time.close();

         #endif

      } 


      #ifdef MPICF
         // Add value of kpoints for all cpus JRH
         if (vmpi::my_rank == 0){
            std::cout << "Reducing processor array to node array" << std::endl;
         }
         MPI_Reduce(&Skx_FFT_array_R[0], &Skx_FFT_array_R_node[time*nk], nk, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
         MPI_Reduce(&Skx_FFT_array_I[0], &Skx_FFT_array_I_node[time*nk], nk, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
         // MPI_Gather(&Skx_FFT_array_R[0], nk, MPI_DOUBLE, &Skx_FFT_array_R_node[time*nk], 100, MPI_DOUBLE, 0, MPI_COMM_WORLD);
         // MPI_Gather(&Skx_FFT_array_I[0], nk, MPI_DOUBLE, &Skx_FFT_array_I_node[time*nk], 100, MPI_DOUBLE, 0, MPI_COMM_WORLD);
         if (vmpi::my_rank == 0){
            std::cout << "Reduction completed." << std::endl;
         } 
      #endif
      
      
      
      // if(vmpi::my_rank == 0){
      //    for (int i = 0; i < nk*internal::numtimepoints; i++){
      //       file_K_time << i <<" "<<  spinwaves::Skx_FFT_array_R_node[i]<<" "<<  spinwaves::Skx_FFT_array_I_node[i] << "\n";
      //    }
      //    file_K_time.close();
      // }
    
	   return;

   }

} // end of sw namespace

