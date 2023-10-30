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
   void spin_wave( const std::vector<double>& atom_coords_x,
                   const std::vector<double>& atom_coords_y,
                   const std::vector<double>& atom_coords_z,
                   const int time ){

      // std::cout<< "Test Spin waves part2...."<<std::endl;
      int Na=spinwaves::internal::kx_FFT_array.size();

      std::fill(Skx_FFT_array_R.begin(), Skx_FFT_array_R.end(), 0.0);
      std::fill(Skx_FFT_array_I.begin(), Skx_FFT_array_I.end(), 0.0);
      
      for(unsigned int uca=0;uca<Na;uca++){
         

         // JRH change filename
         std::ofstream file_K_time;
         std::stringstream sstr;
         sstr << "K_vs_time_" << std::setw(4) << std::setfill('0') << std::to_string(uca) << ".dat";
         file_K_time.open(sstr.str(),std::ios_base::app);

         kx=spinwaves::internal::kx_FFT_array[uca];
         ky=spinwaves::internal::ky_FFT_array[uca];
         kz=spinwaves::internal::kz_FFT_array[uca];
         // spinwaves::Skx_FFT_array_R[uca] = 0.0;
         // spinwaves::Skx_FFT_array_I[uca] = 0.0;

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
               Skx_FFT_array_R[uca] += sx*cosK;
               Skx_FFT_array_I[uca] += sx*sinK;
               
               // file_K_time << vmpi::my_rank << " " << atom << " " << " " << rx << " " << ry << " " << rz << " " << atoms::z_spin_array[atom] << "\n";
               // int mat=atoms::type_array[atom];

               // if (mat < 9) {
               //    skx_R += sx; 
               //    Skx_FFT_array_R[uca] += sx*cosK;
               //    Skx_FFT_array_I[uca] += sx*sinK;
               // }
            }

            // Add value of kpoints for all cpus JRH
            MPI_Allreduce(MPI_IN_PLACE, &Skx_FFT_array_R[0], Na, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &Skx_FFT_array_I[0], Na, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

         #else

            for(int atom=0;atom<vmpi::num_core_atoms+vmpi::num_bdry_atoms;atom++){
               rx=atom_coords_x[atom];
               ry=atom_coords_y[atom];
               rz=atom_coords_z[atom];
               arg=-kx*rx - ky*ry - kz*rz ;  
               cosK= cos(arg);
               sinK= sin(arg);
               sx=atoms::x_spin_array[atom];
               Skx_FFT_array_R[uca] += sx*cosK;
               Skx_FFT_array_I[uca] += sx*sinK;

               // file_K_time << vmpi::my_rank << " " << atom << " " << " " << rx << " " << ry << " " << rz << " " << atoms::z_spin_array[atom] << "\n";
               // int mat=atoms::type_array[atom];
               
               // if (mat < 9) {
               //    // skx_R += sx; 
               //    Skx_FFT_array_R[uca] += sx*cosK;
               //    Skx_FFT_array_I[uca] += sx*sinK;
               // }
            }

         #endif


         if(vmpi::my_rank == 0){
            file_K_time << time <<" "<<  spinwaves::Skx_FFT_array_R[uca]<<" "<<  spinwaves::Skx_FFT_array_I[uca] << "\n";
            file_K_time.close();
         }
      } 
    
	   return;

   }

} // end of sw namespace

