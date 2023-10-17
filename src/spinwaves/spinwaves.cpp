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
//   std::vector <double> Skx_FFT_array;
//   std::vector <double> Sky_FFT_array;
//   std::vector <double> Skz_FFT_array;

   //----------------------------------------------------------------------------
   // Function to initialize sw module
   //----------------------------------------------------------------------------
   void spin_wave( const std::vector<double>& atom_coords_x,
                   const std::vector<double>& atom_coords_y,
                   const std::vector<double>& atom_coords_z,
                   const int time ){

//    std::cout<< "Test Spin waves part2...."<<std::endl;
    int Na=spinwaves::internal::kx_FFT_array.size();
 

    std::fill(Skx_FFT_array_R.begin(), Skx_FFT_array_R.end(), 0.0);
    std::fill(Skx_FFT_array_I.begin(), Skx_FFT_array_I.end(), 0.0);
    for(unsigned int uca=0;uca<Na;uca++){
         std::ofstream file_K_time;
		 // JRH change filename
		 std::stringstream sstr;
		 //sstr << "K_vs_time_" << std::setw(4) << std::setfill('0') << std::to_string(uca) << ".dat";
         sstr << "K_vs_time" << 
			"_kx_" << std::setw(6) << std::to_string(spinwaves::internal::kx_FFT_array[uca]) << 
			"_ky_" << std::setw(6) << std::to_string(spinwaves::internal::ky_FFT_array[uca]) << 
			"_kz_" << std::setw(6) << std::to_string(spinwaves::internal::kz_FFT_array[uca]) << 		
			".dat";

	   file_K_time.open(sstr.str(),std::ios_base::app);
       const double kx=spinwaves::internal::kx_FFT_array[uca];
       const double ky=spinwaves::internal::ky_FFT_array[uca];
       const double kz=spinwaves::internal::kz_FFT_array[uca];
       double skx_R=0.0;
       double skx_I=0.0;
//       spinwaves::Skx_FFT_array_R[uca] = 0.0;
//       spinwaves::Skx_FFT_array_I[uca] = 0.0;
      double Skx_FFT_mat1=0;
      double Skx_FFT_mat2=0;

       for(int atom=0;atom<atoms::num_atoms;atom++){
           const double rx=atom_coords_x[atom];
           const double ry=atom_coords_y[atom];
           const double rz=atom_coords_z[atom];
           double arg=-kx*rx - ky*ry - kz*rz ;  
           double cosK= cos(arg);
           double sinK= sin(arg);
           double sx=atoms::x_spin_array[atom];
           int mat=atoms::type_array[atom];
           if (mat < 9) {
           skx_R += sx; 
           spinwaves::Skx_FFT_array_R[uca] += sx*cosK;
           spinwaves::Skx_FFT_array_I[uca] += sx*sinK;
           }
         //if (uca==0) std::cout<<uca<<"\t"<<atom<<"\t"<<rx<<"\t"<<ry<<"\t"<<rz<<"\t"<<cosK<<"\t"<<std::endl;
 //      if (uca==0) std::cout<<uca<<"\t"<<time<<"\t"<<sx<<"\t"<<skx_R<<"\t"<<Skx_FFT_array_R[uca]<<"\t"<<cosK<<"\t"<<std::endl;

       }
  //     if (uca==0) std::cout<<uca<<"\t"<<time<<"\t"<<"\t"<<skx_R<<"\t"<<Skx_FFT_array_R[uca]<<"\t"<<std::endl;
       file_K_time << uca<<" "<<  spinwaves::Skx_FFT_array_R[uca]<<" "<<  spinwaves::Skx_FFT_array_I[uca] << "\n";
       file_K_time.close();

    } 
      return;

   }

} // end of sw namespace

