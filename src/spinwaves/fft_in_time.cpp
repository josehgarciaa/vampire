//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Joel Hirst 2023. All rights reserved.
//
//   Email: j.r.hirst@shu.ac.uk
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
#include "vmpi.hpp"
#include <cmath>
#include "fstream"
#include "atoms.hpp"
#include <fftw3.h>


namespace spinwaves {



   void fft_in_time(){

      std::cout << "CALCULATING FFT IN TIME " << std::endl;
      
      const int real = 0;
      const int imag = 1;

      if (vmpi::my_rank == 0){
         std::cout << "numtimepoints " << internal::numtimepoints << std::endl;
         std::cout << "numkpoints " << internal::numkpoints << std::endl;
      }

      if (vmpi::my_rank == 0){
         std::cout << "creating fourier stuff." << std::endl;
      }


      std::vector <double> Skx_FFT_array_R_transposed;
      std::vector <double> Skx_FFT_array_I_transposed;
      if (vmpi::my_rank == 0){
         Skx_FFT_array_R_transposed.resize(internal::numtimepoints*internal::numkpoints);
         Skx_FFT_array_I_transposed.resize(internal::numtimepoints*internal::numkpoints);	
      }

      std::vector<fftw_complex> combined_real_imag(internal::numtimepoints);
      std::vector<fftw_complex> combined_real_imag_fftd(internal::numtimepoints);
      fftw_plan fft_in_time;
      fft_in_time = fftw_plan_dft_1d(internal::numtimepoints, &combined_real_imag[0], &combined_real_imag_fftd[0], FFTW_FORWARD, FFTW_MEASURE);
      
      
      if (vmpi::my_rank == 0){
         std::cout << "fourier stuff created." << std::endl;
      }

      // the loop below transposes from:
      //
      // |    t0         t1      t2     |   ..... 
      // |[k0 -> kn][k0 -> kn][k0 -> kn]|
      //
      // to
      //
      // |    k0         k1      k2     |   ..... 
      // |[t0 -> tn][t0 -> tn][t0 -> tn]|
      //
      //
      // This then has to be scattered to the the processors. Example:
      //
      // |          rank 0              ||          rank 1              |
      // |    k0         k1      k2     ||    k0         k1      k2     |   ..... 
      // |[t0 -> tn][t0 -> tn][t0 -> tn]||[t0 -> tn][t0 -> tn][t0 -> tn]|
      //
      // Then loop over the number of k's on each rank
      //
      //
      // Not always the case that Nkpoints/nranks will be an integer. Need to round up and zero pad the array.
      // For example. 100 timesteps, 50 kpoints ran on 4 ranks will mean there is (100*50) / 4 data points per rank.
      // WIll have to round up and zero-pad.


      // rearrange array for scatter - I think this bit needs to be done for serial and parallel
      if (vmpi::my_rank == 0){
         for (int k = 0; k < internal::numkpoints; k++){
            for (int time=0; time < internal::numtimepoints; time++){   
               // Fill FFTW arrays with values from spacial FFT.
               Skx_FFT_array_R_transposed[k*internal::numtimepoints + time] = Skx_FFT_array_R_node[time*internal::numkpoints + k];
               Skx_FFT_array_I_transposed[k*internal::numtimepoints + time] = Skx_FFT_array_I_node[time*internal::numkpoints + k];
            }     
         }


      }

      #ifdef MPICF



         int nk_per_rank = std::ceil(static_cast<double>(internal::numkpoints) / static_cast<double>(vmpi::num_processors));
         int scatterlength = nk_per_rank * internal::numtimepoints;
         std::vector<double> Skx_FFT_array_R_scatter;
         std::vector<double> Skx_FFT_array_I_scatter;
         Skx_FFT_array_R_scatter.resize(scatterlength,0.0);
         Skx_FFT_array_I_scatter.resize(scatterlength,0.0);

         if (vmpi::my_rank == 0){

            // Not always the case that Nkpoints/nranks will be an integer. Need to round up and zero pad the array.
            // For example. 100 timesteps, 50 kpoints ran on 4 ranks will mean there is (100*50) / 4 data points per rank.
            // WIll have to round up and zero-pad.

            std::cout << "nk_per_rank " << nk_per_rank << std::endl;
            Skx_FFT_array_R_transposed.resize(nk_per_rank * vmpi::num_processors * internal::numtimepoints,0.0);
            Skx_FFT_array_I_transposed.resize(nk_per_rank * vmpi::num_processors * internal::numtimepoints,0.0);
            std::cout << "Skx_FFT_array_R_transposed.size " << Skx_FFT_array_R_transposed.size() << std::endl;
            for (int i = 0; i < Skx_FFT_array_R_transposed.size(); i++){
               std::cout << Skx_FFT_array_R_transposed[i] << std::endl;

            }


            // calculate number of elements for the scatter on each rank. Create variable on each rank of that size
            std::cout << "Number of kpoints in the scattered arrays " << scatterlength << std::endl;
            std::cout << "Skx_FFT_array_R_scatter.size " << Skx_FFT_array_R_scatter.size() << std::endl;
            std::cout << "trying to scatter..." << std::endl;

         }



         MPI_Scatter(&Skx_FFT_array_R_transposed[0], scatterlength, MPI_DOUBLE, &Skx_FFT_array_R_scatter[0], scatterlength, MPI_DOUBLE, 0, MPI_COMM_WORLD);
         MPI_Scatter(&Skx_FFT_array_I_transposed[0], scatterlength, MPI_DOUBLE, &Skx_FFT_array_I_scatter[0], scatterlength, MPI_DOUBLE, 0, MPI_COMM_WORLD); 

         if (vmpi::my_rank == 0){
            std::cout << "scattered" << std::endl;
         }
         for (int k = 0; k < nk_per_rank; k++){
               if (vmpi::my_rank == 0){
                  std::cout << k << std::endl;
               }
            for (int time=0; time < internal::numtimepoints; time++){

               // Fill FFTW arrays with values from spacial FFT.
               combined_real_imag[time][real] = Skx_FFT_array_R_scatter[k*internal::numtimepoints + time];
               combined_real_imag[time][imag] = Skx_FFT_array_R_scatter[k*internal::numtimepoints + time];
               // std::cout << k << " " << combined_real_imag[time][real] << " " << combined_real_imag[time][imag] << " ";
               // std::cout << Skx_FFT_array_R[k*internal::numtimepoints + time] << " " <<  Skx_FFT_array_I[k*internal::numtimepoints + time] << std::endl;

            }

            fftw_execute(fft_in_time);

            // =============================================================================================================================================
            // Calculate one sided spectrum ================================================================================================================
            // =============================================================================================================================================
            int j2;
            double os1, os2, os[internal::numtimepoints/2];

            for (int j1 = 0; j1 < internal::numtimepoints/2; j1++){
               j2 = internal::numtimepoints-j1-1;     
               os1 = combined_real_imag_fftd[j1][real] * combined_real_imag_fftd[j1][real] + combined_real_imag_fftd[j1][imag] * combined_real_imag_fftd[j1][imag];
               os2 = combined_real_imag_fftd[j2][real] * combined_real_imag_fftd[j2][real] + combined_real_imag_fftd[j2][imag] * combined_real_imag_fftd[j2][imag];

               os[j1] = os1 + os2;
            }
            // =============================================================================================================================================
            // =============================================================================================================================================
            // =============================================================================================================================================


            // =============================================================================================================================================
            // normalise the amplitudes to the largest value for each kpoint ===============================================================================
            // =============================================================================================================================================
            double largest = os[0];
            double index;

            // Find largest value in k_z array
            for (int j1 = 1; j1 < internal::numtimepoints/2; j1++){
               if (largest < os[j1]){
                  largest = os[j1];
                  index = j1;
               }
            }
            // normlise each value
            for (int j1 = 0; j1 < internal::numtimepoints/2; j1++){
               os[j1] /= largest;
            }
            // =============================================================================================================================================
            // =============================================================================================================================================
            // =============================================================================================================================================


            // =============================================================================================================================================
            // Output values to file =======================================================================================================================
            // =============================================================================================================================================
            std::ofstream file_K_time;
            std::stringstream sstr;
            sstr << "K_vs_time_" << std::setw(4) << std::setfill('0') << std::to_string(k+vmpi::my_rank*nk_per_rank) << ".dat";
            file_K_time.open(sstr.str(),std::ios_base::app);
            for (int time=0; time < internal::numtimepoints/2; time++){
               // file_K_time <<  combined_real_imag_fftd[time][real] << " " <<  combined_real_imag_fftd[time][imag] << "\n";
               file_K_time << os[time] << "\n";
            }
            // file_K_time.close();
            // =============================================================================================================================================
            // =============================================================================================================================================
            // =============================================================================================================================================

         }
         
      #else

      for (int k = 0; k < internal::numkpoints; k++){
         for (int time=0; time < internal::numtimepoints; time++){

            // Fill FFTW arrays with values from spacial FFT.
            combined_real_imag[time][real] = Skx_FFT_array_R_transposed[k*internal::numtimepoints + time];
            combined_real_imag[time][imag] = Skx_FFT_array_R_transposed[k*internal::numtimepoints + time];
            std::cout << k << " " << combined_real_imag[time][real] << " " << combined_real_imag[time][imag] << " ";
            std::cout << Skx_FFT_array_R[k*internal::numtimepoints + time] << " " <<  Skx_FFT_array_I[k*internal::numtimepoints + time] << std::endl;

         }
         fftw_execute(fft_in_time);


         // =============================================================================================================================================
         // Calculate one sided spectrum ================================================================================================================
         // =============================================================================================================================================
         int j2;
         double os1, os2, os[internal::numtimepoints/2];

         for (int j1 = 0; j1 < internal::numtimepoints/2; j1++){
            j2 = internal::numtimepoints-j1-1;     
            os1 = combined_real_imag_fftd[j1][real] * combined_real_imag_fftd[j1][real] + combined_real_imag_fftd[j1][imag] * combined_real_imag_fftd[j1][imag];
            os2 = combined_real_imag_fftd[j2][real] * combined_real_imag_fftd[j2][real] + combined_real_imag_fftd[j2][imag] * combined_real_imag_fftd[j2][imag];

            os[j1] = os1 + os2;
         }
         // =============================================================================================================================================
         // =============================================================================================================================================
         // =============================================================================================================================================






         // =============================================================================================================================================
         // normalise the amplitudes to the largest value for each kpoint ===============================================================================
         // =============================================================================================================================================
         double largest = os[0];
         double index;

         // Find largest value in k_z array
         for (int j1 = 1; j1 < internal::numtimepoints/2; j1++){
            if (largest < os[j1]){
               largest = os[j1];
               index = j1;
            }
         }
         // normlise each value
         for (int j1 = 0; j1 < internal::numtimepoints/2; j1++){
            os[j1] /= largest;
         }
         // =============================================================================================================================================
         // =============================================================================================================================================
         // =============================================================================================================================================




         // =============================================================================================================================================
         // Output values to file =======================================================================================================================
         // =============================================================================================================================================
         std::ofstream file_K_time;
         std::stringstream sstr;
         sstr << "K_vs_time_" << std::setw(4) << std::setfill('0') << std::to_string(k) << ".dat";
         file_K_time.open(sstr.str(),std::ios_base::app);
         for (int time=0; time < internal::numtimepoints/2; time++){
            // file_K_time <<  combined_real_imag_fftd[time][real] << " " <<  combined_real_imag_fftd[time][imag] << "\n";
            file_K_time << os[time] << "\n";
         }
         file_K_time.close();
         // =============================================================================================================================================
         // =============================================================================================================================================
         // =============================================================================================================================================

      }


      fftw_destroy_plan(fft_in_time);

      #endif
      
   }

   
   




} 

