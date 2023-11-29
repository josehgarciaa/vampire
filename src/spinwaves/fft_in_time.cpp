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
#include "stopwatch.hpp"
#include "iostream"
#include <iomanip>
#include <sstream>
#include "vio.hpp"

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



      // Start time for time series fourier transform
      // Set timer for runtime
      std::cout     << "Starting calculation of time series Discrete Fourier Transform." << std::endl;
      zlog << zTs() << "Starting calculation of time series Discrete Fourier Transform." << std::endl;
      stopwatch_t fft_stopwatch;
      fft_stopwatch.start();
   
      const int real = 0;
      const int imag = 1;

      std::vector<double> os;
      os.resize(internal::nt/2,0.0);
      
      std::vector <double> skx_r_node_transposed;
      std::vector <double> skx_i_node_transposed;
      if (vmpi::my_rank == 0){
         skx_r_node_transposed.resize(internal::nt*internal::nk);
         skx_i_node_transposed.resize(internal::nt*internal::nk);	
      }


      std::vector<fftw_complex> combined_real_imag(internal::nt);
      std::vector<fftw_complex> combined_real_imag_fftd(internal::nt);
      fftw_plan fft_in_time;
      fft_in_time = fftw_plan_dft_1d(internal::nt, &combined_real_imag[0], &combined_real_imag_fftd[0], FFTW_FORWARD, FFTW_MEASURE);
      

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
      //
      //
      // Is it possible using scatterv to send to correct node in a single process....
      // To be memory efficient we need to have a
      
      // rearrange array for scatter - I think this bit needs to be done for serial and parallel
      // rearrange array for scatter - I think this bit needs to be done for serial and parallel
      // rearrange array for scatter - I think this bit needs to be done for serial and parallel
      // rearrange array for scatter - I think this bit needs to be done for serial and parallel
      // rearrange array for scatter - I think this bit needs to be done for serial and parallel
      // if (vmpi::my_rank == 0){
      //    for (int k = 0; k < internal::nk; k++){
      //       for (int time=0; time < internal::nt; time++){   
      //          // Fill FFTW arrays with values from spacial FFT.
      //          skx_r_node_transposed[k*internal::nt + time] = skx_r_node[time*internal::nk + k];
      //          skx_i_node_transposed[k*internal::nt + time] = skx_i_node[time*internal::nk + k];
      //       }     
      //    }


      // }
      // rearrange array for scatter - I think this bit needs to be done for serial and parallel
      // rearrange array for scatter - I think this bit needs to be done for serial and parallel
      // rearrange array for scatter - I think this bit needs to be done for serial and parallel
      // rearrange array for scatter - I think this bit needs to be done for serial and parallel

      


      #ifdef MPICF

         // int nk_per_rank = std::ceil(static_cast<double>(internal::nk) / static_cast<double>(vmpi::num_processors));
         // int scatterlength = nk_per_rank * internal::nt;
         // std::vector<double> skx_r_scatter;
         // std::vector<double> skx_i_scatter;
         // skx_r_scatter.resize(scatterlength,0.0);
         // skx_i_scatter.resize(scatterlength,0.0);



         if (vmpi::my_rank == 0){

            // Not always the case that Nkpoints/nranks will be an integer. Need to round up and zero pad the array.
            // For example. 100 timesteps, 50 kpoints ran on 4 ranks will mean there is (100*50) / 4 data points per rank.
            // WIll have to round up and zero-pad.

            std::cout     << "Number of k-points allocated to each rank: " << nk_per_rank << std::endl;
            zlog << zTs() << "Number of k-points allocated to each rank: " << nk_per_rank << std::endl;

            // zero pad the array
            skx_r_node_transposed.resize(nk_per_rank * vmpi::num_processors * internal::nt,0.0);
            skx_i_node_transposed.resize(nk_per_rank * vmpi::num_processors * internal::nt,0.0);
         }



         // make a mask that determines which kpoints will be dft'd. This prevents dft of lots of 0s on last rank
         std::vector<int> kmask;
         kmask.resize(nk_per_rank * vmpi::num_processors,0);
         for (int i = 0; i < internal::nk; i++){
            kmask[i] = 1;
         }
         std::cout     << "Created mask for k-points." << std::endl;
         zlog << zTs() << "Created mask for k-points." << std::endl;

         // // // scatter array to every processor
         // MPI_Scatter(&skx_r_node_transposed[0], scatterlength, MPI_DOUBLE, &skx_r_scatter[0], scatterlength, MPI_DOUBLE, 0, MPI_COMM_WORLD);
         // MPI_Scatter(&skx_i_node_transposed[0], scatterlength, MPI_DOUBLE, &skx_i_scatter[0], scatterlength, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
         // std::cout     << "Scattered points to each rank." << std::endl;
         // zlog << zTs() << "Scattered points to each rank." << std::endl;



         // loop over the k-points on each rank
         for (int k = 0; k < nk_per_rank; k++){

            // if the mask is 1, the calculate the dft.
            if (kmask[k+vmpi::my_rank*nk_per_rank] == 1){

               // populate fftw_complex vector
               for (int time=0; time < internal::nt; time++){
                  combined_real_imag[time][real] = skx_r_scatter[k*internal::nt + time];
                  combined_real_imag[time][imag] = skx_r_scatter[k*internal::nt + time];
               }

               // exexcute the fft
               fftw_execute(fft_in_time);

               // calculate one sided spectrum
               spinwaves::internal::one_sided_spectrum(os, combined_real_imag_fftd);

               // normalise to largest amplitude for each k-value
               spinwaves::internal::normalise_spectrum(os);

               // write each k-value to file
               spinwaves::internal::write_to_file(os, k, nk_per_rank);
            }
         }
         
      #else

         for (int k = 0; k < internal::nk; k++){
            for (int time=0; time < internal::nt; time++){

               // Fill FFTW arrays with values from spacial FFT.
               combined_real_imag[time][real] = skx_r_node_transposed[k*internal::nt + time];
               combined_real_imag[time][imag] = skx_r_node_transposed[k*internal::nt + time];

            }

            // exexcute the fft
            fftw_execute(fft_in_time);

            // calculate one sided spectrum
            spinwaves::internal::one_sided_spectrum(os, combined_real_imag_fftd);            

            // normalise to largest amplitude for each k-value
            spinwaves::internal::normalise_spectrum(os);

            // write each k-value to file
            spinwaves::internal::write_to_file(os, k, 0);
         }

      #endif

      // Start time for time series fourier transform
      // Set timer for runtime
      std::cout     << "Total duration of time series discrete fourier transform [s]: " << fft_stopwatch.elapsed_seconds() << std::endl;
      zlog << zTs() << "Total duration of time series discrete fourier transform [s]: " << fft_stopwatch.elapsed_seconds() << std::endl;

      // destroy fftw3 plan
      fftw_destroy_plan(fft_in_time);
      
   }
} 

