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

    namespace internal {

        std::ofstream file_K_time;
        int j2;
        double os1, os2;
        std::vector<double> os;
        const int real = 0;
        const int imag = 1;

        // =============================================================================================================================================
        // Calculate one sided spectrum ================================================================================================================
        // =============================================================================================================================================
        void one_sided_spectrum(std::vector<double>& os, 
                                std::vector<fftw_complex>& combined_real_imag_fftd){

            for (int j1 = 0; j1 < internal::nt/2; j1++){
                j2 = internal::nt-j1-1;     
                os1 = combined_real_imag_fftd[j1][real] * combined_real_imag_fftd[j1][real] + combined_real_imag_fftd[j1][imag] * combined_real_imag_fftd[j1][imag];
                os2 = combined_real_imag_fftd[j2][real] * combined_real_imag_fftd[j2][real] + combined_real_imag_fftd[j2][imag] * combined_real_imag_fftd[j2][imag];

                os[j1] = os1 + os2;
            }
        }
        

        void normalise_spectrum(std::vector<double>& os){

            double largest = os[0];
            double index;


            // Find largest value in k_z array
            for (int j1 = 1; j1 < internal::nt/2; j1++){
                if (largest < os[j1]){
                    largest = os[j1];
                    index = j1;
                }
                }
                // normlise each value
                for (int j1 = 0; j1 < internal::nt/2; j1++){
                os[j1] /= largest;
            }

        }

        void write_to_file(std::vector<double>& os, int k, int nk_per_rank){

            std::stringstream sstr;
            
            #ifdef MPICF
                sstr << "k" << std::setw(4) << std::setfill('0') << std::to_string(k+vmpi::my_rank*nk_per_rank) << ".dat";
            #else 
                sstr << "k" << std::setw(4) << std::setfill('0') << std::to_string(k) << ".dat";
            #endif


            file_K_time.open(sstr.str(),std::ios_base::app);
            for (int time=0; time < internal::nt/2; time++){
            // file_K_time <<  combined_real_imag_fftd[time][real] << " " <<  combined_real_imag_fftd[time][imag] << "\n";
            file_K_time << os[time] << "\n";
            }
            file_K_time.close();

        }

    }

}