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

}