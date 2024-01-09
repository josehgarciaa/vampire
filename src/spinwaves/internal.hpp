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

#ifndef SW_INTERNAL_H_
#define SW_INTERNAL_H_
//
//---------------------------------------------------------------------
// This header file defines shared internal data structures and
// functions for the sw module. These functions and
// variables should not be accessed outside of this module.
//---------------------------------------------------------------------

// C++ standard library headers

// Vampire headers
#include "spinwaves.hpp"

// sw module headers
#include "internal.hpp"
#include <fftw3.h>
#include "vector"

namespace spinwaves{

   namespace internal{

      //-------------------------------------------------------------------------
      // Internal data type definitions
      //-------------------------------------------------------------------------

      //-----------------------------------------------------------------------------
      // internal materials class for storing material parameters
      //-----------------------------------------------------------------------------
      class mp_t{

          private:

          public:

             //------------------------------
             // material parameter variables
             //------------------------------
             double test;

             // constructor
             mp_t (const unsigned int max_materials = 100):
                test(0.0) // constructor initialisation of test variable
             {
                // constructor body for initialising more complex data/arrays
             }; // end of constructor

       }; // end of internal::mp class

      //-------------------------------------------------------------------------
      // Internal shared variables
      //-------------------------------------------------------------------------

      extern bool enabled; // bool to enable module

      extern std::vector<internal::mp_t> mp; // array of material properties
      extern std::vector <double> kx_list;
      extern std::vector <double> ky_list;
      extern std::vector <double> kz_list;
      extern std::vector <double> structure_factor_array_R, structure_factor_array_I;
      extern std::string kpath_filename;

      // JRH path vectors 
      extern std::vector<double> pathx;
      extern std::vector<double> pathy;
      extern std::vector<double> pathz;

      // JRH fourier prefactors
      extern std::vector<double> cos_k;
      extern std::vector<double> sin_k;

      //JRH time variabls
      extern int nt;
      extern int nk;

      // JRH fft-in-time variables
      extern std::string reduc_ver;
      extern bool oss, cm, normk;

      //-------------------------------------------------------------------------
      // Internal function declarations
      //-------------------------------------------------------------------------

      void path_sc();
      extern void path_bcc();
      extern void path_fcc();
      extern void path_hcp();
      extern void path_mn2au();
      extern void determine_path();
      extern void determine_kpoints(
                     const double system_dimensions_x, 
                     const double system_dimensions_y, 
                     const double system_dimensions_z,	
                     const double unit_cell_size_x, 
                     const double unit_cell_size_y, 
                     const double unit_cell_size_z);
      extern void calculate_fourier_prefactor(const std::vector<double>& rx, const std::vector<double>& ry, const std::vector<double>& rz);

      // post analysis functions
      extern void normalise_each_kpoint(std::vector<fftw_complex>& os);
      extern void write_to_file(std::vector<fftw_complex>& os, int k);
      extern void one_sided_spectrum(std::vector<fftw_complex>& os);
      extern void complex_magnitude(std::vector<fftw_complex>& os);


      extern int gcd(int a, int b);

      // Peak finding at some point?????/
      // extern void write_peaks_to_file();
      // extern void calculate_spectrum_peaks();
      

   } // end of internal namespace

} // end of sw namespace

#endif //SW_INTERNAL_H_
