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
      extern std::vector <double> kx_FFT_array;
      extern std::vector <double> ky_FFT_array;
      extern std::vector <double> kz_FFT_array;
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
      extern int numtimepoints;
      extern int numkpoints;

      //-------------------------------------------------------------------------
      // Internal function declarations
      //-------------------------------------------------------------------------

      void path_sc();
      extern void path_bcc();
      extern void path_fcc();
      extern void path_hcp();
      // extern void path_heusler();
      // extern void path_honeycomb();
      // extern void path_honeycomb_alpha();
      // extern void path_honeycomb_beta();
      // extern void path_kagome();
      extern void path_mn2au();
      // extern void path_NdFeB();
      // extern void path_rock_salt();
      // extern void path_spinel();
      // extern void path_spinel_layered();
      // extern void path_SmFeN();
      extern void determine_path();
      extern void determine_kpoints(
                     const double system_dimensions_x, 
                     const double system_dimensions_y, 
                     const double system_dimensions_z,	
                     const double unit_cell_size_x, 
                     const double unit_cell_size_y, 
                     const double unit_cell_size_z);
      extern void calculate_fourier_prefactor(const std::vector<double>& rx, const std::vector<double>& ry, const std::vector<double>& rz);
      extern int gcd(int a, int b);
      

   } // end of internal namespace

} // end of sw namespace

#endif //SW_INTERNAL_H_
