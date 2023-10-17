//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) JoelHirst 2018. All rights reserved.
//
//   Email: j.r.hirst@hallam.shu.ac.uk
//
//------------------------------------------------------------------------------
//

#ifndef SPINWAVES_H_
#define SPINWAVES_H_

// C++ standard library headers
#include <string>

// Vampire headers
#include "spinwaves.hpp"

//--------------------------------------------------------------------------------
// Namespace for variables and functions for spinwaves module
//--------------------------------------------------------------------------------
namespace spinwaves{

   //-----------------------------------------------------------------------------
   // Function to initialise spinwaves module
   //-----------------------------------------------------------------------------
   void initialize();

   //---------------------------------------------------------------------------
   // Function to process input file parameters for spinwaves module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line);

   //---------------------------------------------------------------------------
   // Function to process material parameters
   //---------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index);

} // end of spinwaves namespace

#endif //SPINWAVES_H_
