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
#include <algorithm>
#include <string>

// Vampire headers
#include "spinwaves.hpp"
#include "errors.hpp"
#include "vio.hpp"

// sw module headers
#include "internal.hpp"

namespace spinwaves{

   //---------------------------------------------------------------------------
   // Function to process input file parameters for sw module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line){

      // Check for valid key, if no match return false
      std::string prefix="spinwaves";
      if(key!=prefix) return false;

      //  ------------------------------------------------------------------
      // check if file has been specified for k path -----------------------
      //  ------------------------------------------------------------------
      std::string test="kfile";
      if(word==test){
         std::string kpath_file=value;
         // strip quotes
         kpath_file.erase(std::remove(kpath_file.begin(), kpath_file.end(), '\"'), kpath_file.end());
         test="";
         // if filename not blank set ucf file name
         if(kpath_file!=test){
            spinwaves::internal::kpath_filename=kpath_file;
            return true;
         }
         else{
            terminaltextcolor(RED);
            std::cerr << "Error - empty filename in control statement \'spinwaves:" << word << "\' on line " << line << " of input file" << std::endl;
            terminaltextcolor(WHITE);
            return false;
         }
      }
      //  ------------------------------------------------------------------


      //  ------------------------------------------------------------------
      // smoothing for each kpoint -----------------------------------------
      //  ------------------------------------------------------------------
      test="one-sided";
      if(word==test){

         std::cout << value << std::endl;
         if (value == "false"){
            internal::oss = false;
            return true;
         }
         else if (value == "true"){
            internal::oss = true;
            return true;
         }
         else {
            terminaltextcolor(RED);
            std::cerr << "Error - Unknown value in control statement \'spinwaves:" << word << " = " << value << "\' on line " << line << " of input file" << std::endl;
            terminaltextcolor(WHITE);
            return false;
         }
      }
      //  ------------------------------------------------------------------

      //  ------------------------------------------------------------------
      // component of magnetisation to use for fourier transform -----------
      //  ------------------------------------------------------------------
      test="complex-magnitude";
      if(word==test){
         if (value == "false"){
            internal::cm = false;
            return true;
         }
         else if (value == "true"){
            internal::cm = true;
            return true;
         }
         else {
            terminaltextcolor(RED);
            std::cerr << "Error - Unknown value in control statement \'spinwaves:" << word << " = " << value << "\' on line " << line << " of input file" << std::endl;
            terminaltextcolor(WHITE);
            return false;
         }
      }
      //  ----------------------------------------------------------------

      //  ------------------------------------------------------------------
      // component of magnetisation to use for fourier transform -----------
      //  ------------------------------------------------------------------
      test="normalise-each-k";
      if(word==test){
         if (value == "false"){
            internal::normk = false;
            return true;
         }
         else if (value == "true"){
            internal::normk = true;
            return true;
         }
         else {
            terminaltextcolor(RED);
            std::cerr << "Error - Unknown value in control statement \'spinwaves:" << word << " = " << value << "\' on line " << line << " of input file" << std::endl;
            terminaltextcolor(WHITE);
            return false;
         }
      }
      //  ------------------------------------------------------------------
      
      //  ------------------------------------------------------------------
      // component of magnetisation to use for fourier transform -----------
      //  ------------------------------------------------------------------
      test="reduction-method";
      if(word==test){
         std::string reduc_string=value;
         reduc_string.erase(std::remove(reduc_string.begin(), reduc_string.end(), '\"'), reduc_string.end());
         internal::reduc_ver=reduc_string;

         if (internal::reduc_ver == "rank0" || internal::reduc_ver == "direct_scatter"){
            return true;
         }
         else {
            terminaltextcolor(RED);
            std::cerr << "Error - Unexpected method for MPI reduction in spinwave model. \'spinwaves:" << word << " = " << value << "\' on line " << line << " of input file" << std::endl;
            terminaltextcolor(WHITE);
            return false;
         }
      }
      //  -----

      //--------------------------------------------------------------------
      // Keyword not found
      //--------------------------------------------------------------------
      return false;

   }

   //---------------------------------------------------------------------------
   // Function to process material parameters
   //---------------------------------------------------------------------------
   bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index){

      // add prefix string
      std::string prefix="material:";

      // Check for material id > current array size and if so dynamically expand mp array
      if((unsigned int) super_index + 1 > internal::mp.size() && super_index + 1 < 101) internal::mp.resize(super_index + 1);

      //--------------------------------------------------------------------
      // Keyword not found
      //--------------------------------------------------------------------
      return false;

   }

} // end of sw namespace

