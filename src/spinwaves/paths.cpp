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
#include "unitcell.hpp"
#include "internal.hpp"
#include "errors.hpp"
#include "atoms.hpp"
#include "vmpi.hpp"
#include "vio.hpp"
#include "sim.hpp"


// sw module headers
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>



namespace spinwaves{

    namespace internal {


        void determine_path(){

            // if kpath_filename has not been found in spiwnaves/interface.cpp use a built in path
            if (spinwaves::internal::kpath_filename==""){
                std::cout << "sw_crystal_structure = " << uc::sw_crystal_structure<< std::endl;
                
                // Determine which path to take depending on crystal structure
                if(uc::sw_crystal_structure      == "sc"                ) spinwaves::internal::path_sc();
                else if(uc::sw_crystal_structure == "bcc"               ) spinwaves::internal::path_bcc();
                else if(uc::sw_crystal_structure == "bcc-110"           ) spinwaves::internal::path_bcc();
                else if(uc::sw_crystal_structure == "fcc"               ) spinwaves::internal::path_fcc();
                else if(uc::sw_crystal_structure == "fcc-111"           ) spinwaves::internal::path_fcc();
                else if(uc::sw_crystal_structure == "hcp"               ) spinwaves::internal::path_hcp();
                else if(uc::sw_crystal_structure == "heusler"           ) spinwaves::internal::path_bcc();
                // else if(uc::sw_crystal_structure == "honeycomb"      ) spinwaves::internal::path_honeycomb();
                // else if(uc::sw_crystal_structure == "alpha-honeycomb") spinwaves::internal::path_honeycomb_alpha();
                // else if(uc::sw_crystal_structure == "beta-honeycomb" ) spinwaves::internal::path_honeycomb_beta();
                // else if(uc::sw_crystal_structure == "kagome"         ) spinwaves::internal::path_kagome();
                else if(uc::sw_crystal_structure == "mn2au"             ) spinwaves::internal::path_mn2au();
                // else if(uc::sw_crystal_structure == "NdFeB"          ) spinwaves::internal::path_NdFeB();
                // else if(uc::sw_crystal_structure == "rocksalt"       ) spinwaves::internal::path_rock_salt();
                // else if(uc::sw_crystal_structure == "spinel"         ) spinwaves::internal::path_spinel();
                // else if(uc::sw_crystal_structure == "spinel-layered" ) spinwaves::internal::path_spinel_layered();
                // else if(uc::sw_crystal_structure == "SmFeN"          ) spinwaves::internal::path_SmFeN();
                else{
                    terminaltextcolor(RED);
                    std::cerr << "Error: Unknown spinwaves crystal_type "<< uc::sw_crystal_structure << " found during spinwave path initialisation. Exiting." << std::endl;
                    terminaltextcolor(WHITE);
                    zlog << zTs() << "Error: Unknown spinwaves crystal_type "<< uc::sw_crystal_structure << " found during spinwave path initialisation. Exiting." << std::endl;
                    err::vexit();
                }
            }
            // else read the path from the file specified in interface.cpp
            else {

                // check file exists
                std::ifstream file_read_K_path;
                file_read_K_path.open(spinwaves::internal::kpath_filename);

                if (file_read_K_path.is_open()) {
                    std::cout << "Spinwaves path file " << spinwaves::internal::kpath_filename << " succesfully opened." << std::endl;
                    zlog << zTs()  << "Spinwaves path file " << spinwaves::internal::kpath_filename << " succesfully opened." << std::endl;
                } 
                else {
                    terminaltextcolor(RED);
                    std::cerr << "Error: Cannot find spinwaves path file " << spinwaves::internal::kpath_filename << ". Exiting." << std::endl;
                    terminaltextcolor(WHITE);
                    zlog << zTs() << "Error: Cannot find spinwaves path file "<< spinwaves::internal::kpath_filename << ". Exiting." << std::endl;
                    err::vexit();
                }

                // read contents of path file
                int linecount = 0;
                std::string line;
                while (std::getline(file_read_K_path, line)) {
                    std::istringstream iss(line);
                    linecount++;
                    double val1, val2, val3;
                    if (iss >> val1 >> val2 >> val3) {
                        spinwaves::internal::pathx.push_back(val1);
                        spinwaves::internal::pathy.push_back(val2);
                        spinwaves::internal::pathz.push_back(val3);
                    } 
                    else {
                        std::cerr << "Error reading line: " << line << std::endl;
                    }
                }

                std::cout << "Spinwaves path file " << spinwaves::internal::kpath_filename << " contains " << linecount << " lines." << std::endl;
                zlog << zTs()  << "Spinwaves path file" << spinwaves::internal::kpath_filename << " contains " << linecount << " lines." << std::endl;

                // close file and check
                file_read_K_path.close();
                std::cout << "Spinwaves path file " << spinwaves::internal::kpath_filename << " has been closed." << std::endl;
                zlog << zTs()  << "Spinwaves path file " << spinwaves::internal::kpath_filename << " has been closed." << std::endl;

            }
        }

        void determine_kpoints(const double system_dimensions_x,
                    const double system_dimensions_y,
                    const double system_dimensions_z,
					const double unit_cell_size_x,
					const double unit_cell_size_y,
					const double unit_cell_size_z){
            
            int len=spinwaves::internal::pathx.size();

            // Get reciprocal lattice vectors from cubic unit cell JRH 26/10/23
            // Based on equations through link: http://lampx.tugraz.at/~hadley/ss1/crystaldiffraction/fourier/reciprocal_lattice.php
            // The assumption is made that the unit cell is cubic
            double b[3];
            b[0] = 2.0 * M_PI * (unit_cell_size_y * unit_cell_size_z) / (unit_cell_size_x * (unit_cell_size_y * unit_cell_size_z));
            b[1] = 2.0 * M_PI * (unit_cell_size_z * unit_cell_size_x) / (unit_cell_size_x * (unit_cell_size_y * unit_cell_size_z));
            b[2] = 2.0 * M_PI * (unit_cell_size_x * unit_cell_size_y) / (unit_cell_size_x * (unit_cell_size_y * unit_cell_size_z));
            
            // Generate k-points
            std::cout << "Determining k-points..." << std::endl;
            // std::cout << "Len " << std::endl;

            for (int row=0; row<len; row+=2){

                double cellx = system_dimensions_x/unit_cell_size_x;
                double celly = system_dimensions_y/unit_cell_size_y;
                double cellz = system_dimensions_z/unit_cell_size_z;
                
                double kx0=spinwaves::internal::pathx[row];
                double ky0=spinwaves::internal::pathy[row];
                double kz0=spinwaves::internal::pathz[row];
                double kx1=spinwaves::internal::pathx[row+1];
                double ky1=spinwaves::internal::pathy[row+1];
                double kz1=spinwaves::internal::pathz[row+1];

                // std::cout << "kx " << kx0 << " " << kx1 << std::endl;
                // std::cout << "ky " << ky0 << " " << ky1 << std::endl;
                // std::cout << "kz " << kz0 << " " << kz1 << std::endl;
                double distancex = cellx * (kx1-kx0);
                double distancey = celly * (ky1-ky0);
                double distancez = cellz * (kz1-kz0);
                

                // std::cout << "distance: " << distancex << " " << distancey << " " << distancez << std::endl;
                int res1 = spinwaves::internal::gcd(static_cast<int>(round(distancex)), static_cast<int>(round(distancey)));
                int common_denom = spinwaves::internal::gcd(res1, static_cast<int>(round(distancez)));

                double kx=kx0;
                double ky=ky0;
                double kz=kz0;

                for (int k = 0; k < common_denom; k++){

                    kx = kx + distancex/static_cast<double>(common_denom)/cellx;
                    ky = ky + distancey/static_cast<double>(common_denom)/celly;
                    kz = kz + distancez/static_cast<double>(common_denom)/cellz;

                    // make sure to convert to units of 2pi/latconst
                    spinwaves::internal::kx_FFT_array.push_back(b[0] * kx);
                    spinwaves::internal::ky_FFT_array.push_back(b[1] * ky);
                    spinwaves::internal::kz_FFT_array.push_back(b[2] * kz);
                    // std::cout << b[0] * kx  << " " << b[1] * ky << " " << b[2] * kz << std::endl;
                    // std::cout << kx  << " " << ky << " " << kz << std::endl;
                }

            }
            
            
            // define number of time and kpoitns 
            internal::numkpoints = spinwaves::internal::kx_FFT_array.size();
            internal::numtimepoints = (sim::total_time / sim::partial_time);
            std::cout << "total number of spinwave kpoints "    << internal::numkpoints << std::endl;
            std::cout << "total number of spinwave timepoints " << internal::numtimepoints << std::endl;
            internal::structure_factor_array_R.resize(internal::numkpoints,0.0);
            internal::structure_factor_array_I.resize(internal::numkpoints,0.0);

            // This array is the calculated at each timepoint
            // spinwaves::Skx_FFT_array_R.resize(length,0.0);
            // spinwaves::Skx_FFT_array_I.resize(length,0.0);
            // const double mem = spinwaves::Skx_FFT_array_R.size()*sizeof(double)*(sim::total_time / sim::partial_time)/ 1.0e6 + 
            //                    spinwaves::Skx_FFT_array_I.size()*sizeof(double)*(sim::total_time / sim::partial_time)/ 1.0e6;


            spinwaves::Skx_FFT_array_R.resize(internal::numtimepoints*internal::numkpoints,0.0);
            spinwaves::Skx_FFT_array_I.resize(internal::numtimepoints*internal::numkpoints,0.0);	
            const double mem = spinwaves::Skx_FFT_array_R.size()*sizeof(double)/ 1.0e6 + spinwaves::Skx_FFT_array_I.size()*sizeof(double)/ 1.0e6;

            // Calculate memory requirements of k_vs_time arrays and inform user. There's 6 arrays x,y,z with real and imag.
            std::cout << "Spinwave module using  " << mem << " MB of RAM for k vs time arrays for a simulation time of " << sim::total_time;
            std::cout << " timesteps with the dynamic structure factor calculated every " << sim::partial_time << " timesteps." << std::endl;
        }


        void calculate_fourier_prefactor(const std::vector<double>& rx, const std::vector<double>& ry, const std::vector<double>& rz){
            
            double arg;
            int nk=spinwaves::internal::kx_FFT_array.size();
            spinwaves::internal::cos_k.resize(nk*atoms::num_atoms,0.0);
            spinwaves::internal::sin_k.resize(nk*atoms::num_atoms,0.0);

            for (int k=0; k < nk; k++){

                // std::cout << "test1" << std::endl;
                double kx = spinwaves::internal::kx_FFT_array[k];
                double ky = spinwaves::internal::ky_FFT_array[k];
                double kz = spinwaves::internal::kz_FFT_array[k];
                // std::cout << "test2" << std::endl;
           
                #ifdef MPICF

                for(int atom=0;atom<vmpi::num_core_atoms+vmpi::num_bdry_atoms;atom++){
                    arg=-kx*rx[atom] - ky*ry[atom] - kz*rz[atom];  
                    spinwaves::internal::cos_k[k*atoms::num_atoms+atom] = cos(arg);
                    spinwaves::internal::sin_k[k*atoms::num_atoms+atom] = sin(arg);                    
                }

                #else 

                for(int atom=0;atom<atoms::num_atoms;atom++){
                    arg=-kx*rx[atom] - ky*ry[atom] - kz*rz[atom];  
                    spinwaves::internal::cos_k[k*atoms::num_atoms+atom] = cos(arg);
                    spinwaves::internal::sin_k[k*atoms::num_atoms+atom] = sin(arg);                      
                }
            
                #endif

            }
            

            // Calculate memory requirements and inform user
            double mem = 0.0;
            mem += spinwaves::internal::cos_k.size()*sizeof(double) / 1.0e6;
            mem += spinwaves::internal::sin_k.size()*sizeof(double) / 1.0e6;

            std::cout << "Spinwave module using " << mem << " MB of RAM for fourier prefactor." << std::endl;

        }
    } 
} 