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
#include <string>
#include <sstream>

//sergiu for SW
#include "unitcell.hpp"
#include "errors.hpp"
#include "sim.hpp"
#include "vector"
#include <cmath>
#include "fstream"
#include "atoms.hpp"



namespace spinwaves{
	//   std::vector <double> Skx_FFT_array;
	//   std::vector <double> Sky_FFT_array;
	//   std::vector <double> Skz_FFT_array;

	//----------------------------------------------------------------------------
	// Function to initialize sw module
	//----------------------------------------------------------------------------
	void initialize(const double system_dimensions_x,
					const double system_dimensions_y,
					const double system_dimensions_z,
					const double total_num_unit_cells_x,
					const double total_num_unit_cells_y,
					const double total_num_unit_cells_z,
					const double unit_cell_size_x,
					const double unit_cell_size_y,
					const double unit_cell_size_z,
					const std::vector <unitcell::atom_t>& atom,
					const std::vector<double>& atom_coords_x,
					const std::vector<double>& atom_coords_y,
					const std::vector<double>& atom_coords_z){


		std::cout<< "Test Spin waves...."<<std::endl;
		int Na=atom.size();
		const double twopi_over_a=2.0*M_PI/unit_cell_size_x;
		const double twopi_over_b=2.0*M_PI/unit_cell_size_y;
		const double twopi_over_c=2.0*M_PI/unit_cell_size_z;
		const double toll_Fikj=1e-5; // tollerance for the structure factor values smaller than tollerance are considered 0

		// number of steps for sw calculations
		int numsteps = (sim::total_time / sim::partial_time);

		
		// std::ofstream file_Reciprocal_lattice;
		// file_Reciprocal_lattice.open("K_space.dat");

		std::cout << " "<<std::endl;
		std::cout << "System dimensions: "<<system_dimensions_x<<"\t"<<system_dimensions_y<<"\t"<<system_dimensions_z<<std::endl;
		std::cout << "Total number of unit cells: "<<total_num_unit_cells_x<<"\t"<<total_num_unit_cells_y<<"\t"<<total_num_unit_cells_z<<std::endl;
		std::cout << "Unit cell size: "<<unit_cell_size_x<<"\t"<<unit_cell_size_y<<"\t"<<unit_cell_size_z<<"\n"<<std::endl;
		std::cout << "atoms in one unit cell: "<<Na<<std::endl;

		for(unsigned int uca=0; uca<atom.size(); uca++){
			std::cout<< "atom index, x,y,z: "<<atom[uca].x<<"\t"<<atom[uca].y<<"\t"<<atom[uca].z<<std::endl;
		}
		std::cout<<" "<<std::endl;


		// JRH determine whether to use path for a predefined crystal or whether path has been specified by user.
		spinwaves::internal::determine_path();
		std::cout << "PATH DETERMINED " << std::endl;

		// JRH calculate the values of $k$ at which to calculate the DSF for the user specified system dimensions.
		spinwaves::internal::determine_kpoints(system_dimensions_x, system_dimensions_y, system_dimensions_z, unit_cell_size_x, unit_cell_size_y, unit_cell_size_z);
		std::cout << "KPOINT DETERMINED " << std::endl;

		// determine prefactor that will be used in fourier transform. sin(k_x*r_x) etc.
		spinwaves::internal::calculate_fourier_prefactor(atom_coords_x, atom_coords_y, atom_coords_x);
		std::cout << "FOURIER PREFACTOR DETERMINED " << std::endl;

		// err::vexit();

		//##################################################################################################
		//#### Compute the structure factor
		for(unsigned int uca=0;uca<spinwaves::internal::kx_FFT_array.size();uca++){
			const double kx=spinwaves::internal::kx_FFT_array[uca];
			const double ky=spinwaves::internal::ky_FFT_array[uca];
			const double kz=spinwaves::internal::kz_FFT_array[uca];
			for(int atom=0;atom<atoms::num_atoms;atom++){
				const double rx=atom_coords_x[atom];
				const double ry=atom_coords_y[atom];
				const double rz=atom_coords_z[atom];
				double arg=-kx*rx - ky*ry - kz*rz ;
				double cosK= cos(arg);
				double sinK= sin(arg);
				double sx=1;
				spinwaves::internal::structure_factor_array_R[uca] += sx*cosK;
				spinwaves::internal::structure_factor_array_I[uca] += sx*sinK;
			}
			// std::cout<<uca<<"\t"<<"\t"<<spinwaves::internal::structure_factor_array_R[uca]<<"\t"<<spinwaves::internal::structure_factor_array_I[uca]<<"\t"<<std::endl;
		}
		//#######################################################################################################################

		// for(unsigned int uca=0;uca<spinwaves::Skx_FFT_array_R.size();uca++){
		// 	std::ofstream file_K_time;
		// 	std::stringstream sstr;
		// 	// JRH change name of files
		// 	sstr << "K_vs_time_" << std::setw(4) << std::setfill('0') << std::to_string(uca) << ".dat";
		// 	//sstr << "K_vs_time" << 
		// 	//"_kx_" << std::setw(6) << std::to_string(spinwaves::internal::kx_FFT_array[uca]) << 
		// 	//"_ky_" << std::setw(6) << std::to_string(spinwaves::internal::ky_FFT_array[uca]) << 
		// 	//"_kz_" << std::setw(6) << std::to_string(spinwaves::internal::kz_FFT_array[uca]) << 		
		// 	//".dat";

		// 	file_K_time.open(sstr.str());
		// 	file_K_time.close();
		// 	double Fijk=sqrt(spinwaves::internal::structure_factor_array_I[uca]*spinwaves::internal::structure_factor_array_I[uca] + spinwaves::internal::structure_factor_array_R[uca]*spinwaves::internal::structure_factor_array_R[uca]);
		// 	file_Reciprocal_lattice << uca<<"\t"<<spinwaves::internal::kx_FFT_array[uca]<<"\t"<< spinwaves::internal::ky_FFT_array[uca]<<"\t"<< spinwaves::internal::kz_FFT_array[uca]<<"\t"<< Fijk<<"\t"<< spinwaves::internal::structure_factor_array_R[uca]<<"\t"<< spinwaves::internal::structure_factor_array_I[uca] <<std::endl;
		// }

		// std::cout<< "There are  " << Nktotal << " k-points in total for "<<N_K_path <<" K-paths" << std::endl;
		std::cout<< "SW initialisation completed."<<std::endl;
		std::cout<<" "<<std::endl;
		// file_Reciprocal_lattice.close();


		return;

	}

} // end of sw namespace
//
