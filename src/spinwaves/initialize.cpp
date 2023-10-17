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
		const double pi=3.14159265358979323846264338327;
		const double twopi_over_a=2*pi/unit_cell_size_x;
		const double twopi_over_b=2*pi/unit_cell_size_y;
		const double twopi_over_c=2*pi/unit_cell_size_z;
		const double toll_Fikj=1e-5; // tollerance for the structure factor values smaller than tollerance are considered 0
		std::ofstream file_Reciprocal_lattice;
		file_Reciprocal_lattice.open("K_space.dat");

		int N_points=total_num_unit_cells_x*5.0; // Number of points for each K-path
		std::cout<<" "<<std::endl;
		std::cout<< "System dimensions: "<<system_dimensions_x<<"\t"<<system_dimensions_y<<"\t"<<system_dimensions_z<<std::endl;
		std::cout<< "Total number of unit cells: "<<total_num_unit_cells_x<<"\t"<<total_num_unit_cells_y<<"\t"<<total_num_unit_cells_z<<std::endl;
		std::cout<< "Unit cell size: "<<unit_cell_size_x<<"\t"<<unit_cell_size_y<<"\t"<<unit_cell_size_z<<"\n"<<std::endl;
		std::cout<< "atoms in one unit cell: "<<Na<<std::endl;

		for(unsigned int uca=0; uca<atom.size(); uca++){
			std::cout<< "atom index, x,y,z: "<<atom[uca].x<<"\t"<<atom[uca].y<<"\t"<<atom[uca].z<<std::endl;
		}
		std::cout<<" "<<std::endl;


		// Read the K-space path
		int N_K_path;
		int ink=0;  // total numer of K-point generated

		double kx1,ky1,kz1,kx2,ky2,kz2; 
		std::ifstream file_read_K_path;
		file_read_K_path.open("K_path.inp");

		// check if file has been read in JRH	
		if (file_read_K_path.is_open()) {
			std::cout << "k path file succesfully opened." << std::endl;
		} 
		else {
			// File opening failed
			std::cerr << "Error: Unable to open file 'K_path.inp'" << std::endl;
		}

		//file_read_K_path >> N_K_path; JRH
		//int Nktotal=N_points*N_K_path; JRH
	

		// read specific kpoints from file JRH
		int Nktotal = 0;
		std::string line;
		while (std::getline(file_read_K_path, line)) {
			std::istringstream iss(line);
			Nktotal++;
			double val01, val02, val1, val2, val3;
			if (iss >> val01 >> val02 >> val1 >> val2 >> val3) {
				spinwaves::internal::kx_FFT_array.push_back(val1*2.0*M_PI);
				spinwaves::internal::ky_FFT_array.push_back(val2*2.0*M_PI);
				spinwaves::internal::kz_FFT_array.push_back(val3*2.0*M_PI);
			} 
			else {
				std::cerr << "Error reading line: " << line << std::endl;
			}
		}

		//std::cout << "N_points: \t" << N_points << std::endl;
		//std::cout << "N_K_path: \t" << N_K_path << std::endl;
		std::cout<<"temp:number of Kpoints is  "<<Nktotal<<"\t"<<N_K_path<<std::endl;

		spinwaves::internal::kx_FFT_array.resize(Nktotal,0.0);
		spinwaves::internal::ky_FFT_array.resize(Nktotal,0.0);
		spinwaves::internal::kz_FFT_array.resize(Nktotal,0.0);
		spinwaves::internal::structure_factor_array_R.resize(Nktotal,0.0);
		spinwaves::internal::structure_factor_array_I.resize(Nktotal,0.0);

		spinwaves::Skx_FFT_array_R.resize(Nktotal,0.0);
		spinwaves::Sky_FFT_array_R.resize(Nktotal,0.0);
		spinwaves::Skz_FFT_array_R.resize(Nktotal,0.0);
		spinwaves::Skx_FFT_array_I.resize(Nktotal,0.0);
		spinwaves::Sky_FFT_array_I.resize(Nktotal,0.0);
		spinwaves::Skz_FFT_array_I.resize(Nktotal,0.0);

		


		// JRH comment out
		//for (int i=0;i<N_K_path;i++){
		//	file_read_K_path >> kx1>>ky1>>kz1>>kx2>>ky2>>kz2;
		//	std::cout<<kx1<<"\t"<<ky1<<"\t"<<kz1<<"\t"<<kx2<<"\t"<<ky2<<"\t"<<kz2<<"\t \n"<<std::endl;
		//	for(unsigned int ik=0;ik<N_points;ik++){
		//		spinwaves::internal::kx_FFT_array[ink] = twopi_over_a*(kx1 + (kx2-kx1)*ik/(N_points+0.0));
		//		spinwaves::internal::ky_FFT_array[ink] = twopi_over_b*(ky1 + (ky2-ky1)*ik/(N_points+0.0)); 
		//		spinwaves::internal::kz_FFT_array[ink] = twopi_over_c*(kz1 + (kz2-kz1)*ik/(N_points+0.0));
		//		//std::cout<<spinwaves::internal::kx_FFT_array[ink]<<"\t"<<spinwaves::internal::ky_FFT_array[ink]<<"\t"<<spinwaves::internal::kz_FFT_array[ink]<<std::endl;
		//		ink=ink+1;
		//	}
		//}


		// JRH comment out
		//std::cout<<"Number of K points "<<Nktotal<<std::endl;
		//if (ink != Nktotal) {
		//	std::cout<<"Errorare "<<std::endl;
		//	exit(1);
		//}

		// Calculate memory requirements and inform user
		const double mem = 1e6*double(Nktotal) * ( 6*sizeof(double)) / 1.0e6;
		std::cout     << "SW required  " << mem << " MB of RAM for 1ns (10^6 time steps)" << std::endl;


		//##################################################################################################
		//#### Compute the structure factor
		for(unsigned int uca=0;uca<Nktotal;uca++){
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
			std::cout<<uca<<"\t"<<"\t"<<spinwaves::internal::structure_factor_array_R[uca]<<"\t"<<spinwaves::internal::structure_factor_array_I[uca]<<"\t"<<std::endl;
		}
		//#######################################################################################################################


		for(unsigned int uca=0;uca<Nktotal;uca++){
			std::ofstream file_K_time;
			std::stringstream sstr;
			// JRH change name of files
			//sstr << "K_vs_time_" << std::setw(4) << std::setfill('0') << std::to_string(uca) << ".dat";
			sstr << "K_vs_time" << 
			"_kx_" << std::setw(6) << std::to_string(spinwaves::internal::kx_FFT_array[uca]) << 
			"_ky_" << std::setw(6) << std::to_string(spinwaves::internal::ky_FFT_array[uca]) << 
			"_kz_" << std::setw(6) << std::to_string(spinwaves::internal::kz_FFT_array[uca]) << 		
			".dat";

			file_K_time.open(sstr.str());
			file_K_time.close();
			double Fijk=sqrt(spinwaves::internal::structure_factor_array_I[uca]*spinwaves::internal::structure_factor_array_I[uca] + spinwaves::internal::structure_factor_array_R[uca]*spinwaves::internal::structure_factor_array_R[uca]);
			file_Reciprocal_lattice << uca<<"\t"<<spinwaves::internal::kx_FFT_array[uca]<<"\t"<< spinwaves::internal::ky_FFT_array[uca]<<"\t"<< spinwaves::internal::kz_FFT_array[uca]<<"\t"<< Fijk<<"\t"<< spinwaves::internal::structure_factor_array_R[uca]<<"\t"<< spinwaves::internal::structure_factor_array_I[uca] <<std::endl;
		}

		std::cout<< "There are  " << Nktotal << " k-points in total for "<<N_K_path <<" K-paths" << std::endl;
		std::cout<< "SW initialisation completed."<<std::endl;
		std::cout<<" "<<std::endl;
		file_Reciprocal_lattice.close();
		return;

	}

} // end of sw namespace
//
