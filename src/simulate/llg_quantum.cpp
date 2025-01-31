//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans 2020. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// Standard Libraries
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>
#include <complex>
#include <iomanip>

// Library for FFT
#include <fftw3.h>

// Vampire Header files
#include "atoms.hpp"
#include "errors.hpp"
#include "LLG.hpp"
#include "material.hpp"
#include "sim.hpp"

namespace LLGQ_arrays{

    class FFTWPlanCache {
    public:
        // Singleton access
        static FFTWPlanCache& instance() {
            static FFTWPlanCache cache;
            return cache;
        }

        // Get or create plans for given size
        std::pair<fftw_plan, fftw_plan> get_plans(int N, double* in, 
                fftw_complex* out, double* result) {
            if (N != cached_N) {
                cleanup();
                forward = fftw_plan_dft_r2c_1d(N, in, out, FFTW_MEASURE);
                backward = fftw_plan_dft_c2r_1d(N, out, result, FFTW_MEASURE);
                cached_N = N;
            }
            return {forward, backward};
        }

        // Cleanup on destruction
        ~FFTWPlanCache() {
            cleanup();
        }

    private:
        // Private constructor for singleton
        FFTWPlanCache() : cached_N(0), forward(nullptr), backward(nullptr) {}
        
        // Prevent copying
        FFTWPlanCache(const FFTWPlanCache&) = delete;
        FFTWPlanCache& operator=(const FFTWPlanCache&) = delete;

        void cleanup() {
            if (cached_N > 0) {
                fftw_destroy_plan(forward);
                fftw_destroy_plan(backward);
                forward = nullptr;
                backward = nullptr;
                cached_N = 0;
            }
        }

        int cached_N;
        fftw_plan forward;
        fftw_plan backward;
    };


    void calculate_random_fields(int realizations, int n, double dt, double T);
    void assign_unique_indices(int realizations);
    void precompute_sqrt_PSD(int n, double dt, double T);
    double PSD(const double& omega, const double& T);

	// Local arrays for LLG integration
	std::vector <double> x_euler_array;
	std::vector <double> y_euler_array;
	std::vector <double> z_euler_array;

	std::vector <double> x_heun_array;
	std::vector <double> y_heun_array;
	std::vector <double> z_heun_array;

	std::vector <double> x_spin_storage_array;
	std::vector <double> y_spin_storage_array;
	std::vector <double> z_spin_storage_array;

	std::vector <double> x_initial_spin_array;
	std::vector <double> y_initial_spin_array;
	std::vector <double> z_initial_spin_array;

	std::vector<double> sqrt_PSD_buffer;

    void assign_unique_indices(int realizations) {
        const int num_atoms = atoms::num_atoms;

        // Ensure we have enough unique indices for all atoms
        if (realizations < num_atoms * 3) {
            throw std::runtime_error("Not enough realizations for unique indices");
        }

        std::vector<int> all_indices(realizations);
        std::iota(all_indices.begin(), all_indices.end(), 0); // Fill with 0, 1, ..., realizations-1

        atoms::atom_idx_x.resize(num_atoms);
        atoms::atom_idx_y.resize(num_atoms);
        atoms::atom_idx_z.resize(num_atoms);

        for (int atom = 0; atom < num_atoms; atom++) {
            int base_index = atom * 3;
            atoms::atom_idx_x[atom] = all_indices[base_index];
            atoms::atom_idx_y[atom] = all_indices[base_index + 1];
            atoms::atom_idx_z[atom] = all_indices[base_index + 2];
        }
    }


    void calculate_random_fields(int realizations, int n, double dt, double T) {
        // Pre-allocate result vectors
        atoms::noise_field.clear();
        const double S0 = mp::material[0].S0;
        const double inv_sqrt_S0 = 1.0 / std::sqrt(S0);
        
        
        // Pre-allocate all buffers
        double* __restrict in = (double*)fftw_malloc(sizeof(double) * n);
        fftw_complex* __restrict out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (n/2 + 1));
        double* __restrict result = (double*)fftw_malloc(sizeof(double) * n);
        
        // Get cached FFTW plans once
        auto& planner = FFTWPlanCache::instance();
        auto [forward, backward] = planner.get_plans(n, in, out, result);
        
        // Pre-calculate constants
        const double df = 1.0 / (n * dt);
        const double two_pi = 2.0 * M_PI;
        const double norm_factor = 1.0 / n;
        
        // Setup random number generation
        static thread_local std::random_device rd;
        static thread_local std::mt19937 gen(rd());
        std::normal_distribution<> dist(0.0, 1.0 / std::sqrt(dt));
        
        atoms::noise_field.resize(realizations, std::vector<double>(n));


        // Progress bar setup
        const int bar_width = 50;
        int last_printed_percent = -1;


        // Precompute square root of PSD
        if (sqrt_PSD_buffer.empty()) {
            precompute_sqrt_PSD(n, dt, T);
        }

        for (int r = 0; r < realizations; ++r) {
            // Generate white noise directly into input buffer
            for (int i = 0; i < n; ++i) {
                in[i] = dist(gen);
            }

            // Forward FFT
            fftw_execute(forward);

            // Apply PSD
            for (int i = 0; i <= n/2; i++) {
                const double magnitude = LLGQ_arrays::sqrt_PSD_buffer[i];
                out[i][0] *= magnitude;
                out[i][1] *= magnitude;
            }

            // Inverse FFT
            fftw_execute(backward);

            // Store normalized result directly in noise_field
            for (int j = 0; j < n; ++j) {
                atoms::noise_field[r][j] = result[j] * norm_factor * inv_sqrt_S0;
            }

            // Progress bar update
            int current_percent = static_cast<int>(std::round((r + 1) * 100.0 / realizations));
            if (current_percent > last_printed_percent) {
                float progress = static_cast<float>(r + 1) / realizations;
                int pos = static_cast<int>(bar_width * progress);

                std::cout << "\r[";
                for (int i = 0; i < bar_width; ++i) {
                    if (i < pos) std::cout << "=";
                    else if (i == pos) std::cout << ">";
                    else std::cout << " ";
                }
                std::cout << "] " << std::setw(3) << current_percent << "%";
                std::cout.flush();

                last_printed_percent = current_percent;
            }
        }

        // Cleanup
        fftw_free(in);
        fftw_free(out);
        fftw_free(result);
        
        std::cout << std::endl;
    }

	double PSD(const double& omega, const double& T) {
		const double A = mp::material[0].A;
		const double Gamma = ::mp::material[0].Gamma;
		const double omega0 = ::mp::material[0].omega0;

		double x =   omega / (2 *  T);  // hbar and kB constants
		double lorentzian = A * Gamma *    omega / ((omega0 * omega0 - omega * omega) * (omega0 * omega0 - omega * omega) + Gamma * Gamma * omega * omega);
		double coth = (x < 1e-10) ? 1.0 / x : 1.0 / tanh(x);  // Stabilize coth calculation near zero

	
		switch (sim::noise_type) {
			case 0: // Classical Noise
				return 2*T* A * Gamma / ((omega0 * omega0 - omega * omega) * (omega0 * omega0 - omega * omega) + Gamma * Gamma * omega * omega);
			case 1: // Quantum Noise
				if (omega==0) return 2*T* A * Gamma / (omega0 * omega0  * omega0 * omega0);
				else return coth * lorentzian;
			case 2: // Semiquantum Noise
				if (omega==0) return 2*T* A * Gamma / (omega0 * omega0  * omega0 * omega0);
				else return (coth-1) * lorentzian;

			default:
				std::cout << "Default Noise: Classical Noise" << std::endl;
				return 1.0 / x * lorentzian;
		}
	}

    void precompute_sqrt_PSD(int n, double dt, double T) {
        sqrt_PSD_buffer.resize(n/2 + 1);
        double df = 1.0 / (n * dt);
        for (int i = 0; i <= n/2; ++i) {
            double omega = 2.0 * M_PI * i * df;
            sqrt_PSD_buffer[i] = std::sqrt(PSD(omega, T));
        }
    }




	bool LLG_set=false; ///< Flag to define state of LLG arrays (initialised/uninitialised)

}

namespace sim{

/// @brief LLG Initialisation function
///
/// @details Resizes arrays used for Heun integration
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2011. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    07/02/2011
///
/// @return EXIT_SUCCESS
///
/// @internal
///	Created:		07/02/2011
///	Revision:	  ---
///=====================================================================================
///



int LLGQinit(){

    //std::cout << "LLG Quantum Init" << std::endl;

	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "sim::LLG_init has been called" << std::endl;}

	using namespace LLGQ_arrays;

	x_spin_storage_array.resize(atoms::num_atoms,0.0);
	y_spin_storage_array.resize(atoms::num_atoms,0.0);
	z_spin_storage_array.resize(atoms::num_atoms,0.0);

	x_initial_spin_array.resize(atoms::num_atoms,0.0);
	y_initial_spin_array.resize(atoms::num_atoms,0.0);
	z_initial_spin_array.resize(atoms::num_atoms,0.0);

	x_euler_array.resize(atoms::num_atoms,0.0);
	y_euler_array.resize(atoms::num_atoms,0.0);
	z_euler_array.resize(atoms::num_atoms,0.0);

	x_heun_array.resize(atoms::num_atoms,0.0);
	y_heun_array.resize(atoms::num_atoms,0.0);
	z_heun_array.resize(atoms::num_atoms,0.0);

    atoms::x_w_array.resize(atoms::num_atoms, 0.0);
    atoms::y_w_array.resize(atoms::num_atoms, 0.0);
    atoms::z_w_array.resize(atoms::num_atoms, 0.0);
    atoms::x_v_array.resize(atoms::num_atoms, 0.0);
    atoms::y_v_array.resize(atoms::num_atoms, 0.0);
    atoms::z_v_array.resize(atoms::num_atoms, 0.0);

	// Set number of realizations (full field)
	const int num_atoms = atoms::num_atoms;
    int realizations = num_atoms * 3 + 4;
	atoms::noise_index = 0;

    // Disable external thermal field calculations
    sim::hamiltonian_simulation_flags[3] = 0;


    //std::cout << "Realizations: " << realizations << std::endl;
	// Assign unique indices for random fields
	assign_unique_indices(realizations);

	// Precompute square root of PSD
	const int n = static_cast<int>(sim::equilibration_time) + 1;
    //std::cout << "Precompute sqrt PSD" << std::endl;
	precompute_sqrt_PSD(n, mp::dt, sim::temperature);
    //std::cout << "Calculate random fields" << std::endl;
    calculate_random_fields(realizations, n, mp::dt, sim::temperature);

	LLG_set=true;

  	return EXIT_SUCCESS;
}

namespace internal{

/// @brief LLG Heun Integrator Corrector
///
/// @callgraph
/// @callergraph
///
/// @details Integrates the system using the LLG and Heun solver
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2011. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    07/02/2011
///
/// @return EXIT_SUCCESS
///
/// @internal
///	Created:		05/02/2011
///	Revision:	  ---
///=====================================================================================
///
// ODE system for spin dynamics
void spinDynamics(const std::vector<double>& y, const std::vector<double>& b, std::vector<double>& dydt) {
    const double A = mp::material[0].A;
    const double Gamma = ::mp::material[0].Gamma;
    const double omega0 = ::mp::material[0].omega0;
    const double gyromagnetic_ratio = 1.;

    // Compute dS/dt
    dydt[0] = gyromagnetic_ratio * (y[1]*(b[2]+y[5]) - y[2]*(b[1]+y[4]));
    dydt[1] = gyromagnetic_ratio * (y[2]*(b[0]+y[3]) - y[0]*(b[2]+y[5]));
    dydt[2] = gyromagnetic_ratio * (y[0]*(b[1]+y[4]) - y[1]*(b[0]+y[3]));

    // Compute dV/dt
    dydt[3] = y[6];
    dydt[4] = y[7];
    dydt[5] = y[8];

	// Compute dW/dt
    dydt[6] = -omega0*omega0 * y[3] - Gamma * y[6] + A * gyromagnetic_ratio * y[0];
    dydt[7] = -omega0*omega0 * y[4] - Gamma * y[7] + A * gyromagnetic_ratio * y[1];
    dydt[8] = -omega0*omega0 * y[5] - Gamma * y[8] + A * gyromagnetic_ratio * y[2];
}

void llg_quantum_step(){
    //std::cout << "LLG Quantum Step" << std::endl;

	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "sim::mLLG has been called" << std::endl;}

	using namespace LLGQ_arrays;

	// Check for initialisation of LLG integration arrays
	if(LLG_set==false) sim::LLGQinit();

    //std::cout << "LLG initialized" << std::endl;
	// Local variables
	const int num_atoms = atoms::num_atoms;
    std::vector<double> y(9);
	std::vector<double> H(3);
    std::vector<double> dydt(9);

    //std::cout << "Storage arrays" << std::endl;
	// Storage arrays
    std::vector<std::vector<double>> k1_storage(num_atoms, std::vector<double>(9));
    std::vector<std::vector<double>> k2_storage(num_atoms, std::vector<double>(9));
    std::vector<std::vector<double>> k3_storage(num_atoms, std::vector<double>(9));
    std::vector<std::vector<double>> k4_storage(num_atoms, std::vector<double>(9));
    std::vector<std::vector<double>> y_pred_storage(num_atoms, std::vector<double>(9));
    std::vector<std::vector<double>> y_in_storage(num_atoms, std::vector<double>(9));
	int idx_x, idx_y, idx_z;



    //std::cout << "Store initial state" << std::endl;
	// Store initial state
    for (int atom = 0; atom < num_atoms; atom++) {
        //std::cout << "Atom: " << atom << std::endl;
        y_in_storage[atom][0] = atoms::x_spin_array[atom];
        y_in_storage[atom][1] = atoms::y_spin_array[atom];
        y_in_storage[atom][2] = atoms::z_spin_array[atom];
        y_in_storage[atom][3] = atoms::x_v_array[atom];
        y_in_storage[atom][4] = atoms::y_v_array[atom];
        y_in_storage[atom][5] = atoms::z_v_array[atom];
        y_in_storage[atom][6] = atoms::x_w_array[atom];
        y_in_storage[atom][7] = atoms::y_w_array[atom];
        y_in_storage[atom][8] = atoms::z_w_array[atom];
    }


    //std::cout << "Calculate fields" << std::endl;
	// Calculate fields
	calculate_spin_fields(0, num_atoms);
    calculate_external_fields(0, num_atoms);


	//std::cout << "K1 Step" << std::endl;


    // K1 Step
    for (int atom = 0; atom < num_atoms; ++atom) {
		
		idx_x = atoms::atom_idx_x[atom];
		idx_y = atoms::atom_idx_y[atom];
		idx_z = atoms::atom_idx_z[atom];

        H[0] = atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom] + atoms::noise_field[idx_x][atoms::noise_index];
        H[1] = atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom] + atoms::noise_field[idx_y][atoms::noise_index];
        H[2] = atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom] + atoms::noise_field[idx_z][atoms::noise_index];

        spinDynamics(y_in_storage[atom], H, k1_storage[atom]);

        for (size_t i = 0; i < 9; ++i) {
            y_pred_storage[atom][i] = y_in_storage[atom][i] + 0.5 * mp::dt * k1_storage[atom][i];
        }

        double S_magnitude = std::sqrt(y_pred_storage[atom][0] * y_pred_storage[atom][0] + y_pred_storage[atom][1] * y_pred_storage[atom][1] + y_pred_storage[atom][2] * y_pred_storage[atom][2]);

        y_pred_storage[atom][0] /= S_magnitude;
        y_pred_storage[atom][1] /= S_magnitude;
        y_pred_storage[atom][2] /= S_magnitude;

        atoms::x_spin_array[atom] = y_pred_storage[atom][0];
        atoms::y_spin_array[atom] = y_pred_storage[atom][1];
        atoms::z_spin_array[atom] = y_pred_storage[atom][2];
        atoms::x_v_array[atom] = y_pred_storage[atom][3];
        atoms::y_v_array[atom] = y_pred_storage[atom][4];
        atoms::z_v_array[atom] = y_pred_storage[atom][5];
        atoms::x_w_array[atom] = y_pred_storage[atom][6];
        atoms::y_w_array[atom] = y_pred_storage[atom][7];
        atoms::z_w_array[atom] = y_pred_storage[atom][8];
    }


	// Update fields for k2
    calculate_spin_fields(0, num_atoms);
    //calculate_thermal_fields_half(0, num_atoms);

    // K2 Step
    for (int atom = 0; atom < num_atoms; ++atom) {
        idx_x = atoms::atom_idx_x[atom];
		idx_y = atoms::atom_idx_y[atom];
		idx_z = atoms::atom_idx_z[atom];

        H[0] = atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom] + (atoms::noise_field[idx_x][atoms::noise_index] + atoms::noise_field[idx_x][atoms::noise_index + 1]) / 2.0;
        H[1] = atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom] + (atoms::noise_field[idx_y][atoms::noise_index] + atoms::noise_field[idx_y][atoms::noise_index + 1]) / 2.0;
        H[2] = atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom] + (atoms::noise_field[idx_z][atoms::noise_index] + atoms::noise_field[idx_z][atoms::noise_index + 1]) / 2.0;

        spinDynamics(y_pred_storage[atom], H, k2_storage[atom]);

        for (size_t i = 0; i < 9; ++i) {
            y_pred_storage[atom][i] = y_in_storage[atom][i] + 0.5 * mp::dt * k2_storage[atom][i];
        }

		// Normalize spin length
		double S_magnitude = std::sqrt(y_pred_storage[atom][0] * y_pred_storage[atom][0] + y_pred_storage[atom][1] * y_pred_storage[atom][1] + y_pred_storage[atom][2] * y_pred_storage[atom][2]);

        y_pred_storage[atom][0] /= S_magnitude;
        y_pred_storage[atom][1] /= S_magnitude;
        y_pred_storage[atom][2] /= S_magnitude;

        atoms::x_spin_array[atom] = y_pred_storage[atom][0];
        atoms::y_spin_array[atom] = y_pred_storage[atom][1];
        atoms::z_spin_array[atom] = y_pred_storage[atom][2];
        atoms::x_v_array[atom] = y_pred_storage[atom][3];
        atoms::y_v_array[atom] = y_pred_storage[atom][4];
        atoms::z_v_array[atom] = y_pred_storage[atom][5];
        atoms::x_w_array[atom] = y_pred_storage[atom][6];
        atoms::y_w_array[atom] = y_pred_storage[atom][7];
        atoms::z_w_array[atom] = y_pred_storage[atom][8];
    }


	// Update fields for k3
    calculate_spin_fields(0, num_atoms);

    // K3 Step
    for (int atom = 0; atom < num_atoms; ++atom) {
        idx_x = atoms::atom_idx_x[atom];
		idx_y = atoms::atom_idx_y[atom];
		idx_z = atoms::atom_idx_z[atom];

        H[0] = atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom] + (atoms::noise_field[idx_x][atoms::noise_index] + atoms::noise_field[idx_x][atoms::noise_index + 1]) / 2.0;
        H[1] = atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom] + (atoms::noise_field[idx_y][atoms::noise_index] + atoms::noise_field[idx_y][atoms::noise_index + 1]) / 2.0;
        H[2] = atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom] + (atoms::noise_field[idx_z][atoms::noise_index] + atoms::noise_field[idx_z][atoms::noise_index + 1]) / 2.0;

        spinDynamics(y_pred_storage[atom], H, k3_storage[atom]);

        for (size_t i = 0; i < 9; ++i) {
            y_pred_storage[atom][i] = y_in_storage[atom][i] + mp::dt * k3_storage[atom][i];
        }


		// Normalize spin length
		double S_magnitude = std::sqrt(y_pred_storage[atom][0] * y_pred_storage[atom][0] + y_pred_storage[atom][1] * y_pred_storage[atom][1] + y_pred_storage[atom][2] * y_pred_storage[atom][2]);

        y_pred_storage[atom][0] /= S_magnitude;
        y_pred_storage[atom][1] /= S_magnitude;
        y_pred_storage[atom][2] /= S_magnitude;

        atoms::x_spin_array[atom] = y_pred_storage[atom][0];
        atoms::y_spin_array[atom] = y_pred_storage[atom][1];
        atoms::z_spin_array[atom] = y_pred_storage[atom][2];
        atoms::x_v_array[atom] = y_pred_storage[atom][3];
        atoms::y_v_array[atom] = y_pred_storage[atom][4];
        atoms::z_v_array[atom] = y_pred_storage[atom][5];
        atoms::x_w_array[atom] = y_pred_storage[atom][6];
        atoms::y_w_array[atom] = y_pred_storage[atom][7];
        atoms::z_w_array[atom] = y_pred_storage[atom][8];
    }


	// Update fields for k4
    calculate_spin_fields(0, num_atoms);
    
    //K4 step
    for (int atom = 0; atom < num_atoms; ++atom) {
        idx_x = atoms::atom_idx_x[atom];
		idx_y = atoms::atom_idx_y[atom];
		idx_z = atoms::atom_idx_z[atom];

        H[0] = atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom] + atoms::noise_field[idx_x][atoms::noise_index+1];
        H[1] = atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom] + atoms::noise_field[idx_y][atoms::noise_index+1];
        H[2] = atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom] + atoms::noise_field[idx_z][atoms::noise_index+1];

        spinDynamics(y_pred_storage[atom], H, k4_storage[atom]);

        for (size_t i = 0; i < 9; ++i) {
            y_pred_storage[atom][i] = y_in_storage[atom][i] + (mp::dt / 6.0) * (k1_storage[atom][i] + 2.0 * k2_storage[atom][i] + 2.0 * k3_storage[atom][i] + k4_storage[atom][i]);
        }

        double S_magnitude = std::sqrt(y_pred_storage[atom][0] * y_pred_storage[atom][0] + y_pred_storage[atom][1] * y_pred_storage[atom][1] + y_pred_storage[atom][2] * y_pred_storage[atom][2]);

        y_pred_storage[atom][0] /= S_magnitude;
        y_pred_storage[atom][1] /= S_magnitude;
        y_pred_storage[atom][2] /= S_magnitude;

        atoms::x_spin_array[atom] = y_pred_storage[atom][0];
        atoms::y_spin_array[atom] = y_pred_storage[atom][1];
        atoms::z_spin_array[atom] = y_pred_storage[atom][2];
        atoms::x_v_array[atom] = y_pred_storage[atom][3];
        atoms::y_v_array[atom] = y_pred_storage[atom][4];
        atoms::z_v_array[atom] = y_pred_storage[atom][5];
        atoms::x_w_array[atom] = y_pred_storage[atom][6];
        atoms::y_w_array[atom] = y_pred_storage[atom][7];
        atoms::z_w_array[atom] = y_pred_storage[atom][8];
    }

	// Increment noise index
	atoms::noise_index += 1;

	return;
}


} // end of intental namespace

} // end of sim namespace
