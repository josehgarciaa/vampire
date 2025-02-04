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
#include "random.hpp"

namespace LLGQ_arrays{

    //---------------------------------------------------------------------------
    // Local arrays and variables
    //---------------------------------------------------------------------------

    // Auxiallary variables
    std::vector <double> x_w_array;
	std::vector <double> y_w_array;
	std::vector <double> z_w_array;
	std::vector <double> x_v_array;
	std::vector <double> y_v_array;
	std::vector <double> z_v_array;

    // Storage arrays for RK4
    std::vector<std::vector<double>> k1_storage;
    std::vector<std::vector<double>> k2_storage;
    std::vector<std::vector<double>> k3_storage;
    std::vector<std::vector<double>> k4_storage;
    std::vector<std::vector<double>> y_pred_storage;
    std::vector<std::vector<double>> y_in_storage;

    // Arrays for noise generation
	std::vector<double> sqrt_PSD_buffer;
    std::vector<std::vector<double>> noise_field;
	double noise_index;

    // Indices for random fields
	std::vector<double> atom_idx_x;
    std::vector<double> atom_idx_y;
	std::vector<double> atom_idx_z;


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

//---------------------------------------------------------------------------
// Class for FFTW plan management
//---------------------------------------------------------------------------
void calculate_random_fields(int realizations, int n, double dt, double T);
void assign_unique_indices(int realizations);
void precompute_sqrt_PSD(int n, double dt, double T);
double PSD(const double& omega, const double& T);



int LLGQinit(){
	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "sim::LLG_init has been called" << std::endl;}

	using namespace LLGQ_arrays;


    x_w_array.resize(atoms::num_atoms, 0.0);
    y_w_array.resize(atoms::num_atoms, 0.0);
    z_w_array.resize(atoms::num_atoms, 0.0);
    x_v_array.resize(atoms::num_atoms, 0.0);
    y_v_array.resize(atoms::num_atoms, 0.0);
    z_v_array.resize(atoms::num_atoms, 0.0);

    k1_storage.resize(atoms::num_atoms, std::vector<double>(9));
    k2_storage.resize(atoms::num_atoms, std::vector<double>(9));
    k3_storage.resize(atoms::num_atoms, std::vector<double>(9));
    k4_storage.resize(atoms::num_atoms, std::vector<double>(9));
    y_pred_storage.resize(atoms::num_atoms, std::vector<double>(9));
    y_in_storage.resize(atoms::num_atoms, std::vector<double>(9));


	// Set number of realizations (full field, for final release allow even smaller number of realizations)
	const int num_atoms = atoms::num_atoms;
    int realizations = num_atoms * 3 + 4;
	LLGQ_arrays::noise_index = 0;

    // Disable external thermal field calculations
    sim::hamiltonian_simulation_flags[3] = 0;

	// Assign unique indices for random fields
	assign_unique_indices(realizations);

	// Precompute square root of PSD
	const int n = static_cast<int>(sim::equilibration_time) + 1;
	precompute_sqrt_PSD(n, mp::dt, sim::temperature);

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

inline void spinDynamics(const double* y,const double* H,double* dydt);

void llg_quantum_step(){

	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "sim::mLLG has been called" << std::endl;}

	using namespace LLGQ_arrays;

	// Check for initialisation of LLG integration arrays
	if(LLG_set==false) sim::LLGQinit();


	// Local variables
	std::vector<double> H(3);
    const int num_atoms = atoms::num_atoms;


    

	// Store initial state
    for (int atom = 0; atom < num_atoms; atom++) {
        y_in_storage[atom][0] = atoms::x_spin_array[atom];
        y_in_storage[atom][1] = atoms::y_spin_array[atom];
        y_in_storage[atom][2] = atoms::z_spin_array[atom];
        y_in_storage[atom][3] = LLGQ_arrays::x_v_array[atom];
        y_in_storage[atom][4] = LLGQ_arrays::y_v_array[atom];
        y_in_storage[atom][5] = LLGQ_arrays::z_v_array[atom];
        y_in_storage[atom][6] = LLGQ_arrays::x_w_array[atom];
        y_in_storage[atom][7] = LLGQ_arrays::y_w_array[atom];
        y_in_storage[atom][8] = LLGQ_arrays::z_w_array[atom];
    }


	// Calculate fields
	calculate_spin_fields(0, num_atoms);
    calculate_external_fields(0, num_atoms);

    const double dt = mp::dt;
    const double half_dt = 0.5 * mp::dt;
    const double dt_over_6 = mp::dt / 6.0;

    // K1 Step
    for (int atom = 0; atom < num_atoms; ++atom) {
        // Update fields for k1
        H[0] = atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom] + LLGQ_arrays::noise_field[atom_idx_x[atom]][LLGQ_arrays::noise_index];
        H[1] = atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom] + LLGQ_arrays::noise_field[atom_idx_y[atom]][LLGQ_arrays::noise_index];
        H[2] = atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom] + LLGQ_arrays::noise_field[atom_idx_z[atom]][LLGQ_arrays::noise_index];

        

        // Calculate K1
        spinDynamics(y_in_storage[atom].data(), H.data(), k1_storage[atom].data());

        // Calculate y_pred for k2
       
        
        for (size_t i = 0; i < 9; ++i) {
            y_pred_storage[atom][i] = y_in_storage[atom][i] + half_dt * k1_storage[atom][i];
        }

        // Normalize spin length
        double S_magnitude = std::sqrt(y_pred_storage[atom][0] * y_pred_storage[atom][0] + y_pred_storage[atom][1] * y_pred_storage[atom][1] + y_pred_storage[atom][2] * y_pred_storage[atom][2]);
        const double inv_mag = 1.0 / S_magnitude;
        y_pred_storage[atom][0] *= inv_mag;
        y_pred_storage[atom][1] *= inv_mag;
        y_pred_storage[atom][2] *= inv_mag;

        // Update spin for field calculation
        atoms::x_spin_array[atom] = y_pred_storage[atom][0];
        atoms::y_spin_array[atom] = y_pred_storage[atom][1];
        atoms::z_spin_array[atom] = y_pred_storage[atom][2];
    }


	// Update spin fields for k2
    calculate_spin_fields(0, num_atoms);

    // K2 Step
    for (int atom = 0; atom < num_atoms; ++atom) {
        // Update fields for k2
        H[0] = atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom] + (LLGQ_arrays::noise_field[atom_idx_x[atom]][LLGQ_arrays::noise_index] + LLGQ_arrays::noise_field[atom_idx_x[atom]][LLGQ_arrays::noise_index + 1]) / 2.0;
        H[1] = atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom] + (LLGQ_arrays::noise_field[atom_idx_y[atom]][LLGQ_arrays::noise_index] + LLGQ_arrays::noise_field[atom_idx_y[atom]][LLGQ_arrays::noise_index + 1]) / 2.0;
        H[2] = atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom] + (LLGQ_arrays::noise_field[atom_idx_z[atom]][LLGQ_arrays::noise_index] + LLGQ_arrays::noise_field[atom_idx_z[atom]][LLGQ_arrays::noise_index + 1]) / 2.0;

        
        // Calculate K2
        spinDynamics(y_pred_storage[atom].data(), H.data(), k2_storage[atom].data());

        // Calculate y_pred for k3
        for (size_t i = 0; i < 9; ++i) {
            y_pred_storage[atom][i] = y_in_storage[atom][i] + half_dt * k2_storage[atom][i];
        }

		// Normalize spin length
		double S_magnitude = std::sqrt(y_pred_storage[atom][0] * y_pred_storage[atom][0] + y_pred_storage[atom][1] * y_pred_storage[atom][1] + y_pred_storage[atom][2] * y_pred_storage[atom][2]);
        const double inv_mag = 1.0 / S_magnitude;
        y_pred_storage[atom][0] *= inv_mag;
        y_pred_storage[atom][1] *= inv_mag;
        y_pred_storage[atom][2] *= inv_mag;

        // Update spin for field calculation
        atoms::x_spin_array[atom] = y_pred_storage[atom][0];
        atoms::y_spin_array[atom] = y_pred_storage[atom][1];
        atoms::z_spin_array[atom] = y_pred_storage[atom][2];
    }


	// Update spin fields for k3
    calculate_spin_fields(0, num_atoms);

    // K3 Step
    for (int atom = 0; atom < num_atoms; ++atom) {
        // Update fields for k2
        H[0] = atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom] + (LLGQ_arrays::noise_field[atom_idx_x[atom]][LLGQ_arrays::noise_index] + LLGQ_arrays::noise_field[atom_idx_x[atom]][LLGQ_arrays::noise_index + 1]) / 2.0;
        H[1] = atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom] + (LLGQ_arrays::noise_field[atom_idx_y[atom]][LLGQ_arrays::noise_index] + LLGQ_arrays::noise_field[atom_idx_y[atom]][LLGQ_arrays::noise_index + 1]) / 2.0;
        H[2] = atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom] + (LLGQ_arrays::noise_field[atom_idx_z[atom]][LLGQ_arrays::noise_index] + LLGQ_arrays::noise_field[atom_idx_z[atom]][LLGQ_arrays::noise_index + 1]) / 2.0;

        // Calculate K3
        spinDynamics(y_pred_storage[atom].data(), H.data(), k3_storage[atom].data());

        // Calculate y_pred for k4
        for (size_t i = 0; i < 9; ++i) {
            y_pred_storage[atom][i] = y_in_storage[atom][i] + dt * k3_storage[atom][i];
        }

		// Normalize spin length
		double S_magnitude = std::sqrt(y_pred_storage[atom][0] * y_pred_storage[atom][0] + y_pred_storage[atom][1] * y_pred_storage[atom][1] + y_pred_storage[atom][2] * y_pred_storage[atom][2]);
        const double inv_mag = 1.0 / S_magnitude;
        y_pred_storage[atom][0] *= inv_mag;
        y_pred_storage[atom][1] *= inv_mag;
        y_pred_storage[atom][2] *= inv_mag;

        // Update spin for field calculation
        atoms::x_spin_array[atom] = y_pred_storage[atom][0];
        atoms::y_spin_array[atom] = y_pred_storage[atom][1];
        atoms::z_spin_array[atom] = y_pred_storage[atom][2];
    }


	// Update fields for k4
    calculate_spin_fields(0, num_atoms);

    //K4 step
    for (int atom = 0; atom < num_atoms; ++atom) {
        // Update fields for k4
        H[0] = atoms::x_total_spin_field_array[atom] + atoms::x_total_external_field_array[atom] + LLGQ_arrays::noise_field[atom_idx_x[atom]][LLGQ_arrays::noise_index+1];
        H[1] = atoms::y_total_spin_field_array[atom] + atoms::y_total_external_field_array[atom] + LLGQ_arrays::noise_field[atom_idx_y[atom]][LLGQ_arrays::noise_index+1];
        H[2] = atoms::z_total_spin_field_array[atom] + atoms::z_total_external_field_array[atom] + LLGQ_arrays::noise_field[atom_idx_z[atom]][LLGQ_arrays::noise_index+1];

        // Calculate K4
        spinDynamics(y_pred_storage[atom].data(), H.data(), k4_storage[atom].data());
        
        // Final update and normalization for spin components
        for (size_t i = 0; i < 3; ++i) {
            y_pred_storage[atom][i] = y_in_storage[atom][i] + dt_over_6 * (k1_storage[atom][i] + 2.0 * k2_storage[atom][i] + 2.0 * k3_storage[atom][i] + k4_storage[atom][i]);
        }

        double S_magnitude = std::sqrt(y_pred_storage[atom][0] * y_pred_storage[atom][0] + y_pred_storage[atom][1] * y_pred_storage[atom][1] + y_pred_storage[atom][2] * y_pred_storage[atom][2]);
        const double inv_mag = 1.0 / S_magnitude;
        y_pred_storage[atom][0] *= inv_mag;
        y_pred_storage[atom][1] *= inv_mag;
        y_pred_storage[atom][2] *= inv_mag;

        atoms::x_spin_array[atom] = y_pred_storage[atom][0];
        atoms::y_spin_array[atom] = y_pred_storage[atom][1];
        atoms::z_spin_array[atom] = y_pred_storage[atom][2];

        // Update auxiliary variables
        LLGQ_arrays::x_v_array[atom] = y_in_storage[atom][3] + (dt_over_6) * (k1_storage[atom][3] + 2.0 * k2_storage[atom][3] + 2.0 * k3_storage[atom][3] + k4_storage[atom][3]);
        LLGQ_arrays::y_v_array[atom] = y_in_storage[atom][4] + (dt_over_6) * (k1_storage[atom][4] + 2.0 * k2_storage[atom][4] + 2.0 * k3_storage[atom][4] + k4_storage[atom][4]);
        LLGQ_arrays::z_v_array[atom] = y_in_storage[atom][5] + (dt_over_6) * (k1_storage[atom][5] + 2.0 * k2_storage[atom][5] + 2.0 * k3_storage[atom][5] + k4_storage[atom][5]);
        LLGQ_arrays::x_w_array[atom] = y_in_storage[atom][6] + (dt_over_6) * (k1_storage[atom][6] + 2.0 * k2_storage[atom][6] + 2.0 * k3_storage[atom][6] + k4_storage[atom][6]);
        LLGQ_arrays::y_w_array[atom] = y_in_storage[atom][7] + (dt_over_6) * (k1_storage[atom][7] + 2.0 * k2_storage[atom][7] + 2.0 * k3_storage[atom][7] + k4_storage[atom][7]);
        LLGQ_arrays::z_w_array[atom] = y_in_storage[atom][8] + (dt_over_6) * (k1_storage[atom][8] + 2.0 * k2_storage[atom][8] + 2.0 * k3_storage[atom][8] + k4_storage[atom][8]);
    }

	// Increment noise index
	LLGQ_arrays::noise_index += 1;

	return;
}

inline void spinDynamics(const double* y, const double* H, double* dydt) {
    const double A = mp::material[0].A;
    const double Gamma = mp::material[0].Gamma;
    const double omega0 = mp::material[0].omega0;

    // dS/dt = S × (H + v)
    dydt[0] =  (y[1]*(H[2]+y[5]) - y[2]*(H[1]+y[4]));
    dydt[1] =  (y[2]*(H[0]+y[3]) - y[0]*(H[2]+y[5]));
    dydt[2] =  (y[0]*(H[1]+y[4]) - y[1]*(H[0]+y[3]));
    
    
    // dv/dt = w
    dydt[3] = y[6];
    dydt[4] = y[7];
    dydt[5] = y[8];
    
    // dw/dt = -ω₀²v - Γw + AS
    dydt[6] = -omega0*omega0*y[3] - Gamma*y[6] + A*y[0];
    dydt[7] = -omega0*omega0*y[4] - Gamma*y[7] + A*y[1];
    dydt[8] = -omega0*omega0*y[5] - Gamma*y[8] + A*y[2];
}


} // end of internal namespace

void assign_unique_indices(int realizations) {
    const int num_atoms = atoms::num_atoms;

    // Ensure we have enough unique indices for all atoms
    if (realizations < num_atoms * 3) {
        throw std::runtime_error("Not enough realizations for unique indices");
    }

    std::vector<int> all_indices(realizations);
    std::iota(all_indices.begin(), all_indices.end(), 0); // Fill with 0, 1, ..., realizations-1

    LLGQ_arrays::atom_idx_x.resize(num_atoms);
    LLGQ_arrays::atom_idx_y.resize(num_atoms);
    LLGQ_arrays::atom_idx_z.resize(num_atoms);

    for (int atom = 0; atom < num_atoms; atom++) {
        int base_index = atom * 3;
        LLGQ_arrays::atom_idx_x[atom] = all_indices[base_index];
        LLGQ_arrays::atom_idx_y[atom] = all_indices[base_index + 1];
        LLGQ_arrays::atom_idx_z[atom] = all_indices[base_index + 2];
    }
}


void calculate_random_fields(int realizations, int n, double dt, double T) {
    // Pre-allocate result vectors
    LLGQ_arrays::noise_field.clear();
    const double S0 = mp::material[0].S0;
    const double inv_sqrt_S0 = 1.0 / std::sqrt(S0);
    
    
    // Pre-allocate all buffers
    double* __restrict in = (double*)fftw_malloc(sizeof(double) * n);
    fftw_complex* __restrict out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (n/2 + 1));
    double* __restrict result = (double*)fftw_malloc(sizeof(double) * n);
    
    // Create FFTW plans
    fftw_plan forward = fftw_plan_dft_r2c_1d(n, in, out, FFTW_MEASURE);
    fftw_plan backward = fftw_plan_dft_c2r_1d(n, out, result, FFTW_MEASURE);
    
    // Pre-calculate constants
    const double df = 1.0 / (n * dt);
    const double two_pi = 2.0 * M_PI;
    const double norm_factor = 1.0 / n;
    
    // Setup random number generation
    static thread_local std::random_device rd;
    static thread_local std::mt19937 gen(rd());
    std::normal_distribution<> dist(0.0, 1.0 / std::sqrt(dt));
    
    LLGQ_arrays::noise_field.resize(realizations, std::vector<double>(n));


    // Progress bar setup
    const int bar_width = 50;
    int last_printed_percent = -1;


    // Precompute square root of PSD
    if (LLGQ_arrays::sqrt_PSD_buffer.empty()) {
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
            LLGQ_arrays::noise_field[r][j] = result[j] * norm_factor * inv_sqrt_S0;
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
    LLGQ_arrays::sqrt_PSD_buffer.resize(n/2 + 1);
    double df = 1.0 / (n * dt);
    for (int i = 0; i <= n/2; ++i) {
        double omega = 2.0 * M_PI * i * df;
        LLGQ_arrays::sqrt_PSD_buffer[i] = std::sqrt(PSD(omega, T));
    }
}

} // end of sim namespace
