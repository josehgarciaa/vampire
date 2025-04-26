Alucard
==============

A fork of Vampire (by Richard Evans)

Alucard is a high performance general purpose code for the atomistic simulation of magnetic materials. It is a fork of the Vampire project, inheriting its scientific capabilities and codebase, and is maintained in compliance with the original GNU GPL v2 license.

By Richard Evans (original Vampire author)

Capabilities
---------------
Vampire is designed to be highly flexible to deal with a wide variety of problems using a diverse set of simulation tools and methods. The capabilities of the code can be summarised broadly in terms of the simulation methods, standard problems, structural properties and features of the code, all of which can be combined to tackle almost any problem.

**Simulation methods**
-Stochastic Landau-Lifshitz-Gilbert equation (spin dynamics)
-Monte Carlo metropolis
-Constrained Monte Carlo metropolis

**Standard calculations**
-Ultrafast spin dynamics
-Hysteresis loops
-Curie temperature
-Temperature dependent anisotropy
-Temperature dependent energy barriers
-Field cooling
-Heat assisted and conventional magnetic recording
-Laser induced spin dynamics

**Structural properties**
-Bulk-like systems
-Thin films
-Nanoparticles - spheres, cubes, truncated octahedra, cylinders
-Voronoi granular structures
-Nanoparticle arrays
-Core-shell nanoparticles
-Multilayer thin films
-Interface roughness and intermixing
-Dilute magnetic systems
-Lithographically defined geometries
-SC, FCC, HCP, and BCC crystal structures
-User-defined atomic structures - for example from Molecular Dynamics simulations

**Magnetic properties**
-Ferromagnets
-Antiferromagnets
-Ferrimagnets
-Spin glass
-Single-ion, 2-ion and cubic anisotropies
-Scalar, vector and tensor forms of exchange including the DM interaction
-User-defined Hamiltonian from ab-initio Density Functional Theory (DFT) calculations
-Demagnetisation fields (macrocell approximation)

**Code features**
-Modular object-oriented C++
-Simple to use textfile input
-High performance code
-Parallelisation using the MPI library
-Variety of geometric decomposition algorithms
-Usable on a laptop to a supercomputer with thousands of cores
-Output to PoVRAY for visualisation and publication quality graphics
-Output to rasmol/jmol for structural inspection
-Minimal dependence on external libraries for portability
-Freely available open source code

Installation
---------------

**Dependencies:**
- CUDA Toolkit (>= 11.0)
- cuSPARSE, cuBLAS (provided by CUDA Toolkit)
- Thrust (provided by CUDA Toolkit)
- C++17 compatible compiler
- GNU Make

**Build and Install:**
```sh
make cuda
sudo make install
```

**Testing:**
After installation, you can run the basic tests:
```sh
cd tests/basic
./run_basic_tests.sh
```

This will run a set of three simple smoke tests to verify your installation.

License
---------------
This project is a fork of Vampire and is distributed under the GNU General Public License v2 (GPL v2), in compliance with the original Vampire license. See the license file for details.

