# GPU Modernization Roadmap

This document outlines the prioritized files and modules in the Vampire codebase that require modernization for GPU (CUDA) support. The goal is to eliminate local (hdr) dependencies and ensure all GPU code is linked directly to official CUDA libraries.

---

## 1. Highest Priority: Core CUDA Implementation

These files contain the main CUDA logic and should be updated first. All local includes (e.g., `#include "cuda.hpp"`, `#include "cusplibrary-0.5.1/..."`, etc.) must be removed and replaced with official CUDA library calls.

- `src/cuda/llg_heun.cu` / `llg_heun.hpp`
- `src/cuda/statistics.cu` / `statistics.hpp`
- `src/cuda/exchange_fields.cu` / `exchange_fields.hpp`
- `src/cuda/dipole.cu`
- `src/cuda/initialize.cu`
- `src/cuda/finalize.cu`
- `src/cuda/data.cu` / `data.hpp`
- `src/cuda/transfer.cu`
- `src/cuda/internal.cu` / `internal.hpp`
- `src/cuda/external_fields.cu`
- `src/cuda/config.cu`
- `src/cuda/cuda_utils.hpp`
- `src/cuda/typedefs.hpp`
- `src/cuda/cuda_timer.h`

---

## 2. Medium Priority: GPU Interface and Utility

These files manage the interface between the main code and the GPU modules. They should be updated to remove any local header dependencies and ensure all GPU calls are routed through official CUDA APIs.

- `src/gpu/interface.cpp`
- `src/gpu/initialize.cpp`
- `src/gpu/finalize.cpp`
- `src/gpu/data.cpp`
- `src/gpu/config.cpp`
- `src/gpu/llg_heun.cpp`
- `src/gpu/statistics.cpp`
- `src/gpu/transfer.cpp`

---

## 3. Lower Priority: Build System and Configuration

Update the build system to remove references to local CUDA headers and libraries, and ensure proper linking to CUDA libraries.

- `src/cuda/makefile`
- `src/gpu/makefile`
- Top-level `makefile` (or migrate to CMake for better CUDA support)

---

## 4. Remove Local/Deprecated Libraries

The following local libraries and headers should be removed from the project:

- `hdr/cuda.hpp`
- `hdr/cusplibrary-0.5.1/` (entire directory)
- `hdr/cub-1.5.2/` (entire directory)
- Any other local CUDA-related headers in `hdr/`

---

## 5. General Guidelines for Modernization

- Replace all local CUDA utility and typedef headers with official CUDA, cuBLAS, cuRAND, and cuSPARSE headers.
- Use modern CUDA memory management (e.g., `cudaMallocManaged`, streams, events).
- Update kernel launches to use modern syntax and error checking.
- Remove all deprecated or custom random number generation in favor of cuRAND.
- Use C++17 or newer features where possible for host code.

---

## 6. Next Steps

1. Refactor and test each file/module in the order above.
2. Ensure all GPU code compiles and links only against official CUDA libraries.
3. Remove all local CUDA/CUSP/CUB headers and libraries.
4. Update documentation to reflect new dependencies and build instructions.

---

**Note:**  
This roadmap focuses on code modernization and dependency cleanup only. The workflow and scientific logic should remain unchanged during this phase.
