## HETP: HETerogeneous-vectorized-or-Parallel 

[![Ubuntu](https://github.com/sjmiller204/HETerogeneous-vectorized-or-Parallel/actions/workflows/ubuntu.yml/badge.svg)](https://github.com/sjmiller204/HETerogeneous-vectorized-or-Parallel/actions/workflows/ubuntu.yml)

HETP is an aerosol thermodynamic equilibrium solver written in modern Fortran based on ISORROPIA II (which is written in FORTRAN 77).  HETP solves only the 'forward' metastable state of the NH4+/Na+/Ca2+/K+/Mg2+/SO42–/NO3–/Cl–/H2O system.  The main publication for HETP is avaiable online at:  https://gmd.copernicus.org/preprints/gmd-2023-159/.

This repository contains the HETP case-by-case implementation, with a simple interface to call HETP for a single set of input conditions.  The file 'hetp_main.F90' is a Fortran script to interact with HETP, the file 'hetp_mod.F90' holds the HETP code and the file 'mach_hetp_mod.F90' holds parameters used in HETP.  Input required are the total gas + aerosol concentrations of eight precursor species with units of mol/m^3 air (i.e., sulfate, ammonium, nitrate, sodium, chloride, calcium, potassium, magnesium), the relative humidity (on a 0-1 scale) and the air temperature (units of K).

## How to build and run HETP test

A HETP standalone test can be built using CMake and a fortran compiler. It runs the program in src/Test/hetp_main.F90. After downloading the model navigate to the main directory and execute the following commands. If your build is successful an executable will be placed in the build/bin directory.

```
cd /path/to/hetp
mkdir build
cd build
cmake ..
make -j
```

If you wish to build HETP with compiler debug flags on simply run the following command in your build folder, and then rebuild and rerun after copying your new executable to your run directory.

```
cmake . -DCMAKE_BUILD_TYPE=Debug
```

To run the simple HETP test do the following from a directory that contains the executable. Inputs and outputs will be printed to the terminal screen.

```
./hetp_test
```

## How to use HETP in an external model

How to connect HETP to an external model is dependent on the external model's build system. In general, all files needed are stored within the src/Core directory. The simplest way to connect HETP is therefore to add files in that directory to the parent model's build list. If the external model uses CMake you can also connect HETP such that its CMake files are used. To do this you must set parameter HETP_EXTERNAL_CONFIG to TRUE within a CMake file in the parent model. This will bypass compiling the src/Test directory files and creating a HETP standalone executable.
