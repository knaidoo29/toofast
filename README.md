# TwoFast

Two and three point correlation function estimator in c++ with OpenMPI parallelisation.

## Usage

Don't use this yet... I'm still testing it and things may not work as they are supposed to.
Full disclosure: this was built as an educational exercise to give myself an excuse to write
in c/c++ and learn how to parallelise using OpenMPI.

## Installing

Edit the Makefile in the top directory and then run `make`. This creates a library
for two point correlation function estimation. In the script folder edit the Makefile
and again run `make`. The scripts take in a parameter file of which examples exist
in the script folder.

## Running

Running the scripts within the script directory.

Non-parallel version:

``./script_two_point_3d.o paramfile_two_point_3d_5000.ini``

Parallel version:

``/path/to/mpirun -n 4 ./script_two_point_3d_omp.o paramfile_two_point_3d_5000.ini``

## What can it do?

Brute force pair distance binning for auto and cross correlation used to calculate
the two point correlation function in 3D (need to add 2D, tomographic and spherical polar).

## To do list

Brute force method:
- Three point correlation function using Szapudi-Szalay estimator.

Fast method:
- KD-Tree method for small scales.
- Pixelised method for large scales.
