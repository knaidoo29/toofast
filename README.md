# TwoFast

Two point correlation function estimator in c++ with OpenMPI parallelisation and python
management class for generating input and parameter files, and then subsequently
running the routines.

## Usage

Don't use this yet... I'm still testing it and things may not work as they are supposed to.

## Installing

Edit the Makefile in the top directory and then run `make`. This creates a library
for two point correlation function estimation. In the script folder edit the Makefile
and again run `make`. The scripts take in a parameter file of which examples exist
in the script folder.

## Running

Run the scripts from the main directory with a paramfile (examples in the paramfile folder).

Non-parallel version:

``./TWOPOINT paramfile.ini``

Parallel version:

``/path/to/mpirun -n 4 ./TWOPOINT_MPI paramfile.ini``

## What can it do?

Brute force pair distance binning for auto and cross correlation used to calculate
the two point correlation function in 2D/3D and tomographic coordinates.

## To do list

Brute force method:
- Spherical polar coordinates, with binning for parallel and perpendicular to the line-of-sight.
- Three point correlation function using Szapudi-Szalay estimator.

Fast method:
- KD-Tree method for small scales.
- Pixelised method for large scales.
