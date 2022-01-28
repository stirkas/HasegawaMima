#pragma once

#include <array>
#include <complex>
#include <cstddef>
#include <math.h>
#include <boost/multi_array.hpp>

namespace HM {

//General parameters.
size_t nt = 20000;              //Num of timesteps.
size_t nx = 256;                //Num of x-steps.
size_t ny = 256;                //Num of y-steps.
size_t saveRate  = 50;          //How often to save data.
size_t numFrames = nt/saveRate; //Total number of saves.
double lx = 2*M_PI/.15;         //Box size in x-dim.
double ly = 2*M_PI/.15;         //Box size in y-dim.
double dx = lx/nx;              //Stepsize in x-dim.
double dy = ly/ny;              //Stepsize in y-dim.
double dt = 2e-1;               //Size of timestep.

//Params specific to model (Haotian's ETG).
double tau      = 0;   //(T_e/T_i)
double eta      = 0;   //(r_n/r_T) - where (grad_x(A_e)/A_e) = -1/r_A
double mRat     = 0;   //(m_e/m_i)
double kappa    = 0.1; //HasegawaMima eq. const.
double rnByRhoI = 0;   //(r_n/rho_i)

//Extra flags.
bool loadGENE = false; //Possibility for loading data from GENE ASCII output.


//Arrays
typedef boost::multi_array<double, 1> Arr1d;
typedef Arr1d::extent_gen extents1d;
Arr1d xGrid;
Arr1d yGrid;
Arr1d kxGrid;
Arr1d kyGrid;
Arr1d kxAliased;
Arr1d kyAliased;
extents1d ext1d;

typedef boost::multi_array<double, 2> Arr2d;
typedef Arr2d::extent_gen extents2d;
Arr2d kConst;
extents2d ext2d;

typedef boost::multi_array<std::complex<double>, 2> Arrc2d;
typedef Arrc2d::extent_gen extentsc2d;
Arrc2d phi;
Arrc2d phik;
extentsc2d extc2d;

typedef boost::multi_array<double, 3> Arr3d;
typedef Arr3d::extent_gen extents3d;
Arr3d phit;
Arr3d phikt;
extents3d extc3d;

} //Close namespace.