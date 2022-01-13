#pragma once

#include <array>
#include <complex>
#include <cstddef>
#include <math.h>
#include <boost/multi_array.hpp>

namespace HM {

//General parameters.
//TODO: Don't compile if nx,ny not even...
constexpr size_t nt = 200000;             //Num of timesteps.
constexpr size_t nx = 256;                //Num of x-steps.
constexpr size_t ny = 256;                //Num of y-steps.
constexpr size_t saveRate  = 250;         //How often to save data.
constexpr size_t numFrames = nt/saveRate; //Total number of saves.
constexpr double lx = 2*M_PI/.15;         //Box size in x-dim.
constexpr double ly = 2*M_PI/.15;         //Box size in y-dim.
constexpr double dx = lx/nx;              //Stepsize in x-dim.
constexpr double dy = ly/ny;              //Stepsize in y-dim.
constexpr double dt = 2e-1;               //Size of timestep.

//Params specific to model (Haotian's ETG).
constexpr double tau      = 0;   //(T_e/T_i)
constexpr double eta      = 0;   //(r_n/r_T) - where (grad_x(A_e)/A_e) = -1/r_A
constexpr double mRat     = 0;   //(m_e/m_i)
constexpr double kappa    = 0.1; //HasegawaMima eq. const.
constexpr double rnByRhoI = 0;   //(r_n/rho_i)

//Extra flags.
constexpr bool loadGENE = false; //Possibility for loading data from GENE ASCII output.


//Arrays
std::array<double, nx> xGrid     = {0}; //Ranges from [0, (nx-1)*dx] w/ nx vals.
std::array<double, ny> yGrid     = {0};
std::array<double, nx> kxGrid    = {0}; //Ranges from [-pi*nx/lx, pi*nx/lx) w/ nx vals.
std::array<double, ny> kyGrid    = {0};
std::array<double, nx> kxAliased = {0};
std::array<double, ny> kyAliased = {0};

double kConst[ny][nx]             = {0};
std::complex<double> phi[ny][nx]  = {0};
std::complex<double> phik[ny][nx] = {0};

std::array<std::complex<double>[nx][ny], numFrames> phit  = {0};
std::array<std::complex<double>[nx][ny], numFrames> phikt = {0};

} //Close namespace.