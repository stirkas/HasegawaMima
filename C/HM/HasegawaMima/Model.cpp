#include "Model.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <iomanip> //TODO: Remove
#include <iostream>

#include <fftw3.h>

#include "Definitions.hpp"
namespace HM {

//Useful things for the actual advancing logic.
using namespace std::complex_literals;
fftw_plan phiForeFFT;
fftw_plan phiBackFFT;
fftw_plan phikxFFT;
fftw_plan phikyFFT;
fftw_plan zetakxFFT;
fftw_plan zetakyFFT;
fftw_plan advTotFFT;
double kxSqGrid[nx] = {0};
double kySqGrid[ny] = {0};
std::complex<double> phix[ny][nx]      = {0};
std::complex<double> phiy[ny][nx]      = {0};
std::complex<double> phikx[ny][nx]     = {0};
std::complex<double> phiky[ny][nx]     = {0};
std::complex<double> zetak[ny][nx]     = {0};
std::complex<double> zetax[ny][nx]     = {0};
std::complex<double> zetay[ny][nx]     = {0};
std::complex<double> zetakx[ny][nx]    = {0};
std::complex<double> zetaky[ny][nx]    = {0};
std::complex<double> phikyTrue[ny][nx] = {0};
std::complex<double> advTot[ny][nx]    = {0};
std::complex<double> advkTot[ny][nx]   = {0};
void Advance(std::complex<double> (&phi)[ny][nx]);

void Initialize()
{
   //Generate FFT plans now because speedy FFTW_MEASURE option involves clearing data...
   phiForeFFT = fftw_plan_dft_2d(nx, ny, reinterpret_cast<fftw_complex*>(&phi[0][0]),   reinterpret_cast<fftw_complex*>(&phik[0][0]),
                                 FFTW_FORWARD,  FFTW_MEASURE);
   phiBackFFT = fftw_plan_dft_2d(nx, ny, reinterpret_cast<fftw_complex*>(&phik[0][0]),  reinterpret_cast<fftw_complex*>(&phi[0][0]),
                                 FFTW_BACKWARD, FFTW_MEASURE);
   phikxFFT   = fftw_plan_dft_2d(nx, ny, reinterpret_cast<fftw_complex*>(&phikx[0][0]), reinterpret_cast<fftw_complex*>(&phix[0][0]),
                                 FFTW_BACKWARD, FFTW_MEASURE);
   phikyFFT   = fftw_plan_dft_2d(nx, ny, reinterpret_cast<fftw_complex*>(&phiky[0][0]), reinterpret_cast<fftw_complex*>(&phiy[0][0]),
                                 FFTW_BACKWARD, FFTW_MEASURE);
   zetakxFFT  = fftw_plan_dft_2d(nx, ny, reinterpret_cast<fftw_complex*>(&zetakx[0][0]),reinterpret_cast<fftw_complex*>(&zetax[0][0]),
                                 FFTW_BACKWARD, FFTW_MEASURE);
   zetakyFFT  = fftw_plan_dft_2d(nx, ny, reinterpret_cast<fftw_complex*>(&zetaky[0][0]),reinterpret_cast<fftw_complex*>(&zetay[0][0]),
                                 FFTW_BACKWARD, FFTW_MEASURE);
   advTotFFT  = fftw_plan_dft_2d(nx, ny, reinterpret_cast<fftw_complex*>(&advTot[0][0]),reinterpret_cast<fftw_complex*>(&advkTot[0][0]),
                                 FFTW_FORWARD,  FFTW_MEASURE);

   //Initialize positional and fourier grids.
   for (size_t i = 0; i < nx; ++i)
   {
      xGrid[i]  = i*dx;
      kxGrid[i] = -M_PI*nx/lx + 2*M_PI*i/lx;
   }
   for (size_t j = 0; j < ny; ++j)
   {
      yGrid[j]  = j*dy;
      kyGrid[j] = -M_PI*ny/ly + 2*M_PI*j/ly;
   }

   //Shift kx and ky to how fftw will index phik data.
   std::rotate(kxGrid.begin(), kxGrid.begin() + nx/2, kxGrid.end());
   std::rotate(kyGrid.begin(), kyGrid.begin() + ny/2, kyGrid.end());

   //Load up aliasing vectors. They remove highest 1/3 of k-modes (on each side) to keep nonlinear vals in our k-range (2/3 rule)
   //Store in nonshifted mode to use for calculations.
   kxAliased.fill(1);
   kyAliased.fill(1);
   for (size_t i = ceil(nx/3); i < 2*floor(nx/3) + nx%3; ++i)
      kxAliased[i] = 0;
   for (size_t j = ceil(ny/3); j < 2*floor(ny/3) + ny%3; ++j)
      kyAliased[j] = 0;

   //Calculate kconst needed for HM algorithm. kConst = 1/(1 + kx**2 + ky**2).
   for (size_t j = 0; j < ny; ++j)
      for (size_t i = 0; i < nx; ++i)
         kConst[j][i] = 1/(1 + kxSqGrid[i] + kySqGrid[j]);

   //Set up ICs.
   //Random strong mode + random weaker modes + random phase shifts in each.
   constexpr size_t mx = 8;
   constexpr size_t my = 8;
   std::complex<double> A[my][mx];
   A[0][2] = exp(2.1i);
   A[0][3] = 0.1*exp(1.7i);
   A[7][3] = 0.05*exp(1.1i);
   A[5][2] = 0.05*exp(0.1i);
   A[3][6] = 0.05*exp(0.7i);
   A[4][7] = 0.05*exp(2.7i);
   for (size_t j = 0; j < ny; ++j)
      for (size_t i = 0; i < nx; ++i)
         for (size_t m2 = 0; m2 < my; ++m2)
            for (size_t m1 = 0; m1 < mx; ++m1)
               phi[j][i] = phi[j][i] + A[m2][m1]*exp(1i*kxGrid[m1]*xGrid[i]+1i*kyGrid[m2]*yGrid[j]);

   //Transpose phi and reset phi to real part - so also do diagonals...
   //TODO: Do all these arrays really need to be kept real until final storage???
   double temp = 0;
   for (size_t j = 0; j < ny; ++j)
   {
      for (size_t i = 0; i < j; ++i)
      {
         temp      = phi[j][i].real();
         phi[j][i] = phi[i][j].real();
         phi[i][j] = temp;
      }
   }

   fftw_execute(phiForeFFT);

   //Store phi and phik into the time arrays.
   for (size_t j = 0; j < ny; ++j)
      for (size_t i = 0; i < nx; ++i)
      {
         phit[0][j][i]  = phi[j][i].real();
         phikt[0][j][i] = std::abs(phik[j][i]);
      }

   std::cout << "Data initialization complete." << std::endl;
}

void Run()
{
   std::cout << "Starting main data loop." << std::endl;
   std::cout << "Running for " << nt << " frames and saving data every " << saveRate << " frames." << std::endl;

   //Set up necessary arrays once for all looping.
   //Runge-Kutta pieces.
   std::complex<double> rk[ny][nx]  = {0};
   std::complex<double> rk1[ny][nx] = {0};
   std::complex<double> rk2[ny][nx] = {0};
   std::complex<double> rk3[ny][nx] = {0};
   std::complex<double> rk4[ny][nx] = {0};

   for (size_t t = 0; t < nt; ++t)
   {
      //Calculate each RK piece.
      for (size_t j = 0; j < ny; ++j)
         for (size_t i = 0; i < nx; ++i)
            rk1[j][i] = phik[j][i];
      Advance(rk1);

      for (size_t j = 0; j < ny; ++j)
         for (size_t i = 0; i < nx; ++i)
            rk2[j][i] = phik[j][i] + 0.5*dt*rk1[j][i];
      Advance(rk2);

      for (size_t j = 0; j < ny; ++j)
         for (size_t i = 0; i < nx; ++i)
            rk3[j][i] = phik[j][i] + 0.5*dt*rk2[j][i];
      Advance(rk3);

      for (size_t j = 0; j < ny; ++j)
         for (size_t i = 0; i < nx; ++i)
            rk4[j][i] = phik[j][i] + dt*rk3[j][i];
      Advance(rk4);

      //Put together all the RK pieces and time-advance phik.
      for (size_t j = 0; j < ny; ++j)
         for (size_t i = 0; i < nx; ++i)
         {
            rk[j][i]    = (rk1[j][i] + 2.0*rk2[j][i] + 2.0*rk3[j][i] + rk4[j][i])/6.0;
            phik[j][i] += dt*rk[j][i];
         }

      fftw_execute(phiBackFFT);
      //Renormalize inverse transforms.
      for (size_t j = 0; j < ny; ++j)
         for (size_t i = 0; i < nx; ++i)
            phi[j][i] /= std::complex<double>(nx*ny,0);

      //Convert and store run data. (First frame already saved above...)
      if (t%saveRate == 0 && t > 0)
      {
         std::cout << "Frame: " << t << "/" << nt << std::endl;
         for (size_t j = 0; j < ny; ++j)
         {
            for (size_t i = 0; i < nx; ++i)
            {
               phit[t/saveRate][j][i]  = phi[j][i].real();
               phikt[t/saveRate][j][i] = std::abs(phik[j][i]);
            }
         }
      }
   }

   fftw_destroy_plan(phiForeFFT);
   fftw_destroy_plan(phiBackFFT);
   fftw_destroy_plan(phikxFFT);
   fftw_destroy_plan(phikyFFT);
   fftw_destroy_plan(zetakxFFT);
   fftw_destroy_plan(zetakyFFT);
   fftw_destroy_plan(advTotFFT);
   fftw_cleanup();

   std::cout << "Finished storing run data." << std::endl;
}

void Advance(std::complex<double> (&phik)[ny][nx])
{
   for (size_t j = 0; j < ny; ++j)
      for (size_t i = 0; i < nx; ++i)
      {
         zetak[j][i]     = -(kxSqGrid[i] + kySqGrid[j])*phik[j][i];

         phikx[j][i]     = 1i*kxGrid[i]*phik[j][i]*kxAliased[i]*kyAliased[j];
         phikyTrue[j][i] = 1i*kyGrid[j]*phik[j][i]; //Need to hold on to this unaliased for eqn below.
         zetakx[j][i]    = 1i*kxGrid[i]*zetak[j][i]*kxAliased[i]*kyAliased[j];
         zetaky[j][i]    = 1i*kyGrid[j]*zetak[j][i]*kxAliased[i]*kyAliased[j];
         phiky[j][i]     = phikyTrue[j][i]*kxAliased[i]*kyAliased[j];
      }

   fftw_execute(phikxFFT);
   fftw_execute(phikyFFT);
   fftw_execute(zetakxFFT);
   fftw_execute(zetakyFFT);
   //Renormalize inverse transforms.
   for (size_t j = 0; j < ny; ++j)
      for (size_t i = 0; i < nx; ++i)
      {
         phix[j][i]  /= std::complex<double>(nx*ny,0);
         phiy[j][i]  /= std::complex<double>(nx*ny,0);
         zetax[j][i] /= std::complex<double>(nx*ny,0);
         zetay[j][i] /= std::complex<double>(nx*ny,0);
      }

   for (size_t j = 0; j < ny; ++j)
      for (size_t i = 0; i < nx; ++i)
         advTot[j][i] = phix[j][i].real()*zetay[j][i].real() - zetax[j][i].real()*phiy[j][i].real();

   fftw_execute(advTotFFT);

   for (size_t j = 0; j < ny; ++j)
      for (size_t i = 0; i < nx; ++i)
         phik[j][i] = kConst[j][i]*(advkTot[j][i] - kappa*phikyTrue[j][i]);
}

} //Close namespace