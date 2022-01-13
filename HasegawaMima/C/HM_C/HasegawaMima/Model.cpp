#include "Model.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>

#include <fftw3.h>

#include "Definitions.hpp"
namespace HM {

//Useful things for the actual advancing logic.
fftw_plan p1;
fftw_plan p2;
fftw_plan p3;
fftw_plan p4;
fftw_plan p5;
bool fftw_init = false;
std::complex<double> phix[ny][nx]      = {0};
std::complex<double> phiy[ny][nx]      = {0};
std::complex<double> phikx[ny][nx]     = {0};
std::complex<double> phiky[ny][nx]     = {0};
std::complex<double> zetak[ny][nx]     = {0};
std::complex<double> zetax[ny][nx]     = {0};
std::complex<double> zetay[ny][nx]     = {0};
std::complex<double> zetakx[ny][nx]    = {0};
std::complex<double> zetaky[ny][nx]    = {0};
std::complex<double> phikyTemp[ny][nx] = {0};
std::complex<double> phiTot[ny][nx]    = {0};
std::complex<double> phikTot[ny][nx]   = {0};
void Advance(std::complex<double> (&phi)[ny][nx]);

void Initialize()
{
   //Initialize positional and fourier grids.
   for (size_t i = 0; i < nx; ++i)
   {
      xGrid[i]     = i*dx;
      kxGrid[i]    = -M_PI*nx/lx + 2*M_PI*i/lx;
      //kxAliased[i] = true ? kxGrid[i] : 0;
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
         kConst[j][i] = 1/(1 + kxGrid[i]*kxGrid[i] + kyGrid[j]*kyGrid[j]);

   //Set up ICs.
   //Random strong mode + random weaker modes + random phase shifts in each.
   using namespace std::complex_literals;
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
            for (size_t m1 = 0; m2 < mx; ++m2)
               phi[j][i] = phi[j][i]+A[m2][m1]*exp(1i*kxGrid[m1]*xGrid[i]+1i*kyGrid[m2]*yGrid[j]);

   //Transpose phi and take real part.
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

   //Load FFT of phi into phik.
   auto plan = fftw_plan_dft_2d(nx, ny, reinterpret_cast<fftw_complex*>(&phi[0]), reinterpret_cast<fftw_complex*>(&phik[0]),
                                FFTW_FORWARD, FFTW_ESTIMATE);

   fftw_execute(plan);
   fftw_destroy_plan(plan);

   //Store phi and phik into the time arrays.
   for (size_t j = 0; j < ny; ++j)
      for (size_t i = 0; i < nx; ++i)
      {
         phit[0][j][i]  = phi[j][i].real();
         phikt[0][j][i] = sqrt(phik[j][i].real()*phik[j][i].real() + phik[j][i].imag()*phik[j][i].imag());
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
            rk2[j][i] = 0.5*dt*rk1[j][i];
      Advance(rk2);

      for (size_t j = 0; j < ny; ++j)
         for (size_t i = 0; i < nx; ++i)
            rk3[j][i] = 0.5*dt*rk2[j][i];
      Advance(rk3);

      for (size_t j = 0; j < ny; ++j)
         for (size_t i = 0; i < nx; ++i)
            rk4[j][i] = 0.5*dt*rk4[j][i];
      Advance(rk4);

      //Put together all the RK pieces and time-advance phik.
      for (size_t j = 0; j < ny; ++j)
         for (size_t i = 0; i < nx; ++i)
         {
            rk[j][i]    = (rk1[j][i] + 2.0*rk2[j][i] + 2.0*rk3[j][i] + rk4[j][i])/6.0;
            phik[j][i] += dt*rk[j][i];
         }

      //Convert and store run data. (First frame already saved above...)
      if (t%saveRate == 0 && t != 0)
         std::cout << "Frame: " << t << "/" << nt << std::endl;
         for (size_t j = 0; j < ny; ++j)
            for (size_t i = 0; i < nx; ++i)
            {
               phit[t/saveRate][j][i]  = phi[j][i].real();
               phikt[t/saveRate][j][i] = sqrt(phik[j][i].real()*phik[j][i].real() + phik[j][i].imag()*phik[j][i].imag());
            }
   }

   fftw_destroy_plan(p1);
   fftw_destroy_plan(p2);
   fftw_destroy_plan(p3);
   fftw_destroy_plan(p4);
   fftw_destroy_plan(p5);

   std::cout << "Finished storing run data." << std::endl;
}

void Advance(std::complex<double> (&phik)[ny][nx])
{
   using namespace std::complex_literals;

   for (size_t j = 0; j < ny; ++j)
      for (size_t i = 0; i < nx; ++i)
      {
         zetak[j][i]     = -(kxGrid[i]*kxGrid[i] + kyGrid[j]*kyGrid[j])*phik[j][i];

         phikx[j][i]     = 1i*kxGrid[i]*phik[j][i]*kxAliased[i]*kyAliased[j];
         phikyTemp[j][i] = 1i*kyGrid[j]*phik[j][i]; //Need to hold on to this for eqn below.
         zetakx[j][i]    = 1i*kxGrid[i]*zetak[j][i]*kxAliased[i]*kyAliased[j];
         zetaky[j][i]    = 1i*kyGrid[j]*zetak[j][i]*kxAliased[i]*kyAliased[j];
         phiky[j][i]     = phiky[j][i]*kyAliased[j]*kxAliased[i]*kyAliased[j];
      }

   //Fourier transform per pseudo-spectral algorithm.
   if (!fftw_init)
   {
      p1 = fftw_plan_dft_2d(nx, ny, reinterpret_cast<fftw_complex*>(&phikx[0]), reinterpret_cast<fftw_complex*>(&phix[0]),
                            FFTW_BACKWARD, FFTW_ESTIMATE);
      p2 = fftw_plan_dft_2d(nx, ny, reinterpret_cast<fftw_complex*>(&phiky[0]), reinterpret_cast<fftw_complex*>(&phiy[0]),
                            FFTW_BACKWARD, FFTW_ESTIMATE);
      p3 = fftw_plan_dft_2d(nx, ny, reinterpret_cast<fftw_complex*>(&zetakx[0]), reinterpret_cast<fftw_complex*>(&zetax[0]),
                            FFTW_BACKWARD, FFTW_ESTIMATE);
      p4 = fftw_plan_dft_2d(nx, ny, reinterpret_cast<fftw_complex*>(&zetaky[0]), reinterpret_cast<fftw_complex*>(&zetay[0]),
                            FFTW_BACKWARD, FFTW_ESTIMATE);
      p5 = fftw_plan_dft_2d(nx, ny, reinterpret_cast<fftw_complex*>(&phiTot[0]), reinterpret_cast<fftw_complex*>(&phikTot[0]),
                            FFTW_FORWARD, FFTW_ESTIMATE);
      fftw_init = true;
   }
   fftw_execute(p1);
   fftw_execute(p2);
   fftw_execute(p3);
   fftw_execute(p4);
   
   for (size_t j = 0; j < ny; ++j)
      for (size_t i = 0; i < nx; ++i)
         phiTot[j][i] = phix[j][i]*zetay[j][i] - zetax[j][i]*phiy[j][i];

   fftw_execute(p5);

   for (size_t j = 0; j < ny; ++j)
      for (size_t i = 0; i < nx; ++i)
         phik[j][i] = kConst[j][i]*phikTot[j][i] - kappa*phikyTemp[j][i];
}

} //Close namespace