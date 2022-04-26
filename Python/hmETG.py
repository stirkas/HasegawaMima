import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy.interpolate as terp
import os
import sys
import random
import time

#Add custom script dir.
#sys.path.append(os.getcwd())
#from geneScripts import plotPhi as genePlot

def setMainVars(nt_, nx_, ny_, lx_, ly_, dt_, mRat_, tau_, eta_, rnByRhoI_):
    global nt, nx, ny, dx, dy, lx, ly, nt, dt, mRat, tau, eta, rnByRhoI, phi
    nt = nt_
    nx = nx_
    ny = ny_
    lx = lx_
    ly = ly_
    dt = dt_
    mRat = mRat_
    tau = tau_
    eta = eta_
    rnByRhoI = rnByRhoI_

    dx = lx/nx
    dy = ly/ny

    phi = np.zeros((nx,ny))

def createGrid(kxOffset = 0, kyOffset = 0):
    global nx, ny, dx, dy #Vars to use.
    global x, y, kx, ky, kxd, kyd, X, Y, KX, KY, KXD, KYD #Vars to change.
    
    #Create initial grids.
    x = np.arange(nx)*dx
    y = np.arange(ny)*dy

    #Create grid vals for fourier space. Note, drop end of linspace to match the fft of the arange in real space.
    #For GENE - requires extra offset. 1/nx, 1/ny
    kx = (2*np.pi*nx/lx)*np.linspace(-1/2+kxOffset, 1/2, nx, endpoint=False)
    ky = (2*np.pi*ny/ly)*np.linspace(-1/2+kyOffset, 1/2, ny, endpoint=False)
    #Shift to match fft output for calculations.
    kx = np.fft.ifftshift(kx)
    ky = np.fft.ifftshift(ky)

    #Gather up alias vectors. They remove outer 1/3 of k-modes (on each side) to keep nonlinear vals in our k-range (2/3 rule)
    #Store in nonshifted mode to use for calculations.
    kxd = np.r_[np.ones(nx//3),np.zeros(nx//3+nx%3),np.ones(nx//3)]
    kyd = np.r_[np.ones(ny//3),np.zeros(ny//3+ny%3),np.ones(ny//3)]

    X,Y = np.meshgrid(x,y)
    KX,KY = np.meshgrid(kx,ky)
    KXD,KYD = np.meshgrid(kxd,kyd)

print("Initializing variables.")

#Init vars changed by each case.
lx = dx = dy = lx = ly = nt = dt = tau = eta = rnByRhoI = 1
nx = ny = 64 #Just a nice base grid size.
mRat = 1/2000
plotSize = waveFreq = 1
x = y = kx = ky = kxd = kyd = np.zeros(1)
X = Y = KX = KY = KXD = KYD = np.zeros(1)
phi = np.zeros(1)
caseString = "Unspecified"

#Set vars shared among cases.
nt = 100000
initialCase = 3
showPlot = False
saveAnim = False
currentlySavingVideo = False #Turn on once files are being saved.
saveData = True

#Set up vars specific to routines.
#Allocate initial conditions.
if (initialCase == 1):
    #2 strong modes.
    caseString = "TwoStrongModes"
    plotSize = 40
    setMainVars(nt, nx, ny, (2*np.pi/.15)*np.sqrt(mRat), (2*np.pi/.15)*np.sqrt(mRat), 10**-5, mRat, 1, 3, 500)
    createGrid()
    phi = np.cos(waveFreq*2*np.pi*Y/ly)*np.cos(waveFreq*2*np.pi*X/lx)
elif (initialCase == 2):
    #Gaussian
    caseString = "Gaussian"
    plotSize = 3/2
    setMainVars(nt, nx, ny, (2*np.pi/.15)*np.sqrt(mRat), (2*np.pi/.15)*np.sqrt(mRat), 10**-5, mRat, 1, 3, 500)
    createGrid()
    phi = np.exp(-((X-lx/2)**2 + (Y-ly/2)**2)/4)
elif (initialCase == 3):
    #Random strong mode + random weaker modes + random phase shifts in each.
    caseString = "StrongModeWithWeakerPerturbations"
    plotSize = 40
    mx = my = 8
    nx = ny = 256
    setMainVars(nt, nx, ny, (2*np.pi/.15)*np.sqrt(mRat), (2*np.pi/.15)*np.sqrt(mRat), 10**-5, mRat, 1, 3, 500)
    createGrid()

    A=np.zeros((mx,my))+0.*1j
    phi=np.zeros((nx,ny))+0.*1j
    A[0,2]=1*np.exp(2.1J)
    A[0,3]=0.1*np.exp(1.7J)
    A[7,3]=0.05*np.exp(1.1J)
    A[5,2]=0.05*np.exp(0.1J)
    A[3,6]=0.05*np.exp(0.7J)
    A[4,7]=0.05*np.exp(2.7J)
    for i in range(nx):
        for j in range(ny):
            for m1 in range(mx):
                for m2 in range(my):
                    phi[i,j]= phi[i,j]+A[m1,m2]*np.exp(1j*kx[m1]*x[i]+1j*ky[m2]*y[j])
    phi=np.transpose(np.real(phi))
    phi = phi * 10e-3 #Scale to work nice in ETG realm.
elif (initialCase == 4): #Fredys ICs
    plotSize = 40
    caseString = 'FredyICs'
    dt = 10**-3 #5
    setMainVars(nt, nx, ny, (2*np.pi/.15)*np.sqrt(mRat), (2*np.pi/.15)*np.sqrt(mRat), dt, mRat, 1, 3, 500)
    createGrid()

    mx = 16
    my = 16
    px=np.linspace(-mx/2,mx/2,mx)*2.*np.pi/lx
    px[8]=0
    px[7]=0
    py=np.arange(my)*2.*np.pi/ly
    A=np.zeros((mx,my))+0.*1j
    phi=np.zeros((nx,ny))+0.*1j
    # had at 0.3-0.5 loop for 64 nx.
    # had 0-8 for both modes

    for k in range(mx):
         rand=random.uniform(0,2*np.pi)
         A[int(random.uniform(4,12)),int(random.uniform(0,8))]=random.uniform(0.08,0.15)*np.exp(1j*rand)

    A[8,int(random.uniform(3*my/8,4*my/8))] = 0.11*np.exp(1j*random.uniform(0,2*np.pi))
    rand=random.uniform(0,2*np.pi)
    A[8,int(my/4)]=1.*np.exp(1j*rand)  

    # Actual IC
    for i in range(nx):
        for j in range(ny):
            for m1 in range(mx):
                for m2 in range(my):
                    phi[i,j]=phi[i,j]+A[m1,m2]*np.exp(1j*px[m1]*x[i]+1j*py[m2]*y[j])
    phi = np.transpose(np.real(phi))
elif (initialCase == 5):
    #Load output from GENE.
    caseString = 'GENE'
    #Setup sim vars.
    plotSize = 40
    rnByRhoI = 213.6 #500 = Haotian's r_n/rho_i #GENE - 213.6
    dt = (1/rnByRhoI)*(3.5*10**-2)
    setMainVars(300000, 512, 512, 5.63558, 2.96377, dt, mRat, 1, 3.135, rnByRhoI)
    createGrid(1/nx, 1/ny)

    #Begin loading data.
    geneTimesteps = 637
    nxGene = 193
    nyGene = 33
    genePhi = np.zeros((geneTimesteps, nxGene, nyGene)) #[t,x,y]
    #genePlot.readPhi(os.getcwd() + '/geneScripts/GoerlerETG/phi_0021_0025_r.dat', genePhi)

    #Normalize phi from GENE to Haotians ETG eqns. Just involves a factor of rho_star.
    #Transpose because in a plot y is rows and x is columns.
    # 136, 485, 601
    phi = genePhi[136,:,:]*.003 #rhostar for my GENE runs is .003.
    phi = np.transpose(phi)

    #Finally, interpolate to get a more symmetric grid so HM code works better.
    dxGene = lx/nxGene
    dyGene = ly/nyGene
    xGene = np.arange(nxGene)*dxGene
    yGene = np.arange(nyGene)*dyGene
    XG, YG = np.meshgrid(xGene, yGene)

    #Plotting useful for taking snapshots.
    #fig = plt.figure(num=None, figsize=(12,6), dpi=100)
    #plt.rcParams.update({'font.size': 30})
    #plt.contourf(XG, YG, np.transpose(genePhi[601,:,:]), 20, cmap='jet')
    #plt.title("$\\phi$", pad=20)
    #plt.xlabel("x/$\\rho_i$", labelpad=14)
    #plt.ylabel("y/$\\rho_i$", labelpad=14)
    #plt.xticks(fontsize=18)
    #plt.yticks(fontsize=18)
    #plt.grid()
    #plt.tight_layout()
    #plt.savefig('./HasegawaMima/genePhiETG.pdf')

    interpPhi = terp.interp2d(xGene, yGene, phi, kind='linear')
    phi = interpPhi(x, y)

else:
    sys.exit("Invalid initial conditions. Current value: " + caseString + ".")

phik = np.fft.fft2(phi)
print("Loaded initial conditions: " + caseString)

#Allocate space for time routine.
saveRate = 50
print("Running for " + str(nt) + " frames and saving data every " + str(saveRate) + " frames.")
numFrames = nt//saveRate #Save images once every twenty five dt's
phit  = np.zeros((numFrames,ny,nx))
phikt = np.zeros((numFrames,ny,nx))
phit[0,:,:]  = np.real(phi)
phikt[0,:,:] = np.abs(np.fft.fftshift(phik))

#Set up kx^2
kx2 = kx**2
ky2 = ky**2
KX2, KY2 = np.meshgrid(kx2,ky2)
kconst = 1/(1+(1+tau)*mRat*(KX2+KY2)/2) #Save off const for later calculations.

def adv(phik):
    zetak = -(KX2+KY2)*phik

    #Define derivs for main equation.
    phikx  = 1j*KX*phik;  phix  = np.real(np.fft.ifft2(phikx*KXD*KYD))
    phiky  = 1j*KY*phik;  phiy  = np.real(np.fft.ifft2(phiky*KXD*KYD))
    zetakx = 1j*KX*zetak; zetax = np.real(np.fft.ifft2(zetakx*KXD*KYD))
    zetaky = 1j*KY*zetak; zetay = np.real(np.fft.ifft2(zetaky*KXD*KYD))

    term1 = ((1+tau)**2)*mRat*rnByRhoI/4
    term2 = (1+eta)/2
    term3 = tau*(1+tau)*(1+eta)*mRat/4
    derivative = kconst*(term1*np.fft.fft2(phix*zetay-zetax*phiy) + term2*phiky + term3*zetaky)

    return derivative

#Main loop
print("Starting main data loop.")
for it in range(1,nt+1):
    rk1 = adv(phik)
    rk2 = adv(phik + 0.5*dt*rk1)
    rk3 = adv(phik + 0.5*dt*rk2)
    rk4 = adv(phik + dt*rk3)

    rk = (rk1 + 2*rk2 + 2*rk3 + rk4)/6

    phik = phik + dt*rk #RK4 advance.
    
    #Store plots as specified.
    if ((it%(saveRate))==0):
        print("Storing frame data: " + str(it) + "/" + str(nt) + ".")
        phit[(it//saveRate)-1,:,:]  = np.real(np.fft.ifft2(phik))
        phikt[(it//saveRate)-1,:,:] = np.abs(np.fft.fftshift(phik))

print("Finished storing run data.")

def update_anim(it):
    fig.clf()    
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    ax1.clear()
    ax2.clear()
    im1 = ax1.contourf(X,Y, phit[it,:,:])
    im2 = ax2.contourf(np.fft.fftshift(KX), np.fft.fftshift(KY), phikt[it,:,:])
    ax1.grid()
    ax2.grid()
    ax1.title.set_text("$\\phi$")
    ax2.title.set_text("$\\phi_k$")
    global plotSize
    ax2.set_xlim(-plotSize, plotSize)
    ax1.set_xlabel("x/$\\rho_i$")
    ax2.set_xlabel("$k_x$$\\rho_i$")
    ax2.set_ylim(-plotSize, plotSize)
    ax1.set_ylabel("y/$\\rho_i$")
    ax2.set_ylabel("$k_y$$\\rho_i$")
    #fig.colorbar(im1, ax=ax1)
    #fig.colorbar(im2, ax=ax2)
    plt.tight_layout()
    plt.rcParams.update({'font.size': 12})

    if (currentlySavingVideo):
        print("Saving frame: " + str(it+1) + "/" + str(numFrames) + ".")
    else:
        print("Displaying frame: " + str(it+1) + "/" + str(numFrames) + ".")

    if ((it+1)==numFrames): #Since it index starts at 0, but numFrames is a count so doesn't include 0.
        if (showPlot):
            plt.close(fig)

# Set up formatting for the movie files
#Writer = animation.writers['ffmpeg'] #Requires ffmpeg package on linux.
#writer = Writer(fps=30, bitrate=-1, codec='h264')

fig = plt.figure(num=None, figsize=(10,5))
#anim=animation.FuncAnimation(fig,update_anim,frames=numFrames,repeat=False)

#if (showPlot):
#    plt.show()

if (saveData):
    print('Saving data array.')
    np.savez_compressed('hmETG_' + str(initialCase) + '.npz', phit)
    np.savez_compressed('hmETG_' + str(initialCase) + '_k.npz', phikt)

#print("Starting animation.")

#if (saveAnim):
#    print("Saving animation. This will probably take a few minutes...")
#    currentlySaving = True
#    anim.save('hmETG_' + str(initialCase) + '.mp4', writer=writer)
sys.exit("Animation complete.")