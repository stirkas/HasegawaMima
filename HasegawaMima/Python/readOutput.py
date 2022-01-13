import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

initialCase = 3
save = False
phit = np.load('./HasegawaMima/Python/hmETG_' + str(initialCase) + '.npz')['arr_0']
phikt = np.load('./HasegawaMima/Python/hmETG_' + str(initialCase) + '_k.npz')['arr_0']

frame = 1090

plotSize = 40
nx = ny = 256
mRat = 1/2000
lx = ly = (2*np.pi/.15)*np.sqrt(mRat) 
dx = lx/nx
dy = ly/ny

#Create initial grids.
x = np.arange(nx)*dx
y = np.arange(ny)*dy
X,Y = np.meshgrid(x,y)

#Create grid vals for fourier space. Note, drop end of linspace to match the fft of the arange in real space.
#For GENE - requires extra offset. 1/nx, 1/ny
kx = (2*np.pi*nx/lx)*np.linspace(-1/2, 1/2, nx, endpoint=False)
ky = (2*np.pi*ny/ly)*np.linspace(-1/2, 1/2, ny, endpoint=False)
#Shift to match fft output for calculations.
kx = np.fft.ifftshift(kx)
ky = np.fft.ifftshift(ky)
KX,KY = np.meshgrid(kx,ky)

fig = plt.figure(num=None, figsize=(12,6), dpi=100)
padding   = 20
bigText   = 30
smallText = 18
plt.rcParams.update({'font.size': bigText})

fig.clf()   
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
ax1.clear()
ax2.clear()
im1 = ax1.contourf(X,Y, phit[frame,:,:]) #, 20, cmap='jet')
im2 = ax2.contourf(np.fft.fftshift(KX), np.fft.fftshift(KY), phikt[frame,:,:]) #, 20, cmap='jet')
ax1.grid()
ax2.grid()
ax1.set_title("$\\phi$",   pad=padding)
ax2.set_title("$\\phi_k$", pad=padding)
ax2.set_xlim(-plotSize, plotSize)
ax1.set_xlabel("x/$\\rho_i$",    labelpad=padding)
ax2.set_xlabel("$k_x$$\\rho_i$", labelpad=padding)
ax2.set_ylim(-plotSize, plotSize)
ax1.set_ylabel("y/$\\rho_i$",    labelpad=padding)
ax2.set_ylabel("$k_y$$\\rho_i$", labelpad=padding)
ax1.xaxis.set_tick_params(labelsize=smallText)
ax1.yaxis.set_tick_params(labelsize=smallText)
ax2.xaxis.set_tick_params(labelsize=smallText)
ax2.yaxis.set_tick_params(labelsize=smallText)
fig.colorbar(im1, ax=ax1)
fig.colorbar(im2, ax=ax2)
plt.tight_layout()

if (save):
   plt.savefig('./HasegawaMima/hmPhiETG_' + str(initialCase) + '.pdf')
else:
   plt.show()