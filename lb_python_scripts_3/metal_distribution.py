import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pyfits
import numpy as np

fig = plt.figure()
ax1 = fig.add_subplot(211)
ax1.set_yscale('log')
ax1.set_xlabel('particle Z')
ax2 = fig.add_subplot(212)
ax2.set_yscale('log')
ax2.set_xlabel('grid Z')
Zrange = (0,0.25)
tmp = np.array([0,7,26,29,70,2])

for i in np.arange(0,300,50):

    sf = pyfits.open('grid_%.3d.fits'%i)
    Z = sf['PARTICLEDATA'].data.field('metallicity')
    print "snap %.3d:\nmax Z: %g (%d > 0.05)"%(i,Z.max(),Z[Z>0.05].size)
    mm_grid = sf['GRIDDATA'].data.field('mass_metals')
    mg_grid = sf['GRIDDATA'].data.field('mass_gas')
    Zgrid = (mm_grid[mg_grid>0]/mg_grid[mg_grid>0])
    print "max Zgrid: %g (%d > 0.04,%d>0.05)"%(Zgrid.max(),Zgrid[Zgrid>0.04].size,Zgrid[Zgrid>0.05].size)
    if i>0: print (np.sort(Zgrid))[-(tmp[i/50])]
    sf.close()
    ax1.hist(Z,range=Zrange,bins=10,histtype='step',lw=0.008*i+0.4)
    ax2.hist(Zgrid,range=Zrange,bins=10,histtype='step',lw=0.008*i+0.4)
    del Z, Zgrid
    
fig.savefig('metallicities.pdf')
plt.close()
