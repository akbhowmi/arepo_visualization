import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pyfits

f='/n/home00/lblecha/patrik_sunrise_stuff/laura1/mcrx_0441.fits'

cmap=plt.cm.gist_heat		
cmap.set_bad(color='black')

cam_idx = 1

cube = pyfits.open(f)[18+cam_idx].data   # nonscatter data for camera 'cam_idx'

print shape(cube) # gives dimensions of array

im = zeros((300,300))
for i in range(0,299):
    for j in range(0,299):
        im[i,j] = np.sum(cube[0:193,i,j])

print shape(im)

p = plt.imshow(np.log10(im),cmap=cmap,extent=[-1,1-1,1]

savefig('test.png')
