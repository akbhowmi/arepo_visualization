import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import scipy.signal as signal
import pyfits, make_color, glob
from copy import copy
import make_blurred_image

def make_broadband(path='/n/hernquistfs2/lblecha/xmm_ao15prop_plots',
		   snaps=[470],filters='ugz',cam=3,scale='auto'):

	f=pyfits.open('%s/broadband_%d.fits'%(path,snaps[0]))
	all_fnames=np.array([ x.strip() for x in f['FILTERS'].data.field('filter') ])
	f.close()
	fnames = ["%s_SDSS.res"%b for b in filters]
	band=np.array([])
	for fn in fnames:
		ix=np.where(fn==all_fnames)[0]
		if len(ix)!=1:
			print "Error: no unique match for filter name: ",fn
			return -1
		band=np.append(band,ix)
	band=copy(band[::-1])
	#band = [np.where(fnames[i]==all_fnames)[0] for i range(len(fnames))]
	#band = np.array(band).reshape(len(fnames))
	#if (min(len(fi) for fi in band),max(len(fi) for fi in band)) != (1,1)
	print "using these filters for rgb colors:"
	print "fnames: ",fnames
	print "band: ",band

	for snap in snaps:
		print "snap=%d"%snap
		filename='%s/broadband_%d.fits'%(path,snap)
		#make_color.write_color(filename,18+cam,
		#		       '%s/sdss_%s_%d_cam%d_scale1.jpg'%(path,filters,snap,cam),
		#		       band=band,overwrite=True)
		#make_color.write_color(filename,18+cam,
		#		       '%s/sdss_%s_%d_cam%d_setscale.jpg'%(path,filters,snap,cam),
		#		       band=band,overwrite=True,scale=(1.5,1,1.5))		
		if scale=='auto':
			make_color.write_color(filename,18+cam,
					       '%s/sdss_%s_%d_cam%d_autolum.jpg'%(path,filters,snap,cam),
					       band=band,overwrite=True,scale='autolum',autopercentile=0.1)
		else:
			make_color.write_color(filename,18+cam,
					       '%s/sdss_%s_%d_cam%d_myscale.jpg'%(path,filters,snap,cam),
					       band=band,overwrite=True,scale=scale)
	

## final choices for myscale:
#sdss_img.make_broadband(snaps=[170],scale=(1.25,0.9,1.25))
#sdss_img.make_broadband(snaps=[320],scale=(1.15,0.8,1.15))
#sdss_img.make_broadband(snaps=[470],scale=(1.85,1.5,1.85))
#sdss_img.make_broadband(snaps=[620],scale=(1.55,1.2,1.55))
