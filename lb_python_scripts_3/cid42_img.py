import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import scipy.signal as signal
import pyfits
import make_color
import glob
import make_blurred_image

def broadband_img():

	fpath_dual = '/n/scratch2/hernquist_lab/lblecha/cid42_sims/dualbhs/q0.5_fg0.4_M0.333_b0.1_orb3_mseed7e-5_nomrg/hisnapres/'
	fpath_lowk = '/n/scratch2/hernquist_lab/lblecha/cid42_sims/q0.5_fg0.5_M0.333_b0.1_orb3_tc_v1428_hsr/vhisnapres/'
	fpath_hik = '/n/scratch2/hernquist_lab/lblecha/cid42_sims/q0.5_fg0.5_M0.333_b0.1_orb3_tc_v2470_hsr/vhisnapres/'

#############################################
#	rtpath = 'rt_fov1.4_camaxes/'
#	rtpath = 'rt_fov1.4_camaxes_nompm/'
#	rtpath = 'rt_fov1.4_camaxes_lowdust/'
	rtpath = 'rt_fov1.4_camaxes_hidust/'
#	dir = fpath_dual+rtpath
#	snap = str(239)
#	dir = fpath_lowk+rtpath
#	snap = str(218)
	dir = fpath_hik+rtpath
	snap = str(217)
#############################################

	f=dir+'broadband_'+snap+'.hdf5.fits'

	cmap=plt.cm.gist_heat		
	cmap.set_bad(color='black')
	filt_idx=24    				# F814W_WFC; should be line # in filter file - 1
	cam_idx = 0 				# camera 
	kpc2pix = 400/63.  			# # pix/FOV in kpc
	arcsec2kpc = 4.989			# for z = 0.359; change as appropriate
	beamsize_arcsec = 0.04		# change for your specific observations
	#	beamsize_arcsec = 0.085		# change for your specific observations
	#	sigma=6.e-3					# sigma for Gaussian noise; current value completely arbitrary;
	#	sigma=6.e-4					# sigma for Gaussian noise; current value completely arbitrary;
	sigma=6.e-3					# sigma for Gaussian noise; current value completely arbitrary;
	# I think you want Poisson anyway
	# expectation value for Poisson noise; no idea what this should actually be
	# you'll need to figure out normalization
	#	lam=1.e-1					
	#	lam=5.e-4					
								
	im=pyfits.open(f)[18+cam_idx].data
	noise=np.random.normal(0,sigma,im[filt_idx,:,:].shape)		# Gaussian noise
	#	noise=np.random.poisson(lam,im[filt_idx,:,:].shape)			# Poisson noise
	bim=make_blurred_image.blur_image(im[filt_idx,:,:]+noise,beamsize_arcsec*arcsec2kpc*kpc2pix)
	
	#	vmax=np.log10(np.nanmax(bim))		# set max of color scale to max of image
	
	print 'max = ',np.nanmax(bim)
	print 'min = ',np.nanmin(bim)
	
#	vmax=1.4		# set max by hand
	vmax=1.2		# set max by hand
	vmin=vmax-4.5						# set min to 10^-x of vmax
#	vmin=vmax-4.7						# set min to 10^-x of vmax
#        vmin=vmax-6						
	print 'vmax = ', vmax
	print 'vmin = ', vmin

	# BH y & z coords for recoil, vk=1428:
#	bhy=[-0.05939]
#	bhz=[-0.05172]
	# BH y & z coords for recoil, vk=2470:
	bhy=[-0.05679]
	bhz=[-0.05561]
	# BH y & z coords for dual BHs:
#	bhy=[-0.04136]
#	bhz=[0.01642]
#	bhy2=[0.03052]
#	bhz2=[-0.01562]
	
	fig=plt.figure(figsize=(8,8))
#	fig=plt.figure(figsize=(10,8))
	ax=fig.add_axes((0,0,1,1))
	p=ax.imshow(np.log10(np.rot90(bim,2)),cmap=cmap,vmin=vmin,vmax=vmax,aspect='auto',extent=[-1,1,-1,1])	# plot image with log scale
#	p=ax.imshow(np.log10(bim),cmap=cmap,vmin=vmin,vmax=vmax,aspect='auto',extent=[-1,1,-1,1])	# plot image with log scale
	ax.set_xlim(-0.9,0.6)
	ax.set_ylim(-0.7,0.8)
	ax.set_autoscale_on(False)

	### big markers
#	ax.scatter([bhy],[bhz],color='black',s=140)
#	ax.scatter([bhy],[bhz],color='w',s=70)
#	ax.scatter([bhy2],[bhz2],color='black',s=140)
#	ax.scatter([bhy2],[bhz2],color='w',s=70)

        ### small markers
	ax.scatter([bhy],[bhz],color='black',s=20)
#	ax.scatter([bhy2],[bhz2],color='black',s=20)

#	ax.plot([-0.05],[-0.05],'wo')
#	cb=fig.colorbar(p)
#	fig.savefig(dir+'f814w_'+snap+'.png',transparent=False,dpi=80)
#	fig.savefig(dir+'f814w_'+snap+'_nomarker.png',transparent=False,dpi=80)
	fig.savefig(dir+'f814w_'+snap+'_nomarker_sig085.png',transparent=False,dpi=80)
