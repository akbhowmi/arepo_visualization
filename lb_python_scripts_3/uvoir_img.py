import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import scipy.signal as signal
import pyfits, make_color_lb, glob
from copy import copy
from PIL import Image
from PIL import ImageFont
from PIL import ImageDraw 

#import make_blurred_image

### New filter names: ###

#WISE: 
#wise/WISE-W[1234].res

#JWST:
#jwst/miri_F1000W.res          
#jwst/miri_F1130W.res          
#jwst/miri_F1280W.res   
#jwst/miri_F1500W.res          
#jwst/miri_F1800W.res    
#jwst/miri_F2100W.res          
#jwst/miri_F2550W.res          
#jwst/miri_F560W.res           
#jwst/miri_F770W.res           
#jwst/nircam_f070w.res         
#jwst/nircam_f090w.res         
#jwst/nircam_f115w.res         
#jwst/nircam_f150w.res         
#jwst/nircam_f200w.res         
#jwst/nircam_f277w.res         
#jwst/nircam_f356w.res         
#jwst/nircam_f444w.res         

#SDSS:
#sdss/[ugriz]_SDSS.res

#HST:
#hst/acs_f435w.res            
#hst/acs_f606w.res             
#hst/acs_f775w.res            
#hst/acs_f814w.res             
#hst/acs_f850lp.res           
#hst/wfc3_f105w.res            
#hst/wfc3_f125w.res           
#hst/wfc3_f140w.res            
#hst/wfc3_f160w.res           
#hst/wfc3_f275w.res            
#hst/wfc3_f336w.res           
#hst/wfc3_f350lp.res           
#hst/wfc3_f438w.res           


def make_broadband(path='/oasis/projects/nsf/hvd115/lblecha/q0.5_fg0.3_sunruns/test_fiducial_hires/',
		   snaps=[],filter_dir='sdss/',fname='gri',filters=('u_SDSS','g_SDSS','z_SDSS'),
		   cam=3,scale='auto',extra='',plotlabel='',plotlabel2='',label_only=False,plotscale=False):

        if label_only and plotlabel=='' and plotlabel2=='':
                print("make_broadband called with label_only=True but no plotlabel[2] defined. Exiting...")
                sys.exit()

	print '%s/broadband_%.3d.fits'%(path,snaps[0])
	f=pyfits.open('%s/broadband_%.3d.fits'%(path,snaps[0]))
	all_fnames=np.array([x.strip() for x in f['FILTERS'].data.field('filter')])
	f.close()
	if filter_dir and filter_dir[-1] !='/': filter_dir=filter_dir+'/'
	fnames = ["%s%s.res"%(filter_dir,b) for b in filters]
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

	if extra != '': extra='_'+extra
	if (not isinstance(snaps,list) and not isinstance(snaps,np.ndarray) and 
	    not isinstance(snaps,tuple)): snaps=[snaps]
	for snap in snaps:
		print "snap=%d"%snap
		filename='%s/broadband_%.3d.fits'%(path,snap)
		if scale=='auto':
                        img_fbase='%s/%s_%s_%.3d_cam%d_autolum%s'%(path,filter_dir.strip('/'),
								   fname,snap,cam,extra)
                        if not label_only:
                                make_color_lb.write_color(filename,f.index_of('CAMERA0-BROADBAND')+cam,img_fbase+'.jpg',
                                                          band=band,overwrite=True,scale='autolum',autopercentile=0.1,
                                                          plotlabel=plotlabel,plotlabel2=plotlabel2,plotscale=plotscale)
                else:
			img_fbase='%s/%s_%s_%.3d_cam%d_myscale%s'%(path,filter_dir.strip('/'),
								   fname,snap,cam,extra)
                        if not label_only:
                                make_color_lb.write_color(filename,f.index_of('CAMERA0-BROADBAND')+cam,img_fbase+'.jpg',
                                                          band=band,overwrite=True,scale=scale,
                                                          plotlabel=plotlabel,plotlabel2=plotlabel2,plotscale=plotscale)
	
		#if plotlabel !='':
		#	img = Image.open(img_fbase+'.jpg')
		#	draw = ImageDraw.Draw(img)
		#	# font = ImageFont.truetype(<font-file>, <font-size>)
		#	font = ImageFont.truetype("/n/home00/lblecha/envs/MY_PYENV/fonts/DejaVuSans.ttf", 18)
		#	#draw.text((x, y),"Sample Text",(r,g,b))
		#	draw.text((15, 10),plotlabel,(255,255,255),font=font)
		#	#draw.text((0.1, 0.1),plotlabel,(255,255,255))
		#	if plotscale:
		#		#draw plotscale
		#		draw.line((25,285,55,285),width=2)
		#		draw.text((15,260),'5 kpc',font=font)
                #       img.save(img_fbase+'_lbl.jpg')
			

## final choices for myscale:
#sdss_img.make_broadband(snaps=[170],scale=(1.25,0.9,1.25))
#sdss_img.make_broadband(snaps=[320],scale=(1.15,0.8,1.15))
#sdss_img.make_broadband(snaps=[470],scale=(1.85,1.5,1.85))
#sdss_img.make_broadband(snaps=[620],scale=(1.55,1.2,1.55))
