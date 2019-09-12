import simread.readsnapHDF5 as ws
import numpy as np
import matplotlib.pyplot as plt
import plotting.images.gas_images     as gasimages
import plotting.images.stellar_images as stellarimages
import visualization.image_maker as pfh_image_maker
import glob
import sys
import os 
import util.hdf5lib as hdf5lib
import matplotlib.patheffects as PathEffects


class frame_maker:
    def __init__(self, path, min_time=None, max_time=None, n_frames=None, 
			fov=15.0, npixels=256,
			tag='',
                        frametag='',
			interp=False,
			zoom_factor=0.0,
			IFU=False,
			plottype='gas',
                        threecolor=True,
                        dynrange=0.6e5, maxden=0.01,
                        plot_stars=False,
                        theta=0,
                        proj_depth=0,
                        arepo=0,
                        set_bh_center_and_leave=0,
                        center_pos=[0,0,0],
			**kwargs):

	print path
	self.sdir = path
	snap_list = np.array(glob.glob(path+'/snap*'))	#_*hdf5'))
	print snap_list
        try:
	    all_snapnrs = np.array([int( file[file.index('snapshot_')+9:file.index('.hdf5')] ) for file in snap_list], dtype=int)
        except:
            try:
                all_snapnrs = np.array([int( file[file.index('snap_')+5:file.index('.hdf5')] ) for file in snap_list], dtype=int)
            except:
                all_snapnrs = np.array([int( file[file.index('snapdir_')+8:] ) for file in snap_list], dtype=int)



	all_snapnrs = np.sort(all_snapnrs)

#	all_snapnrs = all_snapnrs[ all_snapnrs > 0 ]

	print all_snapnrs
	first_snapnr = all_snapnrs[0]

	new_snap_list = [None] * all_snapnrs.shape[0] #(snap_list.shape[0]
	self.all_snapnrs = all_snapnrs

	snap_timing  = np.zeros(all_snapnrs.shape[0])
        for index in np.arange(all_snapnrs.shape[0]):
            try:
                this_file = path+"snapshot_"+str(index+int(first_snapnr)).zfill(3)+".hdf5"
	        new_snap_list[index] = this_file
#	        print this_file
	        head = ws.snapshot_header(this_file[:this_file.index('.hdf5')])
	        snap_timing[index] = head.time
            except:
                this_file = path+"snapdir_"+str(index+int(first_snapnr)).zfill(3)+"/snapshot_"+str(index+int(first_snapnr)).zfill(3)	#+".0.hdf5"
                new_snap_list[index] = this_file
#                print this_file
                head = ws.snapshot_header(this_file)	#[:this_file.index('.hdf5')])
                snap_timing[index] = head.time
#                snap_timing[index] = -1.0
#		print "snapshot failed..."
	

 	self.interp		= interp   
	self.savetag		= tag
        self.frametag           = frametag
	self.snap_list          = np.array(new_snap_list)
	self.snap_timing	= snap_timing
	self.fov		= fov
	self.npixels		= npixels
	self.maxden		= maxden
	self.dynrange		= dynrange
	self.zoom_factor	= zoom_factor
	self.IFU		= IFU			# True / False
	self.plottype		= plottype		# gas, sfr, stellar XYZ band
	self.threecolor		= threecolor		# True / False 
        self.plot_stars         = plot_stars
        self.theta              = theta
        self.proj_depth         = proj_depth
        self.arepo              = arepo

        print "snap timing:"
        print snap_timing

	if (not os.path.exists( "./plots/frames/"+self.savetag ) ):
            os.makedirs("./plots/frames/"+self.savetag )

	if min_time==None:
	    self.min_movie_time = np.min(snap_timing)
	else:
	    self.min_movie_time = min_time
	if max_time==None:
	    self.max_movie_time = np.max(snap_timing)
	else:
	    self.max_movie_time = max_time
	if n_frames != None:
 	    self.n_frames       = n_frames
        else:
	    self.n_frames       = snap_list.shape[0]
	self.dt_frames = (self.max_movie_time - self.min_movie_time) / (1.0*self.n_frames-1.0)

        if center_pos == [0,0,0]:
            print "here!"
            print set_bh_center_and_leave
            if set_bh_center_and_leave==1:
                bhpos = ws.read_block(self.snap_list[0], 'POS ', 5)[0]
                self.center_pos=bhpos
            else:
                self.center_pos=np.array([0,0,0])
        else:
            self.center_pos = center_pos



    def make_one_frame(self, iii, frame_type=0, snapnr_label=False, overwrite=True, show_time_label=False, **kwargs):
        
        frame_name = "./plots/frames/"+self.savetag+"/frame_"+self.frametag+self.savetag+"_"+str(self.all_snapnrs[iii]).zfill(4)+".png"
        if ( overwrite or not os.path.exists( frame_name ) ):
	    this_time = iii * self.dt_frames
	    time_diff = self.snap_timing - this_time - self.snap_timing[0]

            try:
	        index = np.where(np.abs(time_diff) == np.min( np.abs(time_diff) ))[0][0]
            except:
                print np.where(np.abs(time_diff) == np.min( np.abs(time_diff) ))
                sys.exit()

	    if time_diff[index] < 0 and self.interp==True:        # on the negative side, get the positive side
	            this_file1 = self.snap_list[index]
	            this_file2 = self.snap_list[index+1]
	            interp_coef = (this_time - self.snap_timing[index]) / (  self.snap_timing[index+1] - self.snap_timing[index])
	    elif time_diff[index] > 0 and self.interp==True:
	            this_file1 = self.snap_list[index-1]
	            this_file2 = self.snap_list[index]
	            interp_coef = (this_time - self.snap_timing[index-1]) / (  self.snap_timing[index+1] - self.snap_timing[index])
	    else:
		    this_file1 = self.snap_list[index]
		    this_file2 = self.snap_list[index]
		    interp_coef = 0.0


            fov_min = (1.0 - self.zoom_factor) *  self.fov
            fov_max = self.fov
            log_fac = np.log10( fov_max / fov_min )
            fac = (3.1416 /2.0 + np.arctan( (1.0*iii / (1.0*self.n_frames)) * 2 * 3.14159 ) ) / (3.1416)     #0.0 -> 1.0 
            this_fov = fov_max * 10.0**( -1.0 * fac * log_fac   )
            if self.proj_depth > 0:
                this_thick = self.proj_depth
            else:
                this_thick = this_fov

            print " this_fov and thickness: "
            print this_fov, this_thick

	    if frame_type==0:	# single gas projection
                image = gasimages.gas_image( '', 0,
	                        this_file=this_file1, next_file=this_file2, interp_coef=interp_coef,
	                        xrange=[-self.fov,self.fov], yrange=[-self.fov,self.fov], 
				dynrange=self.dynrange,maxden=self.maxden, pixels=self.npixels, 
				interp=self.interp,
				**kwargs )

	        fig = plt.figure(figsize=(5,5))
	        ax = fig.add_subplot(1,1,1)
	        fig.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.0, wspace=0.0, hspace=0.0)
	        imgplot = ax.imshow(image,origin='lower', extent=[-self.fov, self.fov, -self.fov, self.fov])
	        ax.set_xlim([-self.fov, self.fov])
	        ax.set_ylim([-self.fov, self.fov])
	        fig.savefig(frame_name, dpi=self.npixels)
	    elif frame_type==1:
                image1 = gasimages.gas_image( '', 0,
                                this_file=this_file1, next_file=this_file2, interp_coef=interp_coef,
                                xrange=[-self.fov,self.fov], yrange=[-self.fov,self.fov],
                                dynrange=self.dynrange,maxden=self.maxden, pixels=self.npixels,
                                interp=self.interp,
                                **kwargs )

                image2 = gasimages.gas_image( '', 0,
                                this_file=this_file1, next_file=this_file2, interp_coef=interp_coef,
                                xrange=[-self.fov,self.fov], yrange=[-self.fov,self.fov],
                                dynrange=self.dynrange,maxden=self.maxden, pixels=self.npixels,
                                interp=self.interp,
                                **kwargs )

                fig = plt.figure(figsize=(10,5))
                ax = fig.add_subplot(1,2,1)
                fig.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.0, wspace=0.0, hspace=0.0)
                imgplot = ax.imshow(image1,origin='lower', extent=[-self.fov, self.fov, -self.fov, self.fov])
                ax.set_xlim([-self.fov, self.fov])
                ax.set_ylim([-self.fov, self.fov])

                ax = fig.add_subplot(1,2,2)
                fig.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.0, wspace=0.0, hspace=0.0)
                imgplot = ax.imshow(image2,origin='lower', extent=[-self.fov, self.fov, -self.fov, self.fov])
                ax.set_xlim([-self.fov, self.fov])
                ax.set_ylim([-self.fov, self.fov])
                fig.savefig(frame_name, dpi=self.npixels)


	    elif frame_type==2:
                image = stellarimages.stellar_image( '', 0,
                                this_file=this_file1, next_file=this_file2, interp_coef=interp_coef,
                                xrange=[-self.fov,self.fov], yrange=[-self.fov,self.fov],
                                dynrange=self.dynrange,maxden=self.maxden, pixels=self.npixels,
                                interp=self.interp,
                                **kwargs )

#		massmap = np.transpose( massmap )
                pic = np.zeros_like(image)
                for jjj in np.arange(image.shape[2]):
                    pic[:,:,jjj] = np.transpose( image[:,:,jjj] )


                fig = plt.figure(figsize=(1,1))
                ax = fig.add_subplot(1,1,1)
                fig.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.0, wspace=0.0, hspace=0.0)
		extent = [-this_fov, this_fov, -this_fov, this_fov]
                imgplot = ax.imshow(pic,origin='lower', extent=extent)

                fig.savefig(frame_name, dpi=self.npixels)


	    elif frame_type==3:

		print kwargs.get('center_pos',[0,0,0])[0]
                print self.center_pos
                print self.arepo
                print this_fov
                print this_thick
		show_gasstarxray='gas'
		if self.plottype=='gas':  show_gasstarxray='gas'
                if self.plottype=='sfr':  show_gasstarxray='sfr'
		if self.IFU or self.plottype=='sfr':
			self.threecolor=False

                if this_thick>0:
                    zrange=[-this_thick, this_thick]
                else:
                    zrange=0


                image,massmap = pfh_image_maker.image_maker( self.sdir, self.all_snapnrs[iii],
                                                        snapdir_master='',
                                                        outdir_master='',
                                                        filename_set_manually='tmp_remove',
                                                        center_pos=self.center_pos,
                                                        #center_on_bh=0,
                                                        xrange=[-this_fov, this_fov],
                                                        zrange=zrange,
#                                                        dynrange=self.dynrange, maxden=self.maxden,
							threecolor=self.threecolor,
                                                        pixels=self.npixels, 
							show_gasstarxray=show_gasstarxray,
							theta=self.theta,
                                                        arepo=self.arepo,
							**kwargs)


		print "image values"
		print image
		print image.shape
		print massmap
		print " "



		massmap = np.transpose( massmap )
                pic = np.zeros_like(image)
                for jjj in np.arange(image.shape[2]):
                    pic[:,:,jjj] = np.transpose( image[:,:,jjj] )

		if self.IFU:
		  manga_mask = self.create_ifu_mask( pic )

		
		#star_xyz = ws.read_block(this_file1, 'POS ', 4)


                fig = plt.figure(figsize=(1,1))
                ax = fig.add_subplot(1,1,1)
                fig.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.0, wspace=0.0, hspace=0.0)
		extent = [-this_fov, this_fov, -this_fov, this_fov]
                imgplot = ax.imshow(pic,origin='lower', extent=extent)
	        if self.IFU:
		    imgplot = ax.imshow(np.log10(manga_mask), origin='lower', extent=extent, cmap='Greys', vmin=-1, vmax=1, alpha=0.5)

                if self.plot_stars:
                    try:
                        star_xyz = ws.read_block(this_file1, 'POS ', 4)
                        ax.scatter(star_xyz[:,0]-kwargs.get('center_pos',[0,0,0])[0], star_xyz[:,1] - kwargs.get('center_pos',[0,0,0])[1], s=0.1, lw=0, facecolor='r')
                    except:
		        print "no stars to scatter"

                ax.set_xlim([-1.0*this_fov, this_fov])
                ax.set_ylim([-1.0*this_fov, this_fov])
                plt.axis('off')
                if snapnr_label:
		    ax.text( 0.75*this_fov, 0.85*this_fov, str(iii).zfill(3), fontsize=3, path_effects=[PathEffects.withStroke(linewidth=1,foreground="w")] ) #color='w' )
                if show_time_label:
                    ax.text( 0.05, 0.95, "t={:.2f}Gyr".format(iii * self.dt_frames), fontsize=2, 
                                 ha='left',
                                 transform=ax.transAxes,
                                 path_effects=[PathEffects.withStroke(linewidth=0.2,foreground="w")] )



		if self.IFU:
		    png_file = "./plots/frames/frame_ifu_"+self.savetag+"_"+str(self.all_snapnrs[iii]).zfill(4)+".png"
                    hdf5_file = "./plots/frames/frame_ifu_"+self.savetag+"_"+str(self.all_snapnrs[iii]).zfill(4)+".hdf5"

                    fig.savefig(png_file, dpi=self.npixels)
		    
		    f=hdf5lib.OpenFile(hdf5_file, mode = 'w')
		    
		    group = hdf5lib.CreateGroup(f, "GasSurfaceDensity")
		    hdf5lib.CreateArray(f, group, 'data', massmap.flatten() )
		    group = hdf5lib.CreateGroup(f, "MangaMap")
		    hdf5lib.CreateArray(f, group, 'manga_mask', -1.0*manga_mask.flatten() )
		    f.close()


		else:
                    fig.savefig(frame_name, dpi=self.npixels)
            elif frame_type==4:
                image,massmap = pfh_image_maker.image_maker( self.sdir, self.all_snapnrs[iii],
                                                        snapdir_master='',
                                                        outdir_master='',
                                                        show_time_label=0,
                                                        filename_set_manually='tmp_remove',
                                                        center_on_bh=0,
                                                        xrange=[-this_fov, this_fov],
                                                        dynrange=self.dynrange, maxden=self.maxden,
                                                        pixels=self.npixels, **kwargs)

                pic1 = np.zeros_like(image)
                for jjj in np.arange(3):
                    pic1[:,:,jjj] = np.transpose( image[:,:,jjj] )

                image,massmap = pfh_image_maker.image_maker( self.sdir, self.all_snapnrs[iii],
                                                        snapdir_master='',
                                                        outdir_master='',
                                                        show_time_label=0,
                                                        filename_set_manually='tmp_remove',
                                                        center_on_bh=0,
                                                        xrange=[-this_fov, this_fov],
                                                        dynrange=self.dynrange, maxden=self.maxden,
                                                        pixels=self.npixels, 
							theta=90, **kwargs)

                pic2 = np.zeros_like(image)
                for jjj in np.arange(3):
                    pic2[:,:,jjj] = np.transpose( image[:,:,jjj] )

                fig = plt.figure(figsize=(2,1))
                ax = fig.add_subplot(1,2,1)
                fig.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.0, wspace=0.0, hspace=0.0)
                imgplot = ax.imshow(pic1,origin='lower', extent=[-this_fov, this_fov, -this_fov, this_fov])
                ax.set_xlim([-1.0*this_fov, this_fov])
                ax.set_ylim([-1.0*this_fov, this_fov])
                plt.axis('off')

                ax = fig.add_subplot(1,2,2)
                fig.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.0, wspace=0.0, hspace=0.0)
                imgplot = ax.imshow(pic2,origin='lower', extent=[-this_fov, this_fov, -this_fov, this_fov])
                ax.set_xlim([-1.0*this_fov, this_fov])
                ax.set_ylim([-1.0*this_fov, this_fov])
                plt.axis('off')
                fig.savefig(frame_name, dpi=self.npixels)


	    else:
                fig = plt.figure(figsize=(10,5))
                fig.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.0, wspace=0.0, hspace=0.0)

                ax = fig.add_subplot(1,2,1)
                imgplot = ax.imshow(image1,origin='lower', extent=[-self.fov, self.fov, -self.fov, self.fov])
                ax.set_xlim([-self.fov, self.fov])
                ax.set_ylim([-self.fov, self.fov])
                ax.axes.xaxis.set_ticklabels([])
                ax.axes.yaxis.set_ticklabels([])

                ax = fig.add_subplot(1,2,2)
                imgplot = ax.imshow(image2,origin='lower', extent=[-self.fov/10.0, self.fov/10.0, -self.fov/10.0, self.fov/10.0])
                ax.set_xlim([-self.fov/10.0, self.fov/10.0])
                ax.set_ylim([-self.fov/10.0, self.fov/10.0])

		ax.axes.xaxis.set_ticklabels([])
                ax.axes.yaxis.set_ticklabels([])

                fig.savefig(frame_name, dpi=self.npixels)




    def make_all_frames(self, rank=0, size=1, **kwargs):
        print "Number of frames to process: {:d}".format( self.n_frames )
	for index in range(self.n_frames):
            print (index % size), rank, (index % size) == rank
	    if (index % size) == rank:

                    print "processing frame {:d} on rank {:d}".format( index, rank )
#	        try:
                    self.make_one_frame(index, **kwargs)
#                except:
#		    print "one frame failed "+str(index)


    def create_ifu_mask( self, pic, ifu_type='Manga'):
	if ifu_type=='Manga':
            ifu_x_pos = np.array([])	#( 1.5 + np.arange(13) ) / 15.0  * pic.shape[0]
            ifu_y_pos = np.array([])	#zeros( 13 ) + pic.shape[0] / 2.0

            for jjj in np.arange(13):
                        djjj = 6 - jjj
                        n_fibers = 13 - np.abs( djjj )
                        this_x_pos = ( 1.5 + np.abs(djjj)/2.0 + np.arange(n_fibers) ) / 15.0  * pic.shape[0]
                        this_y_pos = np.zeros_like(this_x_pos) + ( 1.5 + jjj ) / 15.0  * pic.shape[0]
                        ifu_x_pos = np.append( ifu_x_pos, this_x_pos )
                        ifu_y_pos = np.append( ifu_y_pos, this_y_pos )

        ifu_mask = np.ones( (pic.shape[0], pic.shape[1]) )
        ifu_fiber_radius = (ifu_x_pos[1] - ifu_x_pos[0]) / 2.0 * 0.95

        pic_x = np.zeros( (pic.shape[0], pic.shape[1]) )
        pic_y = np.zeros( (pic.shape[0], pic.shape[1]) )

        for jjj in np.arange( pic.shape[0] ):  pic_x[jjj,:] = jjj
        for jjj in np.arange( pic.shape[1] ):  pic_y[:,jjj] = jjj

        for index in np.arange( ifu_x_pos.shape[0] ):
            pix_dr = np.sqrt(   (ifu_x_pos[index] - pic_x)**2 + (ifu_y_pos[index] - pic_y)**2 )
            ifu_mask[ pix_dr < ifu_fiber_radius ] = -1.0*(index+1)

        return np.transpose( ifu_mask )
