import simread.readsnapHDF5 as ws
import numpy as np
import matplotlib.pyplot as plt
import plotting.images.gas_images     as gasimages
import plotting.images.stellar_images as stellarimages
import visualization.image_maker as pfh_image_maker
import glob
import sys
import util.hdf5lib as hdf5lib
#from mpi4py import MPI


class ifu_maker:
    def __init__(self, path, 
			min_time=None, max_time=None, n_frames=None, 
			fov=15.0, npixels=256,
			tag='',
			gas_dynrange=3e2,  gas_maxden=1e-2,
			sfr_dynrange=3e2,  sfr_maxden= 1e-5,
                        star_dynrange=3e2, star_maxden= 1e-5,
			**kwargs):

        if path != None:

	    print path
	    self.sdir = path
	    snap_list = np.array(glob.glob(path+'/snap*_???.hdf5'))
	    all_snapnrs = np.array([int( file[file.index('snapshot_')+9:file.index('.hdf5')] ) for file in snap_list], dtype=int)
	    all_snapnrs = np.sort(all_snapnrs)
	    all_snapnrs = all_snapnrs[ all_snapnrs > 0 ]

	    first_snapnr = all_snapnrs[0]

	    new_snap_list = [None] * all_snapnrs.shape[0] #(snap_list.shape[0]
	    self.all_snapnrs = all_snapnrs

	    snap_timing  = np.zeros(all_snapnrs.shape[0])
            for index in np.arange(all_snapnrs.shape[0]):
                    this_file = path+"snapshot_"+str(index+int(first_snapnr)).zfill(3)+".hdf5"
	            new_snap_list[index] = this_file
	            print this_file
	            head = ws.snapshot_header(this_file[:this_file.index('.hdf5')])
	            snap_timing[index] = head.time
	    

	    self.savetag		= tag
	    self.snap_list          = np.array(new_snap_list)
	    self.snap_timing	= snap_timing
	    self.fov		= fov
	    self.npixels		= npixels
	    self.gas_dynrange       = gas_dynrange
	    self.gas_maxden		= gas_maxden
            self.sfr_dynrange       = sfr_dynrange
            self.sfr_maxden         = sfr_maxden
            self.star_dynrange       = star_dynrange
            self.star_maxden         = star_maxden
	    print self.star_dynrange, self.star_maxden 

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


    def make_one_frame(self, iii, frame_type=0, **kwargs):
	    this_time = iii * self.dt_frames
            this_fov = self.fov

	    pic = np.zeros( (self.npixels, self.npixels) )
	    ifu_r, manga_mask = self.create_ifu_mask( pic )

            hdf5_file = "./plots/frames/frame_ifu_"+self.savetag+"_"+str(self.all_snapnrs[iii]).zfill(4)+".hdf5"
            f=hdf5lib.OpenFile(hdf5_file, mode = 'w')
            group = hdf5lib.CreateGroup(f, "MangaMap")
            hdf5lib.CreateArray(f, group, 'manga_mask', -1.0*manga_mask.flatten() )
            hdf5lib.CreateArray(f, group, 'fiber_r', ifu_r )


            group = hdf5lib.CreateGroup(f, "Header")
            hdf5lib.CreateArray(f, group, 'time', [this_time] )

	    # first do gas mass plot
	    show_gasstarxray='gas'
	    threecolor=False

	    if False:
	        image,massmap = pfh_image_maker.image_maker( self.sdir, self.all_snapnrs[iii],
                                                        snapdir_master='',
                                                        outdir_master='',
                                                        show_time_label=0,
                                                        filename_set_manually='tmp_remove',
                                                        center_on_bh=0,
                                                        xrange=[-this_fov, this_fov],
                                                        dynrange=self.gas_dynrange, maxden=self.gas_maxden,
                                                        threecolor=threecolor,
                                                        pixels=self.npixels,
                                                        show_gasstarxray=show_gasstarxray,
#							include_lighting=0,
                                                        **kwargs) 

	        massmap = np.transpose( massmap )
	        pic = np.zeros_like(image)
                for jjj in np.arange(image.shape[2]):
                    pic[:,:,jjj] = np.transpose( image[:,:,jjj] )

	        fig, ax = setup_my_ifu_figure(this_fov)	#
	        extent = [-this_fov, this_fov, -this_fov, this_fov]
                imgplot = ax.imshow(pic,origin='lower', extent=extent)
	        print "pic min max:"
	        print pic.min(), pic.max(), pic.shape
	        imgplot = ax.imshow(np.log10(manga_mask), origin='lower', extent=extent, cmap='Greys', vmin=-1, vmax=1, alpha=0.5)

	        png_file = "./plots/frames/frame_ifu_gas_"+self.savetag+"_"+str(self.all_snapnrs[iii]).zfill(4)+".png"
	        fig.savefig(png_file, dpi=self.npixels)

                group = hdf5lib.CreateGroup(f, "GasSurfaceDensity")
                hdf5lib.CreateArray(f, group, 'map', massmap.flatten() )

	        gas_ifu_data = self.create_ifu_average( massmap )	# aveage surface density in fiber
	        hdf5lib.CreateArray(f, group, 'fiber_data', gas_ifu_data )
	    # =========================================== #
            # second do sfr mass plot
	    if False:
                show_gasstarxray='sfr'
                threecolor=False

                image,massmap = pfh_image_maker.image_maker( self.sdir, self.all_snapnrs[iii],
                                                        snapdir_master='',
                                                        outdir_master='',
                                                        show_time_label=0,
                                                        filename_set_manually='tmp_remove',
                                                        center_on_bh=0,
                                                        xrange=[-this_fov, this_fov],
                                                        dynrange=self.sfr_dynrange, maxden=self.sfr_maxden,
                                                        threecolor=threecolor,
                                                        pixels=self.npixels,
                                                        show_gasstarxray=show_gasstarxray,
                                                        **kwargs)

                massmap = np.transpose( massmap )
                pic = np.zeros_like(image)
                for jjj in np.arange(image.shape[2]):
                    pic[:,:,jjj] = np.transpose( image[:,:,jjj] )

                fig, ax = setup_my_ifu_figure(this_fov)  
                extent = [-this_fov, this_fov, -this_fov, this_fov]
                imgplot = ax.imshow(pic,origin='lower', extent=extent)
                imgplot = ax.imshow(np.log10(manga_mask), origin='lower', extent=extent, cmap='Greys', vmin=-1, vmax=1, alpha=0.5)

                png_file = "./plots/frames/frame_ifu_sfr_"+self.savetag+"_"+str(self.all_snapnrs[iii]).zfill(4)+".png"
                fig.savefig(png_file, dpi=self.npixels)

                group = hdf5lib.CreateGroup(f, "GasSFR")
                hdf5lib.CreateArray(f, group, 'map', massmap.flatten() )
                sfr_ifu_data = self.create_ifu_total( massmap )   # aveage surface density in fiber
                hdf5lib.CreateArray(f, group, 'fiber_data', sfr_ifu_data )
            # =========================================== #


            # =========================================== #
            # third do stellar light plot
	    if False:
                show_gasstarxray='star'
                threecolor=True

                image,massmap = pfh_image_maker.image_maker( self.sdir, self.all_snapnrs[iii],
                                                            snapdir_master='',
                                                            outdir_master='',
                                                            show_time_label=0,
                                                            filename_set_manually='tmp_remove',
                                                            center_on_bh=0,
                                                            xrange=[-this_fov, this_fov],
                                                            dynrange=self.star_dynrange, maxden=self.star_maxden,
                                                            threecolor=threecolor,
                                                            pixels=self.npixels,
                                                            show_gasstarxray=show_gasstarxray,
                                                            include_lighting=0,
                                                            **kwargs)

                massmap = np.transpose( massmap )
	        print massmap.shape
	        print image.shape

                pic = np.zeros_like(image)
                for jjj in np.arange(image.shape[2]):
                     pic[:,:,jjj] = np.transpose( image[:,:,jjj] )

                fig, ax = setup_my_ifu_figure(this_fov)
                extent = [-this_fov, this_fov, -this_fov, this_fov]
                imgplot = ax.imshow(pic,origin='lower', extent=extent, interpolation='bicubic',aspect='normal')
                imgplot = ax.imshow(np.log10(manga_mask), origin='lower', extent=extent, cmap='Greys', vmin=-1, vmax=1, alpha=0.5)

                png_file = "./plots/frames/frame_ifu_star_"+self.savetag+"_"+str(self.all_snapnrs[iii]).zfill(4)+".png"
                fig.savefig(png_file, dpi=self.npixels)

                group = hdf5lib.CreateGroup(f, "StellarImage")
                hdf5lib.CreateArray(f, group, 'map', pic.flatten() )
            # =========================================== #

            if True:
                show_gasstarxray='gas'
                threecolor=False
                image,tot_massmap = pfh_image_maker.image_maker( self.sdir, self.all_snapnrs[iii],
                                                        snapdir_master='',
                                                        outdir_master='',
                                                        show_time_label=0,
                                                        filename_set_manually='tmp_remove',
                                                        center_on_bh=0,
                                                        xrange=[-this_fov, this_fov],
                                                        dynrange=self.gas_dynrange*1e9, maxden=self.gas_maxden*1e3,
                                                        threecolor=threecolor,
                                                        pixels=self.npixels,
                                                        show_gasstarxray=show_gasstarxray,
                                                        include_lighting=0,
							do_with_colors=0,
							h_rescale_factor=3.0,
                                                        **kwargs)
                show_gasstarxray='gas_metallicity'
                threecolor=False
                image,met_massmap = pfh_image_maker.image_maker( self.sdir, self.all_snapnrs[iii],
                                                        snapdir_master='',
                                                        outdir_master='',
                                                        show_time_label=0,
                                                        filename_set_manually='tmp_remove',
                                                        center_on_bh=0,
                                                        xrange=[-this_fov, this_fov],
                                                        dynrange=self.gas_dynrange*1e9, maxden=self.gas_maxden*1e3,
                                                        threecolor=threecolor,
                                                        pixels=self.npixels,
                                                        show_gasstarxray=show_gasstarxray,
                                                        include_lighting=0,
							h_rescale_factor=3.0,
                                                        **kwargs)

	        metallicity_map = met_massmap / tot_massmap # + np.min( tot_massmap[ tot_massmap > 0 ] ) )
		print " "
                print np.min(met_massmap), np.max(met_massmap), np.median(met_massmap), np.mean(met_massmap)
                print np.min(tot_massmap), np.max(tot_massmap), np.median(tot_massmap), np.mean(tot_massmap)
		print np.min(metallicity_map), np.max(metallicity_map), np.median(metallicity_map), np.mean(metallicity_map)
		print " "
		

#                massmap = np.transpose( massmap )
#                pic = np.zeros_like(image)
#                for jjj in np.arange(image.shape[2]):
#                    pic[:,:,jjj] = np.transpose( image[:,:,jjj] )

                fig, ax = setup_my_ifu_figure(this_fov) #
                extent = [-this_fov, this_fov, -this_fov, this_fov]
		min_met =  2e-2	# the starting metallicity
		max_met =  6e-2
		metallicity_map[ metallicity_map < min_met ] = min_met
		metallicity_map[ metallicity_map > max_met ] = max_met
		metallicity_map[ tot_massmap < self.gas_maxden/self.gas_dynrange ] = min_met


		cmap=plt.get_cmap('afmhot_r')
		scaled_met = np.log10( metallicity_map /  min_met) / ( np.log10(max_met/min_met)  )
		image = cmap( scaled_met ) 


		alpha_vals = np.log10( tot_massmap / ((self.gas_maxden/10.0)/self.gas_dynrange) ) / np.log10( self.gas_dynrange )
		alpha_vals[ alpha_vals > 1.0 ] = 1.0
		alpha_vals[ alpha_vals < 0   ] = 0.0
		image[:,:,3] = image[:,:,3] * alpha_vals

		implot = ax.imshow( image, origin='lower', extent=extent ) #, cmap='afmhoat', vmin=np.log10(min_met), vmax=np.log10(max_met))
                print "pic min max:"
                print pic.min(), pic.max(), pic.shape
                imgplot = ax.imshow(np.log10(manga_mask), origin='lower', extent=extent, cmap='Greys', vmin=-1, vmax=1, alpha=0.5)

                png_file = "./plots/frames/frame_ifu_gas_metallicity_"+self.savetag+"_"+str(self.all_snapnrs[iii]).zfill(4)+".png"
                fig.savefig(png_file, dpi=self.npixels)

                #group = hdf5lib.CreateGroup(f, "GasMetallicity")
                #hdf5lib.CreateArray(f, group, 'map', massmap.flatten() )

                #gas_ifu_data = self.create_ifu_average( massmap )       # aveage surface density in fiber
                #hdf5lib.CreateArray(f, group, 'fiber_data', gas_ifu_data )
            # =========================================== #

	    f.close()




    def make_all_frames(self, rank=0, size=1, **kwargs):
	for index in np.arange(self.n_frames):
	    if (index % size) == rank:
		
                self.make_one_frame(index, **kwargs)


    def set_ifu_fibers( self, npix_total, ifu_type='Manga' ):
        if ifu_type=='Manga':
            ifu_x_pos = np.array([])    
            ifu_y_pos = np.array([])    

            for jjj in np.arange(13):
                        djjj = 6 - jjj
                        n_fibers = 13 - np.abs( djjj )
                        this_x_pos = ( 1.5 + np.abs(djjj)/2.0 + np.arange(n_fibers) ) / 15.0  * npix_total
                        this_y_pos = np.zeros_like(this_x_pos) + ( 1.5 + jjj ) / 15.0  * npix_total
                        ifu_x_pos = np.append( ifu_x_pos, this_x_pos )
                        ifu_y_pos = np.append( ifu_y_pos, this_y_pos )

	return ifu_x_pos, ifu_y_pos

    def create_ifu_mask( self, pic, ifu_type='Manga'):
	npix = pic.shape[0]
	ifu_x_pos, ifu_y_pos = self.set_ifu_fibers( npix )
	
	ifu_r_pos = np.sqrt( ((ifu_x_pos - npix)/2.0)**2 + ((ifu_y_pos - npix)/2.0)**2 )

        ifu_mask = np.ones( (npix, npix) )
        ifu_fiber_radius = (ifu_x_pos[1] - ifu_x_pos[0]) / 2.0 * 0.95

        pic_x = np.zeros( (npix, npix) )
        pic_y = np.zeros( (npix, npix) )

        for jjj in np.arange( npix ):  pic_x[jjj,:] = jjj
        for jjj in np.arange( npix ):  pic_y[:,jjj] = jjj

        for index in np.arange( ifu_x_pos.shape[0] ):
            pix_dr = np.sqrt(   (ifu_x_pos[index] - pic_x)**2 + (ifu_y_pos[index] - pic_y)**2 )
            ifu_mask[ pix_dr < ifu_fiber_radius ] = -1.0*(index+1)

        return ifu_r_pos, np.transpose( ifu_mask )


    def create_ifu_average( self, map ):	# = self.create_ifu_average( massmap )   # aveage surface density in fiber
        npix = map.shape[0]
        ifu_x_pos, ifu_y_pos = self.set_ifu_fibers( npix )
	ifu_fiber_radius = (ifu_x_pos[1] - ifu_x_pos[0]) / 2.0 * 0.95
	pic_x = np.zeros( (npix, npix) )
        pic_y = np.zeros( (npix, npix) )
        for jjj in np.arange( npix ):  pic_x[jjj,:] = jjj
        for jjj in np.arange( npix ):  pic_y[:,jjj] = jjj

	result = np.zeros_like( ifu_x_pos )
        for index in np.arange( ifu_x_pos.shape[0] ):
            pix_dr = np.sqrt(   (ifu_x_pos[index] - pic_x)**2 + (ifu_y_pos[index] - pic_y)**2 )
	    result[index] = np.mean( map[ pix_dr < ifu_fiber_radius ] ) 
	return result


    def create_ifu_total( self, map ):
        npix = map.shape[0]
        ifu_x_pos, ifu_y_pos = self.set_ifu_fibers( npix )
        ifu_fiber_radius = (ifu_x_pos[1] - ifu_x_pos[0]) / 2.0 * 0.95
        pic_x = np.zeros( (npix, npix) )
        pic_y = np.zeros( (npix, npix) )
        for jjj in np.arange( npix ):  pic_x[jjj,:] = jjj
        for jjj in np.arange( npix ):  pic_y[:,jjj] = jjj

        result = np.zeros_like( ifu_x_pos )
        for index in np.arange( ifu_x_pos.shape[0] ):
            pix_dr = np.sqrt(   (ifu_x_pos[index] - pic_x)**2 + (ifu_y_pos[index] - pic_y)**2 )
            result[index] = np.sum( map[ pix_dr < ifu_fiber_radius ] )

        return result




def setup_my_ifu_figure( this_fov ):
    fig = plt.figure(figsize=(1,1))
    ax = fig.add_subplot(1,1,1)
    fig.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.0, wspace=0.0, hspace=0.0)
    ax.set_xlim([-1.0*this_fov, this_fov])
    ax.set_ylim([-1.0*this_fov, this_fov])
    plt.axis('off')

    return fig, ax

            





