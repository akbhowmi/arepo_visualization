import simread.readsnapHDF5 as ws
import numpy as np
import matplotlib.pyplot as plt
import plotting.images.gas_images     as gasimages
import plotting.images.stellar_images as stellarimages
from matplotlib.colors import LogNorm
import visualization.image_maker as pfh_image_maker
import glob
import sys
import os 
import util.hdf5lib as hdf5lib
import h5py
import units.springel_units as units

import matplotlib.patheffects as PathEffects

#from mpi4py import MPI


class phase_diagram:
    def __init__(self, path, min_time=None, max_time=None, n_frames=None, 
			npixels=256,
			tag='',
			interp=False,
			log_rho_min=-6.0, log_rho_max=4.0,
                        log_t_min=1.0,    log_t_max=8.0,
                        ignore_central_region=False,
			**kwargs):

	print path
	self.sdir = path
	snap_list = np.array(glob.glob(path+'/snap*_*.hdf5'))
        snap_list = snap_list[:-1]
        if len(snap_list) == 1:
            snap_list = np.array([ snap_list[0], snap_list[0], snap_list[0] ])
	print snap_list
        try:
	    all_snapnrs = np.array([int( file[file.index('snapshot_')+9:file.index('.hdf5')] ) for file in snap_list], dtype=int)
        except:
            all_snapnrs = np.array([int( file[file.index('snap_')+5:file.index('.hdf5')] ) for file in snap_list], dtype=int)


	all_snapnrs = np.sort(all_snapnrs)
       

	print all_snapnrs
	first_snapnr = all_snapnrs[0]

	new_snap_list = [None] * all_snapnrs.shape[0] #(snap_list.shape[0]
	self.all_snapnrs = all_snapnrs

	snap_timing  = np.zeros(all_snapnrs.shape[0])
        for index in np.arange(all_snapnrs.shape[0]):
            try:
                this_file = path+"snapshot_"+str(index+int(first_snapnr)).zfill(3)+".hdf5"

                this_file = path+"snapshot_"+str(all_snapnrs[index]).zfill(3)+".hdf5"

#index+int(first_snapnr)).zfill(3)+".hdf5"

	        new_snap_list[index] = this_file
	        print this_file
	        head = ws.snapshot_header(this_file[:this_file.index('.hdf5')])
	        snap_timing[index] = head.time
            except:
                this_file = path+"snapshot_"+str(all_snapnrs[index]).zfill(4)+".hdf5"
                new_snap_list[index] = this_file
                print this_file
                head = ws.snapshot_header(this_file[:this_file.index('.hdf5')])
                snap_timing[index] = head.time

                snap_timing[index] = -1.0
		print "snapshot failed..."



 	self.interp		= interp   
	self.savetag		= tag
	self.snap_list          = np.array(new_snap_list)
	self.snap_timing	= snap_timing
	self.npixels		= npixels
        self.log_rho_min        = log_rho_min
        self.log_rho_max        = log_rho_max
        self.log_t_min          = log_t_min    
        self.log_t_max          = log_t_max
        self.ignore_central_region = ignore_central_region

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


    def make_one_frame(self, iii, frame_type=0, snapnr_label=False, overwrite=True, just_wind=False, **kwargs):
        
        
        if just_wind:
            frame_name = "./plots/frames/"+self.savetag+"/phase_winds_"+str(self.all_snapnrs[iii]).zfill(4)+".png"
        else:
            frame_name = "./plots/frames/"+self.savetag+"/phase_"+str(self.all_snapnrs[iii]).zfill(4)+".png"

        if ( overwrite or not os.path.exists( frame_name ) ):
	    this_time = iii * self.dt_frames
	    time_diff = self.snap_timing - this_time - self.snap_timing[0]

	    index = np.where(np.abs(time_diff) == np.min( np.abs(time_diff) ))[0][0]
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

            print "processing {:s} and {:s} for time {:.2f}".format( this_file1, this_file2, this_time )

            with h5py.File(this_file1,mode='r') as f:
                rho      = np.array( f["PartType0/Density"][:] )
		ne        = np.array( f["PartType0/ElectronAbundance"][:] )
                u        = np.array( f["PartType0/InternalEnergy"][:] )
		temperature   = units.gas_code_to_temperature( u, ne )
		rho           = units.gas_code_to_cgs_density( rho )
	        
                if just_wind:
                    pos = np.array( f["PartType0/Coordinates"][:] )
                    vel = np.array( f["PartType0/Velocities"][:] )
                    r  =  np.sqrt( pos[:,0] * pos[:,0] + pos[:,1] * pos[:,1] + pos[:,2] * pos[:,2] ) + 0.001
                    vr = (pos[:,0]*vel[:,0] + pos[:,1]*vel[:,1] + pos[:,2]*vel[:,2] ) / r
                    index = vr > 500.0
                    if np.sum(index) > 1:
                        rho = rho[index]
                        ne  = ne[index]
                        u   = u[index]
                        temperature = temperature[index]
                    else:
                        rho = rho[:1]
                        ne  = ne[:1]
                        u   = u[:1]
                        temperature = temperature[:1]

                if self.ignore_central_region:
                    pos = np.array( f["PartType0/Coordinates"][:] )
                    index = pos[:,2] > 10.0
                    if np.sum(index) > 1:
                        rho = rho[index]
                        ne  = ne[index]
                        u   = u[index]
                        temperature = temperature[index]
                    else:
                        rho = rho[:1]
                        ne  = ne[:1]
                        u   = u[:1]
                        temperature = temperature[:1]

            fig = plt.figure(figsize=(5,5))
	    ax = fig.add_subplot(1,1,1)
            ax.hist2d( np.log10(rho), np.log10(temperature), 
			bins=self.npixels/16.0, 
                        range= [[self.log_rho_min,self.log_rho_max], 
                                [self.log_t_min, self.log_t_max]]  ,
			cmap='Reds',
			norm=LogNorm() )
            ax.set_xlabel(r'$\mathrm{Log(}\rho\mathrm{)}$ $\mathrm{(cm^{-3})}$')
            ax.set_ylabel(r'$\mathrm{Log(T)}$ $\mathrm{(K)}$')
	    ax.set_xlim([self.log_rho_min,self.log_rho_max]); 
            ax.set_ylim([self.log_t_min, self.log_t_max])
            fig.subplots_adjust(left=0.11,right=0.98,top=0.98,bottom=0.12)
            fig.savefig(frame_name, dpi=self.npixels/8.0)
            plt.close(fig)	#.close()


    def make_all_frames(self, rank=0, size=1, **kwargs):
	for index in np.arange(self.n_frames):
	    if (index % size) == rank:
                self.make_one_frame(index, **kwargs)


