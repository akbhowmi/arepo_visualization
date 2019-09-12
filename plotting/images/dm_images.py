#!/usr/bin/env python
"""Routines to determine and manipulate the file name formats of various
    output files."""

__author__ = "Paul Torrey and contributing authors"
__copyright__ = "Copyright 2014, The Authors"
__credits__ = ["Paul Torrey and contributing authors"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Paul Torrey"
__email__ = "ptorrey@mit.harvard.edu"
__status__ = "Beta -- forever."


import simread.readsnapHDF5 as ws
import numpy as np
import plotting.images.projection as projection
import plotting.images.contour_makepic as makepic
import util.calc_hsml as calc_hsml
import util.cast
import matplotlib.pyplot as plt
import units.springel_units
import os.path

from astropy import units as u
from astropy.coordinates import SkyCoord

class dm_image_data:
    def __init__(self, this_file, next_file=None, interp_coef=0.0, interp=True, crude=False, ra_range=None, dec_range=None,
                ra_center=None, dec_center=None, zrange=None, convert_to_sg=False, \
                **kwargs):
        if next_file==None or interp==False:	# no interp
            print "have a single file, no interp"
            self.head = ws.snapshot_header(this_file)
            self.xyz     = np.array(ws.read_block(this_file,'POS ', parttype=1))
            self.mass    = np.array(ws.read_block(this_file,'MASS', parttype=1))
            npart = np.shape(self.xyz)[0]*1.0

            if crude:
                self.hsml    = np.zeros_like( self.mass ) + self.head.boxsize / ( npart**0.333 ) * 1.0
            else:
                sys.exit()
                self.hsml    = calc_hsml.get_particle_hsml( self.xyz[:,0], self.xyz[:,1], self.xyz[:,2] , **kwargs )

            if True:    #ra_center is not None and dec_center is not None:        # need to do some rotation

                if True:        # this likely makes the most sense. MW is at the center.
                    tmp_x = self.xyz[:,0] - self.head.boxsize/2.0
                    tmp_y = self.xyz[:,1] - self.head.boxsize/2.0
                    tmp_z = self.xyz[:,2] - self.head.boxsize/2.0
                else:           # this could be possible, but puts the MW at the origin, at the edge of the box.
                    tmp_x = self.xyz[:,0]
                    tmp_y = self.xyz[:,1]
                    tmp_z = self.xyz[:,2]

                    tmp_x[ tmp_x > self.head.boxsize/2.0 ] -= self.head.boxsize
                    tmp_y[ tmp_y > self.head.boxsize/2.0 ] -= self.head.boxsize
                    tmp_z[ tmp_z > self.head.boxsize/2.0 ] -= self.head.boxsize

                tmp_r = np.sqrt( tmp_x**2 + tmp_y**2 + tmp_z**2 )

                dec = np.arcsin( tmp_z / tmp_r ) * 180.0 / 3.14159
                ra  = np.arctan2( tmp_y, tmp_x ) * 180.0 / 3.14159


                c = SkyCoord(ra=ra*u.degree, dec=dec*u.degree)
                if convert_to_sg:
                    new_coords = c.transform_to('supergalactic')
                    self.xyz[:,0] = new_coords.cartesian.x * tmp_r + self.head.boxsize/2.0
                    self.xyz[:,1] = new_coords.cartesian.y * tmp_r + self.head.boxsize/2.0
                    self.xyz[:,2] = new_coords.cartesian.z * tmp_r + self.head.boxsize/2.0
                else:
                    self.xyz[:,0] = c.cartesian.x * tmp_r + self.head.boxsize/2.0
                    self.xyz[:,1] = c.cartesian.y * tmp_r + self.head.boxsize/2.0
                    self.xyz[:,2] = c.cartesian.z * tmp_r + self.head.boxsize/2.0

                new_r = np.sqrt( self.xyz[:,0] **2 + self.xyz[:,1] **2 + self.xyz[:,2] **2 )


            if ra_range is not None and dec_range is not None:
                tmp_x = self.xyz[:,0] - self.head.boxsize/2.0
                tmp_y = self.xyz[:,1] - self.head.boxsize/2.0
                tmp_z = self.xyz[:,2] - self.head.boxsize/2.0
                tmp_r = np.sqrt( tmp_x**2 + tmp_y**2 + tmp_z**2 )

                dec = np.arcsin( tmp_z / tmp_r )
                ra  = np.arctan2( tmp_y, tmp_x )

                in_box_index = ( ra > np.min(ra_range) ) & ( ra < np.max(ra_range) ) & \
                               ( dec> np.min(dec_range)) & (dec < np.max(dec_range))

                print "Found {:d} out of {:d} particles in our projection wedge".format( np.sum(in_box_index), len(in_box_index) )
#                sys.exit() 
                self.xyz  = self.xyz[in_box_index,:]
                self.mass = self.mass[in_box_index]
                self.hsml = self.hsml[in_box_index]



        else:
            # TODO:  FINISH THIS BLOCK

    	    head = ws.snapshot_header(this_file)

            gas_id1      = np.array(ws.read_block(this_file,'ID  ', parttype=0))
	    order       = np.argsort(gas_id1) 
	    gas_id1 = gas_id1[order]

            gas_xyz1     = np.array(ws.read_block(this_file,'POS ', parttype=0))[order,:]
            gas_hsml1    = np.array(ws.read_block(this_file,'HSML', parttype=0))[order]
            gas_mass1    = np.array(ws.read_block(this_file,'MASS', parttype=0))[order]
            gas_u1       = np.array(ws.read_block(this_file,'U   ', parttype=0))[order]
            gas_nume1    = np.array(ws.read_block(this_file,'NE  ', parttype=0))[order]
            try:
                gas_z1       = np.array(ws.read_block(this_file,'Z   ', parttype=0))
		gas_z1       = np.sum( gas_z1[order,2:], axis=1)
            except:
                gas_z1       = np.array(ws.read_block(this_file,'GZ  ', parttype=0))[order]

            gas_id2      = np.array(ws.read_block(next_file,'ID  ', parttype=0))
            order       = np.argsort(gas_id2) 
            gas_id2 = gas_id2[order]

            gas_xyz2     = np.array(ws.read_block(next_file,'POS ', parttype=0))[order,:]
            gas_hsml2    = np.array(ws.read_block(next_file,'HSML', parttype=0))[order]
            gas_mass2    = np.array(ws.read_block(next_file,'MASS', parttype=0))[order]
            gas_u2       = np.array(ws.read_block(next_file,'U   ', parttype=0))[order]
            gas_nume2    = np.array(ws.read_block(next_file,'NE  ', parttype=0))[order]
            try:
                gas_z2       = np.array(ws.read_block(next_file,'Z   ', parttype=0))
		gas_z2       = np.sum( gas_z2[order,2:], axis=1)
            except:
                gas_z2       = np.array(ws.read_block(next_file,'GZ  ', parttype=0))[order]
            gas_id2      = np.array(ws.read_block(next_file,'ID  ', parttype=0))

	    i1 = np.in1d( gas_id1, gas_id2 )
	    gas_xyz1 = gas_xyz1[i1,:]
	    gas_hsml1= gas_hsml1[i1]
	    gas_mass1= gas_mass1[i1]
	    gas_u1   = gas_u1[i1]
	    gas_nume1= gas_nume1[i1]
	    gas_z1   = gas_z1[i1]

            i1 = np.in1d( gas_id2, gas_id1 )
            gas_xyz2 = gas_xyz2[i1,:]
            gas_hsml2= gas_hsml2[i1]
            gas_mass2= gas_mass2[i1]
            gas_u2   = gas_u2[i1]
            gas_nume2= gas_nume2[i1]
            gas_z2   = gas_z2[i1]

	    self.gas_xyz = (1.0 - interp_coef) * gas_xyz1   + interp_coef * gas_xyz2
	    self.gas_hsml= (1.0 - interp_coef) * gas_hsml1  + interp_coef * gas_hsml2
	    self.gas_mass= (1.0 - interp_coef) * gas_mass1  + interp_coef * gas_mass2
	    self.gas_u   = (1.0 - interp_coef) * gas_u1     + interp_coef * gas_u2
	    self.gas_nume= (1.0 - interp_coef) * gas_nume1  + interp_coef * gas_nume2
	    self.gas_z   = (1.0 - interp_coef) * gas_z1     + interp_coef * gas_z2





def dm_image(dir='./', snapnum=0, this_file=None, center=None, cosmo_wrap=False, 
                maxden=None, dynrange=None, **kwargs):
    snap_ext = str(snapnum).zfill(3)
    if this_file==None:
	print "setting snapshot"
        this_file = dir+'/snapshot_'+snap_ext
        if not os.path.exists(this_file+'.hdf5'):
            this_file = dir+'/snapdir_'+snap_ext+'/snapshot_'+snap_ext

#    this_file = this_file+'.hdf5'
    print this_file

    data = dm_image_data(this_file, **kwargs)

    try: 
        bh_xyz      = np.array(ws.read_block(this_file,'POS ', parttype=5))
    except:
        bh_xyz      = np.array( [0, 0, 0] )
    

    if cosmo_wrap:
        boxsize = np.max( data.xyz ) - np.min( data.xyz )

    if center != None:
        data.xyz[:,0] -= np.mean(center[0])
        data.xyz[:,1] -= np.mean(center[1])
        data.xyz[:,2] -= np.mean(center[2])
    else:
        data.xyz[:,0] -= np.mean(bh_xyz[:,0])
        data.xyz[:,1] -= np.mean(bh_xyz[:,1])
        data.xyz[:,2] -= np.mean(bh_xyz[:,2])

    if cosmo_wrap:
        data.xyz[ data.xyz[:,0] > boxsize/2.0 ,0] -= boxsize 
        data.xyz[ data.xyz[:,1] > boxsize/2.0 ,1] -= boxsize 
        data.xyz[ data.xyz[:,2] > boxsize/2.0 ,2] -= boxsize 

        data.xyz[ data.xyz[:,0] < -1.0* boxsize/2.0 ,0] += boxsize 
        data.xyz[ data.xyz[:,1] < -1.0* boxsize/2.0 ,1] += boxsize 
        data.xyz[ data.xyz[:,2] < -1.0* boxsize/2.0 ,2] += boxsize 

    
#    color_weights = units.springel_units.gas_code_to_temperature( gas_data.gas_u, gas_data.gas_nume)
    color_weights = 1e4
    weights = data.mass
    h_all   = data.hsml

    try:
      xrange = kwargs['xrange']
      yrange = kwargs['yrange']
      if xrange != None and yrange != None:
        fov = np.max( [xrange, yrange]  )
        ok = (data.xyz[:,0] > -1.1*fov) & (data.xyz[:,0] < 1.1*fov) & \
             (data.xyz[:,1] > -1.1*fov) & (data.xyz[:,1] < 1.1*fov) & \
             (data.xyz[:,2] > -1.1*fov) & (data.xyz[:,2] < 1.1*fov)  

        data.xyz = data.xyz[ok,:]
#	color_weights = color_weights[ok]
	weights = weights[ok]
	h_all  = h_all[ok]
    except:
	print "failed to clip"

    try:
      zrange = kwargs['zrange']
      if zrange != None:
        ok = (data.xyz[:,2] > np.min(zrange) ) & (data.xyz[:,2] < np.max(zrange) )

        data.xyz = data.xyz[ok,:]
        weights = weights[ok]
        h_all  = h_all[ok]
    except:
        print "failed to clip"
        sys.exit()

    
    ok = np.arange(h_all.shape[0])
    
    #weights = h_all 
    print np.min(weights), np.max(weights) 
    print data.xyz.shape
    print h_all.shape
    print weights.shape

    image, out_u, out_g, out_r = projection.raytrace_projection_compute(
			     data.xyz[:,0], data.xyz[:,1], data.xyz[:,2], h_all, weights, 
			     weights, weights, weights, 0, 0, 0, TRIM_PARTICLES=0,
			     **kwargs)
    print image.shape

    print "image max/min is  {:g}/{:g}".format(np.min(image), np.max(image) )
    if((maxden != None) and (dynrange != None)):
        image[image > maxden] = maxden
        image[image < maxden / dynrange ] = maxden / dynrange
        print "image clipping is {:g}/{:g}".format(maxden, maxden / dynrange )
    else:
        print "no clipping (maxden/dynrange not set)..."


#gas_raytrace_temperature( \
#                             gas_map_temperature_cuts, \
#                             gas_data.gas_xyz[ok,0], gas_data.gas_xyz[ok,1], gas_data.gas_xyz[ok,2], 
#			     color_weights[ok], weights[ok], h_all[ok], \
#                             kernel_width=kernel_widths,\
#                             isosurfaces = 1,\
#                             KAPPA_UNITS = 2.0885*np.array([1.1,2.0,1.5]),\
#                             add_temperature_weights = 0,\
#                             **kwargs)


#    image24, massmap =  makepic.make_threeband_image_process_bandmaps( out_r,out_g,out_u, \
        #                                                              **kwargs)
 
    print np.min(image), np.max(image)
    return np.log10(image)


