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
#import pfh_python.visualization.contour_makepic as makepic
import util.calc_hsml as calc_hsml
import util.cast
import matplotlib.pyplot as plt
import units.springel_units
import os.path

class gas_image_data:
    def __init__(self, this_file, next_file=None, interp_coef=0.0, interp=True, **kwargs):
        if next_file==None or interp==False:	# no interp
	    print "have a single file, no interp"
#	    self.head = ws.snapshot_header(this_file)
            self.gas_xyz     = np.array(ws.read_block(this_file,'POS ', parttype=0))
            self.gas_mass    = np.array(ws.read_block(this_file,'MASS', parttype=0))
            self.gas_u       = np.array(ws.read_block(this_file,'U   ', parttype=0))
	    try:
		vol = np.array(ws.read_block(this_file,'VOL ', parttype=0))
		self.gas_hsml = (vol * 3.0 / (4.0 * 3.14159 ))**0.333
	    except:
                print "failed to load gas hsml via volumes.  Trying hsml values directly."
                try:
                    self.gas_hsml    = np.array(ws.read_block(this_file,'HSML', parttype=0))
                except:
                    print "failed to load gas hsml via hsml (directly).  Loading via mass and density."
                    try:
                        self.gas_hsml = 3.0*(  3.0/(4.0*3.14159) * self.gas_mass / np.array( ws.read_block( this_file, 'RHO ', parttype=0)) ) **0.33333
                    except:
                        print "failed to load gas hsml values for a third time via mass and density.  stopping."
                        sys.exit()

            try:
                self.gas_sfr = np.array( ws.read_block( this_file, 'SFR ', parttype=0) )
            except:
                self.gas_sfr = np.zeros_like( gas_mass )
  

            try:
                self.gas_nume    = np.array(ws.read_block(this_file,'NE  ', parttype=0))
            except:
                print "failed to load gas nume."
                self.gas_nume    = np.ones_like( self.gas_mass )

            try:
                self.gas_z       = np.array(ws.read_block(this_file,'Z   ', parttype=0))
            except:
                try:
                    self.gas_z       = np.array(ws.read_block(this_file,'GZ  ', parttype=0))
                except:
                    print "failed to load metallicities.  Setting to solar."
                    self.gas_z       = np.ones_like( self.gas_mass ) * 0.02
	else:
#    	    head = ws.snapshot_header(this_file)

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




def simple_gas_image(dir='./', snapnum=0, this_file=None, center=None, dynrange=1e3, maxden=0.01, **kwargs):
    snap_ext = str(snapnum).zfill(3)
    if this_file==None:
        this_file = dir+'/snapshot_'+snap_ext
        if not os.path.exists(this_file+'.hdf5'):
            this_file = dir+'/snapdir_'+snap_ext+'/snapshot_'+snap_ext

def old_repo_gas_image(dir='./', snapnum=0, band_ids = [9,10,11], this_file=None, center=None, edge_on=False, cosmo_wrap=False, massmap=False, **kwargs):


    if this_file==None:
        snap_ext = str(snapnum).zfill(3)
        this_file = dir+'/snapshot_'+snap_ext
        if not os.path.exists(this_file+'.hdf5'):
            this_file = dir+'/snapdir_'+snap_ext+'/snapshot_'+snap_ext
    else:
        this_file = this_file+'.hdf5'

    print this_file
    data = gas_image_data(this_file, **kwargs)
    try:
        bh_xyz      = np.array(ws.read_block(this_file,'POS ', parttype=5))
    except:
        bh_xyz      = np.array( [0, 0, 0] )

    if center != None:
        data.xyz[:,0] -= np.mean(bh_xyz[:,0])
        data.xyz[:,1] -= np.mean(bh_xyz[:,1])
        data.xyz[:,2] -= np.mean(bh_xyz[:,2])
    weights = data.gas_mass
    h_all   = data.gas_hsml

    try:
      xrange = kwargs['xrange']
      yrange = kwargs['yrange']
    except:
      xrange = [-1e10, 1e10]
      yrange = [-1e10, 1e10]
    fov = np.max( [xrange, yrange]  )

    ok = (data.gas_xyz[:,0] > -1.1*fov) & (data.gas_xyz[:,0] < 1.1*fov) & \
             (data.gas_xyz[:,1] > -1.1*fov) & (data.gas_xyz[:,1] < 1.1*fov) & \
             (data.gas_xyz[:,2] > -1.1*fov) & (data.gas_xyz[:,2] < 1.1*fov)
    data.gas_xyz = data.gas_xyz[ok,:]
    weights = weights[ok]
    h_all  = h_all[ok]
    
    ok = np.arange(h_all.shape[0])
    image, out_u, out_g, out_r = projection.raytrace_projection_compute(
                             data.gas_xyz[:,0], data.gas_xyz[:,1], data.gas_xyz[:,2], h_all, weights,
                             weights, weights, weights, 0, 0, 0, TRIM_PARTICLES=0
                             )
    image[image > maxden] = maxden
    image[image < maxden / dynrange ] = maxden / dynrange
    return np.log10(image)






def sz_gas_image(dir='./', snapnum=0, this_file=None, center=None, dynrange=None, maxden=None, **kwargs):
    snap_ext = str(snapnum).zfill(3)
    if this_file==None:
        this_file = dir+'/snapshot_'+snap_ext
        if not os.path.exists(this_file+'.hdf5'):
            this_file = dir+'/snapdir_'+snap_ext+'/snapshot_'+snap_ext

    this_file = this_file+'.hdf5'
    print this_file
    data = gas_image_data(this_file, **kwargs)
    try:
        bh_xyz      = np.array(ws.read_block(this_file,'POS ', parttype=5))
    except:
        bh_xyz      = np.array( [0, 0, 0] )

    if center != None:
        data.xyz[:,0] -= np.mean(bh_xyz[:,0])
        data.xyz[:,1] -= np.mean(bh_xyz[:,1])
        data.xyz[:,2] -= np.mean(bh_xyz[:,2])
    
    effective_mass = data.gas_nume * data.gas_mass * data.gas_temperature

    weights = effective_mass
    h_all   = data.gas_hsml

    try:
      xrange = kwargs['xrange']
      yrange = kwargs['yrange']
    except:
      xrange = [-1e10, 1e10]
      yrange = [-1e10, 1e10]
    fov = np.max( [xrange, yrange]  )

    ok = (data.gas_xyz[:,0] > -1.1*fov) & (data.gas_xyz[:,0] < 1.1*fov) & \
             (data.gas_xyz[:,1] > -1.1*fov) & (data.gas_xyz[:,1] < 1.1*fov) & \
             (data.gas_xyz[:,2] > -1.1*fov) & (data.gas_xyz[:,2] < 1.1*fov)
    data.gas_xyz = data.gas_xyz[ok,:]
    weights = weights[ok]
    h_all  = h_all[ok]

    ok = np.arange(h_all.shape[0])
    image, out_u, out_g, out_r = projection.raytrace_projection_compute(
                             data.gas_xyz[:,0], data.gas_xyz[:,1], data.gas_xyz[:,2], h_all, weights,
                             weights, weights, weights, 0, 0, 0, TRIM_PARTICLES=0
                             )
    print "image max/min is  {:g}/{:g}".format(np.min(image), np.max(image) )
    if((maxden != None) and (dynrange != None)):
        image[image > maxden] = maxden
        image[image < maxden / dynrange ] = maxden / dynrange
        print "image clipping is {:g}/{:g}".format(maxden, maxden / dynrange )
    else:
        print "no clipping (maxden/dynrange not set)..."

    return np.log10(image)


def sfr_image(  sfr_image=True, massmap=False, **kwargs):
    image, massmap = gas_image( sfr_image=True, massmap=False, **kwargs )
    return image,massmap
 
def gas_image(dir='./', snapnum=0, band_ids = [9,10,11], this_file=None, center=None, edge_on=False, cosmo_wrap = False, \
                        massmap=False, projaxis=0, sfr_image=False, unit_length_mpc=False, include_lighting=True, **kwargs):
    snap_ext = str(snapnum).zfill(3)
    if this_file==None:
	print "setting snapshot"
        this_file = dir+'/snapshot_'+snap_ext
        if not os.path.exists(this_file+'.hdf5'):
            this_file = dir+'/snapdir_'+snap_ext+'/snapshot_'+snap_ext
    print this_file

    gas_data = gas_image_data(this_file, **kwargs)

    if cosmo_wrap:
        boxsize = np.max( gas_data.gas_xyz ) - np.min( gas_data.gas_xyz )
    
    if not( center is None ):
        gas_data.gas_xyz[:,0] -= np.mean(center[0])
        gas_data.gas_xyz[:,1] -= np.mean(center[1])
        gas_data.gas_xyz[:,2] -= np.mean(center[2])
    else:
        bh_xyz      = np.array(ws.read_block(this_file,'POS ', parttype=5))
        gas_data.gas_xyz[:,0] -= np.mean(bh_xyz[:,0])
        gas_data.gas_xyz[:,1] -= np.mean(bh_xyz[:,1])
        gas_data.gas_xyz[:,2] -= np.mean(bh_xyz[:,2])

    if cosmo_wrap:
        gas_data.gas_xyz[ gas_data.gas_xyz[:,0] > boxsize/2.0 ,0] -= boxsize
        gas_data.gas_xyz[ gas_data.gas_xyz[:,1] > boxsize/2.0 ,1] -= boxsize
        gas_data.gas_xyz[ gas_data.gas_xyz[:,2] > boxsize/2.0 ,2] -= boxsize

        gas_data.gas_xyz[ gas_data.gas_xyz[:,0] < -1.0* boxsize/2.0 ,0] += boxsize
        gas_data.gas_xyz[ gas_data.gas_xyz[:,1] < -1.0* boxsize/2.0 ,1] += boxsize
        gas_data.gas_xyz[ gas_data.gas_xyz[:,2] < -1.0* boxsize/2.0 ,2] += boxsize

    gas_map_temperature_cuts=np.array([300., 2.0e4, 3.0e5 ])
    kernel_widths=np.array([0.8,0.3,0.6])

    if unit_length_mpc:
        KAPPA_UNITS = 2.0885*np.array([1.1,2.0,1.5]) / 1000.0 / 1000.0
    else:
        KAPPA_UNITS = 2.0885*np.array([1.1,2.0,1.5])

    if massmap:
        color_weights = np.ones_like( gas_data.gas_mass )
    else:
        # this just returns temperature so that we can determine colors
        color_weights = units.springel_units.gas_code_to_temperature( gas_data.gas_u, gas_data.gas_nume)        

    if sfr_image:
        weights = gas_data.gas_sfr
    else:
        weights = gas_data.gas_mass

    h_all   = gas_data.gas_hsml

    try:
      xrange = kwargs['xrange']
      yrange = kwargs['yrange']
      if xrange != None and yrange != None:
        fov = np.max( [xrange, yrange]  )
        print np.min(gas_data.gas_xyz[:,0]), np.max(gas_data.gas_xyz[:,0])
        ok = (gas_data.gas_xyz[:,0] > -1.1*fov) & (gas_data.gas_xyz[:,0] < 1.1*fov) & \
             (gas_data.gas_xyz[:,1] > -1.1*fov) & (gas_data.gas_xyz[:,1] < 1.1*fov) & \
             (gas_data.gas_xyz[:,2] > -1.1*fov) & (gas_data.gas_xyz[:,2] < 1.1*fov)  

        print "Found xrange/yrange arguments in gas_images() call.  Clipping down to {:d} particles from {:d}".format( np.sum( ok ), len( ok ) )
        gas_data.gas_xyz = gas_data.gas_xyz[ok,:]
        color_weights = color_weights[ok]
        weights = weights[ok]
        h_all  = h_all[ok]
    except:
	print "failed to clip"

    try:
        ok = weights > 0
        if np.sum(ok) < len(ok):
            print "Found particles with negative weights in gas_images() call.  Clipping down to {:d} particles from {:d}".format( np.sum( ok ), len( ok ) )
            gas_data.gas_xyz = gas_data.gas_xyz[ok,:]
            color_weights = color_weights[ok]
            weights = weights[ok]
            h_all  = h_all[ok]
    except:
        print "failed to clip based on zero-weight values"
        sys.exit()

    if edge_on:
        x_index = 0
        y_index = 2
        z_index = 1
    elif projaxis==0:
        x_index = 0
        y_index = 1
	z_index = 2
    elif projaxis==1:
        x_index = 2
        y_index = 0
	z_index = 1
    elif projaxis==2:
        x_index = 1
        y_index = 2
	z_index = 0
    else:
        x_index = 0
        y_index = 1
	z_index = 2
    
    if massmap or sfr_image:
        image, out_u, out_g, out_r = projection.raytrace_projection_compute( gas_data.gas_xyz[:,x_index], gas_data.gas_xyz[:,y_index], gas_data.gas_xyz[:,z_index], \
                                        h_all[:], weights[:], weights[:], weights[:], weights[:], \
                                        0.0, 0.0, 0.0, TRIM_PARTICLES=0, \
                                        **kwargs )

        image24, massmap =  makepic.make_threeband_image_process_bandmaps( out_r,out_g,out_u, \
                                                                      **kwargs)
        return np.log10(image), np.log10(image)
#def raytrace_projection_compute( x, y, z, hsml, mass, wt1, wt2, wt3, \
#                                kappa_1, kappa_2, kappa_3, xrange=0, yrange=0, zrange=0, pixels=720, \
#                                TRIM_PARTICLES=1 ):
    else:
        print "executing projection routine."
        out_gas,out_u,out_g,out_r = projection.gas_raytrace_temperature( \
                             gas_map_temperature_cuts, \
                             gas_data.gas_xyz[:,x_index], gas_data.gas_xyz[:,y_index], gas_data.gas_xyz[:,z_index], 
                             color_weights[:], weights[:], h_all[:], \
                             kernel_width=kernel_widths,\
                             isosurfaces = 1,\
                             KAPPA_UNITS = KAPPA_UNITS,\
                             add_temperature_weights = 0,\
                             **kwargs)
        print "mapping u/g/r to three band image."
        image24, massmap =  makepic.make_threeband_image_process_bandmaps( out_r,out_g,out_u, **kwargs)
        print "layering three band image."
        image24 = makepic.layer_band_images(  image24, massmap  )   #   layer_band_images...

        if (include_lighting==1 and np.sum(image24[:,:,:2]) > 0 ):
            print "adding lighting"
            image24 = makepic.include_lighting( image24, massmap )


        if False:               #(include_lighting==1 and np.sum(image24[:,:,:2]) > 0 ):
            print "Lighting is being included!!! "
            light = viscolors.CustomLightSource(azdeg=0,altdeg=65)
            if (len(massmap.shape)>2):
                ## do some clipping to regulate the lighting:
                elevation = massmap.sum(axis=2)
                minden = maxden / dynrange

                print " "
                print elevation.min(), elevation.max(), elevation.mean()
                print minden, maxden
                print " "
                elevation = (elevation - minden) / (maxden - minden)
                elevation[elevation < 0.] = 0.
                elevation[elevation > 1.] = 1.
                elevation *= maxden
                grad_max = maxden / 5.
                grad_max = maxden / 6.
                #image24_lit = light.shade_rgb(image24, massmap.sum(axis=2))
                image24_lit = light.shade_rgb(image24, elevation, vmin=-grad_max, vmax=grad_max)
            else:
                image24_lit = light.shade(massmap, matplotlib.cm.get_cmap('hot'))   #reversed args          # ptorrey -- important
            image24 = image24_lit
        else:
            print "Lighting is not being included :( "

        return image24, massmap


