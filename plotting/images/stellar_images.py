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
import plotting.images.contour_makepic as makepic
import plotting.images.projection as projection
import plotting.colors
import util.calc_hsml as calc_hsml
import util.cast
import util.naming as naming
import matplotlib.pyplot as plt
import units.springel_units as units
import os

import simread.readsubfHDF5 as readsubf
class stellar_image_data:
    def __init__(self, this_file, next_file=None, interp_coef=0.0, interp=True, dust=True, **kwargs):
      if next_file==None or interp==False:    # no interp
        print "have a single file, no interp"
	print this_file

        self.head = ws.snapshot_header(this_file)
        if dust:
            print "loading dust data blocks"
            self.gas_xyz     = np.array(ws.read_block(this_file,'POS ', parttype=0))
            try:
                self.gas_hsml    = np.array(ws.read_block(this_file,'HSML', parttype=0))
            except:
                print "failed to load HSML values directly.  Try using cell volumes instead"
                rho  = np.array(ws.read_block(this_file,'RHO ', parttype=0))
                mass = np.array(ws.read_block(this_file,'MASS', parttype=0))
                vol = mass / rho
                self.gas_hsml    = (0.75 * vol / 3.14159 ) ** 0.33333 * 3 
            self.gas_mass    = np.array(ws.read_block(this_file,'MASS', parttype=0))
            self.gas_u       = np.array(ws.read_block(this_file,'U   ', parttype=0))
            self.gas_rho     = np.array(ws.read_block(this_file,'RHO ', parttype=0))
            self.gas_numh    = np.array(ws.read_block(this_file,'NH  ', parttype=0))
            self.gas_nume    = np.array(ws.read_block(this_file,'NE  ', parttype=0))
            try:
                self.gas_z       = np.array(ws.read_block(this_file,'Z   ', parttype=0))
            except:
                self.gas_z       = np.array(ws.read_block(this_file,'GZ  ', parttype=0))
        else:
            self.gas_xyz    = np.zeros( (10,3) )
            self.gas_hsml   = np.zeros( 10 )
            self.gas_mass   = np.zeros( 10 )
            self.gas_u      = np.zeros( 10 )
            self.gas_rho    = np.zeros( 10 )
            self.gas_numh   = np.zeros( 10 )
            self.gas_nume   = np.zeros( 10 )
            self.gas_z      = np.zeros( 10 )

        self.star_xyz    = np.array(ws.read_block(this_file,'POS ', parttype=4))
        self.star_mass   = np.array(ws.read_block(this_file,'MASS', parttype=4))
        try:
            self.star_bt     = np.array(ws.read_block(this_file,'AGE ',parttype=4))
            self.star_z      = np.array(ws.read_block(this_file,'Z   ',parttype=4))
        except:
            try:
                self.star_bt     = np.array(ws.read_block(this_file,'GAGE',parttype=4))
                self.star_z      = np.array(ws.read_block(this_file,'GZ  ',parttype=4))
            except:
                print "failed to load stars properly.  Setting to zero to avoid failure, but this will not produce an image."
                self.star_xyz    = np.zeros( (10,3) )
                self.star_bt     = np.zeros( (10) )
                self.star_z      = np.zeros( (10) )
                self.star_mass   = np.zeros( (10)  )



class tng_stellar_image_data:
    def __init__(self, this_file, next_file=None, interp_coef=0.0, interp=True, dust=True, fof_num=-1, sub_num=-1,  **kwargs):
        self.head = ws.snapshot_header(this_file)
        if dust:
            print "loading dust data blocks"
            self.gas_xyz     = np.array(ws.read_block(this_file,'POS ', parttype=0))
            try:
                self.gas_hsml    = np.array(ws.read_block(this_file,'HSML', parttype=0))
            except:
                print "failed to load HSML values directly.  Try using cell volumes instead"
                rho  = np.array(ws.read_block(this_file,'RHO ', parttype=0))
                mass = np.array(ws.read_block(this_file,'MASS', parttype=0))
                vol = mass / rho
                self.gas_hsml    = (0.75 * vol / 3.14159 ) ** 0.33333 * 3
            self.gas_mass    = np.array(ws.read_block(this_file,'MASS', parttype=0))
            self.gas_u       = np.array(ws.read_block(this_file,'U   ', parttype=0))
            self.gas_rho     = np.array(ws.read_block(this_file,'RHO ', parttype=0))
            self.gas_numh    = np.array(ws.read_block(this_file,'NH  ', parttype=0))
            self.gas_nume    = np.array(ws.read_block(this_file,'NE  ', parttype=0))
            try:
                self.gas_z       = np.array(ws.read_block(this_file,'Z   ', parttype=0))
            except:
                self.gas_z       = np.array(ws.read_block(this_file,'GZ  ', parttype=0))
        else:
            self.gas_xyz    = np.zeros( (10,3) )
            self.gas_hsml   = np.zeros( 10 )
            self.gas_mass   = np.zeros( 10 )
            self.gas_u      = np.zeros( 10 )
            self.gas_rho    = np.zeros( 10 )
            self.gas_numh   = np.zeros( 10 )
            self.gas_nume   = np.zeros( 10 )
            self.gas_z      = np.zeros( 10 )

        import simread.readhaloHDF5 as hr
        print this_file
        snapdir = this_file[:this_file.index('output')+7]

        run=this_file[this_file.index('Runs/')+5:this_file.index('output') ].replace('/','')
        print run
        snapnum = int(this_file[-3:])
        halo_reader = hr.HaloReader( snapdir, snapnum, run=run )
        self.star_xyz    = halo_reader.read( 'POS ', 4, fof_num, sub_num )
        self.star_bt     = halo_reader.read( 'GAGE', 4, fof_num, sub_num )
        self.star_z      = halo_reader.read( 'GZ  ', 4, fof_num, sub_num )
        self.star_mass   = halo_reader.read( 'MASS', 4, fof_num, sub_num )

        print "loaded stellar data from "+this_file




def stellar_image(dir='./', snapnum=0, band_ids = [9,10,11], 
                  cosmo=False, 
                  this_file=None, center=None,
                  snapbase='snapshot',
                  illustris_tng=False,
                  tng_subnr=None,
                  tng_groupnr=None,
                  cosmo_wrap=False,
                  massmap=False,
                  projaxis=0,
                  **kwargs):

    if this_file==None:
        snap_ext = "000"+str(snapnum)
        snap_ext = snap_ext[-3:]
        this_file = dir+'/'+snapbase+'_'+snap_ext
        if not os.path.exists(this_file+'.hdf5'):
            this_file = dir+'/snapdir_'+snap_ext+'/'+snapbase+'_'+snap_ext

        print "setting snapshot to {:s}".format( this_file )

    if illustris_tng==True:
        print "Making a TNG stellar image."
        star_data = tng_stellar_image_data(this_file, fof_num=tng_groupnr, sub_num=tng_subnr, **kwargs)
    else:
        star_data = stellar_image_data(this_file, **kwargs)

    if cosmo_wrap:
        boxsize = star_data.head.boxsize      #np.max( star_data.gas_xyz ) - np.min( star_data.gas_xyz )

    if illustris_tng and isinstance( tng_subnr, int ):
        image_center = readsubf.subfind_catalog(  dir, snapnum, subcat=True, grpcat=False,
                       keysel=['SubhaloPos'] )
        center = image_center.SubhaloPos[tng_subnr, :] 
    elif illustris_tng and isinstance( tng_groupnr, int ):
        image_center = readsubf.subfind_catalog(  dir, snapnum, subcat=True, grpcat=True,
                       keysel=['SubhaloPos', 'GroupFirstSub'] )
        center = image_center.SubhaloPos[ image_center.GroupFirstSub[tng_groupnr], :]
    elif center is not None:
        center = center
    else:
        try:
            bh_id      = np.array(ws.read_block(this_file,'ID  ', parttype=5))
            bh_xyz      = np.array(ws.read_block(this_file,'POS ', parttype=5))
            center = [ np.mean(bh_xyz[0,0]), np.mean(bh_xyz[0,1]), np.mean(bh_xyz[0,2]) ]
        except:
            center = [0.0, 0.0, 0.0 ]

    star_data.gas_xyz[:,0] -= np.mean(center[0])
    star_data.gas_xyz[:,1] -= np.mean(center[1])
    star_data.gas_xyz[:,2] -= np.mean(center[2])
    star_data.star_xyz[:,0] -= np.mean(center[0])
    star_data.star_xyz[:,1] -= np.mean(center[1])
    star_data.star_xyz[:,2] -= np.mean(center[2])

    if cosmo_wrap:
        star_data.star_xyz[ star_data.star_xyz[:,0] > boxsize/2.0 ,0] -= boxsize
        star_data.star_xyz[ star_data.star_xyz[:,1] > boxsize/2.0 ,1] -= boxsize
        star_data.star_xyz[ star_data.star_xyz[:,2] > boxsize/2.0 ,2] -= boxsize
        star_data.star_xyz[ star_data.star_xyz[:,0] < -1.0* boxsize/2.0 ,0] += boxsize
        star_data.star_xyz[ star_data.star_xyz[:,1] < -1.0* boxsize/2.0 ,1] += boxsize
        star_data.star_xyz[ star_data.star_xyz[:,2] < -1.0* boxsize/2.0 ,2] += boxsize

        star_data.gas_xyz[ star_data.gas_xyz[:,0] > boxsize/2.0 ,0] -= boxsize
        star_data.gas_xyz[ star_data.gas_xyz[:,1] > boxsize/2.0 ,1] -= boxsize
        star_data.gas_xyz[ star_data.gas_xyz[:,2] > boxsize/2.0 ,2] -= boxsize
        star_data.gas_xyz[ star_data.gas_xyz[:,0] < -1.0* boxsize/2.0 ,0] += boxsize
        star_data.gas_xyz[ star_data.gas_xyz[:,1] < -1.0* boxsize/2.0 ,1] += boxsize
        star_data.gas_xyz[ star_data.gas_xyz[:,2] < -1.0* boxsize/2.0 ,2] += boxsize


    try:
        xrange = kwargs['xrange']
        yrange = kwargs['yrange']
        if xrange != None and yrange != None:
            fov = np.max( [xrange, yrange]  )
            ok = (star_data.star_xyz[:,0] > -1.1*fov) & (star_data.star_xyz[:,0] < 1.1*fov) & \
                 (star_data.star_xyz[:,1] > -1.1*fov) & (star_data.star_xyz[:,1] < 1.1*fov) & \
                 (star_data.star_xyz[:,2] > -1.1*fov) & (star_data.star_xyz[:,2] < 1.1*fov) & \
                 (star_data.star_bt > 0)

        if np.sum(ok) > 10:
              star_data.star_xyz = star_data.star_xyz[ok,:]
              star_data.star_bt  = star_data.star_bt[ok]
              star_data.star_z   = star_data.star_z[ok]
              star_data.star_mass= star_data.star_mass[ok]
        else:
              print "WARN:  I found 10 or less particles in the FOV.  Expect an error."
    except:
        print "Didnt find xrange/yrange values.  Not clipping.  Using some default range..."


    if cosmo:
	star_age = np.array([ units.age_from_a(star_data.star_bt[index], a0=star_data.head.time) for index in np.arange(star_data.star_bt.shape[0])  ])
    else:
        star_age = star_data.head.time - star_data.star_bt

    star_hsml   = calc_hsml.get_particle_hsml( star_data.star_xyz[:,0], star_data.star_xyz[:,1], star_data.star_xyz[:,2] , **kwargs ) 
    print "min/max of star_hslm : {:f}|{:f}".format(np.min(star_hsml), np.max(star_hsml)) 

    if projaxis==0:
        x_axis_index = 0
        y_axis_index = 1
        z_axis_index = 2
    if projaxis==1:
        x_axis_index = 2
        y_axis_index = 0
        z_axis_index = 1
    if projaxis==2:
        x_axis_index = 1
        y_axis_index = 2
        z_axis_index = 0

    print x_axis_index, y_axis_index, z_axis_index


    # This calls "stellar raytrace" which itself is a wrapper for "raytrace projection compute".
    # out_gas should have units of mass/area; out_u/g/r have funky units, but are ~luminosity/area
    # NO CLIPPING TO THIS POINT!
    out_gas,out_u,out_g,out_r = projection.stellar_raytrace(    star_data.star_xyz[:,x_axis_index], star_data.star_xyz[:,y_axis_index], star_data.star_xyz[:,z_axis_index], \
                                                                star_data.star_mass[:], star_age[:], star_data.star_z[:], star_hsml[:], \
                                                                star_data.gas_xyz[:,x_axis_index], star_data.gas_xyz[:,y_axis_index], star_data.gas_xyz[:,z_axis_index], star_data.gas_mass, \
                                                                star_data.gas_z, star_data.gas_hsml, \
                                                                band_ids=band_ids, \
                                                                **kwargs)

    print "within stellar_image the u,g,r images have sizes and min/max values of:"
    print out_gas.shape
    print out_u.shape
    print out_g.shape
    print out_r.shape
    print np.min(out_u), np.max(out_u)
    print np.min(out_g), np.max(out_g)
    print np.min(out_r), np.max(out_r)


    
    if(np.array(out_gas).size<=1):
        out_gas=out_u=out_g=out_r=np.zeros((pixels,pixels))
        image24=massmap=np.zeros((pixels,pixels,3))
        print "SIZE OF OUTPUT IMAGE IS ZERO.  RETURNING IMAGES WITH JUST ZEROS."
    else:
        ## make the resulting maps into an image
        #  THIS IS WHERE CLIPPING SHOULD HAPPEN USING maxden and dynrange
        image24, massmap = makepic.make_threeband_image_process_bandmaps( out_r,out_g,out_u, **kwargs)
        print "Making a three-band image from the R-G-B maps."

    print "within stellar_image image24 has a size:"
    print image24.shape
    print "and a min/max in the first channel of"
    print np.min(image24[:,:,0]), np.max(image24[:,:,0])



    include_lighting=0
    if (include_lighting==1):
        #light = matplotlib.colors.LightSource(azdeg=0,altdeg=65)
        light = viscolors.CustomLightSource(azdeg=0,altdeg=65)
        if (len(massmap.shape)>2):
            ## do some clipping to regulate the lighting:
            elevation = massmap.sum(axis=2)
            minden = maxden / dynrange
            elevation = (elevation - minden) / (maxden - minden)
            elevation[elevation < 0.] = 0.
            elevation[elevation > 1.] = 1.
            elevation *= maxden
            grad_max = maxden / 5.
            grad_max = maxden / 6.
            #image24_lit = light.shade_rgb(image24, massmap.sum(axis=2))
            image24_lit = light.shade_rgb(image24, elevation, vmin=-grad_max, vmax=grad_max)
        else:
            image24_lit = light.shade(image24, massmap)
        image24 = image24_lit


    return image24

#    fig = plt.figure(figsize=(5,5))
#    ax = fig.add_subplot(1,1,1)
#    ax.imshow(image24,origin='lower',interpolation='bicubic',aspect='normal');
#    fig.savefig("./plots/dummy.png")








#####  this is old legacy extra gas layer code for stellar images:


#    if False:       #((add_gas_layer_to_image=='')==False):
#            gas_wt = 0.*gx
#            maxden_for_layer=1.*maxden; dynrange_for_layer=1.*dynrange
#            set_percent_maxden_layer=0.; set_percent_minden_layer=0.;
#            gas_hsml_for_extra_layer=gas_hsml
#            if (add_gas_layer_to_image=='Halpha'):
#                    gas_wt = gas_mass * gas_rho*gas_nume*(1.-gas_numh)*(gadget.gas_temperature(gas_u,gas_nume)**(-0.75))/(gadget.gas_mu(gas_nume)**2.)
#                    maxden_for_layer *= 3.e4;
#                    #dynrange_for_layer *= 1.5;
#                    dynrange_for_layer *= 0.7;
#                    dynrange_for_layer *= 2.0;
#            #gas_hsml_for_extra_layer *= 1.1;
#            if (add_gas_layer_to_image=='CO'):
#                    gas_tmp = gadget.gas_temperature(gas_u,gas_nume)
#                    n_tmp = gas_rho * 176.2; # assuming mean molec weight of 2.3 for dense gas
#                    gas_wt = gas_mass * np.exp(-(gas_tmp/8000. + 10./n_tmp));
#                    maxden_for_layer *= 3.e4;
#                    dynrange_for_layer *= 0.7e4;
#            if (add_gas_layer_to_image=='SFR'):
#                    gas_wt = gas_sfr
#                    set_percent_maxden_layer=0.9999;
#                    set_percent_minden_layer=0.01;
#            if (add_gas_layer_to_image=='Xray'):
#                    gas_wt = gas_lxray
#                    maxden_for_layer *= 0.6e3;
#                    dynrange_for_layer *= 0.1;
#                    gas_hsml_for_extra_layer *= 1.2;
#            if (add_gas_layer_to_image=='Zmetal'):
#                    gas_wt = gas_mass*gas_metallicity
#                    set_percent_maxden_layer=0.9999;
#                    set_percent_minden_layer=0.01;
#            gas_wt /= np.sum(gas_wt)
#
#
#            massmap_gas_extra_layer,image_singledepth_extra_layer = cmakepic.simple_makepic(gx,gy,weights=gas_wt,hsml=gas_hsml_for_extra_layer,\
#                            xrange=xr,yrange=yr,
#                            set_dynrng=dynrange_for_layer,set_maxden=maxden_for_layer,
#                            set_percent_maxden=set_percent_maxden_layer,set_percent_minden=set_percent_minden_layer,
#                            color_temperature=0,pixels=pixels,invert_colorscale=1-invert_colors);
