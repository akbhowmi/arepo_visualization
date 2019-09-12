import matplotlib as mpl
mpl.use('Agg')
###import pyfits, glob, pylab, numpy, copy as cpy, stats, vecops, types, h5py, re
import pyfits, glob, h5py, re, sys, os, getopt
from copy import copy
import pj_constants as pjc
import astro_constants as ac
#from pylab import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.collections as mcoll
import wise_agn as wise
import readsnapbhs
import calc_merger_phases as cmph
import colormaps as cmaps
from matplotlib.font_manager import FontProperties
import PyCos as pc

abmag_zp = 3.631e-23 ## W/m^2/Hz                  

## Jarrett et al. 2011
W1_Vega_to_AB = 2.699
W2_Vega_to_AB = 3.339
W3_Vega_to_AB = 5.174
W4_Vega_to_AB = 6.620

    

class agn_metadata:

    def __init__(self, *args, **kwargs):
        if len(args)!=1:
            print "Error: class agn_metadata takes exactly 1 argument (filepath)."
            sys.stdout.flush()
            sys.exit()

        fpath = args[0]
        if not os.path.exists(fpath):
            print "Error: file %s not found."%fpath
            sys.stdout.flush()
            sys.exit()

        skip_snap0 = kwargs.get('skip_snap0',True)
        self.ncam = kwargs.get('ncam',7)
        tmin = kwargs.get('tmin',0.0)        
        tmax = kwargs.get('tmax',15.0)
        tmax_postmrg = kwargs.get('tmax_postmrg',-1.0)
        agnx0 = kwargs.get('agnx0',False)
        grid_nh_res = kwargs.get('grid_nh_res',0.0)
        skip_nh = kwargs.get('skip_nh',False)

        ### load main bh data ###
        tmpdata = load_infofile(fpath,fbase='bh_info',skip_snap0=skip_snap0)
        self.snap = tmpdata[0]
        self.time = tmpdata[1]
        self.posbh1 = tmpdata[2:5]
        self.mdotbh1 = tmpdata[5]
        self.lagn1 = tmpdata[6] if not agnx0 else tmpdata[6]*0.0
        self.lagn_sunsed1 = tmpdata[7] if not agnx0 else tmpdata[7]*0.0
        self.mbh1 = tmpdata[8]
        self.posbh2 = tmpdata[9:12]
        self.mdotbh2 = tmpdata[12]
        self.lagn2 = tmpdata[13] if not agnx0 else tmpdata[13]*0.0
        self.lagn_sunsed2 = tmpdata[14] if not agnx0 else tmpdata[14]*0.0
        self.mbh2 = tmpdata[15]
        self.nbh = np.ones((self.time.size))
        self.nbh[self.lagn2==self.lagn2] = 2
        self.tmrg = self.time[self.nbh==1].min() if self.snap[self.nbh==1].size<self.snap.size else -1

        ### load bh sep data ###
        ssnap,stime,self.bhsep,self.bhprojsep = load_infofile(fpath,fbase='bh_sep',
                                                              skip_snap0=skip_snap0)
        self.compare_arr_size(self.snap, ssnap, astr='bh_info', bstr='bh_sep')
        if self.tmrg>0: self.compare_arr_size(self.bhsep[self.bhsep==self.bhsep], 
                                              self.nbh[self.nbh==2], 
                                              astr='bh_info (2bh)', bstr='bh_sep (2bh)')

        ### load star & gas data ###
        tmpdata = load_infofile(fpath,fbase='star_gas_info',skip_snap0=skip_snap0)
        sgsnap,sgtime,self.mstar,self.mgas,self.mmetals,self.sfr = tmpdata
        self.compare_arr_size(self.snap, sgsnap, astr='bh_info', bstr='star_gas_info')

        ### load luminosity data ###
        tmpdata = load_infofile(fpath,fbase='lum_info',skip_snap0=skip_snap0)
        lsnap,ltime,self.lbol_grid,self.lbol_absorbed,self.ltot_out,self.lir = tmpdata
        self.compare_arr_size(self.snap, lsnap, astr='bh_info', bstr='lum_info')
        
        if not skip_nh:
            ### load column densities ###
            tmp=np.zeros((self.ncam,self.time.size))-1
            bh1_xpix = copy(tmp)
            bh1_ypix = copy(tmp)
            bh1_Ngas_tot_aux = copy(tmp)
            bh1_Ngas_tot_grid = copy(tmp)
            bh1_Ngas_grid = copy(tmp)
            bh2_xpix = copy(tmp)
            bh2_ypix = copy(tmp)
            bh2_Ngas_tot_aux = copy(tmp)
            bh2_Ngas_tot_grid = copy(tmp)
            bh2_Ngas_grid = copy(tmp)
            for c in range(self.ncam):
                rstr='_grid_res_%.3d'%(1000*grid_nh_res) if grid_nh_res>0 else ''
                with open(fpath+'/Ngas_bh_cam%d%s.txt'%(c,rstr)) as f:
                    tmpdata = np.loadtxt(f,unpack=True,dtype='float64')
                    if skip_snap0:
                        sn = tmpdata[0,:]
                        snapmask = (sn>0)
                        tmpdata = tmpdata[:,snapmask]
                    if tmpdata.shape[0]==7:
                        #sn,bh1_xpix[c,:],bh1_ypix[c,:],bh1_Ngas[c,:],bh2_xpix[c,:],bh2_ypix[c,:],bh2_Ngas[c,:]=tmpdata
                        sn,bh1_xpix[c,:],bh1_ypix[c,:],bh1_Ngas_tot_aux[c,:],bh2_xpix[c,:],bh2_ypix[c,:],bh2_Ngas_tot_aux[c,:]=tmpdata
                    elif tmpdata.shape[0]==11:
                        sn,bh1_xpix[c,:],bh1_ypix[c,:],bh1_Ngas_tot_aux[c,:],bh1_Ngas_tot_grid[c,:],bh1_Ngas_grid[c,:],bh2_xpix[c,:],bh2_ypix[c,:],bh2_Ngas_tot_aux[c,:],bh2_Ngas_tot_grid[c,:],bh2_Ngas_grid[c,:]=tmpdata
                    else: 
                        print "Error: invalid format in file fpath+'/Ngas_bh_cam%d%s.txt'%(c,rstr)"
                        print tmpdata.shape
                        sys.stdout.flush()
                        sys.exit()

                    self.compare_arr_size(self.snap, sn, astr='bh_info', bstr='Ngas_bh_cam%d%s'%(c,rstr))

            self.pixbh1 = np.array([bh1_xpix, bh1_ypix])
            self.pixbh2 = np.array([bh2_xpix, bh2_ypix])
            #print "pixbh1,pixbh2 shape:",self.pixbh1.shape,self.pixbh2.shape
            self.Ngas_bh1_aux = bh1_Ngas_tot_aux/2.0
            self.Ngas_bh1_aux[(bh1_xpix<0)|(bh1_ypix<0)]=np.nan
            self.Ngas_bh2_aux = bh2_Ngas_tot_aux/2.0
            self.Ngas_bh2_aux[(bh2_xpix<0)|(bh2_ypix<0)]=np.nan
            self.Ngas_bh1_grid = bh1_Ngas_grid
            self.Ngas_bh2_grid = bh2_Ngas_grid

        ### load wise colors ###
        wtime, wsnap, self.w1w2, self.w2w3 = load_wise(maindir=fpath,subdir='',ncam=self.ncam,
                                                       skip_snap0=skip_snap0)
        self.compare_arr_size(self.snap, wsnap, astr='bh_info', bstr='wise data')

        ### load wise lums ###
        w1snap, w1time, self.w1lum, self.w2lum = load_wise_lums(path=fpath,ncam=self.ncam)
        ### NOTE: for now, using monochromatic fluxes/luminosities, and leaving in SUNRISE units: W/m
        ##self.w1lum = self.w1lum * 1.0e7 ## convert W to erg/s
        ##self.w2lum = self.w2lum * 1.0e7 ## convert W to erg/s
        
        ### apply time cuts ###
        if tmax_postmrg >= 0 and self.tmrg >= 0:
            mask = ((self.time>=tmin)&(self.time<=tmax)&(self.time<=self.tmrg+tmax_postmrg))
        else:
            mask = ((self.time>=tmin)&(self.time<=tmax))
            if tmax_postmrg >=0: print "Warning: No BH merger found. Ignoring keyword tmax_postmrg."
        if self.snap[mask].size < self.snap.size:
            self.snap = self.snap[mask]
            self.time = self.time[mask]
            self.posbh1 = self.posbh1[:,mask]
            self.mdotbh1 = self.mdotbh1[mask]
            self.lagn1 = self.lagn1[mask]
            self.lagn_sunsed1 = self.lagn_sunsed1[mask]
            self.mbh1 = self.mbh1[mask]
            self.posbh2 = self.posbh2[:,mask]
            self.mdotbh2 = self.mdotbh2[mask]
            self.lagn2 = self.lagn2[mask]
            self.lagn_sunsed2 = self.lagn_sunsed2[mask]
            self.mbh2 = self.mbh2[mask]
            self.nbh = self.nbh[mask]
            self.bhsep = self.bhsep[mask]
            self.bhprojsep = self.bhprojsep[:,mask]
            self.mstar = self.mstar[mask]
            self.mgas = self.mgas[mask]
            self.mmetals = self.mmetals[mask]
            self.sfr = self.sfr[mask]
            self.lbol_grid = self.lbol_grid[mask]
            self.lbol_absorbed = self.lbol_absorbed[mask]
            self.ltot_out = self.ltot_out[:,mask]
            self.lir = self.lir[:,mask]
            if not skip_nh:
                self.pixbh1 = self.pixbh1[:,:,mask]
                self.pixbh2 = self.pixbh2[:,:,mask]
                self.Ngas_bh1_aux = self.Ngas_bh1_aux[:,mask]
                self.Ngas_bh2_aux = self.Ngas_bh2_aux[:,mask]
                self.Ngas_bh1_grid = self.Ngas_bh1_grid[:,mask]
                self.Ngas_bh2_grid = self.Ngas_bh2_grid[:,mask]
            self.w1w2 = self.w1w2[:,mask]
            self.w2w3 = self.w2w3[:,mask]
            self.w1lum = self.w1lum[:,mask]
            self.w2lum = self.w2lum[:,mask]
            
        self.dt = np.zeros(self.time.size)
        self.dt[:-1] = self.time[1:]-self.time[:-1]
        self.dt[-1] = 0.0

        self.lagn = copy(self.lagn1)
        self.lagn[self.nbh==2] = self.lagn1[self.nbh==2] + self.lagn2[self.nbh==2]
        self.lagn_sun = copy(self.lagn_sunsed1)
        self.lagn_sun[self.nbh==2] = self.lagn_sunsed1[self.nbh==2] + self.lagn_sunsed2[self.nbh==2]
        self.mbh = copy(self.mbh1)
        self.mbh[self.nbh==2] = self.mbh1[self.nbh==2] + self.mbh2[self.nbh==2]
    
        self.ledd = calc_lbol_edd(self.mbh)
        self.fedd = np.zeros(self.lagn.size)
        self.fedd[self.ledd>0] = self.lagn[self.ledd>0]/self.ledd[self.ledd>0]
        
        self.ledd1 = calc_lbol_edd(self.mbh1)
        self.ledd2 = calc_lbol_edd(self.mbh2)
        self.fedd1 = self.lagn1/self.ledd1
        self.fedd2 = np.zeros(self.lagn.size)
        self.fedd2[self.nbh==2] = self.lagn2[self.nbh==2]/self.ledd2[self.nbh==2]

        self.fagn_lbol_out = self.lagn_sun/self.ltot_out
        self.fagn1_lbol_out = self.lagn_sunsed1/self.ltot_out
        self.fagn2_lbol_out = self.lagn_sunsed2/self.ltot_out
        self.fagn_lbol = self.lagn_sun/self.lbol_grid
        self.fagn1_lbol = self.lagn_sunsed1/self.lbol_grid
        self.fagn2_lbol = self.lagn_sunsed2/self.lbol_grid

        self.bhprojsep[self.bhprojsep!=self.bhprojsep] = -1
        self.bhsep[self.bhsep!=self.bhsep] = -1
        self.lagn2[self.lagn2!=self.lagn2] = -1
        self.lagn_sunsed2[self.lagn_sunsed2!=self.lagn_sunsed2] = -1
        self.fedd2[self.fedd2!=self.fedd2] = -1
        self.fagn2_lbol[self.fagn2_lbol!=self.fagn2_lbol] = -1 

        self.nsnaps = self.snap.size

    def compare_arr_size(self, a, b, astr='a', bstr='b'):        
        if a.size != b.size:
            print "Error: array size mismatch: %s, %s:"%(astr,bstr)
            print "  %d, %d"%(a.size,b.size)
            sys.stdout.flush()
            sys.exit()

    def mask_snaps(self,mask,skip_nh=False):

        assert (isinstance(mask,list) or isinstance(mask,np.ndarray)),\
            'Error: mask snaps requires a list or numpy array as input.'
        self.snap = self.snap[mask]
        self.time = self.time[mask]
        self.posbh1 = self.posbh1[:,mask]
        self.mdotbh1 = self.mdotbh1[mask]
        self.lagn1 = self.lagn1[mask]
        self.lagn_sunsed1 = self.lagn_sunsed1[mask]
        self.mbh1 = self.mbh1[mask]
        self.posbh2 = self.posbh2[:,mask]
        self.mdotbh2 = self.mdotbh2[mask]
        self.lagn2 = self.lagn2[mask]
        self.lagn_sunsed2 = self.lagn_sunsed2[mask]
        self.mbh2 = self.mbh2[mask]
        self.nbh = self.nbh[mask]
        self.bhsep = self.bhsep[mask]
        self.bhprojsep = self.bhprojsep[:,mask]
        self.mstar = self.mstar[mask]
        self.mgas = self.mgas[mask]
        self.mmetals = self.mmetals[mask]
        self.sfr = self.sfr[mask]
        self.lbol_grid = self.lbol_grid[mask]
        self.lbol_absorbed = self.lbol_absorbed[mask]
        self.ltot_out = self.ltot_out[:,mask]
        self.lir = self.lir[:,mask]
        if not skip_nh:
            self.pixbh1 = self.pixbh1[:,:,mask]
            self.pixbh2 = self.pixbh2[:,:,mask]
            self.Ngas_bh1_aux = self.Ngas_bh1_aux[:,mask]
            self.Ngas_bh2_aux = self.Ngas_bh2_aux[:,mask]
            self.Ngas_bh1_grid = self.Ngas_bh1_grid[:,mask]
            self.Ngas_bh2_grid = self.Ngas_bh2_grid[:,mask]
        self.w1w2 = self.w1w2[:,mask]
        self.w2w3 = self.w2w3[:,mask]
        self.w1lum = self.w1lum[:,mask]
        self.w2lum = self.w2lum[:,mask]

        self.dt = self.dt[mask]
        self.lagn = self.lagn[mask]
        self.lagn_sun = self.lagn_sun[mask]
        self.mbh = self.mbh[mask]    
        self.ledd = self.ledd[mask]
        self.fedd = self.fedd[mask]
        self.ledd1 = self.ledd1[mask]
        self.ledd2 = self.ledd2[mask]
        self.fedd1 = self.fedd1[mask]
        self.fedd2 = self.fedd2[mask]
        self.fagn_lbol_out = self.fagn_lbol_out[:,mask]
        self.fagn1_lbol_out = self.fagn1_lbol_out[:,mask]
        self.fagn2_lbol_out = self.fagn2_lbol_out[:,mask]
        self.fagn_lbol = self.fagn_lbol[mask]
        self.fagn1_lbol = self.fagn1_lbol[mask]
        self.fagn2_lbol = self.fagn2_lbol[mask]


class agn_metadata_multisim:

    def __init__(self, *args, **kwargs):
        if len(args)!=2:
            print "Error: class agn_metadata takes exactly 2 arguments (maindir, subdir_arr)."
            sys.stdout.flush()
            sys.exit()

        maindir = args[0]
        subdir_arr = args[1]
        nsim = len(subdir_arr)

        self.ncam = kwargs.get('ncam',7)        
        tmin = kwargs.get('tmin',0.0)        
        tmax = kwargs.get('tmax',15.0)
        tmax_postmrg = kwargs.get('tmax_postmrg',-1.0)
        agnx0 = kwargs.get('agnx0',False)
        flat = kwargs.get('flat',True)

        #self.posbh2 = np.empty((3,0))
        #self.w2w3 = np.empty((self.ncam,0))
        #self.dt = np.empty((0))

        self.snap = self.init_arr(nsim=nsim, flat=flat)
        self.time = self.init_arr(nsim=nsim, flat=flat)
        self.posbh1 = self.init_arr(len0=3, nsim=nsim, flat=flat)
        self.mdotbh1 = self.init_arr(nsim=nsim, flat=flat)
        self.lagn1 = self.init_arr(nsim=nsim, flat=flat)
        self.lagn_sunsed1 = self.init_arr(nsim=nsim, flat=flat)
        self.mbh1 = self.init_arr(nsim=nsim, flat=flat)
        self.posbh2 = self.init_arr(len0=3, nsim=nsim, flat=flat)
        self.mdotbh2 = self.init_arr(nsim=nsim, flat=flat)
        self.lagn2 = self.init_arr(nsim=nsim, flat=flat)
        self.lagn_sunsed2 = self.init_arr(nsim=nsim, flat=flat)
        self.mbh2 = self.init_arr(nsim=nsim, flat=flat)
        self.nbh = self.init_arr(nsim=nsim, flat=flat)
        self.bhsep = self.init_arr(nsim=nsim, flat=flat)
        self.bhprojsep = self.init_arr(len0=self.ncam, nsim=nsim, flat=flat)

        self.mstar = self.init_arr(nsim=nsim, flat=flat)
        self.mgas = self.init_arr(nsim=nsim, flat=flat)
        self.mmetals = self.init_arr(nsim=nsim, flat=flat)
        self.sfr = self.init_arr(nsim=nsim, flat=flat)
        self.lbol_grid = self.init_arr(nsim=nsim, flat=flat)
        self.lbol_absorbed = self.init_arr(nsim=nsim, flat=flat)
        self.ltot_out = self.init_arr(len0=self.ncam, nsim=nsim, flat=flat)
        self.lir = self.init_arr(len0=self.ncam, nsim=nsim, flat=flat)
        self.w1w2 = self.init_arr(len0=self.ncam, nsim=nsim, flat=flat)
        self.w2w3 = self.init_arr(len0=self.ncam, nsim=nsim, flat=flat)
        self.w1lum = self.init_arr(len0=self.ncam, nsim=nsim, flat=flat)
        self.w2lum = self.init_arr(len0=self.ncam, nsim=nsim, flat=flat)
        self.dt = self.init_arr(nsim=nsim, flat=flat)

        self.lagn = self.init_arr(nsim=nsim, flat=flat)
        self.lagn_sun = self.init_arr(nsim=nsim, flat=flat)
        self.mbh = self.init_arr(nsim=nsim, flat=flat)
        self.ledd = self.init_arr(nsim=nsim, flat=flat)
        self.fedd = self.init_arr(nsim=nsim, flat=flat)
        self.ledd1 = self.init_arr(nsim=nsim, flat=flat)
        self.ledd2 = self.init_arr(nsim=nsim, flat=flat)
        self.fedd1 = self.init_arr(nsim=nsim, flat=flat)
        self.fedd2 = self.init_arr(nsim=nsim, flat=flat)

        self.fagn_lbol_out = self.init_arr(len0=self.ncam, nsim=nsim, flat=flat)
        self.fagn1_lbol_out = self.init_arr(len0=self.ncam, nsim=nsim, flat=flat)
        self.fagn2_lbol_out = self.init_arr(len0=self.ncam, nsim=nsim, flat=flat)
        self.fagn_lbol = self.init_arr(nsim=nsim, flat=flat)
        self.fagn1_lbol = self.init_arr(nsim=nsim, flat=flat)
        self.fagn2_lbol = self.init_arr(nsim=nsim, flat=flat)

        self.nsnaps = np.empty((0))
        self.tmrg = np.empty((0))

        for i_sim, subdir in enumerate(subdir_arr):
            path = "%s/%s"%(maindir,subdir)
            if not os.path.exists(path):
                print "Error: path %s not found."%path
                sys.stdout.flush()
                sys.exit()

            #print "i_sim= %d"%i_sim

            data = agn_metadata(path, ncam=self.ncam, tmin=tmin, tmax=tmax,
                                tmax_postmrg=tmax_postmrg, agnx0=agnx0)

            self.nsnaps = np.append(self.nsnaps, len(data.snap))
            self.tmrg = np.append(self.tmrg,data.tmrg)

            #self.snap = np.append(self.snap,data.snap)

            self.snap = self.concat_data(self.snap,data.snap,index=i_sim,flat=flat)
            self.time = self.concat_data(self.time,data.snap,index=i_sim,flat=flat)
            self.posbh1 = self.concat_data(self.posbh1,data.posbh1,index=i_sim,flat=flat)
            self.mdotbh1 = self.concat_data(self.mdotbh1,data.mdotbh1,index=i_sim,flat=flat)
            self.lagn1 = self.concat_data(self.lagn1,data.lagn1,index=i_sim,flat=flat)
            self.lagn_sunsed1 = self.concat_data(self.lagn_sunsed1,data.lagn_sunsed1,
                                                 index=i_sim,flat=flat)
            self.mbh1 = self.concat_data(self.mbh1,data.mbh1,index=i_sim,flat=flat)
            self.posbh2 = self.concat_data(self.posbh2,data.posbh2,index=i_sim,flat=flat)
            self.mdotbh2 = self.concat_data(self.mdotbh2,data.mdotbh2,index=i_sim,flat=flat)
            self.lagn2 = self.concat_data(self.lagn2,data.lagn2,index=i_sim,flat=flat)
            self.lagn_sunsed2 = self.concat_data(self.lagn_sunsed2,data.lagn_sunsed2,
                                                 index=i_sim,flat=flat)
            self.mbh2 = self.concat_data(self.mbh2,data.mbh2,index=i_sim,flat=flat)
            self.nbh = self.concat_data(self.nbh,data.nbh,index=i_sim,flat=flat)
            self.bhsep = self.concat_data(self.bhsep,data.bhsep,index=i_sim,flat=flat)
            self.bhprojsep = self.concat_data(self.bhprojsep,data.bhprojsep,index=i_sim,flat=flat)
        
            self.mstar = self.concat_data(self.mstar,data.mstar,index=i_sim,flat=flat)
            self.mgas = self.concat_data(self.mgas,data.mgas,index=i_sim,flat=flat)
            self.mmetals = self.concat_data(self.mmetals,data.mmetals,index=i_sim,flat=flat)
            self.sfr = self.concat_data(self.sfr,data.sfr,index=i_sim,flat=flat)
            self.lbol_grid = self.concat_data(self.lbol_grid,data.lbol_grid,index=i_sim,flat=flat)
            self.lbol_absorbed = self.concat_data(self.lbol_absorbed,data.lbol_absorbed,
                                                  index=i_sim,flat=flat)
            self.ltot_out = self.concat_data(self.ltot_out,data.ltot_out,index=i_sim,flat=flat)
            self.lir = self.concat_data(self.lir,data.lir,index=i_sim,flat=flat)
            self.w1w2 = self.concat_data(self.w1w2,data.w1w2,index=i_sim,flat=flat)
            self.w2w3 = self.concat_data(self.w2w3,data.w2w3,index=i_sim,flat=flat)
            self.w1lum = self.concat_data(self.w1lum,data.w1lum,index=i_sim,flat=flat)
            self.w2lum = self.concat_data(self.w2lum,data.w2lum,index=i_sim,flat=flat)
            self.dt = self.concat_data(self.dt,data.dt,index=i_sim,flat=flat)

            self.lagn = self.concat_data(self.lagn,data.lagn,index=i_sim,flat=flat)
            self.lagn_sun = self.concat_data(self.lagn_sun,data.lagn_sun,index=i_sim,flat=flat)
            self.mbh = self.concat_data(self.mbh,data.mbh,index=i_sim,flat=flat)
            self.ledd = self.concat_data(self.ledd,data.ledd,index=i_sim,flat=flat)
            self.fedd = self.concat_data(self.fedd,data.fedd,index=i_sim,flat=flat)
            self.ledd1 = self.concat_data(self.ledd1,data.ledd1,index=i_sim,flat=flat)
            self.ledd2 = self.concat_data(self.ledd2,data.ledd2,index=i_sim,flat=flat)
            self.fedd1 = self.concat_data(self.fedd1,data.fedd1,index=i_sim,flat=flat)
            self.fedd2 = self.concat_data(self.fedd2,data.fedd2,index=i_sim,flat=flat)
            
            self.fagn_lbol_out = self.concat_data(self.fagn_lbol_out,data.fagn_lbol_out,index=i_sim,flat=flat)
            self.fagn1_lbol_out = self.concat_data(self.fagn1_lbol_out,data.fagn1_lbol_out,index=i_sim,flat=flat)
            self.fagn2_lbol_out = self.concat_data(self.fagn2_lbol_out,data.fagn2_lbol_out,index=i_sim,flat=flat)
            self.fagn_lbol = self.concat_data(self.fagn_lbol,data.fagn_lbol,index=i_sim,flat=flat)
            self.fagn1_lbol = self.concat_data(self.fagn1_lbol,data.fagn1_lbol,index=i_sim,flat=flat)
            self.fagn2_lbol = self.concat_data(self.fagn2_lbol,data.fagn2_lbol,index=i_sim,flat=flat)


    def init_arr(self, len0=0,len1_max=1200,nsim=1,flat=True):
        if flat:
            arr = np.empty((len0,0)) if len0>0 else np.empty((0))
        else:
            arr = np.zeros((len0,len1_max,nsim)) if len0>0 else np.zeros((len1_max,nsim))
            
        return arr

    def concat_data(self, arr,data,index=0,flat=True):
        test = [1,2] if flat else [2,3]
        if len(arr.shape) not in test:
            print "Error: invalid array shape in concat_data(): ",arr.shape
            sys.exit()

        if flat:
            axis=len(arr.shape)-1
            #print arr.shape, data.shape
            arr = np.append(arr,data,axis=axis)
        else:
            if len(arr.shape)==2:
                arr[:len(data),index] = data
            elif len(arr.shape)==3:
                arr[:,:len(data[0,:]),index] = data

        return arr

def get_rundir(subdir,z=0.0,agnx0=False):

    if subdir=='q0.5_fg0.3_sunruns' and z==0:
        return 'new_fiducial_all_agnx0' if agnx0 else 'new_fiducial_all'
    else:
        if z==0:
            return 'all_agnx0' if agnx0 else 'all'
        else:
            return 'z%.1f_agnx0'%z if agnx0 else 'z%.1f'%z


def define_sunrise_path_arr_z(basepath='/oasis/projects/nsf/hvd115/lblecha',name='fid',
                              z=0,agnx0=False,return_full_path=True):

    fid_arr = ['q1_fg0.3_sunruns',
               'q0.5_fg0.3_sunruns', 
               'q0.5_fg0.3_BT0.1_sunruns',
               'q0.5_fg0.3_BT0.2_sunruns',
               'q0.5_fg0.3fg0.1_BT0BT0.2_sunruns',
               'q0.2_fg0.3_sunruns',
               'q0.5_fg0.1_sunruns',
               'q0.5_fg0.1_BT0.2_sunruns']
    fid_compare_to_hires_arr = [fid_arr[0],fid_arr[1],fid_arr[3],fid_arr[4]]
    fid_major_arr = fid_arr[:5]+fid_arr[6:]
    hires_arr = ['q1_fg0.3_allrx10_sunruns',
                 'q0.5_fg0.3_allrx10_sunruns', 
                 'q0.5_fg0.3_BT0.2_allrx10_sunruns',
                 'q0.5_fg0.3fg0.1_BT0BT0.2_allrx10_sunruns']

    sim_sets = dict(fid=fid_arr,fid_major=fid_major_arr,
                    fid_compare_to_hires=fid_compare_to_hires_arr,
                    hires=hires_arr)

    if name in sim_sets.keys():
        path_arr = ["%s/%s"%(s,get_rundir(s,agnx0=agnx0,z=z)) for s in sim_sets[name]]
        print 'found name %s in sim_sets.keys().'%name
        print path_arr
    else:
        print 'didnt find name %s in sim_sets.keys().'%name
        print sim_sets.keys()
        path_arr = name[0]+"/%s"%get_rundir(name,agnx0=agnx0,z=z)
        print path_arr
        if len([path_arr])>1: 
            name = '%dsims'%len([path_arr])
            if not isinstance(name,list): path_arr = [path_arr]
        else:
            name=name[0]
            path_arr = [name+"/%s"%get_rundir(name,agnx0=agnx0,z=z)]


    if len(path_arr)>1:
        full_path_arr = ['%s/%s'%(basepath,s) for s in path_arr]
        print full_path_arr
        for path in full_path_arr: 
            assert os.path.exists(path),"Error: path %s in sim set %s not found.\n(Defined sets include: %s)."%(path,name," ".join(k for k in sim_sets.keys()))

    else:
        full_path_arr = ["%s/"%basepath+path_arr[0]]
        print full_path_arr
        assert os.path.exists(full_path_arr[0]),"Error: path %s in sim set %s not found.\n(Defined sets include: %s)."%(full_path_arr,name," ".join(k for k in sim_sets.keys()))



    print "Returning list of sunrise sim paths for sim set '%s'."%name

    return full_path_arr if return_full_path else path_arr


def define_sunrise_path_arr(basepath='/oasis/projects/nsf/hvd115/lblecha',name='fid',z=0,
                            sfrhist_w_grid=False,agnx0_compare=False,return_full_path=True):

    fid_arr = ['q1_fg0.3_sunruns/all',
               'q0.5_fg0.3_sunruns/new_fiducial_all', 
               'q0.5_fg0.3_BT0.1_sunruns/all',
               'q0.5_fg0.3_BT0.2_sunruns/all',
               'q0.5_fg0.3fg0.1_BT0BT0.2_sunruns/all',
               'q0.2_fg0.3_sunruns/all',
               'q0.5_fg0.1_sunruns/all',
               'q0.5_fg0.1_BT0.2_sunruns/all']
    test_arr = [fid_arr[0],fid_arr[6]]
    gasrich_arr = fid_arr[:6]
    gaspoor_arr = fid_arr[6:]
    fid_compare_to_hires_arr = [fid_arr[0],fid_arr[1],fid_arr[3],fid_arr[4]]
    fid_major_arr = fid_arr[:5]+fid_arr[6:]
    gasrich_major_arr = fid_arr[:5]
    fid_agnx0_arr = [s+'_agnx0' for s in fid_arr]
    test_agnx0_arr = [fid_agnx0_arr[0],fid_agnx0_arr[6]]
    gasrich_agnx0_arr = fid_agnx0_arr[:6]
    fid_compare_to_hires_agnx0_arr = [fid_agnx0_arr[0],fid_agnx0_arr[1],
                                      fid_agnx0_arr[3],fid_agnx0_arr[4]]
    
    hires_arr = ['q1_fg0.3_allrx10_sunruns/all',
                 'q0.5_fg0.3_allrx10_sunruns/all', 
                 'q0.5_fg0.3_BT0.2_allrx10_sunruns/all',
                 'q0.5_fg0.3fg0.1_BT0BT0.2_allrx10_sunruns/all']
    hires_agnx0_arr = [s+'_agnx0' for s in hires_arr]
    
    lowres_arr = ['q1_fg0.3_allrx0.1_sunruns/all',
                  'q0.5_fg0.3_allrx0.1_sunruns/all', 
                  'q0.5_fg0.3_BT0.2_allrx0.1_sunruns/all',
                  'q0.5_fg0.3fg0.1_BT0BT0.2_allrx0.1_sunruns/all']

    #allres_d1d0_arr = ['q0.5_fg0.3_BT0.2_allrx0.1_sunruns/all',
    #                   'q0.5_fg0.3_BT0.2_sunruns/all',
    #                   'q0.5_fg0.3_BT0.2_allrx5_sunruns/all',
    #                   'q0.5_fg0.3_BT0.2_allrx10_sunruns/all']
    allres_d1d0_arr = ['q0.5_fg0.3_BT0.2_sunruns/all',
                       'q0.5_fg0.3_BT0.2_allrx5_sunruns/all',
                       'q0.5_fg0.3_BT0.2_allrx10_sunruns/all']


    old_fid_arr = ['q1_fg0.3_sunruns/all',
                   'q0.5_fg0.3_sunruns/new_fiducial_all', 
                   'q0.5_fg0.3_BT0.2_sunruns/all',
                   'q0.5_fg0.3fg0.1_BT0BT0.2_sunruns/all',
                   'q0.5_fg0.1_sunruns/all',
                   'q0.5_fg0.1_BT0.2_sunruns/all',
                   'q1_fg0.3_BT0.2_fb0.025_M0.2_largeorb_sunruns/all']
    old_fidmass_arr = old_fid_arr[:-1]
    old_gasrich_arr = old_fid_arr[:4]+[old_fid_arr[-1]]        
    old_gasrich_fidmass_arr = old_fid_arr[:4]
    old_fid_compare_to_agnx0_arr = ['q1_fg0.3_sunruns/all', 
                                    'q0.5_fg0.3_BT0.2_sunruns/all',
                                    'q1_fg0.3_BT0.2_fb0.025_M0.2_largeorb_sunruns/all']
    old_agnx0_arr = [s+'_agnx0' for s in old_fid_compare_to_agnx0_arr]
    
    #defined_subdir_arrs = dict(test=test_arr,fid=fid_arr,fidmass=fidmass_arr,gasrich=gasrich_arr,
    #                           gasrich_fidmass=gasrich_fidmass_arr,hires=hires_arr,
    #                           fid_compare_to_agnx0=fid_compare_to_agnx0_arr,all_agnx0=all_agnx0_arr)
    sim_sets = dict(test=test_arr,test_agnx0=test_agnx0_arr,fid=fid_arr,fid_agnx0=fid_agnx0_arr,
                    gasrich=gasrich_arr,gasrich_agnx0=gasrich_agnx0_arr,
                    gaspoor=gaspoor_arr,fid_compare_to_hires=fid_compare_to_hires_arr,
                    fid_compare_to_hires_agnx0=fid_compare_to_hires_agnx0_arr,
                    fid_major=fid_major_arr,gasrich_major=gasrich_major_arr,
                    hires=hires_arr,hires_agnx0=hires_agnx0_arr,lowres=lowres_arr,
                    allres_d1d0=allres_d1d0_arr,old_fid=old_fid_arr,
                    old_fidmass=old_fidmass_arr,old_gasrich_fidmass=old_gasrich_fidmass_arr,
                    old_fid_compare_to_agnx0=old_fid_compare_to_agnx0_arr,old_agnx0=old_agnx0_arr)
    
    if not agnx0_compare:
        if name in sim_sets.keys():
            path_arr = sim_sets[name]
            print 'found name %s in sim_sets.keys().'%name
            print path_arr
        else:
            print 'didnt find name %s in sim_sets.keys().'%name
            print sim_sets.keys()
            path_arr = name
            if not isinstance(name,list): path_arr = [path_arr]
            if len(path_arr)>1: name = '%dsims'%len(path_arr)

    else:
        path_arr = name
        if not isinstance(name,list): path_arr = [path_arr]
        assert len(path_arr)==1, 'Only 1 sim can be defined for agnx0_compare.'
        path_arr = [path_arr[0],path_arr[0]+'_agnx0'] 

    if sfrhist_w_grid:
        #path_arr = [s+'/sfrhist_w_grid' if '/sfrhist_w_grid' not in s else s for s in path_arr]
        print 'chekcing for sfrhist_w_grid...'
        print path_arr
        for ip,s in enumerate(path_arr):
            if '/sfrhist_w_grid' not in s and os.path.exists('%s/%s/sfrhist_w_grid'%(basepath,s)):
                path_arr[ip] = s+'/sfrhist_w_grid'

    full_path_arr = ['%s/%s'%(basepath,s) for s in path_arr]
    print full_path_arr

    for path in full_path_arr: 
        assert os.path.exists(path),"Error: path %s in sim set %s not found.\n(Defined sets include: %s)."%(path,name," ".join(k for k in sim_sets.keys()))

    print "Returning list of sunrise sim paths for sim set '%s'."%name

    return full_path_arr if return_full_path else path_arr


def get_sim_name(subdir):

    names = {'q1_fg0.3_sunruns': 'A0A0',
             'q0.5_fg0.3_sunruns': 'A1A0',
             'q0.5_fg0.3_BT0.1_sunruns': 'C1C0',
             'q0.5_fg0.3_BT0.2_sunruns': 'D1D0',
             'q0.5_fg0.3fg0.1_BT0BT0.2_sunruns': 'A1E0',
             'q0.2_fg0.3_sunruns': 'A2A0',
             'q0.5_fg0.1_sunruns': 'B1B0',
             'q0.5_fg0.1_BT0.2_sunruns': 'E1E0',
             'q1_fg0.3_allrx10_sunruns': 'A0A0rx10',
             'q0.5_fg0.3_allrx10_sunruns': 'A1A0rx10',
             'q0.5_fg0.3_BT0.2_allrx10_sunruns': 'D1D0rx10',
             'q0.5_fg0.3fg0.1_BT0BT0.2_allrx10_sunruns': 'A1E0rx10',
             'q1_fg0.3_allrx0.1_sunruns': 'A0A0rx0.1',
             'q0.5_fg0.3_allrx0.1_sunruns': 'A1A0rx0.1',
             'q0.5_fg0.3_BT0.2_allrx0.1_sunruns': 'D1D0rx0.1',
             'q0.5_fg0.3fg0.1_BT0BT0.2_allrx0.1_sunruns': 'A1E0rx0.1',
             'q0.5_fg0.3_BT0.2_allrx5_sunruns': 'D1D0rx5'}

    
    #sd = subdir.split('/sfrhist_w_grid')[0]
    sd = subdir.split('/')[0]
    if sd or '/'+sd in names.keys():
        return names[sd]
    else:
        #print "No match for sim name for subdir %s."%sd
        #return ''
        return subdir
        

def rotate(v,theta,phi):

    if len(v) != 3:
        print "Error: function 'rotate' requires a 3-d vector."
        print "v = ",v
        return -1

    ## THIS IS THE WAY TO DO IT. FINALLY.
    ## perform a counterclockwise rotation by phi about the z axis
    ## followed by a clockwise rotation by theta about the y axis
    xnew = v[0]*np.cos(theta)*np.cos(phi) - v[2]*np.sin(theta) + v[1]*np.cos(theta)*np.sin(phi)
    ynew = v[1]*np.cos(phi) - v[0]*np.sin(phi)
    znew = v[2]*np.cos(theta) + v[0]*np.cos(phi)*np.sin(theta) + v[1]*np.sin(theta)*np.sin(phi)


    return np.array([xnew,ynew,znew])

def rotate_arr(v,theta,phi):

    shape = v.shape
    assert len(shape)==2 and 3 in shape, "Error: function 'rotate' requires an array of 3-d vectors."
    if shape[1]==3: v=v.transpose()

    ## THIS IS THE WAY TO DO IT. FINALLY.
    ## perform a counterclockwise rotation by phi about the z axis
    ## followed by a clockwise rotation by theta about the y axis
    xnew = v[0,:]*np.cos(theta)*np.cos(phi) - v[2,:]*np.sin(theta) + v[1,:]*np.cos(theta)*np.sin(phi)
    ynew = v[1,:]*np.cos(phi) - v[0,:]*np.sin(phi)
    znew = v[2,:]*np.cos(theta) + v[0,:]*np.cos(phi)*np.sin(theta) + v[1,:]*np.sin(theta)*np.sin(phi)


    return np.array([xnew,ynew,znew]).transpose() if shape[1]==3 else np.array([xnew,ynew,znew])


def print_arr_info(arr,nanvals=False):

    if nanvals:
        print "min/max/mean/med: %g %g %g %g"%(np.nanmin(arr),np.nanmax(arr),np.nanmean(arr),np.nanmedian(arr))
    else:
        print "min/max/mean/med: %g %g %g %g"%(arr.min(),arr.max(),arr.mean(),np.median(arr))


def calc_lbol_edd(mbh):
    return mbh*ac.MSUN * 4*np.pi*ac.G*ac.MP*ac.C / (ac.THOMSON)


def get_particle_data(gadgetdir,path,snap,verbose=False,riaf=False):

    sf=pyfits.open(path+'/grid_%.3d.fits'%snap,memmap=True)
    age=sf['PARTICLEDATA'].data.field('age')
    snaptime = sf['SFRHIST'].header['SNAPTIME']/1.0e9
    ix_bh = np.where(age!=age)[0]
    mstar_tot = sf['PARTICLEDATA'].data['mass'][age==age].sum()
    del age

    nbh = len(ix_bh)
    if nbh < 1 or nbh > 2:
        print "Error: %d particles have age=nan."%nbh
        return None
    #print "%d BHs found in snap %d."%(nbh,snap)

    id_bh=np.repeat(0,2).astype('uint64')
    id_bh[:nbh]=sf['PARTICLEDATA'].data.field('ID')[ix_bh] 
    #ix_sunidsort = np.argsort(id_bh) if nbh>1 else np.array([0])
    #id_bh = id_bh[ix_sunidsort]

    pos_bh=np.repeat(np.nan,6).reshape(2,3)
    pos_bh[:nbh,:]=sf['PARTICLEDATA'].data.field('position')[ix_bh,:] # [kpc]
    #pos_bh = pos_bh[ix_sunidsort,:]
    #llam_bh=pdata.field('L_lambda')[ix_bh,:] # [W/m]    
    ## this is lbol *only* in the Sunrise SED wavelength range; not using.
    lbol_sunsed_bh = np.repeat(np.nan,2)
    lbol_sunsed_bh[:nbh]=sf['PARTICLEDATA'].data.field('L_bol')[ix_bh] * 1e7 # [erg/s]
    ### this is the dynamical mass; not using.
    #m_bh=np.repeat(np.nan,2)
    #m_bh[:nbh]=sf['PARTICLEDATA'].data.field('mass')[ix_bh] 
    
    snapbh = readsnapbhs.snapshot_bhs(gadgetdir+'/snapshot_%.3d.hdf5'%snap,riaf=riaf)

    m_bh = np.repeat(np.nan,2)
    mdot_bh = np.repeat(np.nan,2)
    #m_bh[:nbh] = snapbh.mass[ix_gadidsort]
    #m_bh[:nbh] = snapbh.mass
    lbol_bh = np.repeat(np.nan,2)
    #lbol_bh[:nbh] = snapbh.lbol
    lbol_edd_bh = np.repeat(np.nan,2)
    #lbol_edd_bh[:nbh] = snapbh.lbol_edd
    reff_gad = np.repeat(np.nan,2)

    if not np.array_equal(snapbh.ids,id_bh[:nbh]): 
        ix_gadidsort = np.argsort(snapbh.ids)
        assert np.array_equal(snapbh.ids[ix_gadidsort],id_bh[:nbh]), \
            "Error: (%s, snap%d) mismatch in gadget and sunrise id arrays."%(path,snap)
    else: ix_gadidsort = np.arange(nbh).astype('int')

    m_bh[:nbh] = snapbh.mass[ix_gadidsort]
    mdot_bh[:nbh] = snapbh.mdot[ix_gadidsort]
    lbol_bh[:nbh] = snapbh.lbol[ix_gadidsort]
    lbol_edd_bh[:nbh] = snapbh.lbol_edd[:nbh][ix_gadidsort]
    #ix_gadidsort = np.argsort(snapbh.ids)
    func_lbol_edd_bh = calc_lbol_edd(m_bh)
    reff = snapbh.rad_eff[ix_gadidsort]

    sfrh_hdr=sf['SFRHIST'].header
    origin = np.array([ sfrh_hdr['translate_origin%s'%s] for s in ['X','Y','Z'] ])
    g1_cen = np.array([ sfrh_hdr['galcenter1%s'%s] for s in ['X','Y','Z'] ])
    if 'galcenter2X' in sfrh_hdr:
        g2_cen = np.array([ sfrh_hdr['galcenter2%s'%s] for s in ['X','Y','Z'] ])
    else: g2_cen = np.repeat(np.nan,3)
    ## pos_bh coords are related to bh pos in snapshot file by: 
    ##    pos_bh = snap_pos - origin
                          
    sf.close()

    if verbose:
        print "Loaded %d galaxies and %d BHs from sunrise snap %d."%((g2_cen[0]==g2_cen[0])+1, nbh, snap)
        print "sanity checks:"
        print "mbh from gadget snap: ",m_bh
        print "mdot from gadget snap: ", mdot_bh
        print "mdot_edd from gadget snap: ", snapbh.mdot_edd[ix_gadidsort]
        print "rad_eff from gadget snap: ", reff     
        print "lbol_sunsed in sunrise: ",lbol_sunsed_bh
        print "lbol from gadget snap: ",lbol_bh
        print "lbol_edd from gadget snap: ", lbol_edd_bh
        print "lbol_edd from func: ", func_lbol_edd_bh
        print "pos from gadget snap (translated to sunrise grid origin): ", snapbh.pos[ix_gadidsort,:]-origin
        print "pos_bh from sunrise: ",pos_bh
    assert np.allclose(lbol_edd_bh[m_bh==m_bh],func_lbol_edd_bh[m_bh==m_bh],rtol=1.0e-4), \
        "Error: (%s, snap%d) mismatch in lbol_edd: %g %g %g %g"%(path,snap,lbol_edd_bh[0],lbol_edd_bh[1],
                                                                 func_lbol_edd_bh[0],func_lbol_edd_bh[1])
    assert np.allclose(snapbh.pos[ix_gadidsort,:]-origin,pos_bh[m_bh==m_bh,:],rtol=0.01,atol=1.0e-4), \
        "Error: (%s, snap%d) mismatch in pos_bh from sunrise  & gadget."%(path,snap)

    #return snaptime, pos_bh, lbol_bh, g1_cen, g2_cen
    #return snaptime, pos_bh, lbol_bh, m_bh, g1_cen, g2_cen, mstar_tot
    return snaptime, pos_bh, mdot_bh, lbol_bh, lbol_sunsed_bh, m_bh, reff, g1_cen, g2_cen, mstar_tot

def overlap_fac(xbh,ybh,zbh,xg,yg,zg,dr,dA):

    assert xg.size==yg.size==zg.size==dr.size,'Error: x,y,z,dr vectors must have equal length.'

    ### this doesn't account for the (many!) cases where dA>dr
    ### need to do something like this:
    ### if (overlap): ovrlp_x = np.min(xg+dr,xbh+dA)-np.max(xg,xbh-dA)
    ### actually this quantity is negative if no overlap, so can just 
    #tmp_x = xbh+dA - xg
    #ovrlp_x = np.zeros((tmp_x.size))
    #ovrlp_x[(tmp_x>0)&(tmp_x<dr)] = tmp_x[(tmp_x>0)&(tmp_x<dr)]
    #ovrlp_x[(tmp_x>=dr)&(tmp_x<2*dA)] = dr[(tmp_x>=dr)&(tmp_x<2*dA)]
    #ovrlp_x[(tmp_x>=2*dA)&(tmp_x<2*dA+dr)] = (dr[(tmp_x>=2*dA)&(tmp_x<2*dA+dr)] - 
    #                                          (tmp_x[(tmp_x>=2*dA)&(tmp_x<2*dA+dr)]-2*dA))
    ovrlp_x = np.zeros((xg.size))
    tmp_x = np.array([np.minimum(xg[i]+dr[i],xbh+dA)-np.maximum(xg[i],xbh-dA) for i in range(xg.size)])
    ovrlp_x[tmp_x>0] = tmp_x[tmp_x>0]
    fx = ovrlp_x / dr

    #tmp_y = ybh+dA - yg
    #ovrlp_y = np.zeros((tmp_y.size))
    #ovrlp_y[(tmp_y>0)&(tmp_y<dr)] = tmp_y[(tmp_y>0)&(tmp_y<dr)]
    #ovrlp_y[(tmp_y>=dr)&(tmp_y<2*dA)] = dr[(tmp_y>=dr)&(tmp_y<2*dA)]
    #ovrlp_y[(tmp_y>=2*dA)&(tmp_y<2*dA+dr)] = (dr[(tmp_y>=2*dA)&(tmp_y<2*dA+dr)] - 
    #                                          (tmp_y[(tmp_y>=2*dA)&(tmp_y<2*dA+dr)]-2*dA))
    ovrlp_y = np.zeros((yg.size))
    tmp_y = np.array([np.minimum(yg[i]+dr[i],ybh+dA)-np.maximum(yg[i],ybh-dA) for i in range(yg.size)])
    ovrlp_y[tmp_y>0] = tmp_y[tmp_y>0]
    fy = ovrlp_y / dr

    tmp_z = zg+dr - zbh
    ovrlp_z = np.zeros((tmp_z.size))
    ovrlp_z[(tmp_z>0)&(tmp_z<dr)] = tmp_z[(tmp_z>0)&(tmp_z<dr)]
    ovrlp_z[(tmp_z>=dr)] = dr[(tmp_z>=dr)]
    fz = ovrlp_z / dr

    return fx,fy,fz

def calc_grid_nh(path,snap,verbose=False,riaf=False,res=1.0,rot_pos_bh=[],theta=[],phi=[],Ngas_tot_aux=[]):

    assert len(rot_pos_bh.shape)==3 and (rot_pos_bh.shape[0] in (1,2)) and rot_pos_bh.shape[2]==3,\
        'Error: keyword rot_pos_bh must have shape (nbh,ncam,3)'
    nbh=rot_pos_bh.shape[0]
    ncam=rot_pos_bh.shape[1]
    assert len(theta)==ncam and len(phi)==ncam, 'theta & phi must have length ncam (%d)'%ncam

    print "\n retrieving grid data and calculating Ngas_bh..."

    #Ngas = np.zeros((nbh,ncam)) -1
    Ngas = np.zeros((2,ncam)) -1
    Ngas_tot = np.zeros((2,ncam)) -1

    sf=pyfits.open(path+'/grid_%.3d.fits'%snap,memmap=True)
    gstruct=sf['GRIDSTRUCTURE']

    if 'qpos' not in gstruct.data.columns.names:
        print "WARNING: No grid structure data in file %s/grid_%.3d.fits"%(path,snap)
        return Ngas, Ngas_tot

    isRefined=gstruct.data['structure']
    qpos=gstruct.data['qpos']
    level=gstruct.data['level']
    gscale=np.array([gstruct.header['MAX%s'%v]-gstruct.header['MIN%s'%v] for v in ('X','Y','Z')])
    gmin=np.array([gstruct.header['MIN%s'%v] for v in ('X','Y','Z')])
    ## bitwise left-shift operator returns inf values when we attempt array operations
    pos = np.zeros((level.size,3))-1.0
    for j in np.arange(3):
        for l in np.arange(level.size):
            pos[l,j] = qpos[l,j]/(1<<level[l])
        pos[:,j] = pos[:,j]*gscale[j] + gmin[j]
    if verbose: print qpos.shape, level.shape, pos.shape
    level_f = level[~isRefined]
    pos_f = pos[~isRefined]
    if verbose: print level_f.shape, pos_f.shape

    mg=sf['GRIDDATA'].data['mass_gas']
    #mm=sf['GRIDDATA'].data['mass_metals']
    cvol=sf['GRIDDATA'].data['cell_volume']
    rho=mg/cvol
    dr=cvol**(1.0/3)
    #print mg.shape, mm.shape, cvol.shape
    #print pos.min(),pos.max()
    #print pos_f.min(),pos_f.max()
    #print (pos_f[:,0]+dr).min(),(pos_f[:,0]+dr).max()
    #print mg.min(),mg.max()
    #print cvol.min(),cvol.max()
    #print rho.min(),rho.max()
    #print dr.min(),dr.max()
    #print mg[cvol**(1/3.)>50.]

    sf.close()
    del gstruct

    assert (res>0 and res<np.max(pos_f+np.tile(cvol,(3,1)).transpose())),\
        'Error: invalid resolution for nh calc: %g.'%res
                                 

    dA=0.5*res
    if verbose:
        print "dA, min/max dr: %g %g %g"%(dA,dr.min(),dr.max())
        print dr.size,dr[dr>dA].size
    for c in range(ncam):
        rot_pos_grid = rotate_arr(pos_f,theta[c],phi[c])
        xg_rot = rot_pos_grid[:,0]
        yg_rot = rot_pos_grid[:,1]
        zg_rot = rot_pos_grid[:,2]
        if verbose: print "\ncamera %d..."%c

        for i in range(nbh):
            if verbose: print "BH %d..."%i
            xbh_rot,ybh_rot,zbh_rot = rot_pos_bh[i,c,0],rot_pos_bh[i,c,1],rot_pos_bh[i,c,2]
            if verbose: print "rotated BH pos: %g %g %g"%(xbh_rot,ybh_rot,zbh_rot)

            fx,fy,fz = overlap_fac(xbh_rot,ybh_rot,zbh_rot,xg_rot,yg_rot,zg_rot,dr,dA)
            ix_tot = np.where((fx>0)&(fy>0))[0]
            ix_full_tot = np.where((fx==1)&(fy==1))[0]
            ix_partial_tot = np.where((fx>0)&(fy>0)&((fx<1)|(fy<1)))[0]
            ix = np.where((fx>0)&(fy>0)&(fz>0))[0]
            ix_full = np.where((fx==1)&(fy==1)&(fz==1))[0]
            ix_partial = np.where((fx>0)&(fy>0)&(fz>0)&((fx<1)|(fy<1)|(fz<1)))[0]

            if verbose:
                print "min/max fx:",fx.min(),fx.max()
                print "min/max fy:",fy.min(),fy.max()
                print "min/max fz:",fz.min(),fz.max()
                print ix.size, ix_full.size, ix_partial.size
                print ix_tot.size, ix_full_tot.size+ix_partial_tot.size,ix_full_tot.size, ix_partial_tot.size
                print fz.size
                print np.where(fz==0)[0].size
                print np.where(fz>0)[0].size
                print np.where((fz>0)&(fz<1))[0].size + np.where(fz==1)[0].size


            #cvol_in_col = cvol[ix] * fx[ix]*fy[ix]*fz[ix]
            cvol_in_col = cvol * fx*fy*fz
            mg_in_col = mg * fx*fy*fz
            assert np.isclose(mg_in_col.sum(),(mg[ix]*fx[ix]*fy[ix]*fz[ix]).sum()),'mg_in_col incorrect: %g %g'%(mg_in_col.sum(),(mg[ix]*fx[ix]*fy[ix]*fz[ix]).sum())
            mg_in_tot_col = mg[ix_tot] * fx[ix_tot]*fy[ix_tot]
            #print fx[ix].size,fx[ix].size-fx[ix_full].size
            #print fx[(fx<1)&(fx>0)].size,fx[fx>0].size
            #print fy[(fy<1)&(fy>0)].size,fy[fy>0].size
            #print fz[(fz<1)&(fz>0)].size,fz[fz>0].size
            #print fx[(fx>0)&(fy>0)&(fz>0)].size

            #Ngas[i,c] = mg_in_col.sum()*ac.MSUN/ac.MP/(res*ac.KPC)**2
            #Ngas_tot[i,c] = mg_in_tot_col.sum()*ac.MSUN/ac.MP/(res*ac.KPC)**2
            Ngas[i,c] = mg_in_col.sum()*pjc.MSUN_CGS/pjc.MP_CGS/(res*pjc.KPC_CGS)**2
            Ngas_tot[i,c] = mg_in_tot_col.sum()*pjc.MSUN_CGS/pjc.MP_CGS/(res*pjc.KPC_CGS)**2

            if verbose:
                print "cvol[ix].sum()=",cvol[ix].sum()
                print "cvol_in_col.sum()=",cvol_in_col.sum()
                print "mg[ix].sum()/mg.sum()=",mg[ix].sum()/mg.sum()
                print "mg_in_col.sum()/mg.sum()=",mg_in_col.sum()/mg.sum()
                print "mg_in_tot_col.sum()=",mg_in_tot_col.sum()/mg.sum()
                
                print "total unrefined gas cells in grid: %d"%mg.shape
                print "total intersecting with column (res %g kpc): %d"%(res,mg[ix].size)
                print "# fully intersecting with column (res %g kpc): %d"%(res,mg[ix_full].size)
                print "# partially intersecting with column (res %g kpc): %d"%(res,mg[ix_partial].size)
                print "total intersecting with column - no z cut (res %g kpc): %d"%(res,mg[ix_tot].size)

                print "Ngas: ",Ngas[i,c]
                print "0.5*Ngas_tot: ",0.5*Ngas_tot[i,c]
                print "Ngas/(0.5Ngas_tot):",Ngas[i,c]/(0.5*Ngas_tot[i,c])
                print "Ngas_tot: ",Ngas_tot[i,c]
                if len(Ngas_tot_aux)>0:
                    print "Ngas_tot_aux: ",Ngas_tot_aux[i,c]
                    print "Ngas_tot/Ngas_tot_aux:",Ngas_tot[i,c]/Ngas_tot_aux[i,c]
            
    return Ngas, Ngas_tot
        
def write_bhfile(gadgetdir,path,snap,riaf=False):

    time,pos_bh,mdot_bh,lbol_bh,lbol_sunsed_bh,m_bh,reff,g1_cen,g2_cen,mstar_tot=get_particle_data(gadgetdir,path,snap,riaf=riaf)

    ### bh_info.txt file format (columns): ###
    ### snapshot num, time, x pos bh1, y pos bh1, z pos bh1, mdot bh1, lbol bh1, lbol sunsed bh1, mass bh1, x pos bh2, y pos bh2, z pos bh2, mdot bh2, lbol bh2, lbol sunsed bh2, mass bh2
    with open(path+'/bh_info.txt','a') as fbh:        
        pstr1=" ".join([str(p) for p in pos_bh[0,:]])
        pstr2=" ".join([str(p) for p in pos_bh[1,:]])
        fbh.write(("%d %g "+2*"%s %g %g %g %g "+"\n")%(snap,time,pstr1,mdot_bh[0],lbol_bh[0],lbol_sunsed_bh[0],
                                                       m_bh[0],pstr2,mdot_bh[1],lbol_bh[1],lbol_sunsed_bh[1],m_bh[1]))


def calc(gadgetdir,path,snap,cam=range(7),quickimage=False,global_info_only=False,
         grid_nh_res=0,riaf=False,write_all_files=True,verbose=False):

    #print matplotlib.get_backend()
    #time, pos_bh, lbol_bh, g1_cen, g2_cen = get_particle_data(path,snap)
    #time, pos_bh, lbol_bh, g1_cen, g2_cen, mstar_tot = get_particle_data(path,snap)
    time,pos_bh,mdot_bh,lbol_bh,lbol_sunsed_bh,m_bh,reff,g1_cen,g2_cen,mstar_tot=get_particle_data(gadgetdir,path,snap,riaf=riaf)

    mf=pyfits.open(path+'/mcrx_%.3d.fits'%snap,memmap=True)
    ##iq=mf['INTEGRATED_QUANTITIES']
    print 'opened file %s/mcrx_%.3d.fits'%(path,snap)

    mgas_tot = mf['MAKEGRID'].header['M_G_TOT']
    mmetals_tot = mf['MAKEGRID'].header['M_metals_tot']
    sfr_tot = mf['MAKEGRID'].header['SFR_TOT']

    ## load luminosities and convert to erg/s
    lbol_grid = mf['INTEGRATED_QUANTITIES'].header['L_bol_grid'] * 1e7
    lbol_absorbed = mf['INTEGRATED_QUANTITIES'].header['L_bol_absorbed'] * 1e7 
    lout = np.array([ mf['INTEGRATED_QUANTITIES'].header['L_OUT%d'%n] for n in cam ]) * 1e7
    lir = np.array([ mf['INTEGRATED_QUANTITIES'].header['L_IR%d'%n] for n in cam ]) * 1e7

    if verbose: 
        print "lbol_grid = %g"%lbol_grid
        print "lbol_absorbed = %g"%lbol_absorbed
        print "lout = ",lout
        print "lir = ",lir

    if write_all_files:
        with open(path+"/star_gas_info.txt",'a') as fsg:
            fsg.write(("%d "+5*"%g "+"\n")%(snap,time,mstar_tot,mgas_tot,mmetals_tot,sfr_tot))
            
            with open(path+'/lum_info.txt','a') as fplum:
                fplum.write(("%d "+3*"%g "+"%s %s\n")%(snap,time,lbol_grid,lbol_absorbed,
                                                       " ".join([str(l) for l in lout]),
                                                       " ".join([str(l) for l in lir])))

    nbh = len(pos_bh[pos_bh[:,0]==pos_bh[:,0],0])
    #print "pos_bh:",pos_bh
    #print "lbol_bh:",lbol_bh
    #print "nbh:",nbh
    #print "g1_cen:",g1_cen
    #print "g2_cen:",g2_cen

    bh_3d_sep = np.nan
    if nbh==2: bh_3d_sep = np.sqrt(np.sum((pos_bh[1,j]-pos_bh[0,j])**2 for j in range(3)))

    if write_all_files:
        with open(path+'/bh_info.txt','a') as fbh:        
            pstr1=" ".join([str(p) for p in pos_bh[0,:]])
            pstr2=" ".join([str(p) for p in pos_bh[1,:]])
            fbh.write(("%d %g "+2*"%s %g %g %g %g "+"\n")%(snap,time,pstr1,mdot_bh[0],lbol_bh[0],
                                                           lbol_sunsed_bh[0],m_bh[0],pstr2,mdot_bh[1],
                                                           lbol_bh[1],lbol_sunsed_bh[1],m_bh[1]))

    if global_info_only:   
        return

    if write_all_files:
        fp = open(path+'/column_density.out','a')
        fp.write("\n\nprocessing snap %d (t=%g Gyr): nbh=%d\n\n"%(snap,time,nbh))
        fp.write("lbol bh1,bh2: %g %g\n"%(lbol_bh[0],lbol_bh[1]))
        fp.write("bh1 pos: %s\n"%(" ".join([str(p) for p in pos_bh[0,:]])))
        if nbh==2:
            fp.write("bh2 pos: %s\n"%(" ".join([str(p) for p in pos_bh[1,:]])))
            fp.write("3d bh sep: %g\n"%(bh_3d_sep))
            
        fp.write("g1_cen: %s\n"%(" ".join([str(g1) for g1 in g1_cen])))
        fp.write("g2_cen: %s\n"%(" ".join([str(g2) for g2 in g2_cen])))


    rotated_pos_bh = np.zeros((nbh,len(cam),3))
    bh_los_pix = np.zeros((2,len(cam),2)) -1
    bh_proj_sep = np.repeat(np.nan,(len(cam)))
    Ngas_bh_tot_aux = np.zeros((2,len(cam))) -1
    Ngas_tot_aux_max = np.zeros((len(cam))) -1
    pix_Ngas_tot_aux_max = np.zeros((len(cam),2)) -1
    theta = np.zeros((len(cam)))
    phi = np.zeros((len(cam)))
    if quickimage:
        fig=plt.figure(figsize=(6,6))
        mpl.rcParams.update({'font.size': 10})
        subp0=331

    for c in cam:
        print "processing camera %d"%c
        if write_all_files: fp.write("\ncamera %d:\n"%c)

        campar=mf['camera%d-parameters'%c].header
        theta[c]=campar['theta']
        phi[c]=campar['phi']
        xpix = campar['xsize']
        ypix = campar['ysize']
        x_kpc_per_pix = mf['camera%d'%c].header['CD2_2']
        y_kpc_per_pix = mf['camera%d'%c].header['CD1_1']
        aux=mf['CAMERA%d-AUX'%c]
        Ngas_tot_aux=aux.data[0,:,:] * pjc.MSUN_CGS / (pjc.KPC_CGS)**2 / pjc.MP_CGS ## Msun/kpc^2 to cm^-2/mu
        Ngas_tot_aux_max[c] = Ngas_tot_aux.max()
        pix_Ngas_tot_aux_max[c,:] = np.where(Ngas_tot_aux==Ngas_tot_aux.max())

        if write_all_files:
            fp.write("x,y kpc_per_pix: %g, %g\n"%(x_kpc_per_pix,y_kpc_per_pix))
            fp.write("theta=%g, phi=%g\n"%(theta[c],phi[c]))
            fp.write("campos: %g %g %g\n"%(campar['camposx'],campar['camposy'],campar['camposz']))
            tmp=rotate((0,0,1),theta[c],phi[c])*10000
            fp.write("should match this: %g %g %g\n"%(tmp[0],tmp[1],tmp[2]))
            fp.write("max Ngas_tot_aux: %g (pos=[%d,%d])\n"%(Ngas_tot_aux.max(),pix_Ngas_tot_aux_max[c,0],
                                                             pix_Ngas_tot_aux_max[c,1]))

        if quickimage:
            plt.clf()
            plt.close()
            im = mf['camera%d'%c].data.sum(axis=0)
            im[im==0] = 1.0e-20
            ax = fig.add_subplot(subp0+c)
            #plt.xlim(0,xpix)
            #plt.ylim(0,ypix)
            ax.set_xlim(0,xpix)
            ax.set_ylim(0,ypix)
            ax.imshow(np.log10(im),vmin=-1,vmax=1.1*np.log10(im).max())

        for i in range(nbh):

            rotated_pos_bh[i,c,:] = rotate(pos_bh[i,:],theta[c],phi[c])
            if i==1:
                bh_proj_sep[c] = np.sqrt(np.sum((rotated_pos_bh[1,c,j]-
                                                 rotated_pos_bh[0,c,j])**2 for j in range(2)))
                if write_all_files:
                    fp.write("3d sep should match this: %g\n"%(np.sqrt(np.sum((rotated_pos_bh[1,c,j]-
                                                                               rotated_pos_bh[0,c,j])**2 
                                                                              for j in range(3)))))
                    fp.write("proj. bh sep: %g\n"%(bh_proj_sep[c]))


            #bh_los_pix[i,c,:] = (np.int(rotated_pos_bh[i,c,0]/x_kpc_per_pix)+xpix/2,
            #                     np.int(rotated_pos_bh[i,c,1]/y_kpc_per_pix)+ypix/2)
            bh_los_pix[i,c,:] = (round(rotated_pos_bh[i,c,0]/x_kpc_per_pix)+xpix/2,
                                 round(rotated_pos_bh[i,c,1]/y_kpc_per_pix)+ypix/2)

            if (bh_los_pix[i,c,0]<0 or bh_los_pix[i,c,0]>=xpix or
                bh_los_pix[i,c,1]<0 or bh_los_pix[i,c,1]>=ypix ):
                print "BH %d outside camera %d FOV for snap %d."%(i,c,snap)
                print "rotated_pos_bh: ",rotated_pos_bh[i,c,:]
                print "bh_los_pix: ",bh_los_pix[i,c,:]
                print "xpix,ypix: ",xpix,ypix
                print "x,y size [kpc]: ",xpix*x_kpc_per_pix,ypix*y_kpc_per_pix
                continue

            if verbose:
                print "rotated_pos_bh:",rotated_pos_bh[i,c,:]
                print "x_kpc_per_pix:",x_kpc_per_pix
                print "y_kpc_per_pix:",y_kpc_per_pix
                print "bh_los_pix:",bh_los_pix[i,c,:]
                print "Ngas_tot_aux.shape=",Ngas_tot_aux.shape
            Ngas_bh_tot_aux[i,c] = Ngas_tot_aux[bh_los_pix[i,c,0],bh_los_pix[i,c,1]] ## Msun/kpc^2

            if quickimage: ax.plot(bh_los_pix[i,c,1],bh_los_pix[i,c,0],'rx',markersize=8)

            if write_all_files:
                fp.write("bh %d:\n"%(i+1))
                fp.write("pos: %s\n"%(" ".join([str(p) for p in pos_bh[i,:]])))
                fp.write("rotated pos: %s\n"%(" ".join([str(p) for p in rotated_pos_bh[i,c,:]])))
                fp.write("bh_los_pix: %s\n"%(" ".join([str(p) for p in bh_los_pix[i,c,:]])))
                fp.write("Ngas_bh_tot_aux [cm^-2/mu]: %g\n"%Ngas_bh_tot_aux[i,c])
            


    mf.close()
    if write_all_files: fp.close()

    res=grid_nh_res if grid_nh_res>0 else x_kpc_per_pix
    Ngas_bh_grid,Ngas_bh_tot_grid=calc_grid_nh(path,snap,verbose=False,riaf=False,res=res,
                                               rot_pos_bh=rotated_pos_bh,theta=theta,phi=phi,
                                               Ngas_tot_aux=[])
    ### Ngas_bh_cam[N]_*txt file format (columns): ###
    ### snapshot num, bh1 x pixel, bh1 y pixel, bh1 tot NH (aux calculation), bh1 tot NH (grid calculation), bh1 NH (grid), bh2 x pixel, bh2 y pixel, bh2 tot NH (aux calculation), bh2 tot NH (grid calculation), bh2 NH (grid)
    ##### Notes: 
    ##### 'aux' calculation is from SUNRISE auxiliary output, at the pixel scale; includes total column through the whole galaxy (not just to the BH)
    ##### 'grid' calculation uses the actual gas grid from SUNRISE
    #####       This is scaled to an arbitrary resolution as specified in the filename ('grid_res_[X]' where X is the resolution in pc).
    #####       For best results this should be scaled to the resolution of the Gadget simulation, ~ 2.8*rsoft.
    #####       'Ngas_bh_tot_grid' variables are the total column through the whole galaxy (not just to the BH).
    #####       'Ngas_bh_grid' variables are the column density integrated only from the position of the BH.
    for c in cam:
        rstr='_grid_res_%.3d'%(1000*grid_nh_res) if grid_nh_res>0 else ''
        ngfp = open(path+"/Ngas_bh_cam%d%s.txt"%(c,rstr),"a")
        ngfp.write("%s\n"%( ("%d "%snap+("%g "*10))%(bh_los_pix[0,c,0],bh_los_pix[0,c,1],
                                                     Ngas_bh_tot_aux[0,c],Ngas_bh_tot_grid[0,c],
                                                     Ngas_bh_grid[0,c],bh_los_pix[1,c,0],
                                                     bh_los_pix[1,c,1],Ngas_bh_tot_aux[1,c],
                                                     Ngas_bh_tot_grid[1,c],Ngas_bh_grid[1,c]) ))
        ngfp.close()
    
    compfp = open(path+"/compare_ngas%s.out"%(rstr),"a")
    compfp.write("\nsnap %d\n"%snap)
    for i in range(nbh):
        compfp.write("BH %d\n"%i)
        compfp.write("Ngas_bh_grid: %s\n"%(" ".join("%g"%n for n in Ngas_bh_grid[i,:])))
        compfp.write("Ngas_bh_tot_grid: %s\n"%(" ".join("%g"%n for n in Ngas_bh_tot_grid[i,:])))
        compfp.write("Ngas_bh_tot_aux: %s\n"%(" ".join("%g"%n for n in Ngas_bh_tot_aux[i,:])))
        compfp.write("Ngas_bh_tot_grid/Ngas_bh_tot_aux: %s\n"%(" ".join("%g"%n 
                                                                        for n in (Ngas_bh_tot_grid[i,:]/
                                                                                  Ngas_bh_tot_aux[i,:]))))
        compfp.write("<Ngas_bh_grid>: %g\n"%Ngas_bh_grid[i,:].mean())
        compfp.write("<Ngas_bh_tot_grid>: %g\n"%Ngas_bh_tot_grid[i,:].mean())
        compfp.write("<Ngas_bh_tot_aux>: %g\n"%Ngas_bh_tot_aux[i,:].mean())
        compfp.write("<Ngas_bh_tot_grid>/<Ngas_bh_tot_aux>: %g\n"%(Ngas_bh_tot_grid[i,:].mean()/
                                                                   Ngas_bh_tot_aux[i,:].mean()))
    compfp.close()

    if write_all_files:
        bhsfp = open(path+"/bh_sep.txt","a")
        psstr=" ".join("%g"%s for s in bh_proj_sep)
        bhsfp.write(("%d %g %g "+psstr+"\n")%(snap,time,bh_3d_sep))
        bhsfp.close()
        
    if quickimage and Ngas_bh_tot.max()>0:
        fig.subplots_adjust(wspace=0.3,hspace=0.2)
        fig.savefig(path+'/quickimage_snap%d.eps'%snap,format='eps',dpi=1200)
        plt.clf()
        plt.cla()
        plt.close()


def process_snaps(path='/oasis/projects/nsf/hvd115/lblecha/',
                  gadgetdir='p3new_merger_sims/q0.5_fg0.3_nomrg',
                  subdir='q0.5_fg0.3_sunruns/new_fiducial_all',
                  startsnap=-1,global_info_only=False,write_bhfile_only=False,
                  write_ngasfiles_only=False,quickimage=False,append=False,
                  overwrite_all=False, riaf=False, grid_nh_res=0):
    
    if write_bhfile_only and write_ngasfiles_only:
        print "Error: cannot simultaneously define write_bhfile_only and write_ngasfiles_only."
        return
    if global_info_only and write_ngasfiles_only:
        print "Error: cannot simultaneously define global_info_only and write_ngasfiles_only."
        return
    
    mcrxfiles = glob.glob(path+subdir+'/mcrx*fits')
    snaps = np.sort(np.array([np.int(re.search('(?<=mcrx_)\d+',f).group(0)) 
                              for f in mcrxfiles]))
    print path+'/'+subdir
    print startsnap
    print snaps
    print snaps.size
    if startsnap>=0:
        nsnaps = snaps.size
        print nsnaps
        print isinstance(snaps,np.ndarray)
        mask = (snaps>=startsnap)
        print isinstance(mask,np.ndarray)
        print isinstance((snaps>startsnap),np.ndarray)
        print mask
        #snaps = snaps[(snaps>=startsnap)]
        snaps = snaps[mask]
        print snaps
        print "processing %d of %d snaps in %s"%(snaps.size,nsnaps,subdir)
    else:
        print "processing %d snaps (%d-%d) in %s"%(snaps.size,snaps.min(),snaps.max(),subdir)

    rstr='_grid_res_%.3d'%(1000*grid_nh_res) if grid_nh_res>0 else ''

    write_all_files=True
    if write_bhfile_only:
        outfiles = (glob.glob(path+subdir+'/bh_info.txt'))
    elif write_ngasfiles_only:
        write_all_files=False
        outfiles = (glob.glob(path+subdir+'/Ngas_bh_cam*%s.txt'%rstr)+
                    glob.glob(path+subdir+'/compare_ngas%s.out'%rstr))        
        all_Ngas_files = (glob.glob(path+subdir+'/Ngas*txt'))
        other_Ngas_files = [t for t in all_Ngas_files if t not in outfiles]
    else:
        outfiles = (glob.glob(path+subdir+'/*_info.txt'))
        if not global_info_only:
            outfiles = outfiles + (glob.glob(path+subdir+'/bh_sep.txt')+
                                   glob.glob(path+subdir+'/Ngas_bh_cam*%s.txt'%rstr)+
                                   glob.glob(path+subdir+'/column_density.out')+
                                   glob.glob(path+subdir+'/compare_ngas%s.out'%rstr))
            all_Ngas_files = (glob.glob(path+subdir+'/Ngas*txt'))
            other_Ngas_files = [t for t in all_Ngas_files if t not in outfiles]

    if not write_bhfile_only:
        print "outfiles:",len(outfiles),outfiles
        print "all_Ngas_files:",len(all_Ngas_files),all_Ngas_files
        print "other_Ngas_files:",len(other_Ngas_files),other_Ngas_files
        if other_Ngas_files:
            print "Warning: found existing Ngas files at different res:"
            print other_Ngas_files
            print "Generating new set of Ngas files: %s"%('Ngas_bh_cam*%s.txt'%rstr)

    if outfiles and not append:
        if overwrite_all or raw_input("Overwrite existing output files? (y/N) ")=="y":
            print "WARNING: Overwriting all existing output files."
            for of in outfiles: rm(of)

    for snap in snaps:
        print "\n\nprocessing snap %d..."%snap
        if write_bhfile_only:
            write_bhfile(path+gadgetdir, path+subdir, snap, riaf=riaf)
        else:
            calc(path+gadgetdir, path+subdir, snap, global_info_only=global_info_only, 
                 quickimage=quickimage,riaf=riaf, grid_nh_res=grid_nh_res,write_all_files=write_all_files)


def rm(filename):
    try:
        os.remove(filename)
    except:
        pass


def ABmag(Leff,ewidth_nu):

    return -2.5 * np.log10(Leff/(ewidth_nu*abmag_zp)) + 5 + 85.19





def img_mags_bhlos(path,snaps,bh1_xpix,bh1_ypix,bh2_xpix,bh2_ypix,
                   cam=range(7),make_image=False):

    if np.isscalar(snaps): 
        snaps=[snaps]
        bh1_xpix=bh1_xpix.reshape(len(bh1_xpix),1)
        bh1_ypix=bh1_ypix.reshape(len(bh1_ypix),1)
        bh2_xpix=bh2_xpix.reshape(len(bh2_xpix),1)
        bh2_ypix=bh2_ypix.reshape(len(bh2_ypix),1)
        if len(cam)>len(bh1_xpix[:,0]):
            print "Error: too many cameras specified for BH LOS data."
            print bh1_xpix.shape,len(cam)
            sys.exit()
    w1w2_bhpix = np.zeros((2,len(cam),len(snaps)))*np.nan
    w2w3_bhpix = np.zeros((2,len(cam),len(snaps)))*np.nan
    Av_bhpix = np.zeros((2,len(cam),len(snaps)))*np.nan

    for i,snap in enumerate(snaps):
        bb_fname=(path+'/broadband_%.3d.fits'%snap)
        print "\nopening file %s..."%bb_fname
        try:
            fbb = pyfits.open(bb_fname)
        except IOError:
            "WARNING: Could not open file %s."%(bb_fname)
            #sys.exit()
            continue
        filterdata=fbb['FILTERS'].data
        filternames = np.array([ x.strip() for x in filterdata['filter'] ])
        lambda_eff = filterdata['lambda_eff'].astype('float64')
        ewidth_lambda = filterdata['ewidth_lambda'].astype('float64')
        ewidth_nu = filterdata['ewidth_nu'].astype('float64')
        ixw1 = filternames.tolist().index('wise/WISE-W1.res')
        ixw2 = filternames.tolist().index('wise/WISE-W2.res')
        ixw3 = filternames.tolist().index('wise/WISE-W3.res')
        ixV = filternames.tolist().index('johnson/V_Johnson.res')
        for c in cam:
            print "\ncamera %d:"%c
            campar_hdr=fbb['CAMERA%d-PARAMETERS'%c].header
            xpixels = campar_hdr['XSIZE']
            ypixels = campar_hdr['YSIZE']
            sr_per_pix = np.float64(campar_hdr['solid_angle'])/(1.0*xpixels*ypixels)
            cameradist = np.float64(campar_hdr['cameradist'])
            print "sr_per_pix: ",sr_per_pix
            print "cameradist: ",cameradist
            codeSB_to_W = sr_per_pix * 4*np.pi * (cameradist*pjc.KPC_SI)**2
            print "codeSB_to_W: ",codeSB_to_W

            im_W1 = fbb['CAMERA%d-BROADBAND'%c].data[ixw1,:,:].astype('float64')
            print "min/max/shape im_W1:",im_W1.min(),im_W1.max(),im_W1.shape
            im_Leff_W1 = im_W1 * codeSB_to_W * ewidth_lambda[ixw1]
            #im_AB_W1 = -2.5 * np.log10(im_Leff_W1/(ewidth_nu[ixw1]*abmag_zp)) + 5 + 85.19
            im_AB_W1 = ABmag(im_Leff_W1,ewidth_nu[ixw1])
            print "min/max/tot im_AB_W1:",im_AB_W1.min(),im_AB_W1.max(),ABmag(im_Leff_W1.sum(),ewidth_nu[ixw1])
            print "real tot AB:",fbb['FILTERS'].data.field('AB_mag%d'%c)[ixw1]            
            im_Vega_W1 = im_AB_W1 - W1_Vega_to_AB

            im_W2 = fbb['CAMERA%d-BROADBAND'%c].data[ixw2,:,:].astype('float64')
            im_Leff_W2 = im_W2 * codeSB_to_W * ewidth_lambda[ixw2]
            im_AB_W2 = ABmag(im_Leff_W2,ewidth_nu[ixw2])
            print "min/max/tot im_AB_W2:",im_AB_W2.min(),im_AB_W2.max(),ABmag(im_Leff_W2.sum(),ewidth_nu[ixw2])
            print "real tot AB:",fbb['FILTERS'].data.field('AB_mag%d'%c)[ixw2]   
            im_Vega_W2 = im_AB_W2 - W2_Vega_to_AB

            im_W3 = fbb['CAMERA%d-BROADBAND'%c].data[ixw3,:,:].astype('float64')
            im_Leff_W3 = im_W3 * codeSB_to_W * ewidth_lambda[ixw3]
            im_AB_W3 = ABmag(im_Leff_W3,ewidth_nu[ixw3])
            print "min/max/tot im_AB_W3:",im_AB_W3.min(),im_AB_W3.max(),ABmag(im_Leff_W3.sum(),ewidth_nu[ixw3])
            print "real tot AB:",fbb['FILTERS'].data.field('AB_mag%d'%c)[ixw3]            
            im_Vega_W3 = im_AB_W3 - W3_Vega_to_AB

            im_V = fbb['CAMERA%d-BROADBAND'%c].data[ixV,:,:].astype('float64')
            print "min/max/shape im_V:",im_V.min(),im_V.max(),im_V.shape
            im_Leff_V = im_V * codeSB_to_W * ewidth_lambda[ixV]
            im_AB_V = ABmag(im_Leff_V,ewidth_nu[ixV])
            print "min/max/tot im_AB_V:",im_AB_V.min(),im_AB_V.max(),ABmag(im_Leff_V.sum(),ewidth_nu[ixV])
            print "real tot AB:",fbb['FILTERS'].data.field('AB_mag%d'%c)[ixV]     

            im_V_ns = fbb['CAMERA%d-BROADBAND-NONSCATTER'%c].data[ixV,:,:].astype('float64')
            im_Leff_V_ns = im_V_ns * codeSB_to_W * ewidth_lambda[ixV]
            im_AB_V_ns = ABmag(im_Leff_V_ns,ewidth_nu[ixV])
            print "min/max/tot im_AB_V_ns:",im_AB_V_ns.min(),im_AB_V_ns.max(),ABmag(im_Leff_V_ns.sum(),ewidth_nu[ixV])
            print "real tot AB:",fbb['FILTERS'].data.field('AB_mag_nonscatter%d'%c)[ixV]     

            im_Vega_W1minusW2 = im_Vega_W1 - im_Vega_W2 
            im_Vega_W2minusW3 = im_Vega_W2 - im_Vega_W3 
            im_Av = im_AB_V - im_AB_V_ns

            print "min/max/tot im_Leff_W1:",im_Leff_W1.min(),im_Leff_W1.max(),im_Leff_W1.sum()
            print "real tot=",fbb['FILTERS'].data.field('L_lambda_eff%d'%c)[ixw1]*ewidth_lambda[ixw1]
            print "min/max im_Vega_W1:",im_Vega_W1.min(),im_Vega_W1.max()
            print "min/max/tot im_Leff_W2:",im_Leff_W2.min(),im_Leff_W2.max(),im_Leff_W2.sum()
            print "real tot=",fbb['FILTERS'].data.field('L_lambda_eff%d'%c)[ixw2]*ewidth_lambda[ixw2]
            print "min/max im_Vega_W2:",im_Vega_W2.min(),im_Vega_W2.max()
            print "min/max/tot im_Leff_W3:",im_Leff_W3.min(),im_Leff_W3.max(),im_Leff_W3.sum()
            print "real tot=",fbb['FILTERS'].data.field('L_lambda_eff%d'%c)[ixw3]*ewidth_lambda[ixw3]
            print "min/max im_Vega_W3:",im_Vega_W3.min(),im_Vega_W3.max()
            print "min/max/tot im_Leff_V:",im_Leff_V.min(),im_Leff_V.max(),im_Leff_V.sum()
            print "real tot=",fbb['FILTERS'].data.field('L_lambda_eff%d'%c)[ixV]*ewidth_lambda[ixV]
            print "min/max im_AB_V:",im_AB_V.min(),im_AB_V.max()
            print "min/max/tot im_Leff_V_ns:",im_Leff_V_ns.min(),im_Leff_V_ns.max(),im_Leff_V_ns.sum()
            print "real tot=",fbb['FILTERS'].data.field('L_lambda_eff_nonscatter%d'%c)[ixV]*ewidth_lambda[ixV]
            print "min/max im_AB_V_ns:",im_AB_V_ns.min(),im_AB_V_ns.max()
            print "\nmin/max/median im_Vega_W1minusW2:",im_Vega_W1minusW2.min(),im_Vega_W1minusW2.max(),np.median(im_Vega_W1minusW2)
            print "min/max/median im_Vega_W2minusW3:",im_Vega_W2minusW3.min(),im_Vega_W2minusW3.max(),np.median(im_Vega_W1minusW2)
            print "min/max/median im_Av:",im_Av.min(),im_Av.max(),np.median(im_Av)

            if bh1_xpix[c,i]>=0 and bh1_ypix[c,i]>=0:
                w1w2_bhpix[0,c,i] = im_Vega_W1minusW2[bh1_ypix[c,i],bh1_xpix[c,i]]
                w2w3_bhpix[0,c,i] = im_Vega_W2minusW3[bh1_ypix[c,i],bh1_xpix[c,i]]
                Av_bhpix[0,c,i] = im_Av[bh1_ypix[c,i],bh1_xpix[c,i]]
                print "bh1_xpix w1w2: ",w1w2_bhpix[0,c,i]
                print "bh1_xpix w2w3: ",w2w3_bhpix[0,c,i]
                print "bh1_xpix Av: ",Av_bhpix[0,c,i]
            if bh2_xpix[c,i]>=0 and bh2_ypix[c,i]>=0:
                w1w2_bhpix[1,c,i] = im_Vega_W1minusW2[bh2_ypix[c,i],bh2_xpix[c,i]]
                w2w3_bhpix[1,c,i] = im_Vega_W2minusW3[bh2_ypix[c,i],bh2_xpix[c,i]]
                Av_bhpix[1,c,i] = im_Av[bh2_ypix[c,i],bh2_xpix[c,i]]
                print "bh2_xpix w1w2: ",w1w2_bhpix[1,c,i]
                print "bh2_xpix w2w3: ",w2w3_bhpix[1,c,i]
                print "bh2_xpix Av: ",Av_bhpix[1,c,i]

    if make_image:
        print "make_image not coded yet"

        plt.clf()
        plt.close()
        im[im==0] = 1.0e-20
        ax = fig.add_subplot(subp0+c)
        ax.set_xlim(0,xpix)
        ax.set_ylim(0,ypix)
        ax.imshow(np.log10(im),vmin=-1,vmax=1.1*np.log10(im).max())
        if quickimage: ax.plot(bh_los_pix[i,c,1],bh_los_pix[i,c,0],'rx',markersize=8)
        fig.subplots_adjust(wspace=0.3,hspace=0.2)
        fig.savefig(path+'/quickimage_snap%d.eps'%snap,format='eps',dpi=1200)
        plt.clf()
        plt.cla()
        plt.close()


    return w1w2_bhpix,w2w3_bhpix,Av_bhpix

def plot_mag_img(im,xpix,ypix,bh1_xpix,bh2):
    return

def load_infofile(path,fbase='bh_info',fext='',skip_snap0=True,tmax_postmrg=-1):

    assert fbase in ['bh_info','bh_sep','star_gas_info','lum_info'], \
        "invalid filebase: %s.txt"%fbase
    fname='%s_%s.txt'%(fbase,fext) if len(fext)>0 else '%s.txt'%fbase

    with open('%s/%s'%(path,fname)) as f:
        tmpdata = np.loadtxt(f,unpack=True,dtype='float64')
        if skip_snap0:
            snap = tmpdata[0,:]
            snapmask = (snap>0)
            tmpdata = tmpdata[:,snapmask]
        if tmax_postmrg>=0:
            time = tmpdata[1,:]
        if fbase=='bh_info' or fbase=='star_gas_info':
            return tmpdata
        elif fbase=='bh_sep':
            sepsnap,septime,bh_3d_sep = tmpdata[:3]
            bh_proj_sep = np.array([col for col in tmpdata[3:]])
            return sepsnap,septime,bh_3d_sep,bh_proj_sep
        else:
            lsnap,ltime,lbol_grid,lbol_absorbed = tmpdata[:4]
            ltot_out = np.array([col for col in tmpdata[4:11]])
            lir = np.array([col for col in tmpdata[11:]])
            return lsnap,ltime,lbol_grid,lbol_absorbed,ltot_out,lir

def load_wise(maindir='/oasis/projects/nsf/hvd115/lblecha/',
              subdir='q1_fg0.3_sunruns/all',fext='',ncam=7,skip_snap0=True): 
    
    path=maindir+'/'+subdir
    fname12='w1w2_color_%s.txt'%fext if len(fext)>0 else 'w1w2_color.txt'
    fname23='w2w3_color_%s.txt'%fext if len(fext)>0 else 'w2w3_color.txt'
    w1w2_file = glob.glob("%s/%s"%(path,fname12))
    w2w3_file = glob.glob("%s/%s"%(path,fname23))
    if len(w1w2_file)==1 and len(w2w3_file)==1:
        tmpdata = np.loadtxt(w1w2_file[0],unpack=True)
        if skip_snap0:
            snapnum_wise = tmpdata[0,:]
            snapmask = (snapnum_wise>0)
            tmpdata = tmpdata[:,snapmask]
        snapnum_wise = tmpdata[0,:]
        time_wise = tmpdata[1,:]   
        w1w2 = tmpdata[2:,:]   

        tmpdata_b = np.loadtxt(w2w3_file[0],unpack=True)
        if skip_snap0:
            snapnum_w2w3 = tmpdata_b[0,:]
            snapmask = (snapnum_w2w3>0)
            tmpdata_b = tmpdata_b[:,snapmask]
        snapnum_w2w3 = tmpdata_b[0,:]
        time_w2w3 = tmpdata_b[1,:]   
        w2w3 = tmpdata_b[2:,:]   

        assert snapnum_wise.size==snapnum_w2w3.size, \
            "array size mismatch: w1w2_color.txt, w2w3_color.txt (%d %d)"%(snapnum_wise.size,snapnum_w2w3.size)
    else:
        print "Could not read WISE colors from file. Calculating..."
        print 'path=',path
        tmpdata = wise.wise_colors(path=path,subdir='',ncam=ncam,write_to_txt=True)
        snapnum_wise,time_wise,w1w2,w2w3 = tmpdata
        if skip_snap0:
            print snapnum_wise
            snapmask = (snapnum_wise>0)
            snapnum_wise = snapnum_wise[snapmask]
            time_wise = time_wise[snapmask]
            w1w2 = w1w2[:,snapmask]
            w2w3 = w2w3[:,snapmask]

    return time_wise, snapnum_wise, w1w2, w2w3


def load_wise_lums(path='/oasis/projects/nsf/hvd115/lblecha/q1_fg0.3_sunruns/all',ncam=7): 
    
    fname1='W1_Llambda.txt' 
    fname2='W2_Llambda.txt' 
    w1_file = glob.glob("%s/%s"%(path,fname1))
    w2_file = glob.glob("%s/%s"%(path,fname2))
    if len(w1_file)==1 and len(w2_file)==1:
        tmpdata = np.loadtxt(w1_file[0],unpack=True)
        snapnum_w1 = tmpdata[0,:]
        time_w1 = tmpdata[1,:]   
        w1 = tmpdata[2:,:]   

        tmpdata_b = np.loadtxt(w2_file[0],unpack=True)
        snapnum_w2 = tmpdata_b[0,:]
        time_w2 = tmpdata_b[1,:]   
        w2 = tmpdata_b[2:,:]   

    return time_w1, snapnum_w1, w1, w2

def plot(path='/oasis/projects/nsf/hvd115/lblecha/',
         subdir='q0.5_fg0.3_allrx10_sunruns/whole_sed_lowres_newcode',
         xlim=(),skip_snap0=True, nerrskip=2, cam=range(7),
         squish=True, dust_to_metal_ratio=0.4,grid_nh_res=0,pubstyle=False,
         panel_labels=False):

    ## q0.5_fg0.3_allrx10 (old): tmrg = 1.3531/0.7 = 1.93 Gyr
    ## q0.5_fg0.3: tmrg = 1.56181/0.7 = 2.23 Gyr

    mpl.rcParams.update({'font.size': 13})
    carr = ['m','darkblue','b','darkcyan','cyan','darkgreen','green']
    ncam=len(cam)

    #nerrskip2=np.int(np.floor(nerrskip/2.0))
    nerrskip2=1

    ### load txt files:
    tmpdata = load_infofile(path+'/'+subdir,fbase='bh_info',skip_snap0=skip_snap0)
    #snap,time,bh1x,bh1y,bh1z,bh1_lbol,bh2x,bh2y,bh2z,bh2_lbol = tmpdata
    snap,time,bh1x,bh1y,bh1z,bh1_mdot,bh1_lbol,bh1_lbol_sunsed,bh1_mass,bh2x,bh2y,bh2z,bh2_mdot,bh2_lbol,bh2_lbol_sunsed,bh2_mass = tmpdata

    dt = np.append(0.0,time[1:]-time[:-1])
    #dt = np.append(time[1:]-time[:-1],0.0)

    nsnaps=len(snap)
    ix_2bh = np.where(bh2_lbol==bh2_lbol)[0]
    nsnaps_2bh=len(snap[ix_2bh])
    tmrg = (time[bh2_lbol!=bh2_lbol])[0]
    print "tmrg = ",tmrg
    lagn = bh1_lbol
    lagn[ix_2bh] = bh1_lbol[ix_2bh] + bh2_lbol[ix_2bh]
    lagn_tiled = np.tile(lagn,(ncam,1))
    lagn_sun = bh1_lbol_sunsed
    lagn_sun[ix_2bh] = bh1_lbol_sunsed[ix_2bh] + bh2_lbol_sunsed[ix_2bh] 
    lagn_sun_tiled = np.tile(lagn_sun,(ncam,1))
    bh_mass = bh1_mass
    bh_mass[ix_2bh] = bh1_mass[ix_2bh] + bh2_mass[ix_2bh]
    ledd = calc_lbol_edd(bh_mass)
    fedd = np.zeros(lagn.size)
    fedd[ledd>0] = lagn[ledd>0]/ledd[ledd>0]
    print "min/max/med lagn:",lagn.min(),lagn.max(),np.median(lagn)
    print "min/max/med bh1_mass:",bh1_mass.min(),bh1_mass.max(),np.median(bh1_mass)
    print "min/max/med bh2_mass:",bh2_mass.min(),bh2_mass.max(),np.median(bh2_mass)
    print "min/max/med ledd:",ledd.min(),ledd.max(),np.median(ledd)
    print "min/max/med fedd:",fedd.min(),fedd.max(),np.median(fedd)

    sepsnap,septime,bh_3d_sep,bh_proj_sep = load_infofile(path+'/'+subdir,fbase='bh_sep',skip_snap0=skip_snap0)
    sep_ix_2bh = np.where(bh_3d_sep==bh_3d_sep)[0]
    assert nsnaps_2bh==0 or len(sep_ix_2bh) == nsnaps_2bh, \
        "bh array size mismatch: ncam=%d, nsnaps_2bh=%d (=%d in bh_sep)"%(ncam,nsnaps_2bh,len(sep_ix_2bh))

    tmpdata = load_infofile(path+'/'+subdir,fbase='star_gas_info',skip_snap0=skip_snap0)
    sgsnap,sgtime,mstar,mgas,mmetals,sfr = tmpdata
    fgas = 1.0*mgas / (mgas + mstar)
    fdust = 1.0*mmetals*dust_to_metal_ratio / mgas
    ssfr = 1.0*sfr/mstar
    assert len(snap) == len(sgsnap), \
        "array size mismatch in bh_info.txt, star_gas_info.txt: %d, %d"%(snap.size,sgsnap.size)

    tmpdata = load_infofile(path+'/'+subdir,fbase='lum_info',skip_snap0=skip_snap0)
    lsnap,ltime,lbol_grid,lbol_absorbed,ltot_out,lir = tmpdata
    assert lsnap.size == snap.size, \
        "array size mismatch in bh_info.txt, lum_info.txt: %d, %d"%(snap.size,lsnap.size)

    if (lagn_sun/lbol_grid).max() > 1:
        print "Error: snaps with lagn>ltot:"
        print lsnap[lagn_sun>lbol_grid]
        return

    time_wise, snapnum_wise, w1w2, w2w3 = load_wise(maindir=path,subdir=subdir,ncam=ncam,skip_snap0=skip_snap0)
    wmask = np.in1d(snapnum_wise,snap)
    print "retaining %d of %d elements of snapnum_wise also in snap."%(snapnum_wise[wmask].size,
                                                                       snapnum_wise.size)
    nowise = np.in1d(snap,snapnum_wise,invert=True)
    print "%d of %d elements of snap are not in snapnum_wise."%(snap[nowise].size,snap.size)    
    print snap[nowise]

    snapnum_wise = snapnum_wise[wmask]
    time_wise = time_wise[wmask]
    w1w2 = w1w2[:,wmask]
    w2w3 = w2w3[:,wmask]


    if len(xlim)==0:
        #xlim=(time.min()-0.01,time.max()+0.01)
        xlim=(time.min()+0.02,time.max())
        extra=''
    else: extra = '_zoom'


    initarr=np.zeros((ncam,nsnaps))-1
    bh1_xpix = copy(initarr)
    bh1_ypix = copy(initarr)
    bh1_Ngas_grid = copy(initarr)
    bh1_Ngas_tot_grid = copy(initarr)
    bh1_Ngas_tot_aux = copy(initarr)
    bh2_xpix = copy(initarr)
    bh2_ypix = copy(initarr)
    bh2_Ngas_grid = copy(initarr)
    bh2_Ngas_tot_grid = copy(initarr)
    bh2_Ngas_tot_aux = copy(initarr)

    for c in cam:
        rstr='_grid_res_%.3d'%(1000*grid_nh_res) if grid_nh_res>0 else ''
        with open(path+subdir+'/Ngas_bh_cam%d%s.txt'%(c,rstr)) as f:
            tmpdata = np.loadtxt(f,unpack=True,dtype='float64')
            if skip_snap0:
                sn = tmpdata[0,:]
                snapmask = (sn>0)
                tmpdata = tmpdata[:,snapmask]
            sn,bh1_xpix[c,:],bh1_ypix[c,:],bh1_Ngas_tot_aux[c,:],bh1_Ngas_tot_grid[c,:],bh1_Ngas_grid[c,:],bh2_xpix[c,:],bh2_ypix[c,:],bh2_Ngas_tot_aux[c,:],bh2_Ngas_tot_grid[c,:],bh2_Ngas_grid[c,:]=tmpdata

    bh1_Ngas_aux = bh1_Ngas_tot_aux/2.0
    ## don't include BHs that are outside the FOV:
    bh1_Ngas_aux[(bh1_xpix<0)|(bh1_ypix<0)]=np.nan
    #ix = [np.where((bh1_xpix[:,n]>=0)&(bh1_ypix[:,n]>=0))[0] for n in range(nsnaps)]
    #med_bh1_Ngas_aux = np.median(bh1_Ngas_aux[ix[n],:], axis=0)
    med_bh1_Ngas_aux = np.nanmedian(bh1_Ngas_aux, axis=0)
    min_bh1_Ngas_aux = np.nanmin(bh1_Ngas_aux, axis=0)
    max_bh1_Ngas_aux = np.nanmax(bh1_Ngas_aux, axis=0)
    yerr_bh1_aux=[-min_bh1_Ngas_aux+med_bh1_Ngas_aux,max_bh1_Ngas_aux-med_bh1_Ngas_aux]

    med_bh1_Ngas_grid = np.nanmedian(bh1_Ngas_grid, axis=0)
    min_bh1_Ngas_grid = np.nanmin(bh1_Ngas_grid, axis=0)
    max_bh1_Ngas_grid = np.nanmax(bh1_Ngas_grid, axis=0)
    yerr_bh1_grid=[-min_bh1_Ngas_grid+med_bh1_Ngas_grid,max_bh1_Ngas_grid-med_bh1_Ngas_grid]

    if nsnaps_2bh > 0:
        bh2_Ngas_aux = bh2_Ngas_tot_aux/2.0
        bh2_Ngas_aux[(bh2_xpix<0)|(bh2_ypix<0)]=np.nan
        #ix = [np.where((bh2_xpix[:,n]>=0)&(bh2_ypix[:,n]>=0))[0] for n in range(nsnaps)]
        ##med_bh2_Ngas_aux = np.array([ np.median(bh2_Ngas_aux[:,n]) for n in range(nsnaps) ])
        med_bh2_Ngas_aux = np.nanmedian(bh2_Ngas_aux, axis=0)
        min_bh2_Ngas_aux = np.nanmin(bh2_Ngas_aux, axis=0)
        max_bh2_Ngas_aux = np.nanmax(bh2_Ngas_aux, axis=0)
        yerr_bh2_aux=[-min_bh2_Ngas_aux+med_bh2_Ngas_aux,max_bh2_Ngas_aux-med_bh2_Ngas_aux]

        med_bh2_Ngas_grid = np.nanmedian(bh2_Ngas_grid, axis=0)
        min_bh2_Ngas_grid = np.nanmin(bh2_Ngas_grid, axis=0)
        max_bh2_Ngas_grid = np.nanmax(bh2_Ngas_grid, axis=0)
        yerr_bh2_grid=[-min_bh2_Ngas_grid+med_bh2_Ngas_grid,max_bh2_Ngas_grid-med_bh2_Ngas_grid]

    meanw12 = np.mean(w1w2,axis=0)
    minw12 = np.min(w1w2,axis=0)
    maxw12 = np.max(w1w2,axis=0)
    meanw23 = np.mean(w2w3,axis=0)

    bh_proj_sep_negvals = copy(bh_proj_sep)
    bh_proj_sep_negvals[bh_proj_sep!=bh_proj_sep] = -1

    flbol_grid = lagn_sun/lbol_grid
    flbol_grid_tiled = np.tile(flbol_grid,(ncam,1))
    flbol_out = lagn_sun_tiled/ltot_out
    print "flbol_grid (shape, max): ", flbol_grid.shape, flbol_grid.max()
    print "flbol_out (shape, max): ", flbol_out.shape, flbol_out.max()

    if flbol_out.max()>0:
        print "\nMin fLbol with W1W2>0.8: %g"%flbol_out[w1w2>0.8].min() if w1w2.max()>0.8 else "\n(no snaps with W1W2>0.8)"
        print "Min fLbol with W1W2>0.5: %g"%flbol_out[w1w2>0.5].min() if w1w2.max()>0.5 else "\n(no snaps with W1W2>0.8)"
        print "Min fLbol_grid with mean W1W2>0.8: %g"%flbol_grid[meanw12>0.8].min() if w1w2.max()>0.8 else "\n(no snaps with W1W2>0.8)"
        print "Min fLbol_grid with mean W1W2>0.5: %g"%flbol_grid[meanw12>0.5].min() if w1w2.max()>0.5 else "\n(no snaps with W1W2>0.8)"
        print "Max flbol_out each cam:",flbol_out.max(axis=1)
        mean_flbol_out = flbol_out.mean(axis=0)
        print "\n(Note: for small number of cams, ok if mean flbol_out is slightly > 1.)"
        print "Max flbol_out averaged over cams:",mean_flbol_out.max()
        print "snaps with mean flbol_out > 1:",snap[mean_flbol_out>1]
        print "Max flbol_grid: %g\n"%flbol_grid.max()


    plt.ioff()
    fig = plt.figure(figsize=(5,5))
    plt.xlim(-2,135)
    plt.ylim(21.0,24.5)
    #plt.plot(bh_proj_sep_negvals,np.log10(bh1_Ngas_aux),'bo',ms=2,mec='b')
    #plt.plot(bh_proj_sep_negvals,np.log10(bh2_Ngas_aux),'mo',ms=2,mec='m')
    plt.plot(bh_proj_sep_negvals,np.log10(bh1_Ngas_grid),'bo',ms=2,mec='b')
    plt.plot(bh_proj_sep_negvals,np.log10(bh2_Ngas_grid),'mo',ms=2,mec='m')

    fig.savefig(path+subdir+'/bhsep_vs_Ngas_quickplot%s.pdf'%extra,format='pdf')
    #plt.close('all')
    plt.cla()
    plt.clf()
    

    ##ylim1=(2e42,3e46)
    #ylim1=(2e42,5e46)
    ## best for A1A0rx10 sim:
    ylim1=(2e42,1.5e46)

    ### plot BH sep and lbol ratio
    print "Plotting BH sep and lbol ratio..."
    fig = plt.figure(figsize=(5,6))
    ax1=fig.add_subplot(311)
    plt.xlim(time.min()-0.01,time.max()+0.01)
    plt.yscale('log')
    plt.ylim(0.001,120)
    #plt.xlabel(r'time [Gyr]')
    plt.ylabel(r'BH sep [kpc]',fontsize=14)
    if tmrg>0:
        ax1.plot([tmrg,tmrg],[0.001,100],'k',linewidth=0.75)
            
    if nsnaps_2bh > 0:
        min_psep = np.nanmin(bh_proj_sep[:,ix_2bh], axis=0)
        max_psep = np.nanmax(bh_proj_sep[:,ix_2bh], axis=0)
        ##med_psep = np.array([ np.median(bh_proj_sep[:,n]) for n in range(nsnaps_2bh) ])

        ax1.plot(time[ix_2bh], bh_3d_sep[ix_2bh],'k',linewidth=2)
        ax1.plot(time[ix_2bh], max_psep,'b',linewidth=0.8)
        ax1.plot(time[ix_2bh], min_psep,'b',linewidth=0.8)
        ##ax1.errorbar(time[ix_2bh], med_psep, yerr=(med_psep-min_psep,max_psep-med_psep),color='gray')

    ax2=fig.add_subplot(312)
    plt.xlim(time.min()-0.01,time.max()+0.01)
    plt.yscale('linear')
    ylim=(-1.0,2.8)
    plt.ylim(ylim)
    #plt.xlabel(r'time [Gyr]')
    plt.ylabel(r'SFR [log M$_{\odot}$ yr$^{-1}$]',fontsize=14)
    #print "ylim = ",ylim
    if tmrg>0:
        ax2.plot([tmrg,tmrg],ylim,'k',linewidth=0.75)
    ax2.plot(time, np.log10(sfr),'k',linewidth=2)

    ax3=fig.add_subplot(313)
    plt.xlim(time.min()-0.01,time.max()+0.01)
    plt.yscale('linear')
    ylim=(-1.8,1.8)
    plt.ylim(ylim)
    plt.xlabel(r'time [Gyr]')
    plt.ylabel(r'log L$_{\rm bol1}$/L$_{\rm bol2}$',fontsize=14)
    #print "ylim = ",ylim
    if tmrg>0:
        ax3.plot([tmrg,tmrg],ylim,'k',linewidth=0.75)

    if nsnaps_2bh > 0:    
        lbol_ratio = bh1_lbol[ix_2bh]/bh2_lbol[ix_2bh]
        ax3.plot(time[ix_2bh], np.log10(lbol_ratio),'k',linewidth=2)

    fig.suptitle(subdir)
    fig.subplots_adjust(bottom=0.1,left=0.18,right=0.95,hspace=0.26,top=0.94)
    fig.savefig(path+subdir+'/bhsep_sfr_lbolratio%s.pdf'%extra,format='pdf',dpi=300)
    plt.cla()
    plt.clf()


    ### plot Lbol_bh, Ngas_bh, and W1-W2
    print "Making plot 'obscured_agn'..."
    fig = plt.figure(figsize=(5,5.5))
    ax1 = fig.add_subplot(311)
    plt.yscale('log')
    plt.xlim(xlim)
    plt.ylim(ylim1)
    #plt.xlabel(r'time [Gyr]')
    #plt.ylabel(r'L$_{\rm bol}$ [erg s$^{-1}$]',fontsize=14)
    plt.ylabel(r'L$_{\rm AGN}$ [erg s$^{-1}$]',fontsize=14)
    if tmrg>0:
        ax1.plot([tmrg,tmrg],ylim1,'k',linewidth=0.75)
        
    #lw1=1.5
    #lw2=0.75
    lw1=2
    lw2=1
    ax1.plot(time, bh1_lbol,'b',linewidth=lw1,markersize=2,markeredgecolor='b')

    if nsnaps_2bh > 0:
        #ax1.plot([xlim[1]-17*0.01,xlim[1]-14*0.01],[0.8e46,0.8e46],'b-',linewidth=lw1)
        #ax1.text(xlim[1]-12*0.01,0.55e46,'BH1',color='b')

        ax1.plot(time[ix_2bh], bh2_lbol[ix_2bh],'m',
                 linewidth=lw2,markersize=2,markeredgecolor='m')
        #ax1.plot([xlim[1]-17*0.01,xlim[1]-14*0.01],[2e45,2e45],'m-',linewidth=lw2)
        #ax1.text(xlim[1]-12*0.01,1.3e45,'BH2',color='m')

    if panel_labels:
        font0=FontProperties()
        font=font0.copy()
        font.set_weight('bold')
        ## for pub version of A1A0rx10 only:
        ## snap 176
        plx=2.6e42
        off=0.03
        #plx=4.5e45
        ax1.text(1.721-off,plx,'a',fontsize=12,fontproperties=font)
        ## snap 191
        ax1.text(1.868-off,plx,'b',fontsize=12,fontproperties=font)
        ### snap 206
        #ax1.text(2.014-off,plx,'c',fontsize=12,fontproperties=font)
        ## snap 205
        ax1.text(2.004-off,plx,'c',fontsize=12,fontproperties=font)
        ### snap 220
        #ax1.text(2.15-off,plx,'d',fontsize=12,fontproperties=font)
        ## snap 219
        ax1.text(2.14-off,plx,'d',fontsize=12,fontproperties=font)
        font.set_weight('normal')


    ##ylim2=(2e21,5e24)
    ## best for A1A0rx10 sim with grid_nh_res=0.136:
    #ylim2=(2e21,2.1e24)
    ## best for A1A0rx10 sim with grid_nh_res=0.064:
    ylim2=(2e21,7e24)
    ##ylim2=(1e20,3e23)
    ax2 = fig.add_subplot(312)    
    plt.yscale('log')
    #plt.xlabel(r'time [Gyr]')
    plt.ylabel(r'N$_{\rm H}$ [cm$^{-2}$]',fontsize=14)
    plt.ylim(ylim2)
    plt.xlim(xlim)

    if tmrg>0:
        ax2.plot([tmrg,tmrg],ylim2,'k',linewidth=0.75)

    #ax2.errorbar(time,med_bh1_Ngas_aux,yerr=yerr_bh1, color='b',linewidth=lw2)
    #ax2.errorbar(time[::nerrskip],med_bh1_Ngas_aux[::nerrskip],
    #             yerr=(yerr_bh1[0][::nerrskip],yerr_bh1[1][::nerrskip]), 
    #             color='b',linestyle='None')
    ax2.errorbar(time[dt>0.09][::nerrskip2],med_bh1_Ngas_grid[dt>0.09][::nerrskip2],
                 yerr=(yerr_bh1_grid[0][dt>0.09][::nerrskip2],
                       yerr_bh1_grid[1][dt>0.09][::nerrskip2]), 
                 color='b',linestyle='None',linewidth=0.5)
    ax2.errorbar(time[dt<0.09][::nerrskip],med_bh1_Ngas_grid[dt<0.09][::nerrskip],
                 yerr=(yerr_bh1_grid[0][dt<0.09][::nerrskip],
                       yerr_bh1_grid[1][dt<0.09][::nerrskip]), 
                 color='b',linestyle='None',linewidth=0.5)

    if nsnaps_2bh > 0:
        #ax2.errorbar(time,med_bh2_Ngas_grid,yerr=yerr_bh2, color='m',linewidth=lw1)
        ax2.errorbar(time[dt>0.09][::nerrskip2],med_bh2_Ngas_grid[dt>0.09][::nerrskip2],
                     yerr=(yerr_bh2_grid[0][dt>0.09][::nerrskip2],
                           yerr_bh2_grid[1][dt>0.09][::nerrskip2]), 
                     color='m',linestyle='None',linewidth=0.5)
        ax2.errorbar(time[dt<0.09][::nerrskip],med_bh2_Ngas_grid[dt<0.09][::nerrskip],
                     yerr=(yerr_bh2_grid[0][dt<0.09][::nerrskip],
                           yerr_bh2_grid[1][dt<0.09][::nerrskip]), 
                     color='m',linestyle='None',linewidth=0.5)

    if not pubstyle: ax2.plot(time,med_bh1_Ngas_aux, color='g',linewidth=lw2)
    ax2.plot(time,med_bh1_Ngas_grid, color='b',linewidth=lw1)

    if nsnaps_2bh > 0:
        if not pubstyle: ax2.plot(time,med_bh2_Ngas_aux, color='darkorange',linewidth=lw1)
        ax2.plot(time[time<tmrg],med_bh2_Ngas_grid[time<tmrg], color='m',linewidth=lw2)

        #ax2.plot([xlim[1]-0.025,xlim[1]-0.022],[0.32e24,0.32e24],'m-',linewidth=lw2)
        #ax2.text(xlim[1]-0.22,0.25e24,'BH2 LOS',color='m')

    
    ax3 = fig.add_subplot(313)
    ylim3 = (0.1,1.7)
    plt.xlabel(r'time [Gyr]')
    plt.ylabel(r'W1-W2')
    plt.ylim(ylim3)
    plt.xlim(xlim)

    if tmrg>0:
        ax3.plot([tmrg,tmrg],ylim3,'k',linewidth=0.75)
    ax3.plot(xlim,[0.5,0.5],'k:',linewidth=1.0)
    ax3.plot(xlim,[0.8,0.8],'k-.',linewidth=1.0)
    ax3.errorbar(time_wise, meanw12, yerr=(meanw12-minw12,maxw12-meanw12),color='gray')
    ax3.plot(time_wise, meanw12, 'k', linewidth=1.5)
   
    if not pubstyle: fig.suptitle(subdir)
    fig.subplots_adjust(bottom=0.1,left=0.18,right=0.95,hspace=0.26,top=0.94)
    #fig.savefig(path+subdir+'/Ngas_bh.eps')
    pub_ex=extra+'_pubstyle' if pubstyle else extra
    fig.savefig(path+subdir+'/obscured_agn%s.pdf'%pub_ex,format='pdf')
    plt.cla()
    plt.clf()


    print "Making SFR vs L_AGN plot..."
    fig = plt.figure(figsize=(6,6))
    ax=fig.add_subplot(221)

    #xlim=(1.5,6.1)
    #ylim=(0,2.0)
    #plt.xlim(xlim)
    #plt.ylim(ylim)
    ax.set_xlim(41.5,46.5)
    ax.set_ylim(-1,2.8)
    ax.set_xlabel(r'log($L_{\rm AGN}$)')
    ax.set_ylabel(r'log(SFR)')

    ax.plot(np.log10(lagn),np.log10(sfr),'ko',ms=2)
    ax.plot(np.log10(lagn[meanw12>0.5]),np.log10(sfr[meanw12>0.5]),'g^',ms=4,mec='g')
    ax.plot(np.log10(lagn[meanw12>0.8]),np.log10(sfr[meanw12>0.8]),'rs',ms=4,mec='r')
    
    ax=fig.add_subplot(222)

    ax.set_xlim(-3.5,0.1)
    #ax.set_xticklabels(np.arange(-3.5,0.1,dtype='int'))
    ax.set_xticks(np.arange(-3.5,0.1,dtype='int'))
    ax.set_ylim(-1,2.8)
    ax.set_xlabel(r'log($L_{\rm AGN}/L_{\rm tot}$)')
    #ax.set_ylabel(r'log(SFR)')

    #ax.plot(np.log10(lagn/lbol_grid),np.log10(sfr),'ko',ms=2)
    ax.plot(np.log10(flbol_grid),np.log10(sfr),'ko',ms=2)
    ax.plot(np.log10(lagn[meanw12>0.5]/lbol_grid[meanw12>0.5]),
            np.log10(sfr[meanw12>0.5]),'g^',ms=4,mec='g')
    ax.plot(np.log10(lagn[meanw12>0.8]/lbol_grid[meanw12>0.8]),
            np.log10(sfr[meanw12>0.8]),'rs',ms=4,mec='r')   

    j11mask = j11_wedge_cut(meanw12,meanw23,
                            w1lum_lim=w1lum_lim,w2lum_lim=w2lum_lim,w1lum=w1lum,w2lum=w2lum)

    ax=fig.add_subplot(223)

    #ax.set_xlim(40.5,46.5)
    ax.set_xlim(41.5,46.5)
    ax.set_ylim(-1,2.8)
    ax.set_xlabel(r'log($L_{\rm AGN}$)')
    ax.set_ylabel(r'log(SFR)')

    ax.plot(np.log10(lagn),np.log10(sfr),'ko',ms=2)
    ax.plot(np.log10(lagn[j11mask]),np.log10(sfr[j11mask]),'bD',ms=4,mec='b')
    
    ax=fig.add_subplot(224)

    #ax.set_xlim(-5,0.0)
    ax.set_xlim(-3.5,0.1)
    ax.set_xticks(np.arange(-3.5,0.1,dtype='int'))
    ax.set_ylim(-1,2.8)
    ax.set_xlabel(r'log($L_{\rm AGN}/L_{\rm tot}$)')
    #ax.set_ylabel(r'log(SFR)')

    ax.plot(np.log10(flbol_grid),np.log10(sfr),'ko',ms=2)
    ax.plot(np.log10(lagn[j11mask]/lbol_grid[j11mask]),
            np.log10(sfr[j11mask]),'bD',ms=4,mec='b')

    fig.suptitle(subdir)
    fig.subplots_adjust(bottom=0.1,left=0.12,right=0.96,wspace=0.25,hspace=0.25,top=0.94)
    fig.savefig(path+subdir+'/sfr_vs_lagn%s.pdf'%extra,format='pdf')
    plt.clf()
    plt.cla()


    ### plot evolution in w1-w2 vs w2-w3 color space
    print "Making WISE color-color plot..."
    fs=(5,3.5) if squish else (5,5)
    fig = plt.figure(figsize=fs)
    ax=fig.add_subplot(111)

    xlim=(1.5,6.1)
    #ylim=(0,2.0)
    ylim=(-0.05,1.8)
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.xlabel(r'W2-W3')
    plt.ylabel(r'W1-W2')

    #print np.mean(w1w2,axis=0).shape
    #print np.mean(w2w3,axis=0).shape
    #print np.mean(w1w2,axis=0)
    #print np.mean(w2w3,axis=0)

    plot_wise_agn_cuts(ax,xlim=xlim)

    ### Jarrett et al. 2011 criteria:
    #ax.plot([2.2,4.2],[1.7,1.7],color='gray',ls='--')
    #ax.plot([2.2,2.2],[0.5,1.7],color='gray',ls='--')
    #ax.plot([4.2,4.2],[0.8,1.7],color='gray',ls='--')
    #ax.plot([2.2,4.2],[0.5,0.8],color='gray',ls='--')
    ### Stern et al. criteria:
    #ax.plot(xlim,[0.5,0.5],color='gray',ls=':')
    #ax.plot(xlim,[0.8,0.8],color='gray')

    #ax.plot(np.mean(w2w3,axis=0), np.mean(w1w2,axis=0),'k',linewidth=1)
    ax.plot(np.mean(w2w3,axis=0), np.mean(w1w2,axis=0),'darkblue',linewidth=1.75)
    ax.plot(np.mean(w2w3[:,0]), np.mean(w1w2[:,0]),'g',marker='o')

    ix = np.where(np.abs(time-tmrg)==np.min(np.abs(time-tmrg)))[0]
    #print ix
    #print time[ix]
    #print np.mean(w2w3[:,ix[0]])
    #print np.mean(w1w2[:,ix[0]])
    ax.plot(np.mean(w2w3[:,ix[0]]), np.mean(w1w2[:,ix[0]]),'orange',marker='^')
    ax.plot(np.mean(w2w3[:,-1]), np.mean(w1w2[:,-1]),'r',marker='s')

    fig.suptitle(subdir)
    fig.subplots_adjust(bottom=0.18,left=0.18,right=0.95,hspace=0.3,top=0.92)
    #fig.savefig(path+subdir+'/w1w2_vs_w2w3.eps',format='eps',dpi=300)
    fig.savefig(path+subdir+'/w1w2_vs_w2w3%s.pdf'%extra,format='pdf')
    plt.clf()
    plt.cla()


    ### make scatter plots of quantities vs w1-w2 color:
    print "Making w1w2 scatter plots..."
    ## need these quantities: 
    ## ltot, lir, sfr, mstar, mgas, mmetals, dust-to-metal ratio

    ## from grid.fits [hdu 3 - GRIDDATA, hdu 1 - MAKEGRID, hdu 6 - PARTICLEDATA]: 
    ## gas mass: f[3].header['M_G_TOT'] or f[1].header['M_G_TOT']
    ## sfr: f[3].header['SFR_TOT'] or f[1].header['SFR_TOT']
    ## mmetals: f[3].data['mass_metals'].sum() or f[1].header['M_metals_tot']
    ## cell volume: f[3].data['cell_volume']
    ## age: f[6].data['age']
    ## mstar: f[6].data['mass'][age==age].sum() 

    ## from mcrx.fits [hdu 48 - INTEGRATED_QUANTITIES, hdu 6 - MAKEGRID]:
    ## ltot: f[48].header['L_bol_grid']-f[48].header['L_bol_absorbed']
    ##    OR camera-angle dependent: f[48].header['L_OUTn'] for n in icam
    ## lir: f[48].header['L_IRn'] for n in icam
    ## gas mass: f[6].header['M_G_TOT'] 
    ## sfr: f[6].header['SFR_TOT'] 
    ## mmetals: f[6].header['M_metals_tot']

    #compare_arr = np.dstack((np.log10(lagn/ltot_out).flatten(),w1w2.flatten()))[0]
    #print compare_arr.shape
    #print compare_arr[:,1].shape
    #print snapnum_wise.shape
    #ix = np.where( (compare_arr[:,0] > -1.0)&(compare_arr[:,1] < 0.2) )[0]
    #compare_snap = np.tile(snapnum_wise,ncam)
    #print len(ix)
    #if len(ix)>0:
    #    print "Found w1w2 outliers."
    #    print "snaps: ",compare_snap[ix]
    #    print "min w1w2: ",compare_arr[ix,1]
    #    print "max lagn/ltot_out: ",compare_arr[ix,0]
    #    #raise ValueError()

    fig = plt.figure(figsize=(8,8))

    w1w2_scatter_plot(fig.add_subplot(321),w1w2,np.log10(lagn), xlim=(41.5,47),
                      xlabel=r'L$_{AGN}$ [log erg s$^{-1}$]')
    #w1w2_scatter_plot(fig.add_subplot(322),w1w2,np.log10(flbol_out), 
    #                  altvals=(np.log10(flbol_grid),),
    #                  xlim=(-5,0.5), xlabel=r'log (L$_{\rm AGN}$/L$_{\rm tot}$)')
    w1w2_scatter_plot(fig.add_subplot(322),w1w2,np.log10(flbol_grid_tiled), 
                      xlim=(-5,0.5), xlabel=r'log (L$_{\rm AGN}$/L$_{\rm tot}$)')
    w1w2_scatter_plot(fig.add_subplot(323),w1w2,np.log10(lir), xlim=(41.5,47),
                      xlabel=r'L$_{IR}$ [log erg s$^{-1}$]')
    w1w2_scatter_plot(fig.add_subplot(324),w1w2,np.log10(lir/ltot_out), xlim=(-5,0.5),
                      xlabel=r'log (L$_{IR}$/L$_{tot}$)')
    w1w2_scatter_plot(fig.add_subplot(325),w1w2,np.log10(ltot_out), altvals=(np.log10(lbol_grid),),
                      xlim=(41.5,47), xlabel=r'L$_{tot}$ [log erg s$^{-1}$]')
    w1w2_scatter_plot(fig.add_subplot(326),w1w2,np.log10(fedd),
                      xlim=(-5,0.05), xlabel=r'f$_{\rm Edd}$',ms0=1.5)

    fig.suptitle(subdir)
    fig.subplots_adjust(bottom=0.1,left=0.12,right=0.94,hspace=0.35,wspace=0.3,top=0.92)
    fig.savefig(path+subdir+'/w1w2_lum_scatter_plots%s.pdf'%extra,format='pdf')
    plt.clf()
    plt.cla()

    fig = plt.figure(figsize=(8,8))

    w1w2_scatter_plot(fig.add_subplot(321),w1w2,np.log10(sfr), xlim=(-1,2.8),
                      xlabel=r'SFR [log M$_{\odot}$ yr$^{-1}$]')
    w1w2_scatter_plot(fig.add_subplot(322),w1w2,np.log10(ssfr), xlim=(-11.9,-8),
                      xlabel=r'sSFR [log yr$^{-1}$]')
    w1w2_scatter_plot(fig.add_subplot(323),w1w2,fgas, xlim=(0,0.3),
                      xlabel=r'Mgas/(Mgas+M*)')
    #w1w2_scatter_plot(fig.add_subplot(324),w1w2,fdust, xlim=(0.004,0.015),
    #                  xlabel=r'Mdust/Mgas')
    w1w2_scatter_plot(fig.add_subplot(324),w1w2,np.log10(mmetals*dust_to_metal_ratio), 
                      xlim=(7,9), xlabel=r'Mdust [log M$_{\odot}$]')
    if nsnaps_2bh <= 0:
        w1w2_scatter_plot(fig.add_subplot(325),w1w2,np.log10(med_bh1_Ngas_grid),
                          xlim=(21.1,24.5), xlabel=r'$N_{\rm H} [log cm^{-2}$]',ms0=1.5)
    else:
        w1w2_scatter_plot(fig.add_subplot(325),w1w2,np.log10(med_bh1_Ngas_grid),
                          altvals=(np.log10(med_bh2_Ngas_grid),), xlim=(21,24.5), 
                          xlabel=r'$N_{\rm H}$ [log cm$^{-2}$]',ms0=1.5,altms=(2.5,))
        w1w2_scatter_plot(fig.add_subplot(326),w1w2[:,ix_2bh],np.log10(bh_proj_sep[:,ix_2bh]),
                          altvals=(np.log10(bh_3d_sep[ix_2bh]),), xlim=(-2,2.2),
                          xlabel=r'BH sep [log kpc]',plot_val_last=True)

    fig.suptitle(subdir)
    fig.subplots_adjust(bottom=0.1,left=0.12,right=0.94,hspace=0.35,wspace=0.3,top=0.90)
    fig.savefig(path+subdir+'/w1w2_other_scatter_plots%s.pdf'%extra,format='pdf')
    plt.clf()
    plt.cla()


    return

    ### plot quantities extracted along the LOS to each BH:
    #ix=snapnum_wise.tolist().index(212)
    #w1w2_bhpix,w2w3_bhpix,Av_bhpix = img_mags_bhlos(path+subdir,212,bh1_xpix[:,ix],bh1_ypix[:,ix],
    #                                                bh2_xpix[:,ix],bh2_ypix[:,ix])

    #print "w1w2_bhpix:",w1w2_bhpix.shape
    #print "w1w2_bhpix:",w1w2_bhpix
    #print "w2w3_bhpix:",w2w3_bhpix.shape
    #print "w2w3_bhpix:",w2w3_bhpix
    #print "Av_bhpix:",Av_bhpix.shape
    #print "Av_bhpix:",Av_bhpix


    w1w2_bhpix,w2w3_bhpix,Av_bhpix = img_mags_bhlos(path+subdir,snap,bh1_xpix,bh1_ypix,bh2_xpix,bh2_ypix)

    mean_w1w2_bh1 = np.mean(w1w2_bhpix[0,:,:],axis=0)
    mean_w2w3_bh1 = np.mean(w2w3_bhpix[0,:,:],axis=0)
    mean_Av_bh1 = np.mean(Av_bhpix[0,:,:],axis=0)

    mean_w1w2_bh2 = np.mean(w1w2_bhpix[1,:,:],axis=0)
    mean_w2w3_bh2 = np.mean(w2w3_bhpix[1,:,:],axis=0)
    mean_Av_bh2 = np.mean(Av_bhpix[1,:,:],axis=0)

    wise.plot_LOSerr_vs_t(path,subdir,'bh_los_W1W2',time,mean_w1w2_bh1,
                          np.min(w1w2_bhpix[0,:,:],axis=0),np.max(w1w2_bhpix[0,:,:],axis=0),
                          meanval2=mean_w1w2_bh2, 
                          minval2=np.min(w1w2_bhpix[1,:,:],axis=0),
                          maxval2=np.max(w1w2_bhpix[1,:,:],axis=0),
                          lcolor='b',ecolor='b',lcolor2='m',ecolor2='m',ylabel='W1-W2 (BH LOS)')

    wise.plot_LOSerr_vs_t(path,subdir,'bh_los_W2W3',time,mean_w2w3_bh1,
                          np.min(w2w3_bhpix[0,:,:],axis=0),np.max(w2w3_bhpix[0,:,:],axis=0),
                          meanval2=mean_w2w3_bh2, 
                          minval2=np.min(w2w3_bhpix[1,:,:],axis=0),
                          maxval2=np.max(w2w3_bhpix[1,:,:],axis=0),
                          lcolor='b',ecolor='b',lcolor2='m',ecolor2='m',ylabel='W2-W3 (BH LOS)')

    wise.plot_LOSerr_vs_t(path,subdir,'bh_los_Av',time,mean_Av_bh1,
                          np.min(Av_bhpix[0,:,:],axis=0),np.max(Av_bhpix[0,:,:],axis=0),
                          meanval2=mean_Av_bh2, 
                          minval2=np.min(Av_bhpix[1,:,:],axis=0),
                          maxval2=np.max(Av_bhpix[1,:,:],axis=0),
                          lcolor='b',ecolor='b',lcolor2='m',ecolor2='m',ylabel='Av (BH LOS)')


    ### plot evolution along BH LOS in w1-w2 vs w2-w3 color space
    fig = plt.figure(figsize=(5,5))
    ax=fig.add_subplot(111)

    xlim=(1.5,6.1)
    ylim=(0,2.0)
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.xlabel(r'W2-W3 (BH LOS)')
    plt.ylabel(r'W1-W2 (BH_LOS)')

    plot_wise_agn_cuts(ax,xlim=xlim)

    ax.plot(mean_w2w3_bh1, mean_w1w2_bh1,'b',linewidth=1)
    ax.plot(mean_w2w3_bh2, mean_w1w2_bh2,'k',linewidth=1)

    ax.plot(mean_w2w3_bh1[0],mean_w1w2_bh1[0],'g',marker='o')
    ax.plot(mean_w2w3_bh2[0],mean_w1w2_bh2[0],'g',marker='o')

    ix = np.where(np.abs(time-tmrg)==np.min(np.abs(time-tmrg)))[0]
    ax.plot(mean_w2w3_bh1[ix[0]], mean_w1w2_bh1[ix[0]],'c',marker='^')
    ax.plot(mean_w2w3_bh2[ix[0]], mean_w1w2_bh2[ix[0]],'c',marker='^')

    ax.plot(mean_w2w3_bh1[-1], mean_w1w2_bh1[-1],'r',marker='s')

    fig.suptitle(subdir)
    fig.subplots_adjust(bottom=0.18,left=0.18,right=0.95,hspace=0.3,top=0.92)
    fig.savefig(path+subdir+'/bh_los_w1w2_vs_w2w3%s.pdf'%extra,format='pdf')
    plt.clf()
    plt.cla()
    plt.close()




def color_color_compare(path='/oasis/projects/nsf/hvd115/lblecha/',
                        subdirA='q1_fg0.3_sunruns/all',
                        subdirB='q1_fg0.3_sunruns/all_agnx0',alt_wedge=True,
                        xlim=(),skip_snap0=True, nerrskip=2, fname_ext='',sim_name=True,
                        title='',pubstyle=True,cam=range(7),dust_to_metal_ratio=0.4,
                        squish=True,tmax=-1,agnscale=True):


    ### get tmrg:
    ncam=len(cam)
    tmpdata = load_infofile(path+'/'+subdirA,fbase='bh_info',skip_snap0=skip_snap0)
    snap,time,bh1x,bh1y,bh1z,bh1_mdot,bh1_lbol,bh1_lbol_sunsed,bh1_mass,bh2x,bh2y,bh2z,bh2_mdot,bh2_lbol,bh2_lbol_sunsed,bh2_mass = tmpdata
    ix_2bh = np.where(bh2_lbol==bh2_lbol)[0]
    tmrg = (time[bh2_lbol!=bh2_lbol])[0]
    print "tmrg = ",tmrg

    lagnA = bh1_lbol
    lagnA[ix_2bh] = bh1_lbol[ix_2bh] + bh2_lbol[ix_2bh]


    timeA, snapnumA, w1w2A, w2w3A = load_wise(maindir=path,subdir=subdirA,ncam=ncam,skip_snap0=skip_snap0)

    timeB, snapnumB, w1w2B, w2w3B = load_wise(maindir=path,subdir=subdirB,ncam=ncam,skip_snap0=skip_snap0)

    #wmask = np.in1d(snapnum_wise,snap)
    mask = np.in1d(snapnumA,snapnumB)
    print "retaining %d of %d elements of A also in B."%(snapnumA[mask].size, snapnumA.size)

    nomatch = np.in1d(snapnumB,snapnumA,invert=True)
    print "%d of %d elements of snapnumB are not in snapnumA."%(snapnumB[nomatch].size,snapnumB.size)    
    print snapnumB[nomatch]

    snapnumA = snapnumA[mask]
    timeA = timeA[mask]
    w1w2A = w1w2A[:,mask]
    w2w3A = w2w3A[:,mask]
    lagnA = lagnA[mask]
    

    if tmax<=0: tmax = timeA.max()

    ### plot evolution in w1-w2 vs w2-w3 color space
    print "Making WISE color-color comparison plot..."

    mpl.rcParams.update({'font.size': 13})

    #fs=(5,3.5) if squish else (5,5)
    fs=(5.5,3.75) if squish else (5.5,5.5)
    fig = plt.figure(figsize=fs)
    ax=fig.add_subplot(111)

    xlim=(1.5,6.1)
    ylim=(-0.05,1.8)
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.xlabel(r'W2-W3')
    plt.ylabel(r'W1-W2')

    #plot_wise_agn_cuts(ax,xlim=xlim,c='lightgray')
    plot_wise_agn_cuts(ax,xlim=xlim,include_new=True,alt_wedge=True)

    msA=7
    msB=7

    ixmrg = np.where(np.abs(timeA-tmrg)==np.min(np.abs(timeA-tmrg)))[0]

    xB = np.mean(w2w3B,axis=0)
    yB = np.mean(w1w2B,axis=0)

    xA = np.mean(w2w3A,axis=0)
    yA = np.mean(w1w2A,axis=0)

    ### orig version:
    ##ax.plot(np.mean(w2w3B,axis=0), np.mean(w1w2B,axis=0),'k',linewidth=1)
    ##ax.plot(np.mean(w2w3A,axis=0), np.mean(w1w2A,axis=0),'darkblue',linewidth=1.75)
    ###

    ### time-scale only version:
    ##lc = colorline(x, y, cmap='Greys_r')
    #lcB = colorline(x, y, z=timeB, cmap='gist_yarg_r',linewidth=2,
    #                norm=plt.Normalize(timeB.min(),tmax+0.2))
    ##plt.colorbar(lc)

    ##print timeA.min(),tmax,tmrg
    ##lcA = colorline(x, y, z=timeA, cmap=cmaps.plasma,norm=plt.Normalize(timeA.min(),timeA.max()))
    #lcA = colorline(x, y, z=timeA, cmap=cmaps.viridis,norm=plt.Normalize(timeA.min(),tmax))
    ##lc = colorline(x, y, z=timeA, cmap=cmaps.viridis,norm=plt.Normalize(timeA.min(),timeA.max()))
    ##lc = colorline(x, y, cmap=cmaps.viridis)
    ###lc = colorline(x, y, cmap=cmaps.inferno)
    ###lc = colorline(x, y, cmap=cmaps.magma)
    #plt.colorbar(lcA,label='time [Gyr]')

    ### choose time or lagn as color scale: ###
    if agnscale:
        fname_ext = fname_ext+'_agnscale'
        lagn_min=41.5
        lagn_max=46.5
        ax.plot(xB, yB,'k',linewidth=1)
        lcA = colorline(xA, yA, z=np.log10(lagnA), cmap=cmaps.viridis,norm=plt.Normalize(lagn_min,lagn_max))
        #plt.colorbar(lcA,label=r'log L$_{\rm AGN}$ [erg s$^{-1}$]')
        plt.colorbar(lcA,label=r'log L$_{\rm AGN}$ [erg s$^{-1}$]',ticks=[42,43,44,45,46])
    else:
        lcB = colorline(xB, yB, z=timeB, cmap='gist_yarg_r',linewidth=2,
                        norm=plt.Normalize(timeB.min(),tmax+0.2))
        lcA = colorline(xA, yA, z=timeA, cmap=cmaps.viridis,norm=plt.Normalize(timeA.min(),tmax))
        plt.colorbar(lcA,label='time [Gyr]')

        

    ##p1,=ax.plot(np.mean(w2w3B[:,0]), np.mean(w1w2B[:,0]),'g',marker='o',ms=msB,linestyle='None')
    ##ax.plot(np.mean(w2w3A[:,0]), np.mean(w1w2A[:,0]),'g',marker='o',ms=msA)
    p1,=ax.plot(np.mean(w2w3B[:,0]), np.mean(w1w2B[:,0]),'c',marker='o',ms=msB,linestyle='None')
    ax.plot(np.mean(w2w3A[:,0]), np.mean(w1w2A[:,0]),'c',marker='o',ms=msA)

    p2,=ax.plot(np.mean(w2w3B[:,ixmrg[0]]), np.mean(w1w2B[:,ixmrg[0]]),'orange',marker='^',ms=msB+1,linestyle='None')
    ax.plot(np.mean(w2w3A[:,ixmrg[0]]), np.mean(w1w2A[:,ixmrg[0]]),'orange',marker='^',ms=msA+1)

    p3,=ax.plot(np.mean(w2w3B[:,-1]), np.mean(w1w2B[:,-1]),'r',marker='s',ms=msB,linestyle='None')
    ax.plot(np.mean(w2w3A[:,-1]), np.mean(w1w2A[:,-1]),'r',marker='s',ms=msA)

    ##ax.legend((p1,p2,p3),(r't$_0$',r't$_{\rm mrg}$',r't$_{\rm final}$'),
    ##         loc='upper right',numpoints=1,fontsize=14,handletextpad=0.05)
    #ax.legend((p1,p2,p3),(r't$_0$',r't$_{\rm mrg}$',r't$_{\rm final}$'),
    #          loc='upper right',numpoints=1,fontsize=14,handletextpad=0.05,
    #          labelspacing=0.25,handlelength=1.2)
    ax.legend((p1,p2,p3),(r't$_0$',r't$_{\rm mrg}$',r't$_{\rm final}$'),
              loc='upper left',numpoints=1,fontsize=14,handletextpad=0.05,
              labelspacing=0.25,handlelength=1.2)
    top=0.92
    if title: 
        fig.suptitle(title)
    else:
        #fig.suptitle("%s + %s"%(subdirA, subdirB))
        if pubstyle:
            fig.suptitle(get_sim_name(subdirA))
        else:
            fig.suptitle(subdirA)
    #plt.text(5,1.5,get_sim_name(subdirA))
    fig.subplots_adjust(bottom=0.14,left=0.14,right=0.96,top=top)
    fig.savefig(path+subdirA+'/w1w2_vs_w2w3_compare%s.pdf'%fname_ext,format='pdf')
    plt.clf()
    plt.cla()
    plt.close()

def multicolored_lines():
    """
    http://nbviewer.ipython.org/github/dpsanders/matplotlib-examples/blob/master/colorline.ipynb
    http://matplotlib.org/examples/pylab_examples/multicolored_line.html
    """

    x = np.linspace(0, 4. * np.pi, 100)
    y = np.sin(x)
    fig, ax = plt.subplots()
    lc = colorline(x, y, cmap='hsv')
    plt.colorbar(lc)
    plt.xlim(x.min(), x.max())
    plt.ylim(-1.0, 1.0)
    plt.show()


def colorline(x, y, z=None, cmap='copper', norm=plt.Normalize(0.0, 1.0),
              linewidth=3, alpha=1.0):
    """
    http://nbviewer.ipython.org/github/dpsanders/matplotlib-examples/blob/master/colorline.ipynb
    http://matplotlib.org/examples/pylab_examples/multicolored_line.html
    Plot a colored line with coordinates x and y
    Optionally specify colors in the array z
    Optionally specify a colormap, a norm function and a line width
    """

    # Default colors equally spaced on [0,1]:
    if z is None:
        z = np.linspace(0.0, 1.0, len(x))

    # Special case if a single number:
    # to check for numerical input -- this is a hack
    if not hasattr(z, "__iter__"):
        z = np.array([z])

    z = np.asarray(z)

    segments = make_line_segments(x, y)
    lc = mcoll.LineCollection(segments, array=z, cmap=cmap, norm=norm,
                              linewidth=linewidth, alpha=alpha)

    ax = plt.gca()
    ax.add_collection(lc)

    return lc


def make_line_segments(x,y):
    """
    Create list of line segments from x and y coordinates, in the correct format
    for LineCollection: an array of the form numlines x (points per line) x 2 (x
    and y) array
    """

    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    return segments


def agnx0_compare(path='/oasis/projects/nsf/hvd115/lblecha/',
                  subdir='q1_fg0.3_sunruns/all',
                  subdir0='q1_fg0.3_sunruns/all_agnx0',
                  xlim=(),skip_snap0=True, nerrskip=2, fname_ext='',
                  title=True,cam=range(7),dust_to_metal_ratio=0.4,squish=False):

    ### get tmrg: ###
    ncam=len(cam)
    tmpdata_bh = load_infofile(path+'/'+subdir,fbase='bh_info',skip_snap0=skip_snap0)
    snap,time,bh1x,bh1y,bh1z,bh1_mdot,bh1_lbol,bh1_lbol_sunsed,bh1_mass,bh2x,bh2y,bh2z,bh2_mdot,bh2_lbol,bh2_lbol_sunsed,bh2_mass = tmpdata
    ix_2bh = np.where(bh2_lbol==bh2_lbol)[0]
    tmrg = (time_bh[bh2_lbol!=bh2_lbol])[0]
    print "tmrg = ",tmrg

    ### load fiducial data: ###

    time, snap, w1w2, w2w3 = load_wise(maindir=path,subdir=subdir,ncam=ncam,skip_snap0=skip_snap0)
    assert time_bh.size==time.size,"mismatch in bh_info vs wise for agnx1: %g %g"%(time_bh.size,time.size)
    #if len(time_bh) != len(time):
    #    wmask = np.in1d(snap,snap_bh)
    #    print "retaining %d of %d elements of snap also in snap_bh."%(snap[wmask].size,snap.size)
    #    nowise = np.in1d(snap_bh,snap,invert=True)
    #    print "%d of %d elements of snap_bh are not in snap."%(snap_bh[nowise].size,snap_bh.size)    
    #    print snap_bh[nowise]
    #    snap = snap[wmask]
    #    time = time[wmask]
    #    w1w2 = w1w2[:,wmask]
    #    w2w3 = w2w3[:,wmask]

    tmpdata_lum = load_infofile(path+'/'+subdir,fbase='lum_info',skip_snap0=skip_snap0)
    lsnap,ltime,lbol_grid,lbol_absorbed,ltot_out,lir = tmpdata_lum
    #print [x.shape for x in tmpdata_lum]
    assert lsnap.size == snap.size, \
        "array size mismatch in bh_info.txt, lum_info.txt for agnx1: %d, %d"%(snap.size,lsnap.size)


    ### load agnx0 data: ###

    time0, snap0, w1w20, w2w30 = load_wise(maindir=path,subdir=subdir0,ncam=ncam,skip_snap0=skip_snap0)

    tmpdata_lum0 = load_infofile(path+'/'+subdir0,fbase='lum_info',skip_snap0=skip_snap0)
    lsnap0,ltime0,lbol_grid0,lbol_absorbed0,ltot_out0,lir0 = tmpdata_lum0
    assert lsnap0.size == snap0.size, \
        "array size mismatch in lum_info vs wise for agnx0: %d, %d"%(lsnap.size,snap0.size)

    if snap0.size != snap.size:
        mask = np.in1d(snap,snap0)
        print "retaining %d of %d elements of snap also in snap0."%(snap[mask].size, snap.size)
        nofid = np.in1d(snap0,snap,invert=True)
        print "%d of %d elements of snap0 are not in snap."%(snap0[nofid].size,snap0.size)    
        time = time[mask]
        snap = snap[mask]
        w1w2 = w1w2[:,mask]
        w2w3 = w2w3[:,mask]
        tmpdata_bh = tmpdata_bh[:,mask]
        snap_bh,time_bh,bh1x,bh1y,bh1z,bh1_lbol,bh1_mass,bh2x,bh2y,bh2z,bh2_lbol,bh2_mass=tmpdata_bh
        lsnap,ltime,lbol_grid,lbol_absorbed = (x[mask] for x in tmpdata_lum[:4])
        ltot_out,lir = (x[:,mask] for x in tmpdata_lum[4:])
        #print "lum data, after cuts:",[x.shape for x in lsnap,ltime,lbol_grid,lbol_absorbed,ltot_out,lir]


    ### calculate a few numbers ###

    meanw12 = np.mean(w1w2,axis=0)
    meanw23 = np.mean(w2w3,axis=0)

    meanw120 = np.mean(w1w20,axis=0)
    meanw230 = np.mean(w2w30,axis=0)
    meanlir = np.mean(lir,axis=0)
    meanlir0 = np.mean(lir0,axis=0)

    j11mask = j11_wedge_cut(w1w2,w2w3,
                            w1lum_lim=w1lum_lim,w2lum_lim=w2lum_lim,w1lum=w1lum,w2lum=w2lum)
    j11mask_mean = j11_wedge_cut(meanw12,meanw23,
                                 w1lum_lim=w1lum_lim,w2lum_lim=w2lum_lim,w1lum=w1lum,w2lum=w2lum)
    j11mask0 = j11_wedge_cut(w1w20,w2w30,
                             w1lum_lim=w1lum_lim,w2lum_lim=w2lum_lim,w1lum=w1lum,w2lum=w2lum)
    j11mask_mean0 = j11_wedge_cut(meanw120,meanw230,
                                  w1lum_lim=w1lum_lim,w2lum_lim=w2lum_lim,w1lum=w1lum,w2lum=w2lum)


    fagn_lir = (lir - lir0)/lir
    fagn_meanlir = (meanlir - meanlir0)/meanlir
    fagn_ltot_out = (ltot_out - ltot_out0)/ltot_out
    fagn_lbol = (lbol_grid - lbol_grid0)/lbol_grid
    fagn_w1w2 = (w1w2 - w1w20)/w1w2
    fagn_w2w3 = (w2w3 - w2w30)/w2w3
    fagn_meanw12 = (meanw12 - meanw120)/meanw12
    fagn_meanw23 = (meanw23 - meanw230)/meanw23

    lagn = bh1_lbol
    lagn[ix_2bh] = bh1_lbol[ix_2bh] + bh2_lbol[ix_2bh]
    bh_mass = bh1_mass
    bh_mass[ix_2bh] = bh1_mass[ix_2bh] + bh2_mass[ix_2bh]
    ledd = calc_lbol_edd(bh_mass)
    fedd = np.zeros(lagn.size)
    fedd[ledd>0] = lagn[ledd>0]/ledd[ledd>0]


    ###fagn_lw1 = (lw1 - lw10)/lw1
    ###fagn_lw2 = (lw2 - lw20)/lw2
    ###fagn_lw3 = (lw3 - lw30)/lw3

    ## plots:
    ## w1w20 vs w1w2
    ## w2w30 vs w2w3
    ## fagn_lir vs fagn_w1w2
    ## fagn_lir vs fagn_w2w3
    ## fagn_lir vs w1w2
    ## fagn_lbol vs w1w2
    ### fagn_lw1 vs w1w2
    ## fagn_lir vs sfr (similar to current plot)

    plt.close()
    plt.clf()

    mpl.rcParams.update({'font.size': 11})

    fig = plt.figure(figsize=(8,8))
    kwargs = dict(color='k',marker='o',ms=2,linestyle='None')
    j11kwargs = dict(color='m',marker='D',ms=3,linestyle='None',mec='m')

    ax = fig.add_subplot(331)
    ax.set_xlim(-0.1,1.7)
    ax.set_ylim(-0.1,1.7)
    ax.set_xlabel(r'W1W2$_{\rm noAGN}$')
    ax.set_ylabel(r'W1W2$_{\rm AGN}$')
    ax.plot([-0.2,2],[-0.2,2],'k')
    ax.plot(w1w20,w1w2,**kwargs)
    ax.plot(w1w20[j11mask],w1w2[j11mask],**j11kwargs)

    ax = fig.add_subplot(332)
    ax.set_xlim(1.5,7)
    ax.set_ylim(1.5,7)
    ax.set_xlabel(r'W2W3$_{\rm noAGN}$')
    ax.set_ylabel(r'W2W3$_{\rm AGN}$')
    ax.plot([1,8],[1,8],'k')
    ax.plot(w2w30,w2w3,**kwargs)
    ax.plot(w2w30[j11mask],w2w3[j11mask],**j11kwargs)

    ax = fig.add_subplot(333)
    ax.set_xlim(-2.9,0)
    ax.set_ylim(-0.1,1.7)
    ax.set_xlabel(r'L$_{\rm IR,AGN}$/L$_{\rm IR,tot}$')
    ax.set_ylabel(r'W1W2$_{\rm AGN}$ - W1W2$_{\rm noAGN}$')
    #ax.plot(np.log10(fagn_lir),fagn_w1w2,**kwargs)
    ax.plot(np.log10(fagn_lir),w1w2-w1w20,**kwargs)
    ax.plot(np.log10(fagn_lir[j11mask]),w1w2[j11mask]-w1w20[j11mask],**j11kwargs)

    ax = fig.add_subplot(334)
    ax.set_xlim(-0.1,1.7)
    ax.set_ylim(-0.1,1.7)
    ax.set_xlabel(r'W1W2$_{\rm noAGN}$')
    ax.set_ylabel(r'W1W2$_{\rm AGN}$ - W1W2$_{\rm noAGN}$')
    ax.plot(w1w20,w1w2-w1w20,**kwargs)
    ax.plot(w1w20[j11mask0],w1w2[j11mask0]-w1w20[j11mask0],'bo',ms=3,ls='None',mec='b')
    ax.plot(w1w20[j11mask],w1w2[j11mask]-w1w20[j11mask],**j11kwargs)

    ax = fig.add_subplot(335)
    ax.set_xlim(-2.9,0)
    ax.set_ylim(-2,2)
    ax.set_xlabel(r'log L$_{\rm IR,AGN}$/L$_{\rm IR,tot}$')
    ax.set_ylabel(r'W2W3$_{\rm AGN}$ - W2W3$_{\rm noAGN}$')
    #ax.plot(np.log10(fagn_lir),fagn_w2w3,**kwargs)
    ax.plot(np.log10(fagn_lir),w2w3-w2w30,**kwargs)
    ax.plot(np.log10(fagn_lir[j11mask]),w2w3[j11mask]-w2w30[j11mask],**j11kwargs)

    ax = fig.add_subplot(336)
    ax.set_xlim(-2.0,0)
    ax.set_ylim(-0.1,1.7)
    ax.set_xlabel(r'log L$_{\rm IR,AGN}$/L$_{\rm IR,tot}$')
    ax.set_ylabel(r'W1W2$_{\rm AGN}$')
    ax.plot(np.log10(fagn_lir),w1w2,**kwargs)
    ax.plot(np.log10(fagn_lir[j11mask]),w1w2[j11mask],**j11kwargs)

    ax = fig.add_subplot(337)
    ax.set_xlim(-2.9,0)
    ax.set_ylim(-0.1,1.7)
    ax.set_xlabel(r'log L$_{\rm bol,AGN}$/L$_{\rm bol,tot}$')
    ax.set_ylabel(r'W1W2$_{\rm AGN}$')
    #ax.plot(np.log10(fagn_lbol),meanw12,**kwargs)
    ax.plot(np.log10(fagn_ltot_out),w1w2,**kwargs)
    ax.plot(np.log10(fagn_ltot_out[j11mask]),w1w2[j11mask],**j11kwargs)
    
    #ax = fig.add_subplot(338)
    #ax.set_xlim(-3,0)
    #ax.set_ylim(-3,0)
    #ax.plot([-3,0],[-3,0],'k')
    ##ax.plot(np.log10(fagn_lbol),np.log10(fagn_meanlir),**kwargs)
    #ax.plot(np.log10(fagn_ltot_out),np.log10(fagn_lir),**kwargs)
    #ax.plot(np.log10(fagn_ltot_out[j11mask]),np.log10(fagn_lir[j11mask]),**j11kwargs)
    

    ax = fig.add_subplot(338)
    ax.set_xlim(-2.9,0)
    ax.set_xlabel(r'log f$_{\rm Edd}$')
    ax.set_ylabel(r'log f$_{\rm IR,AGN}$')
    ax.plot(np.log10(fedd),np.log10(fagn_meanlir),**kwargs)
    ax.plot(np.log10(fedd[j11mask_mean]),np.log10(fagn_meanlir[j11mask_mean]),**j11kwargs)

    ax = fig.add_subplot(339)
    ax.set_xlim(-2.9,0)
    ax.set_xlabel(r'log f$_{\rm Edd}$')
    ax.set_ylabel(r'log f$_{\rm bol,AGN}$')
    ax.plot(np.log10(fedd),np.log10(fagn_lbol),**kwargs)
    ax.plot(np.log10(fedd[j11mask_mean]),np.log10(fagn_lbol[j11mask_mean]),**j11kwargs)

    fig.suptitle(subdir)
    fig.subplots_adjust(left=0.08,right=0.95,top=0.95,wspace=0.38,hspace=0.3)
    #fig.subplots_adjust(bottom=0.1,left=0.12,right=0.94,hspace=0.35,wspace=0.3,top=0.92)
    fig.savefig(path+'/'+subdir+'/agnx0_compare_plots.pdf',format='pdf')
    plt.clf()
    plt.cla()
    plt.close()

def set_agn_maskarr_list(var,var1,var2,lims):

    agn_maskarr_list,dual_agn_maskarr_list = (),()
    for lim in lims:
        #agn_maskarr =( (var>=lim),(var<lim) ) 
        #dual_agn_maskarr =( ((var1>=lim)&(var2>=lim)),
        #                     ( ((var1>=lim)&(var2<lim))|
        #                       ((var1<lim)&(var2>=lim)) ),
        #                     ((var1<lim)&(var2<lim)) )
        agn_maskarr =(var>=lim)
        dual_agn_maskarr =( ((var1>=lim)&(var2>=lim)),
                             ( ((var1>=lim)&(var2<lim))|
                               ((var1<lim)&(var2>=lim)) ) )
        #dual_agn_maskarr =( ((var1>=lim)&(var2>=lim)),
        #                     ( ((var1>=lim)&(var2<lim))|
        #                       ((var1<lim)&(var2>=lim)) ),
        #                    ((var1<lim)&(var2<lim)) )
        agn_maskarr_list = agn_maskarr_list+(agn_maskarr,)
        dual_agn_maskarr_list = dual_agn_maskarr_list+(dual_agn_maskarr,)

    return agn_maskarr_list, dual_agn_maskarr_list

def set_agn_wmaskarr_list(var,var1,var2,lims,wmask_list,w1lum_lim=0.0,w2lum_lim=0.0,
                          w1lum=np.array([]),w2lum=np.array([]),old_fmt=False):

    agn_wmaskarr_list,dual_agn_wmaskarr_list = (),()
    #tmplbl = ()
    #wlbls=['w05', 'w08', 'j11']

    for iwise,wmask in enumerate(wmask_list):

        l=0
        ncam=wmask.shape[0]
        vtiled = np.tile(var,(ncam,1))
        v1tiled = np.tile(var1,(ncam,1))
        v2tiled = np.tile(var2,(ncam,1))
        for lim in lims:
            if w1lum_lim>0.0 or w2lum_lim>0.0:
                agn_wmaskarr =( ((vtiled>=lim)&(w1lum>=w1lum_lim)&(w2lum>=w2lum_lim)&(wmask)),
                                ((vtiled>=lim)&((w1lum<w1lum_lim)|(w2lum<w2lum_lim)|(~wmask))),
                                ((vtiled<lim)&(w1lum>=w1lum_lim)&(w2lum>=w2lum_lim)&(wmask)),
                                ((vtiled<lim)&((w1lum<w1lum_lim)|(w2lum<w2lum_lim)|(~wmask))) ) 
            else:
                agn_wmaskarr =( ((vtiled>=lim)&(wmask)),((vtiled>=lim)&(~wmask)),
                                ((vtiled<lim)&(wmask)),((vtiled<lim)&(~wmask)) ) 
            if old_fmt:
                dual_agn_wmaskarr =( ((v1tiled>=lim)&(v2tiled>=lim)&(wmask)),
                                     ((v1tiled>=lim)&(v2tiled>=lim)&(~wmask)),
                                     ( ( ((v1tiled>=lim)&(v2tiled<lim))|
                                         ((v1tiled<lim)&(v2tiled>=lim)) )&(wmask)) )
                dual_agn_wmaskarr = dual_agn_wmaskarr+((~dual_agn_wmaskarr[0])&
                                                       (~dual_agn_wmaskarr[1])&
                                                       (~dual_agn_wmaskarr[2]),)
            else:
                if w1lum_lim>0.0 or w2lum_lim>0.0:
                    dual_agn_wmaskarr =( ((v1tiled>=lim)&(v2tiled>=lim)&(w1lum>=w1lum_lim)&(w2lum>=w2lum_lim)&(wmask)),
                                         ((v1tiled>=lim)&(v2tiled>=lim)&((w1lum<w1lum_lim)|(w2lum<w2lum_lim)|(~wmask))),
                                         ( ( ((v1tiled>=lim)&(v2tiled<lim))|
                                             ((v1tiled<lim)&(v2tiled>=lim)) )&(w1lum>=w1lum_lim)&(w2lum>=w2lum_lim)&(wmask) ),
                                         ( ( ((v1tiled>=lim)&(v2tiled<lim))|
                                             ((v1tiled<lim)&(v2tiled>=lim)) )&((w1lum<w1lum_lim)|(w2lum<w2lum_lim)|(~wmask)) ),
                                         ((v1tiled<lim)&(v2tiled<lim)&(w1lum>=w1lum_lim)&(w2lum>=w2lum_lim)&(wmask)),
                                         ((v1tiled<lim)&(v2tiled<lim)&((w1lum<w1lum_lim)|(w2lum<w2lum_lim)|(~wmask))) )
                else:
                    dual_agn_wmaskarr =( ((v1tiled>=lim)&(v2tiled>=lim)&(wmask)),
                                         ((v1tiled>=lim)&(v2tiled>=lim)&(~wmask)),
                                         ( ( ((v1tiled>=lim)&(v2tiled<lim))|
                                             ((v1tiled<lim)&(v2tiled>=lim)) )&(wmask) ),
                                         ( ( ((v1tiled>=lim)&(v2tiled<lim))|
                                             ((v1tiled<lim)&(v2tiled>=lim)) )&(~wmask) ),
                                         ((v1tiled<lim)&(v2tiled<lim)&(wmask)),
                                         ((v1tiled<lim)&(v2tiled<lim)&(~wmask)) )


            agn_wmaskarr_list = agn_wmaskarr_list+(agn_wmaskarr,)
            dual_agn_wmaskarr_list = dual_agn_wmaskarr_list+(dual_agn_wmaskarr,)
            #tmplbl = tmplbl+('%s_lim%d'%(wlbls[iwise],l),)
            l=l+1

    #print tmplbl
    return agn_wmaskarr_list, dual_agn_wmaskarr_list


def assert_isclose(a,b,err):
    assert np.isclose(a,b), 'Error: %s (%g %g)'%(err,a,b)

def plot_dt_hist_onesim(sep_bins,dsep,dt_tot_hist,dt_hists,dtfrac_hists,
                        leg_lbls=('AGN','AGN (L1+2)','Dual'),eq_bins=True,ncam=7,
                        plot_titles=[],path='/oasis/projects/nsf/hvd115/lblecha/',
                        subdir=None,binstr=None,tmax_postmrg=0.1,fbase='agn',title=None,
                        leg_eachline=False, leg_eachplot=False,xlbl_eachline=True,
                        ylbl_eachline=True,carr=['b','darkblue','m'],lsarr=['-',':','-'],lwarr=[2,2,1]):

    if not isinstance(leg_lbls,tuple): leg_lbls=tuple(leg_lbls)
    if not isinstance(dt_hists,tuple): dt_hists=tuple(dt_hists)
    if not isinstance(dtfrac_hists,tuple): dt_hists=tuple(dtfrac_hists)
    assert len(leg_lbls)==len(dt_hists) and len(leg_lbls)==len(dtfrac_hists),\
        'Error: dt_hists, dtfrac_hists, and leg_lbls must have same length.'
    assert dt_tot_hist.shape == (ncam,len(sep_bins)-1),\
        'Error: dt_tot_hist must have shape (%d,%d)'%(ncam,len(sep_bins)-1)
    assert np.array([h.shape==(9,len(sep_bins)-1) for h in dt_hists+dtfrac_hists]).all(),\
        'Error: dt_hists,&dtfrac_hists elems must have shape (9,%d)'%(len(sep_bins)-1)
    assert len(carr)==len(leg_lbls), 'Please specify a unique color for each histogram line.'
    assert len(lsarr)==len(leg_lbls), 'Please specify a unique linestyle for each histogram line.'

    if len(plot_titles)==1: plot_titles=plot_titles*9
    if len(plot_titles)==0: plot_titles=['']*9

    mpl.rcParams.update({'font.size': 11})
    xlabel = 'proj. sep. [kpc]'
    if eq_bins:
        xvar = 0.5*dsep+sep_bins[:-1] 
        xlim = (sep_bins[0],130) if len(xvar)<=20 else (sep_bins[0]-2*dsep,130)
    else:
        xvar = 0.5+np.arange(len(sep_bins)-1)
        #xlim = (0,len(sep_bins)-1)
        xlim = (0.5,len(sep_bins)-1.5)

    if not eq_bins: xticklbl = [" "]+["%g"%(sep_bins[j]) for j in np.arange(1,len(sep_bins))]
    #tlim = (0, 0.4) if eq_bins else (0, 1.5)
    #hkwargs = dict(drawstyle='steps-post')
    hkwargs = dict(marker='o',mec='none',ms=4,drawstyle='steps-mid')
    chkwargs = dict(marker='o',mec='none',ms=4) if len(xvar)<=20 else dict(marker=None,ms=4)
    #hkwargs = dict(marker='o',ls='-',mec='none',ms=4)
    ptype_arr = ('t','tfrac','tcum') if len(xvar)<=20 else ('tcum',)

    if leg_eachline:
        leg_plots=[i for i in np.arange(9) if i%3==2] 
        leg_eachplot=False
    elif leg_eachplot:
        leg_eachline=False
        leg_plots=range(9)
    else: leg_plots=[8]

    if not title: title=subdir

    for ptype in ptype_arr:
        if ptype == 't': ylim = (0.001, 0.8) if eq_bins else (0.001, 1.5)
        if ptype == 'tfrac': ylim = (0, 1.01)
        if ptype == 'tcum': ylim = (0.001,2.2)
        cstr = '(>r)' if 'cum' in ptype else ''
        ylabel = r'$\Delta$t/t$_{\rm bin}$ %s'%cstr if 'frac' in ptype else r'$\Delta$t %s [Gyr]'%cstr
        fig = plt.figure(figsize=(8,8))
        for i in range(9):
            kwargs = hkwargs
            yarr = ()
            if not 'frac' in ptype: 
                y0 = np.nansum(dt_tot_hist,axis=0)/(1.0*ncam)
                if ptype == 'tcum': 
                    y0 = np.cumsum(y0[::-1])[::-1]
                    kwargs = chkwargs
            for y_j in range(len(dt_hists)):
                if ptype == 'tfrac':
                    yarr = yarr+(np.nan_to_num(dtfrac_hists[y_j][i,:]),)
                elif ptype=='t':
                    yarr = yarr+(dt_hists[y_j][i,:],)
                else:
                    yarr = yarr+(np.cumsum((dt_hists[y_j][i,:])[::-1])[::-1],)
                
            ax = fig.add_subplot(331+i)
            if i>5 or xlbl_eachline: ax.set_xlabel(xlabel)
            if i%3==0 or ylbl_eachline: ax.set_ylabel(ylabel)
            ax.set_xlim(xlim)
            #if not eq_bins: ax.set_xticklabels(xticklbl,fontsize=9)
            if not eq_bins: ax.set_xticklabels(xticklbl)
            ax.set_ylim(ylim)
            lh=None
            if not 'frac' in ptype:
                ax.set_yscale('log')
                p0,=ax.plot(xvar,y0,'lightgray',lw=2,**kwargs)   
                lh = (p0,)
            for y_j in range(len(yarr)):
                p,=ax.plot(xvar,yarr[y_j],carr[y_j],ls=lsarr[y_j],lw=lwarr[y_j],**kwargs)
                lh = (p,) if lh==None else lh+(p,)
            ll=leg_lbls if 'frac' in ptype else ('Total',)+leg_lbls
            if i in leg_plots: 
                ax.legend(lh, ll, fontsize=9, loc='upper right', numpoints=1, handletextpad=0.08)
            ax.set_title(plot_titles[i],fontsize=10)
        fig.subplots_adjust(wspace=0.33,hspace=0.35,left=0.08,right=0.96,bottom=0.08,top=0.94)
        fig.suptitle(title,fontsize=9) 
        plotname='%s/%s_%s_vs_projsep%s_tpost%g.pdf'%(path,ptype,fbase,binstr,tmax_postmrg)
        fig.savefig(plotname,format='pdf')
        plt.cla()
        plt.clf()


def plot_multisim_totals(frac_arr,leg_lbls=('Total','Early','Late'),ncam=7,nplots=9,xvar='',xarr=np.empty([0]),
                         ylbl='',ylim=(-0.02,1.02),plot_titles=[],outdir='/oasis/projects/nsf/hvd115/lblecha/',
                         fbase='f_agn',extra='',tmax_postmrg=0.1,latesep=10,carr=['k','g','m'],
                         msarr=['o','^','s'],mewarr=[1.5,1,1],mszarr=[6,6,6],title='',leg_eachline=False,
                         leg_eachplot=False,xlbl_eachline=True,leg_loc='upper right'):

    mpl.rcParams.update({'font.size': 12})
    assert (nplots>0 & nplots<10), 'Error: nplots must be a number 1-9.'
    nsim = frac_arr[0].shape[0]
    if len(frac_arr[0].shape)==1: 
        tmp=()
        for i in range(len(frac_arr)):
            tmp=tmp+(np.transpose(np.tile(np.array(frac_arr[i]),(nplots,1))),)
        frac_arr = tmp

    if nplots>6:
        figsz=(8,8)
        subp0=331
        lastrow0=6
        hsp=0.3+0.08*xlbl_eachline
        subp_adj=dict(wspace=0.34,hspace=hsp,left=0.08,right=0.96,
                      bottom=0.08,top=0.93)
    elif nplots>3 & nplots<=6:
        figsz=(8,5.5)
        subp0=231
        lastrow0=3
        hsp=0.2+0.1*xlbl_eachline
        subp_adj=dict(wspace=0.34,hspace=hsp,left=0.08,right=0.97,
                      bottom=0.1,top=0.92)
    elif nplots==3:
        figsz=(8,3)
        subp0=131
        lastrow0=0
        subp_adj=dict(wspace=0.36,left=0.08,right=0.97,bottom=0.2,top=0.85)
    elif nplots==2:
        figsz=(5.5,3)
        subp0=121
        lastrow0=0
        subp_adj=dict(wspace=0.36,left=0.12,right=0.96,bottom=0.16,top=0.84)
    else:
        figsz=(3,3)
        subp0=111
        lastrow0=0
        subp_adj=dict(left=0.2,right=0.94,bottom=0.18,top=0.88)

    assert len(leg_lbls)==len(frac_arr), 'Error: leg_lbls must have same len as frac_arr.'
    assert len(carr)==len(frac_arr), 'Error: carr must have same len as frac_arr.'
    assert len(msarr)==len(frac_arr), 'Error: msarr must have same len as frac_arr.'
    assert len(mewarr)==len(frac_arr), 'Error: mewarr must have same len as frac_arr.'
    assert len(mszarr)==len(frac_arr), 'Error: mszarr must have same len as frac_arr.'
    #assert np.array([h.shape==(nsim,9) for h in frac_arr]).all(),\
    #    'Error: frac_arr elems must have shape (%d,9), '%nsim
    #for h in frac_arr: print len(h.shape)
    assert np.array([h.shape==(nsim,nplots) for h in frac_arr]).all(),\
        'Error: frac_arr elems must have shape (%d,%d).'%(nsim,nplots)

    if leg_eachline and nplots>=3:
        leg_plots=[i for i in np.arange(9) if i%3==2] 
        leg_eachplot=False
    elif leg_eachplot:
        leg_eachline=False
        leg_plots=range(9)
    else: leg_plots=[nplots-1]

    #print xvar
    xarr=np.array(xarr)
    #print xarr.shape
    xlbl=np.empty(0)
    xtypes = np.array(['','lgssfr','lgsfr','lglagn','lgltot','flbol'])   
    xlbls = np.array([r'Sim #',r'log max(sSFR) [yr$^{-1}$]',
                      r'log max(SFR) [M$_{\rm sun}$ yr$^{-1}$]',
                      r'log max(L$_{\rm AGN}$) [erg s$^{-1}$]',
                      r'log max(L$_{\rm tot}$) [erg s$^{-1}$]',
                      r'max(L$_{\rm AGN}$ / L$_{\rm tot}$)'])
    if (not isinstance(xvar,list) and not isinstance(xvar,tuple) and 
        not isinstance(xvar,np.ndarray)): 
        if xvar=='':
            xarr = np.tile(np.arange(nsim),(nplots,1))
            xlbl = np.repeat('Sim #',nplots)
            xvar=np.repeat(xvar,nplots)
        else:
            xarr=np.tile(np.array(xarr),(nplots,1))
            xvar=np.repeat(xvar,nplots)
            xlbl=np.empty(0)
            for xv in xvar:
                assert xv in xtypes, 'Error: invalid xvar %s; choose one of %s'%(xv,', '.join(xtypes))
                xlbl = np.append(xlbl,xlbls[np.where(xtypes==xv)[0]])
    else:
        for xv in xvar:
            assert xv in xtypes, 'Error: invalid xvar %s; choose one of %s'%(xv,', '.join(xtypes))
            #print xlbl
            xlbl = np.append(xlbl,xlbls[np.where(xtypes==xv)[0]])
    #print xvar
    #print xarr.shape
    assert len(xvar)==nplots,'Error: xvar must have length %d, not %d.'%(nplots,len(xvar),)
    assert xarr.shape==(nplots,nsim),'Error: xarr must have shape (%d,%d).'%(nplots,nsim)

    fig = plt.figure(figsize=figsz)
    for i in range(nplots):
        dx = np.max(xarr[i,:])-np.min(xarr[i,:])
        xlim = (np.min(xarr[i,:])-0.1*dx,np.max(xarr[i,:])+0.1*dx)
        #print "xlim=",xlim
        ax = fig.add_subplot(subp0+i)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        if xlbl_eachline: 
            ax.set_xlabel(xlbl[i])
        else:
            if i>=lastrow0: ax.set_xlabel(xlbl[i])
        ax.set_ylabel(ylbl)
        if (ylim[1]-ylim[0])/1.0*len(ax.get_yticks())<0.1: ax.set_yticks(np.arange(ylim[0],ylim[1],0.1))
        if xvar[i]=='': ax.set_xticks(xarr[i,:])
        if xvar[i]=='lgssfr': 
            ax.set_xticks((-10.5,-10,-9.5,-9,-8.5))
            ax.set_xticklabels(('-10.5','','-9.5','','-8.5'))
        lh=None
        for y_j in range(len(frac_arr)):
            #p,=ax.plot(frac_arr[y_j][:,i],msarr[y_j],fillstyle='none',
            #           mec=carr[y_j],mew=mewarr[y_j],ms=mszarr[y_j])
            p,=ax.plot(xarr[i,:],frac_arr[y_j][:,i],msarr[y_j],fillstyle='none',
                       mec=carr[y_j],mew=mewarr[y_j],ms=mszarr[y_j])
            lh = (p,) if lh==None else lh+(p,)
        if len(plot_titles)>0: ax.set_title(plot_titles[i],fontsize=10)
        if i in leg_plots: 
            ax.legend(lh, leg_lbls, fontsize=9, loc=leg_loc, numpoints=1, handletextpad=0.08)
    #fig.subplots_adjust(wspace=0.34,hspace=hsp,left=0.08,right=0.96,
    #                    bottom=0.08,top=0.92)
    fig.subplots_adjust(**subp_adj)
    if nplots>1: fig.suptitle(title)
    if tmax_postmrg>0:
        plotname='%s/%s%s_lsep%g_tpost%g.pdf'%(outdir,fbase,extra,latesep,tmax_postmrg)
    else:
        plotname='%s/%s%s_lsep%g.pdf'%(outdir,fbase,extra,latesep)
    fig.savefig(plotname,format='pdf')
    plt.cla()
    plt.clf()


def errbars(hist,errtype='mad',valtype='median'):

    assert valtype in ('median','mean'),'Please choose valtype: "mean" or "median".'
    if valtype=='median': 
        assert errtype in ('full','iqr','mad','mad_sqrtn'),\
            'For valtype="median", please choose errtype: "full","iqr","mad","mad_sqrtn".'
    else: 
        assert errtype in ('full','std'),\
            'For valtype="mean", please choose errtype: "full","std".'

    if valtype=='median':
        height = np.nanmedian(hist,axis=0)
    else: 
        height = np.nanmean(hist,axis=0)

    if errtype=='full':
        ## error bars give full range of data (max,min)                                          
        ymin = np.nanmin(hist,axis=0)
        ymax = np.nanmax(hist,axis=0)
    elif errtype=='iqr':
        ## error bars give interquartile range of data (0.25-0.75)                               
        ymin = np.nanpercentile(hist,25,axis=0)
        ymax = np.nanpercentile(hist,75,axis=0)
    elif errtype=='mad':
        ## error bars give median absolute deviation                                             
        mad = (np.nanmedian( np.abs( hist - np.nanmedian(hist,axis=0) ), axis=0 ))
        ymin = height - mad
        ymax = height + mad
    elif errtype=='mad_sqrtn':
        mad_sqrtn = (np.nanmedian( np.abs( hist - np.nanmedian(hist,axis=0) ), axis=0 ) / 
                     np.sqrt(hist.shape[0]))
        ymin = height - mad_sqrtn
        ymax = height + mad_sqrtn
    elif errtype=='std':
        std = np.nanstd(hist,axis=0)
        ymin = height - std
        ymax = height + std

    print "\n in errbars:"
    print "height:",height
    print "min hist:",ymin
    print "max hist:",ymax
    print "errs[0]:",height-ymin
    print "errs[1]:",ymax-height
    print "\n"

    if np.isscalar(height):
        height= np.array([height])
    return (height-ymin, ymax-height)


def merger_phases(maindir='/oasis/projects/nsf/hvd115/lblecha/',subdir_arr='test', 
                  agnx0=False,obsc_lum_only=False,ngas_only=False,ngas_title='',
                  z=0, wise_fluxlims=[0.0, 0.0, 0.0, 0.0],wise_fluxlim_type='',
                  tmax_postmrg=0.1, latesep=10, xlim=(), skip_snap0=True, equal_bins=True,
                  dsep=10,cam=range(7),dust_to_metal_ratio=0.4, plot_3dsep=False,
                  use_logbins=False, ylog=False,extra='',errtype='mad',valtype_allsim='median',
                  plot_wedge_allsim=False,verbose=False,allsim_only=False,alt_wedge=True,
                  single_sim_plots=True, multisim_plots=True,sfrhist_w_grid=False,grid_nh_res=0,
                  omega_m=0.308, h=0.678):



    wflim_str = ''
    if z>0:
        omega_l = 1.0-omega_m
        cosmo = pc.Cosmology(omega_m, omega_l, 0.0, -1.0, h)
        dL = cosmo.lum_dist(z) * 1000.0*ac.KPC ## Mpc to cm

        if wise_fluxlim_type in ('SN3','SN10'):
            print 'Imposing cuts for WISE sensitivity limits for W1 & W2, assuming %s'%wise_fluxlim_type
            if wise_fluxlim_type=='SN3':
                ### total band fluxes (haven't coded for redshifted SEDs yet)
                #W1fluxlim = 2.157e-15 ## erg/s/cm^2
                #W2fluxlim = 3.605e-15 ## erg/s/cm^2
                ### monochromatic fluxes F_lambda 
                W1fluxlim = 3.256e-16 ## W/m/cm^2
                W2fluxlim = 3.45873e-16 ## W/m/cm^2
            elif wise_fluxlim_type=='SN10':
                ### total band fluxes (haven't coded for redshifted SEDs yet)
                #W1fluxlim = 7.690e-15 ## erg/s/cm^2
                #W2fluxlim = 1.285e-14 ## erg/s/cm^2
                ### monochromatic fluxes F_lambda 
                W1fluxlim = 1.1606e-15 ## W/m/cm^2
                W2fluxlim = 1.23287e-15 ## W/m/cm^2
            w1lum_lim=W1fluxlim * 4*np.pi*dL*dL ##W/m (sunrise units)
            w2lum_lim=W2fluxlim * 4*np.pi*dL*dL 
            wflim_str = wise_fluxlim_type+'_'
            
        elif not wise_fluxlim_type and np.max(wise_fluxlims)>0.0:
            print 'Imposing cuts based on input minimum WISE band fluxes:'
            print wise_fluxlims
            w1lum_lim = wise_fluxlims[0] * 4*np.pi*dL*dL ## W/m (sunrise units)
            w2lum_lim = wise_fluxlims[1] * 4*np.pi*dL*dL
            wflim_str = 'man_wflims_'

        else:
            print "No cuts imposed on minimum WISE band fluxes."
            w1lum_lim=0.0
            w2lum_lim=0.0

        print "z = %g"%z
        print "dL = %g cm"%dL
        print "w1lum_lim = %g"%w1lum_lim
        print "w2lum_lim = %g"%w2lum_lim
        #return
    else:
        if wise_fluxlim_type or np.max(wise_fluxlims)>0.0:
            print "Error: cannot set WISE fluxlims for rest-frame (z=0) calculations."
            return
        w1lum_lim=0.0
        w2lum_lim=0.0


    ncam=len(cam)
    print "plot_3dsep?", plot_3dsep

    name = copy(subdir_arr)
    if z==0:
        subdir_arr = define_sunrise_path_arr(basepath=maindir,name=name,return_full_path=False,
                                             sfrhist_w_grid=sfrhist_w_grid)
        zstr=''
    else:
        assert allsim_only, "Error: only 'allsim' plots have nonzero redshift functionality thus far."
        subdir_arr = define_sunrise_path_arr_z(basepath=maindir,name=name,return_full_path=False,
                                               z=z,agnx0=agnx0)
        zstr="_z%.1f"%z
    path_arr = ['%s/%s'%(maindir,s) for s in subdir_arr]

    nsim = len(path_arr)
    if (isinstance(name,list) or isinstance(name,tuple) or isinstance(name,np.ndarray)):
        extra = '_%dsims'%nsim
    elif nsim==1:
        extra = ''
    else: extra = '_'+name

    rstr='_grid_res_%.3d'%(1000*grid_nh_res) if grid_nh_res>0 else ''

    for s in subdir_arr: 
        if agnx0: 
            assert 'agnx0' in s,"Error: keyword 'agnx0' requires subdirs for agnx0 runs."
        else:
            assert 'agnx0' not in s,"Error: found agnx0 subdir but 'agnx0' keyword not set."

    assert valtype_allsim in ('median','mean'), 'Please enter "median" or "mean" for valtype_allsim.'

    plottitle=" "
    outdir = maindir if nsim>1 else path_arr[0]

    if equal_bins:
        ## dsep now set as keyword (default 10 kpc)
        sep_bins = np.arange(-dsep,200,dsep)
        dlgsep = 0.5    
        lgsep_bins = np.arange(-4,3,dlgsep)
    else:
        print 'Warning: ignoring any set value of keyword dsep for unequal bins.'
        sep_bins = np.array([-3.0, 0, 3, 10, 30, 100, 300])
        ##sep_bins = np.array([-3.0, 0, 1, 3, 10, 30, 100, 300])
        dsep = sep_bins[1:]-sep_bins[:-1]
        lgsep_bins = np.array([-4.0, -2, -1, 0, 0.5, 1, 1.5, 2, 2.5]) 
        dlgsep = lgsep_bins[1:]-lgsep_bins[:-1]

    dtbin=0.04
    #dtbin=0.05
    tpostmrg_bins = np.arange(0,0.25,dtbin)
    ### this puts everything in one bin (for testing)
    #dtbin=0.25
    #tpostmrg_bins = np.arange(0,0.5,dtbin)

    ### set limits for agn masks ###
    fedd_lims = [0.01,0.05,0.1]
    lagn_lims = [1.0e43, 1.0e44, 1.0e45]
    flbol_lims = [0.1, 0.3, 0.5]
    nagnmask = len(fedd_lims)+len(lagn_lims)+(len(flbol_lims))
    nwagnmask = nagnmask * 3
          
    w12_lims_hires = [0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
    w23_lims_hires = [2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5]
    lagn_lims_hires = [10**x for x in [42.5,43.0,43.5,44.0,44.5,45.0,45.5,46.0]]

          
    sepstr = '_bhsep' if plot_3dsep else '_bhprojsep'
    lgsepstr = '_lgbhsep' if plot_3dsep else '_lgbhprojsep'
    sepstring = lgsepstr if use_logbins else sepstr
    binstr = '_eqbins%d'%dsep if equal_bins else ''

    
    ### initialize histogram arrays ###
    sfr_max = np.zeros((nsim))
    ssfr_max = np.zeros((nsim))
    fedd_max = np.zeros((nsim))
    lagn_max = np.zeros((nsim))
    ltot_max = np.zeros((nsim))
    flbol_max = np.zeros((nsim))
    nh1_max = np.zeros((nsim))
    nh2_max = np.zeros((nsim))
    nh1_max_med = np.zeros((nsim))
    nh2_max_med = np.zeros((nsim))
    mstar_init = np.zeros((nsim))
    mstar_final = np.zeros((nsim))
    tlate0 = np.zeros((nsim))
    tmrg = np.zeros((nsim))
    mbh_init = np.zeros((nsim))
    mbh_tlate = np.zeros((nsim))
    mbh_tmrg = np.zeros((nsim))
    mbh_final = np.zeros((nsim))
    t_late_arr = np.zeros((nsim,ncam))
    t_post_arr = np.zeros((nsim,ncam))

    hist3d_dt = np.zeros((nsim,len(sep_bins)-1))
    histproj_dt = np.zeros((nsim,ncam,(len(sep_bins)-1)))

    agn_h0 = np.zeros((nsim,nagnmask,len(sep_bins)-1))
    agn_hist3d_dt = copy(agn_h0)
    noagn_hist3d_dt = copy(agn_h0)
    agn_hist3d_dt_frac = copy(agn_h0)
    noagn_hist3d_dt_frac = copy(agn_h0)
    dualagn_hist3d_dt = copy(agn_h0)
    dualagn_hist3d_dt_frac = copy(agn_h0)
    oneagn_hist3d_dt = copy(agn_h0)
    oneagn_hist3d_dt_frac = copy(agn_h0)
    agn_avghistproj_dt = copy(agn_h0)
    agn_avghistproj_dt_frac = copy(agn_h0)
    dualagn_avghistproj_dt = copy(agn_h0)
    dualagn_avghistproj_dt_frac = copy(agn_h0)
    oneagn_avghistproj_dt = copy(agn_h0)
    oneagn_avghistproj_dt_frac = copy(agn_h0)
    oneortwoagn_avghistproj_dt = copy(agn_h0)
    oneortwoagn_avghistproj_dt_frac = copy(agn_h0)

    agn_hp0 = np.zeros((nsim,nagnmask,ncam,(len(sep_bins)-1)))
    agn_histproj_dt = copy(agn_hp0)
    noagn_histproj_dt = copy(agn_hp0)
    agn_histproj_dt_frac = copy(agn_hp0)
    noagn_histproj_dt_frac = copy(agn_hp0)
    dualagn_histproj_dt = copy(agn_hp0)
    dualagn_histproj_dt_frac = copy(agn_hp0)
    oneagn_histproj_dt = copy(agn_hp0)
    oneagn_histproj_dt_frac = copy(agn_hp0)

    Ngas_bh1_agn_histproj = copy(agn_hp0)
    Ngas_bh2_agn_histproj = copy(agn_hp0)
    Ngas_bh1_dualagn_histproj = copy(agn_hp0)
    Ngas_bh2_dualagn_histproj = copy(agn_hp0)
    Ngas_agn_hist_postmrg = np.zeros((nsim,nagnmask,ncam,(len(tpostmrg_bins)-1)))

    Ngas_bh1_wiseagn_histproj = np.zeros((nsim,nwagnmask,ncam,len(sep_bins)-1))
    Ngas_bh2_wiseagn_histproj = np.zeros((nsim,nwagnmask,ncam,len(sep_bins)-1))
    Ngas_bh1_wisedualagn_histproj = np.zeros((nsim,nwagnmask,ncam,len(sep_bins)-1))
    Ngas_bh2_wisedualagn_histproj = np.zeros((nsim,nwagnmask,ncam,len(sep_bins)-1))
    Ngas_wiseagn_hist_postmrg = np.zeros((nsim,nwagnmask,ncam,(len(tpostmrg_bins)-1)))

    Ngas_bh1_wiseagn_mycut_histproj = np.zeros((nsim,2,ncam,len(sep_bins)-1))
    Ngas_bh2_wiseagn_mycut_histproj = np.zeros((nsim,2,ncam,len(sep_bins)-1))
    Ngas_bh1_wisedualagn_mycut_histproj = np.zeros((nsim,2,ncam,len(sep_bins)-1))
    Ngas_bh2_wisedualagn_mycut_histproj = np.zeros((nsim,2,ncam,len(sep_bins)-1))
    Ngas_wiseagn_mycut_hist_postmrg = np.zeros((nsim,2,ncam,(len(tpostmrg_bins)-1)))

    Ngas_bh1_histproj = np.zeros((nsim,ncam,(len(sep_bins)-1)))
    Ngas_bh2_histproj = np.zeros((nsim,ncam,(len(sep_bins)-1)))
    Ngas_hist_postmrg = np.zeros((nsim,ncam,(len(tpostmrg_bins)-1)))

    agn_global0 = np.zeros((nsim,nagnmask))
    fdualagn_tot = copy(agn_global0)
    fdualagn_early = copy(agn_global0)
    fdualagn_late = copy(agn_global0)
    f1or2agn_tot = copy(agn_global0)
    f1or2agn_early = copy(agn_global0)
    f1or2agn_late = copy(agn_global0)
    f1agn_post = copy(agn_global0)
    ftotagn_tot = copy(agn_global0)
    ftotagn_early = copy(agn_global0)
    ftotagn_late = copy(agn_global0)
    t_totagn_obsc = copy(agn_global0)
    t_totagn_obsc_late = copy(agn_global0)
    t_1agn_obsc_post = copy(agn_global0)
    dM_agn_obsc = copy(agn_global0)
    dM_agn_obsc_adv = copy(agn_global0)

    wagn_hp0 = np.zeros((nsim,nwagnmask,ncam,(len(sep_bins)-1)))
    wise_agn_histproj_dt = copy(wagn_hp0)
    wise_agn_histproj_dt_frac = copy(wagn_hp0)
    wise_nototagn_histproj_dt = copy(wagn_hp0)
    wise_nototagn_histproj_dt_frac = copy(wagn_hp0)
    agn_nowise_histproj_dt = copy(wagn_hp0)
    agn_nowise_histproj_dt_frac = copy(wagn_hp0)
    wise_dualagn_histproj_dt = copy(wagn_hp0)
    wise_dualagn_histproj_dt_frac = copy(wagn_hp0)
    dualagn_nowise_histproj_dt = copy(wagn_hp0)
    dualagn_nowise_histproj_dt_frac = copy(wagn_hp0)
    wise_oneagn_histproj_dt = copy(wagn_hp0)
    wise_oneagn_histproj_dt_frac = copy(wagn_hp0)
    wise_oneortwoagn_histproj_dt = copy(wagn_hp0)
    wise_oneortwoagn_histproj_dt_frac = copy(wagn_hp0)
    wise_oneortwoagn_histproj_dtagn_frac = copy(wagn_hp0)
    oneortwoagn_nowise_histproj_dt = copy(wagn_hp0)
    oneortwoagn_nowise_histproj_dt_frac = copy(wagn_hp0)
    oneortwoagn_nowise_histproj_dtagn_frac = copy(wagn_hp0)

    wagn_h0 = np.zeros((nsim,nwagnmask,(len(sep_bins)-1)))
    wise_agn_avghistproj_dt = copy(wagn_h0)
    wise_agn_avghistproj_dt_frac = copy(wagn_h0)
    wise_nototagn_avghistproj_dt = copy(wagn_h0)
    wise_nototagn_avghistproj_dt_frac = copy(wagn_h0)
    agn_nowise_avghistproj_dt = copy(wagn_h0)
    agn_nowise_avghistproj_dt_frac = copy(wagn_h0)

    wise_dualagn_avghistproj_dt = copy(wagn_h0)
    wise_dualagn_avghistproj_dt_frac = copy(wagn_h0)
    dualagn_nowise_avghistproj_dt = copy(wagn_h0)
    dualagn_nowise_avghistproj_dt_frac = copy(wagn_h0)
    wise_oneagn_avghistproj_dt = copy(wagn_h0)
    wise_oneagn_avghistproj_dt_frac = copy(wagn_h0)
    wise_oneortwoagn_avghistproj_dt = copy(wagn_h0)
    wise_oneortwoagn_avghistproj_dt_frac = copy(wagn_h0)
    wise_oneortwoagn_avghistproj_dtagn_frac = copy(wagn_h0)
    oneortwoagn_nowise_avghistproj_dt = copy(wagn_h0)
    oneortwoagn_nowise_avghistproj_dt_frac = copy(wagn_h0)
    oneortwoagn_nowise_avghistproj_dtagn_frac = copy(wagn_h0)

    wagn_global0 = np.zeros((nsim,nwagnmask))
    ftotagn_wise_tot = copy(wagn_global0)
    ftotagn_wise_early = copy(wagn_global0)
    ftotagn_wise_late = copy(wagn_global0)
    fnototagn_wise_tot = copy(wagn_global0)
    fnototagn_wise_early = copy(wagn_global0)
    fnototagn_wise_late = copy(wagn_global0)
    fnototagn_wise_post = copy(wagn_global0)

    fdualagn_wise_tot = copy(wagn_global0)
    fdualagn_wise_early = copy(wagn_global0)
    fdualagn_wise_late = copy(wagn_global0)
    f1or2agn_wise_tot = copy(wagn_global0)
    f1or2agn_wise_early = copy(wagn_global0)
    f1or2agn_wise_late = copy(wagn_global0)
    f1agn_wise_tot = copy(wagn_global0)
    f1agn_wise_early = copy(wagn_global0)
    f1agn_wise_late = copy(wagn_global0)
    f1agn_wise_post = copy(wagn_global0)
    fnoagn_wise_tot = copy(wagn_global0)
    fnoagn_wise_early = copy(wagn_global0)
    fnoagn_wise_late = copy(wagn_global0)
    fnoagn_wise_post = copy(wagn_global0)

    ftotagn_nowise_tot = copy(wagn_global0)
    ftotagn_nowise_early = copy(wagn_global0)
    ftotagn_nowise_late = copy(wagn_global0)

    fdualagn_nowise_tot = copy(wagn_global0)
    fdualagn_nowise_early = copy(wagn_global0)
    fdualagn_nowise_late = copy(wagn_global0)
    f1or2agn_nowise_tot = copy(wagn_global0)
    f1or2agn_nowise_early = copy(wagn_global0)
    f1or2agn_nowise_late = copy(wagn_global0)
    f1agn_nowise_tot = copy(wagn_global0)
    f1agn_nowise_early = copy(wagn_global0)
    f1agn_nowise_late = copy(wagn_global0)
    f1agn_nowise_post = copy(wagn_global0)

    t_1or2agn_tot_hires = np.zeros((nsim,len(lagn_lims_hires)))
    t_dualagn_tot_hires = np.zeros((nsim,len(lagn_lims_hires)))
    t_1or2agn_late_post_hires = np.zeros((nsim,len(lagn_lims_hires),ncam))
    t_dualagn_late_hires = np.zeros((nsim,len(lagn_lims_hires),ncam))

    t_1or2agn_wise_hires_l44_tot = np.zeros((nsim,len(w12_lims_hires),ncam))
    t_1or2agn_wise_hires_l44_late_post = np.zeros((nsim,len(w12_lims_hires),ncam))
    t_dualagn_wise_hires_l44_tot = np.zeros((nsim,len(w12_lims_hires),ncam))
    t_dualagn_wise_hires_l44_late = np.zeros((nsim,len(w12_lims_hires),ncam))
    t_wise_hires_l44_tot = np.zeros((nsim,len(w12_lims_hires),ncam))
    t_wise_hires_l44_late_post = np.zeros((nsim,len(w12_lims_hires),ncam))

    t_1or2agn_wise_hires_l45_tot = np.zeros((nsim,len(w12_lims_hires),ncam))
    t_1or2agn_wise_hires_l45_late_post = np.zeros((nsim,len(w12_lims_hires),ncam))
    t_dualagn_wise_hires_l45_tot = np.zeros((nsim,len(w12_lims_hires),ncam))
    t_dualagn_wise_hires_l45_late = np.zeros((nsim,len(w12_lims_hires),ncam))
    t_wise_hires_l45_tot = np.zeros((nsim,len(w12_lims_hires),ncam))
    t_wise_hires_l45_late_post = np.zeros((nsim,len(w12_lims_hires),ncam))

    t_j11_tot = np.zeros((nsim,ncam))
    t_j11_late_post = np.zeros((nsim,ncam))
    t_1or2agn_j11_l44_tot = np.zeros((nsim,ncam))
    t_1or2agn_j11_l44_late_post = np.zeros((nsim,ncam))
    t_1or2agn_j11_l45_tot = np.zeros((nsim,ncam))
    t_1or2agn_j11_l45_late_post = np.zeros((nsim,ncam))
    t_dualagn_j11_l44_tot = np.zeros((nsim,ncam))
    t_dualagn_j11_l44_late = np.zeros((nsim,ncam))
    t_dualagn_j11_l45_tot = np.zeros((nsim,ncam))
    t_dualagn_j11_l45_late = np.zeros((nsim,ncam))

    t_mycut_tot = np.zeros((nsim,ncam))
    t_mycut_late_post = np.zeros((nsim,ncam))
    t_1or2agn_mycut_l44_tot = np.zeros((nsim,ncam))
    t_1or2agn_mycut_l44_late_post = np.zeros((nsim,ncam))
    t_1or2agn_mycut_l45_tot = np.zeros((nsim,ncam))
    t_1or2agn_mycut_l45_late_post = np.zeros((nsim,ncam))
    t_dualagn_mycut_l44_tot = np.zeros((nsim,ncam))
    t_dualagn_mycut_l44_late = np.zeros((nsim,ncam))
    t_dualagn_mycut_l45_tot = np.zeros((nsim,ncam))
    t_dualagn_mycut_l45_late = np.zeros((nsim,ncam))

    ttot_arr = np.zeros((nsim))
    t_1or2agn_tot_arr = np.zeros((nsim,nagnmask))
    t_1or2agn_late_post_arr = np.zeros((nsim,nagnmask,ncam))

    t_dualagn_tot_arr = np.zeros((nsim,nagnmask))
    t_dualagn_late_arr = np.zeros((nsim,nagnmask,ncam))

    t_1or2agn_wise_arr = np.zeros((nsim,nwagnmask,ncam))
    t_1or2agn_wise_late_post_arr = np.zeros((nsim,nwagnmask,ncam))


    fedd_plot_titles = [r'f$_{\rm Edd}$>%g'%(lim) for lim in fedd_lims]
    lagn_plot_titles = [r'L$_{\rm AGN}$>10$^{%g}$ [erg s$^{-1}$]'%(np.log10(lim))  for lim in lagn_lims]
    flbol_plot_titles = [r'L$_{\rm AGN}$/L$_{\rm tot}$>%g'%(lim) for lim in flbol_lims]
    agn_plot_titles = plot_titles = fedd_plot_titles + lagn_plot_titles + flbol_plot_titles
    fedd_labels = ['fEdd>%g'%(lim) for lim in fedd_lims]
    lagn_labels = ['Lagn>%g'%(np.log10(lim))  for lim in lagn_lims]
    flbol_labels = ['fLbol>%g'%(lim) for lim in flbol_lims]
    agn_labels = fedd_labels + lagn_labels + flbol_labels
    ixagn_l44 = np.where(np.array(agn_labels)=='Lagn>44')[0]
    ixagn_l45 = np.where(np.array(agn_labels)=='Lagn>45')[0]
    ixagn_fl01 = np.where(np.array(agn_labels)=='fLbol>0.1')[0]
    
    wisemask_labels = ('W1W2>0.5','W1W2>0.8','J11 wedge')
    fedd_wise_titles = ['%s, fEdd>%g'%(lbl,lim) for lbl in wisemask_labels for lim in fedd_lims]
    lagn_wise_titles = ['%s, Lagn>%g'%(lbl,lim) for lbl in wisemask_labels for lim in lagn_lims]
    flbol_wise_titles = ['%s, fLbol>%g'%(lbl,lim) for lbl in wisemask_labels for lim in flbol_lims]
    wise_agn_titles = fedd_wise_titles + lagn_wise_titles + flbol_wise_titles
    wise_selection_labels = ['MIR AGN', 'MIR false neg', 'MIR false pos','']
    wise_duals_labels = ['MIR dual AGN', 'dual AGN, not MIR', 'MIR AGN, not dual','']
    print "wise_agn_titles:",wise_agn_titles
    
    ix_fedd = np.array([k for k in range(nwagnmask) if 'fEdd' in wise_agn_titles[k]])
    ix_lagn = np.array([k for k in range(nwagnmask) if 'Lagn' in wise_agn_titles[k]])
    ix_flbol = np.array([k for k in range(nwagnmask) if 'fLbol' in wise_agn_titles[k]])
    ix_w05 = np.array([k for k in range(nwagnmask) if 'W1W2>0.5' in wise_agn_titles[k]])
    ix_w08 = np.array([k for k in range(nwagnmask) if 'W1W2>0.8' in wise_agn_titles[k]])
    ix_j11 = np.array([k for k in range(nwagnmask) if 'J11' in wise_agn_titles[k]])

    print Ngas_bh1_wiseagn_histproj.shape
    print Ngas_wiseagn_hist_postmrg.shape
    print nwagnmask
    print ix_w05.shape
    print ix_w05
    print ix_w08.shape
    print ix_lagn
    print np.intersect1d(ix_w05,ix_lagn).shape
    print np.intersect1d(ix_w05,ix_lagn)

    wfbase = ('w05','w08','j11')


    ### load & process agn metadata ###
    for i_sim,path in enumerate(path_arr):

        print "\nloading agn metadata for sim %d (%s)..."%(i_sim,path)
        d = agn_metadata(path, ncam=ncam, skip_snap0=skip_snap0,
                         tmax_postmrg = tmax_postmrg, agnx0=agnx0, grid_nh_res=grid_nh_res)
    
        print d.dt.shape, d.posbh1.shape, d.lbol_grid.shape, d.bhprojsep.shape, d.w2w3.shape

        has2bh = (d.nbh==2) 
        ttot = d.dt.sum()
        ttot_arr[i_sim] = ttot
        meanw12 = np.mean(d.w1w2,axis=0)
        meanw23 = np.mean(d.w2w3,axis=0)

        j11mask = j11_wedge_cut(d.w1w2,d.w2w3,
                                w1lum_lim=w1lum_lim,w2lum_lim=w2lum_lim,w1lum=d.w1lum,w2lum=d.w2lum)
        j11mask_mean = j11_wedge_cut(meanw12,meanw23,
                                     w1lum_lim=w1lum_lim,w2lum_lim=w2lum_lim,w1lum=d.w1lum,w2lum=d.w2lum)

        mycutmask = my_wedge_cut(d.w1w2,d.w2w3,alt=alt_wedge,
                                 w1lum_lim=w1lum_lim,w2lum_lim=w2lum_lim,w1lum=d.w1lum,w2lum=d.w2lum)
        mycutmask_mean = my_wedge_cut(meanw12,meanw23,alt=alt_wedge,
                                      w1lum_lim=w1lum_lim,w2lum_lim=w2lum_lim,w1lum=d.w1lum,w2lum=d.w2lum)

        w05mask = (d.w1w2>0.5)
        w05mask_mean = (meanw12>0.5)
        w08mask = (d.w1w2>0.8)
        w08mask_mean = (meanw12>0.8)

        lagn_ratio = np.zeros(d.lagn1.size)
        lagn_ratio[has2bh] = d.lagn1[has2bh]/d.lagn2[has2bh]
        lagn_ratio[lagn_ratio>1.0] = 1.0/lagn_ratio[lagn_ratio>1.0]
        lg_lagn_ratio = np.log10(lagn_ratio)
        
        mdot = d.mdotbh1 + d.mdotbh2

        sfr_max[i_sim] = np.max(d.sfr)
        ssfr_max[i_sim] = np.max(d.sfr/d.mstar)
        fedd_max[i_sim] = np.max(d.fedd)
        lagn_max[i_sim] = np.max(d.lagn)
        ltot_max[i_sim] = np.max(d.lbol_grid)
        flbol_max[i_sim] = np.max(d.fagn_lbol)
        nh1_max[i_sim] = np.max(d.Ngas_bh1_grid)
        nh2_max[i_sim] = np.max(d.Ngas_bh2_grid)
        nh1_max_med[i_sim] = np.max(np.nanmedian(d.Ngas_bh1_grid,axis=0))
        nh2_max_med[i_sim] = np.max(np.nanmedian(d.Ngas_bh2_grid,axis=0))
        mstar_init[i_sim] = d.mstar[0]
        mstar_final[i_sim] = d.mstar[-1]
        tmrg[i_sim] = d.tmrg
        mbh_init[i_sim] = d.mbh[0]
        mbh_final[i_sim] = d.mbh[-1]

        ### calc histograms for each simulation ###
        print "calculating histograms for sim %d..."%i_sim

        dt_tiled = d.dt if plot_3dsep else np.tile(d.dt,(ncam,1))
        lagn_tiled = d.lagn if plot_3dsep else np.tile(d.lagn,(ncam,1))
        Ngas_bh1 = np.nanmedian(d.Ngas_bh1_grid,axis=0) if plot_3dsep else d.Ngas_bh1_grid
        Ngas_bh2 = np.nanmedian(d.Ngas_bh2_grid,axis=0) if plot_3dsep else d.Ngas_bh2_grid
        Ngas_avg = np.nanmean(np.array([Ngas_bh1,Ngas_bh2]),axis=0)

        sepvar = np.log10(d.bhsep) if use_logbins else d.bhsep
        projsepvar = np.log10(d.bhprojsep) if use_logbins else d.bhprojsep

        if len(sep_bins[sep_bins<0])>1 or sep_bins.min()>d.bhprojsep.min():
            print "Error: invalid definition of sep_bins:",sep_bins
            sys.exit()                                           

        ttot_sep_bins = np.histogram(sepvar,weights=d.dt,bins=sep_bins)[0]
        ttot_projsep_bins = np.array([ np.histogram(projsepvar[i,:],weights=d.dt,
                                                    bins=sep_bins)[0] for i in range(ncam) ])
        hist3d_dt[i_sim,:] = ttot_sep_bins
        histproj_dt[i_sim,:,:] = ttot_projsep_bins

        assert_isclose(ttot,hist3d_dt[i_sim,:].sum(),'hist3d_dt sum != ttot.')
        assert np.allclose(histproj_dt[i_sim,:,:].sum(axis=1),ttot),'Error: histproj_dt sum != ttot'

        idx_sep_bins = np.digitize(d.bhsep,bins=sep_bins,right=True)
        #print "idx_sep_bins:",idx_sep_bins.shape
        #print sep_bins.shape, d.bhsep.shape                
        idx_projsep_bins = np.array([ np.digitize(d.bhprojsep[i,:],bins=sep_bins,
                                                  right=True) for i in range(ncam) ])
        #print "idx_projsep_bins:",idx_projsep_bins.shape
        #print sep_bins.shape, d.bhprojsep.shape
        #print sep_bins[idx_projsep_bins].shape

        nproj_tnonzero_eachbin = np.array([len(ttot_projsep_bins[ttot_projsep_bins[:,k]>0,k]) 
                                           for k in range(len(sep_bins)-1)])
        avghistproj_dt = np.nansum(ttot_projsep_bins,axis=0) / (1.0*nproj_tnonzero_eachbin)
        #for j in range(ttot_projsep_bins.shape[1]):
        #    print "\nnproj(t>0): %d"%nproj_tnonzero_eachbin[j]
        #    print ttot_projsep_bins[:,j]

        idx_lgsep_bins = np.digitize(np.log10(d.bhsep),bins=lgsep_bins,right=True)
        idx_lgsep_bins[idx_lgsep_bins==len(lgsep_bins)] = 1
        idx_lgprojsep_bins = np.array([ np.digitize(np.log10(d.bhprojsep[i,:]),bins=lgsep_bins,
                                                    right=True) for i in range(ncam) ])
        idx_lgprojsep_bins[idx_lgprojsep_bins==len(lgsep_bins)] = 1


        ttot_lgsep_bins = np.histogram(np.log10(d.bhsep),weights=d.dt,bins=lgsep_bins)[0]
        ttot_lgprojsep_bins = np.array([ np.histogram(np.log10(d.bhprojsep[i,:]),weights=d.dt,
                                                      bins=lgsep_bins)[0] for i in range(ncam) ])
        nproj_tnonzero_eachlgbin = np.array([len(ttot_lgprojsep_bins[ttot_lgprojsep_bins[:,k]>0,k]) 
                                             for k in range(len(lgsep_bins)-1)])



        ### Time-weighted column density ###
        Ngas_bh1_histproj[i_sim,:,:] = np.array([np.histogram(projsepvar[j,:],bins=sep_bins, 
                                                              weights=d.dt*Ngas_bh1[j,:])[0]
                                                 for j in range(ncam)])/ttot_projsep_bins
        Ngas_bh1_histproj[Ngas_bh1_histproj==0] = np.nan
        wts = np.zeros(Ngas_bh2.shape)
        wts[Ngas_bh2==Ngas_bh2] = dt_tiled[Ngas_bh2==Ngas_bh2]*Ngas_bh2[Ngas_bh2==Ngas_bh2]
        Ngas_bh2_histproj[i_sim,:,:] = np.array([np.histogram(projsepvar[j,:],bins=sep_bins,weights=wts[j,:])[0]
                                                 for j in range(ncam)])/ttot_projsep_bins
        Ngas_bh2_histproj[Ngas_bh2_histproj==0] = np.nan
        
        assert (d.time[~has2bh]-d.tmrg).all() >= 0.0, \
            "Found negative tpostmrg values! (tmrg=%g, min tpostmrg=%g)"%(d.tmrg,d.time[~has2bh].min())
        Ngas_hist_postmrg[i_sim,:,:] = np.array([np.histogram(d.time[~has2bh]-d.tmrg,bins=tpostmrg_bins, 
                                                              weights=d.dt[~has2bh]*Ngas_bh1[j,~has2bh])[0] 
                                                 for j in range(ncam)])/dtbin
        Ngas_hist_postmrg[Ngas_hist_postmrg==0.0] = np.nan


        #dt_frac_sep_bins = d.dt/ttot_sep_bins[idx_sep_bins-1]
        #dt_frac_projsep_bins = dt_tiled/np.array([ttot_projsep_bins[i,idx_projsep_bins[i,:]-1]
        #                                          for i in range(ncam)])
        #print dt_frac_projsep_bins.shape
        #dt_frac_lgsep_bins = np.zeros(d.dt.shape)
        #dt_frac_lgsep_bins[has2bh] = d.dt[has2bh]/ttot_lgsep_bins[idx_lgsep_bins[has2bh]-1]
        #dt_frac_lgprojsep_bins = np.zeros(dt_tiled.shape)
        #dt_frac_lgprojsep_bins[:,has2bh] = dt_tiled[:,has2bh]/np.array([ttot_lgprojsep_bins[i,idx_lgprojsep_bins[i,has2bh]-1] for i in range(ncam)])

        #fedd_maskarr = [(d.fedd>=lim) for lim in fedd_lims]
        fedd_maskarr,dual_fedd_maskarr = set_agn_maskarr_list(d.fedd,d.fedd1,d.fedd2,fedd_lims)
        lagn_maskarr,dual_lagn_maskarr = set_agn_maskarr_list(d.lagn,d.lagn1,d.lagn2,lagn_lims)
        flbol_maskarr,dual_flbol_maskarr = set_agn_maskarr_list(d.fagn_lbol,d.fagn1_lbol,
                                                                d.fagn2_lbol,flbol_lims)
        agn_maskarr_list = fedd_maskarr + lagn_maskarr + flbol_maskarr
        dualagn_maskarr_list = dual_fedd_maskarr + dual_lagn_maskarr + dual_flbol_maskarr

        lagn_maskarr_hires,dual_lagn_maskarr_hires = set_agn_maskarr_list(d.lagn,d.lagn1,d.lagn2,lagn_lims_hires)
        

        avg_wisemask_list = [w05mask_mean, w08mask_mean, j11mask_mean]
        wisemask_list = [w05mask, w08mask, j11mask]
        assert len(wisemask_list)==3, 'Error: wise agn arrays need wisemask_list of len 3, not %d'%len(wisemask_list)

        if w1lum_lim>0.0 or w2lum_lim>0.0:
            print "w1lum_lim=%g"%w1lum_lim
            print "max w1lum=%g"%d.w1lum.max()
            print "min w1lum=%g"%d.w1lum.min()
            #print "max w1lum[L1>44|L2>44]=%g"%d.w1lum[:,((d.lagn1>1.0e44)|(d.lagn2>1.0e44))].max()
            #print "min w1lum[L1>44|L2>44]=%g"%d.w1lum[:,((d.lagn1>1.0e44)|(d.lagn2>1.0e44))].min()                                                      

            print "\nw2lum_lim=%g"%w2lum_lim
            print "max w2lum=%g"%d.w2lum.max()
            print "min w2lum=%g"%d.w2lum.min()
            #print "max w2lum[L1>44|L2>44]=%g"%d.w2lum[:,((d.lagn1>1.0e44)|(d.lagn2>1.0e44))].max()
            #print "min w1lum[L1>44|L2>44]=%g"%d.w2lum[:,((d.lagn1>1.0e44)|(d.lagn2>1.0e44))].min()                                                      
            
            fedd_wmaskarr_list,dual_fedd_wmaskarr_list = set_agn_wmaskarr_list(d.fedd,d.fedd1,d.fedd2,fedd_lims,wisemask_list,
                                                                               w1lum=d.w1lum, w2lum=d.w2lum, w1lum_lim=w1lum_lim,
                                                                               w2lum_lim=w2lum_lim )
            flbol_wmaskarr_list,dual_flbol_wmaskarr_list = set_agn_wmaskarr_list(d.fagn_lbol,d.fagn1_lbol,d.fagn2_lbol,
                                                                                 flbol_lims, wisemask_list,w1lum=d.w1lum, w2lum=d.w2lum, 
                                                                                 w1lum_lim=w1lum_lim,w2lum_lim=w2lum_lim)
            lagn_wmaskarr_list,dual_lagn_wmaskarr_list = set_agn_wmaskarr_list(d.lagn,d.lagn1,d.lagn2,lagn_lims, wisemask_list,
                                                                               w1lum=d.w1lum, w2lum=d.w2lum, w1lum_lim=w1lum_lim,
                                                                               w2lum_lim=w2lum_lim)
        else:
            fedd_wmaskarr_list,dual_fedd_wmaskarr_list = set_agn_wmaskarr_list(d.fedd,d.fedd1,d.fedd2,fedd_lims,wisemask_list)
            flbol_wmaskarr_list,dual_flbol_wmaskarr_list = set_agn_wmaskarr_list(d.fagn_lbol,d.fagn1_lbol,d.fagn2_lbol,
                                                                                 flbol_lims,wisemask_list)
            lagn_wmaskarr_list,dual_lagn_wmaskarr_list = set_agn_wmaskarr_list(d.lagn,d.lagn1,d.lagn2,lagn_lims,wisemask_list)
        
        agn_wmaskarr_list = fedd_wmaskarr_list + lagn_wmaskarr_list + flbol_wmaskarr_list
        dualagn_wmaskarr_list = dual_fedd_wmaskarr_list + dual_lagn_wmaskarr_list + dual_flbol_wmaskarr_list

        wisemask_hires_list = [(d.w1w2>wlim) for wlim in w12_lims_hires]
        print len(wisemask_hires_list)
        print len(wisemask_hires_list[0])
        print len(wisemask_hires_list[0][0])

        
        #fedd_wmaskarr_hires,dual_fedd_wmaskarr_hires = set_agn_wmaskarr_list(d.fedd,d.fedd1,d.fedd2,fedd_lims,wisemask_hires_list)
        #flbol_wmaskarr_hires,dual_flbol_wmaskarr_hires = set_agn_wmaskarr_list(d.fagn_lbol,d.fagn1_lbol,d.fagn2_lbol,
        #                                                                     flbol_lims,wisemask_hires_list)
        #lagn_wmaskarr_hires,dual_lagn_wmaskarr_hires = set_agn_wmaskarr_list(d.lagn,d.lagn1,d.lagn2,lagn_lims,wisemask_hires_list)

        #agn_wmaskarr_hires = fedd_wmaskarr_hires + lagn_wmaskarr_hires + flbol_wmaskarr_hires
        #dualagn_wmaskarr_hires = dual_fedd_wmaskarr_hires + dual_lagn_wmaskarr_hires + dual_flbol_wmaskarr_hires


        if w1lum_lim>0.0 or w2lum_lim>0.0:
            print "w1lum_lim=%g"%w1lum_lim
            print "max w1lum=%g"%d.w1lum.max()
            print "min w1lum=%g"%d.w1lum.min()
            #print "max w1lum[L1>44|L2>44]=%g"%d.w1lum[:,((d.lagn1>1.0e44)|(d.lagn2>1.0e44))].max()
            #print "min w1lum[L1>44|L2>44]=%g"%d.w1lum[:,((d.lagn1>1.0e44)|(d.lagn2>1.0e44))].min()                                                      

            print "\nw2lum_lim=%g"%w2lum_lim
            print "max w2lum=%g"%d.w2lum.max()
            print "min w2lum=%g"%d.w2lum.min()
            #print "max w2lum[L1>44|L2>44]=%g"%d.w2lum[:,((d.lagn1>1.0e44)|(d.lagn2>1.0e44))].max()
            #print "min w1lum[L1>44|L2>44]=%g"%d.w2lum[:,((d.lagn1>1.0e44)|(d.lagn2>1.0e44))].min()                                                      
            l44_wmaskarr_hires,dual_l44_wmaskarr_hires = set_agn_wmaskarr_list(d.lagn,d.lagn1,d.lagn2,[1.0e44],wisemask_hires_list,
                                                                               w1lum=d.w1lum, w2lum=d.w2lum, w1lum_lim=w1lum_lim,
                                                                               w2lum_lim=w2lum_lim )
            fl01_wmaskarr_hires,dual_fl01_wmaskarr_hires = set_agn_wmaskarr_list(d.fagn_lbol,d.fagn1_lbol,d.fagn2_lbol,[0.1],
                                                                                 wisemask_hires_list, w1lum=d.w1lum, w2lum=d.w2lum, 
                                                                                 w1lum_lim=w1lum_lim, w2lum_lim=w2lum_lim )
            l45_wmaskarr_hires,dual_l45_wmaskarr_hires = set_agn_wmaskarr_list(d.lagn,d.lagn1,d.lagn2,[1.0e45],wisemask_hires_list,
                                                                               w1lum=d.w1lum, w2lum=d.w2lum, w1lum_lim=w1lum_lim,
                                                                               w2lum_lim=w2lum_lim )
        else:
            l44_wmaskarr_hires,dual_l44_wmaskarr_hires = set_agn_wmaskarr_list(d.lagn,d.lagn1,d.lagn2,[1.0e44],wisemask_hires_list)
            fl01_wmaskarr_hires,dual_fl01_wmaskarr_hires = set_agn_wmaskarr_list(d.fagn_lbol,d.fagn1_lbol,d.fagn2_lbol,
                                                                                 [0.1],wisemask_hires_list)

            l45_wmaskarr_hires,dual_l45_wmaskarr_hires = set_agn_wmaskarr_list(d.lagn,d.lagn1,d.lagn2,[1.0e45],wisemask_hires_list)
        #print len(l44_wmaskarr_hires)
        #print len(l44_wmaskarr_hires[0])
        #print len(l44_wmaskarr_hires[0][0])
        #print len(l44_wmaskarr_hires[0][0][0])
        #print len(dual_l44_wmaskarr_hires)
        #print len(dual_l44_wmaskarr_hires[0])
        #print len(dual_l44_wmaskarr_hires[0][0])
        #print len(dual_l44_wmaskarr_hires[0][0][0])

        #print d.fedd.shape, d.lagn1.shape,d.fagn1_lbol.shape
        #print len(fedd_maskarr),len(dual_fedd_maskarr)
        #print len(fedd_maskarr[0]),len(dual_fedd_maskarr[0])
        #print len(fedd_wmaskarr_list),len(dual_fedd_wmaskarr_list)
        #print len(fedd_wmaskarr_list[0]),len(dual_fedd_wmaskarr_list[0])
        #print len(fedd_wmaskarr_list[0][1]),len(dual_fedd_wmaskarr_list[0][1])
        #print len(lagn_wmaskarr_list[0]),len(dual_lagn_wmaskarr_list[0])
        #print (lagn_wmaskarr_list[0][0].shape),(dual_lagn_wmaskarr_list[0][0].shape)
        #print len(agn_wmaskarr_list)
        #print len(agn_wmaskarr_list[0])
        #print len(agn_wmaskarr_list[0][0])
        #print len(agn_wmaskarr_list[0][0][0])
        #print len(dualagn_wmaskarr_list)
        #print len(dualagn_wmaskarr_list[0])
        #print len(dualagn_wmaskarr_list[0][0])
        #print len(dualagn_wmaskarr_list[0][0][0])


        idx_sep = np.where(sep_bins[:-1]<latesep)[0][-1]

        t_early = np.array([d.dt[d.bhprojsep[i,:]>=latesep].sum() for i in range(ncam)])
        t_late = np.array([d.dt[(d.bhprojsep[i,:]<latesep)&
                                (d.bhprojsep[i,:]>0)].sum() for i in range(ncam)])
        t_post = np.array([d.dt[d.bhprojsep[i,:]<=0].sum() for i in range(ncam)])
        t_late_arr[i_sim] = t_late
        t_post_arr[i_sim] = t_post
        

        tlate0[i_sim] = np.median([d.time[d.bhprojsep[i,:]<latesep].min() for i in range(ncam)])
        print '\n tlate0: ',tlate0
        idx_tlate0 = np.where(np.abs(d.time-tlate0[i_sim])==np.min(np.abs(d.time-tlate0[i_sim])))[0]
        print '\n t_mrg: ',d.tmrg
        mbh_tlate[i_sim] = d.mbh[idx_tlate0]
        idx_tmrg = np.where(np.abs(d.time-d.tmrg)==np.min(np.abs(d.time-d.tmrg)))[0]
        mbh_tmrg[i_sim] = d.mbh[idx_tmrg]
        

        print "\nt_early:",t_early
        print np.nansum(ttot_projsep_bins[:,idx_sep:],axis=1)
        print "\nt_late:",t_late
        print np.nansum(ttot_projsep_bins[:,1:idx_sep+1],axis=1)
        print "\nt_post:",t_post
        print ttot_projsep_bins[:,0]
        print "\nttot:",ttot
        print "max proj. sep for each cam:",[d.bhprojsep[i,:].max() for i in range(ncam)]
        assert t_early.all()>0 and t_late.all()>0 and t_post.all()>0,'Error: at least one projection has a merger phase with t=0. Consider a different value of latesep (%g) or tmax_postmrg (%g).'%(latesep,tmax_postmrg)

        assert np.allclose(t_early+t_late+t_post,ttot),'Error: time in merger phases !=ttot'
        nproj_tnonzero_early = t_early[t_early>0].size
        nproj_tnonzero_late = t_late[t_late>0].size
        nproj_tnonzero_post = t_post[t_post>0].size
        print "nproj_tnonzero_[early,late,post]: %d %d %d"%(nproj_tnonzero_early,nproj_tnonzero_late,nproj_tnonzero_post)

        if not agnx0:
            print "time when fedd1>fedd: %g"%d.dt[d.fedd1>d.fedd].sum()
            print "time when fedd2>fedd: %g"%d.dt[d.fedd2>d.fedd].sum()
            print d.fedd1.min(),d.fedd1.max(),d.fedd1.mean()
            print d.fedd2[d.fedd2>0].min(),d.fedd2.max(),d.fedd2[d.fedd2>0].mean()
            print d.fedd.min(),d.fedd.max(),d.fedd.mean()
            assert d.dt[d.lagn1>d.lagn].size==0, "wtf!"
            assert d.dt[d.lagn2>d.lagn].size==0, "wtf!"
            assert d.dt[d.fagn1_lbol>d.fagn_lbol].size==0, "wtf!"
            assert d.dt[d.fagn2_lbol>d.fagn_lbol].size==0, "wtf!"


        
        for i_mask_hires in np.arange(len(lagn_maskarr_hires)):

            totmask_hires = lagn_maskarr_hires[i_mask_hires]
            dualmask_hires = dual_lagn_maskarr_hires[i_mask_hires][0]
            onemask_hires = dual_lagn_maskarr_hires[i_mask_hires][1]

            t_dualagn_tot_hires[i_sim,i_mask_hires] = d.dt[dualmask_hires].sum()
            t_oneagn_hires = d.dt[onemask_hires].sum()
            t_1or2agn_tot_hires[i_sim,i_mask_hires] = t_dualagn_tot_hires[i_sim,i_mask_hires] + t_oneagn_hires

            t_dualagn_late_hires[i_sim,i_mask_hires,:] = np.array([d.dt[(d.bhprojsep[i,:]<latesep)&
                                                                        (d.bhprojsep[i,:]>0)&(dualmask_hires)].sum() 
                                                                   for i in range(ncam)])
            t_oneagn_late_hires = np.array([d.dt[(d.bhprojsep[i,:]<latesep)&
                                           (d.bhprojsep[i,:]>0)&(onemask_hires)].sum() for i in range(ncam)])
            t_1or2agn_late_hires = t_oneagn_late_hires + t_dualagn_late_hires[i_sim,i_mask_hires,:]

            t_1agn_post_hires = np.array([d.dt[(d.bhprojsep[i,:]<=0)&(onemask_hires)].sum() for i in range(ncam)])

            t_1or2agn_late_post_hires[i_sim,i_mask_hires,:] = t_1or2agn_late_hires + t_1agn_post_hires
            #f1or2agn_tot_hires[i_sim,i_mask_hires] = t_1or2agn_tot_hires / ttot
            #f1or2agn_late_post_hires[i_sim,i_mask_hires] = ( np.sum( (t_1or2agn_late_hires+t_1agn_post_hires) / (t_late+t_post) ) / 
            #                                                 (1.0*np.maximum(nproj_tnonzero_late,nproj_tnonzero_post)) )


        #for i_mask,mask in enumerate(agn_with_hilolagn_maskarr_list):
        for i_mask in np.arange(len(agn_maskarr_list)):
            #if agnx0: continue

            print "\n\n**%s**:\n"%agn_labels[i_mask]

            totmask = agn_maskarr_list[i_mask]
            dualmask = dualagn_maskarr_list[i_mask][0]
            onemask = dualagn_maskarr_list[i_mask][1]

            ### Time when total AGN luminosity exceeds threshold ###
            t_totagn = d.dt[totmask].sum()
            agn_hist3d_dt[i_sim,i_mask,:] = np.histogram(sepvar[totmask],bins=sep_bins,
                                                         weights=d.dt[totmask])[0]
            t_totagn_obsc[i_sim,i_mask] = d.dt[(totmask)&((np.nanmedian(Ngas_bh1,axis=0)>1.0e23)|
                                                          (np.nanmedian(Ngas_bh2,axis=0)>1.0e23))].sum()

            t_nototagn = d.dt[~totmask].sum()
            agn_hist3d_dt_frac[i_sim,i_mask,:] = agn_hist3d_dt[i_sim,i_mask,:]/ttot_sep_bins

            agn_histproj_dt[i_sim,i_mask,:,:]=np.array([np.histogram(projsepvar[j,totmask],bins=sep_bins, 
                                                                     weights=d.dt[totmask])[0] 
                                                        for j in range(ncam)])
            nproj_tagnnonzero_eachbin=np.array([len(agn_histproj_dt[i_sim,i_mask,agn_histproj_dt[i_sim,i_mask,:,k]>0,k]) 
                                                 for k in range(len(sep_bins)-1)])
            agn_avghistproj_dt[i_sim,i_mask,:] = np.nansum(agn_histproj_dt[i_sim,i_mask,:,:],axis=0) / (1.0*ncam)
            agn_histproj_dt_frac[i_sim,i_mask,:,:] = agn_histproj_dt[i_sim,i_mask,:,:]/ttot_projsep_bins
            agn_avghistproj_dt_frac[i_sim,i_mask,:] = np.nansum(agn_histproj_dt_frac[i_sim,i_mask,:,:],axis=0)/(1.0*nproj_tnonzero_eachbin)
            assert_isclose(agn_hist3d_dt[i_sim,i_mask,:].sum(),np.nansum(agn_avghistproj_dt[i_sim,i_mask,:]),'ttot')

            ftotagn_tot[i_sim,i_mask] = t_totagn / ttot
            print "\ni_mask=%d"%i_mask
            print "t_totagn:",t_totagn
            print "agn_hist3d_dt:",agn_hist3d_dt[i_sim,i_mask].sum()
            print "agn_histproj_dt:",agn_histproj_dt[i_sim,i_mask,:,:].sum(axis=1)
            print "agn_avghistproj_dt:",agn_avghistproj_dt[i_sim,i_mask,:].sum()

            ### Time-weighted AGN column density ###
            #Ngas_bh1_agn_hist3d[i_sim,i_mask,:,:] = np.array([np.histogram(sepvar[totmask],bins=sep_bins, 
            #                                                               weights=d.dt[totmask]*Ngas_bh1[j,totmask])[0] 
            #                                                for j in range(ncam)])/ttot_sep_bins
            Ngas_bh1_agn_histproj[i_sim,i_mask,:,:] = np.array([np.histogram(projsepvar[j,totmask],bins=sep_bins, 
                                                                             weights=d.dt[totmask]*Ngas_bh1[j,totmask])[0]
                                                                for j in range(ncam)])/ttot_projsep_bins
            Ngas_bh1_agn_histproj[Ngas_bh1_agn_histproj==0] = np.nan
            wts = np.zeros(Ngas_bh2.shape)
            wts[Ngas_bh2==Ngas_bh2] = dt_tiled[Ngas_bh2==Ngas_bh2]*Ngas_bh2[Ngas_bh2==Ngas_bh2]
            Ngas_bh2_agn_histproj[i_sim,i_mask,:,:] = np.array([np.histogram(projsepvar[j,totmask],bins=sep_bins, 
                                                                             weights=wts[j,totmask])[0]
                                                                for j in range(ncam)])/ttot_projsep_bins
            Ngas_bh2_agn_histproj[Ngas_bh2_agn_histproj==0] = np.nan
            
            Ngas_agn_hist_postmrg[i_sim,i_mask,:,:] = np.array([np.histogram(d.time[~has2bh][totmask[~has2bh]]-d.tmrg,bins=tpostmrg_bins, weights=d.dt[~has2bh][totmask[~has2bh]]*Ngas_bh1[j,~has2bh][totmask[~has2bh]])[0] for j in range(ncam)])/dtbin
            Ngas_agn_hist_postmrg[Ngas_agn_hist_postmrg==0.0] = np.nan


            ### Time-weighted dual AGN column density ###
            Ngas_bh1_dualagn_histproj[i_sim,i_mask,:,:] = np.array([np.histogram(projsepvar[j,dualmask],bins=sep_bins, 
                                                                             weights=d.dt[dualmask]*Ngas_bh1[j,dualmask])[0]
                                                                for j in range(ncam)])/ttot_projsep_bins
            Ngas_bh1_dualagn_histproj[Ngas_bh1_dualagn_histproj==0] = np.nan
            wts = np.zeros(Ngas_bh2.shape)
            wts[Ngas_bh2==Ngas_bh2] = dt_tiled[Ngas_bh2==Ngas_bh2]*Ngas_bh2[Ngas_bh2==Ngas_bh2]
            Ngas_bh2_dualagn_histproj[i_sim,i_mask,:,:] = np.array([np.histogram(projsepvar[j,dualmask],bins=sep_bins, 
                                                                             weights=wts[j,dualmask])[0]
                                                                for j in range(ncam)])/ttot_projsep_bins
            Ngas_bh2_dualagn_histproj[Ngas_bh2_dualagn_histproj==0] = np.nan

            ### Time with two AGN simultaneously active ###
            t_dualagn = d.dt[dualmask].sum()
            dualagn_hist3d_dt[i_sim,i_mask,:] = np.histogram(sepvar[dualmask],bins=sep_bins,
                                                             weights=d.dt[dualmask])[0]
            dualagn_hist3d_dt_frac[i_sim,i_mask,:] = dualagn_hist3d_dt[i_sim,i_mask,:]/ttot_sep_bins
            dualagn_histproj_dt[i_sim,i_mask,:,:]=np.array([np.histogram(projsepvar[j,dualmask],bins=sep_bins, 
                                                                         weights=d.dt[dualmask])[0] 
                                                            for j in range(ncam)])
            nproj_t2agnnonzero_eachbin=np.array([len(dualagn_histproj_dt[i_sim,i_mask,dualagn_histproj_dt[i_sim,i_mask,:,k]>0,k]) 
                                                 for k in range(len(sep_bins)-1)])
            dualagn_avghistproj_dt[i_sim,i_mask,:] = np.nansum(dualagn_histproj_dt[i_sim,i_mask,:,:],axis=0) / (1.0*ncam)
            dualagn_histproj_dt_frac[i_sim,i_mask,:,:] = dualagn_histproj_dt[i_sim,i_mask,:,:]/ttot_projsep_bins
            dualagn_avghistproj_dt_frac[i_sim,i_mask,:] = np.nansum(dualagn_histproj_dt_frac[i_sim,i_mask,:,:],axis=0)/(1.0*nproj_tnonzero_eachbin)

            ### Time with exactly one AGN active ###
            t_oneagn = d.dt[onemask].sum()
            oneagn_hist3d_dt[i_sim,i_mask,:] = np.histogram(sepvar[onemask],bins=sep_bins,
                                                             weights=d.dt[onemask])[0]
            oneagn_hist3d_dt_frac[i_sim,i_mask,:] = oneagn_hist3d_dt[i_sim,i_mask,:]/ttot_sep_bins
            oneagn_histproj_dt[i_sim,i_mask,:,:]=np.array([np.histogram(projsepvar[j,onemask],bins=sep_bins, 
                                                                        weights=d.dt[onemask])[0] 
                                                           for j in range(ncam)])
            nproj_t1agnnonzero_eachbin=np.array([len(oneagn_histproj_dt[i_sim,i_mask,oneagn_histproj_dt[i_sim,i_mask,:,k]>0,k]) 
                                                 for k in range(len(sep_bins)-1)])
            #print np.nansum(oneagn_histproj_dt[i_sim,i_mask,:,:],axis=0)/(1.0*ncam)
            #print oneagn_hist3d_dt[i_sim,i_mask,:]
            oneagn_avghistproj_dt = np.nansum(oneagn_histproj_dt[i_sim,i_mask,:,:],axis=0) / (1.0*ncam)
            oneagn_histproj_dt_frac[i_sim,i_mask,:,:]=oneagn_histproj_dt[i_sim,i_mask,:,:]/ttot_projsep_bins
            oneagn_avghistproj_dt_frac = np.nansum(oneagn_histproj_dt_frac[i_sim,i_mask,:,:],axis=0)/(1.0*nproj_tnonzero_eachbin)

            ### Time with no AGN active ###
            t_noagn = d.dt[((~dualmask)&(~onemask))].sum()
            noagn_hist3d_dt[i_sim,i_mask,:] = np.histogram(sepvar[((~onemask)&(~dualmask))],bins=sep_bins,
                                                             weights=d.dt[((~onemask)&(~dualmask))])[0]
            noagn_hist3d_dt_frac[i_sim,i_mask,:] = noagn_hist3d_dt[i_sim,i_mask,:]/ttot_sep_bins
            noagn_histproj_dt[i_sim,i_mask,:,:]=np.array([np.histogram(projsepvar[j,((~onemask)&(~dualmask))],bins=sep_bins, 
                                                                       weights=d.dt[((~onemask)&(~dualmask))])[0] 
                                                          for j in range(ncam)])
            nproj_tnoagnnonzero_eachbin=np.array([len(noagn_histproj_dt[i_sim,i_mask,noagn_histproj_dt[i_sim,i_mask,:,k]>0,k]) 
                                                 for k in range(len(sep_bins)-1)])

            noagn_avghistproj_dt = np.nansum(noagn_histproj_dt[i_sim,i_mask,:,:],axis=0) / (1.0*ncam)
            noagn_histproj_dt_frac[i_sim,i_mask,:,:]=noagn_histproj_dt[i_sim,i_mask,:,:]/ttot_projsep_bins
            noagn_avghistproj_dt_frac = np.nansum(noagn_histproj_dt[i_sim,i_mask,:,:]/
                                                  (1.0*ttot_projsep_bins),axis=0) / (1.0*nproj_tnonzero_eachbin)

            print "t_oneagn:",t_oneagn
            print "t_dualagn:",t_dualagn
            print "t_totagn:",t_totagn

            t_1or2agn_tot = t_dualagn + t_oneagn
            t_1or2agn_tot_arr[i_sim,i_mask] = t_1or2agn_tot
            t_dualagn_tot_arr[i_sim,i_mask] = t_dualagn
            ### changing fdual definitions to be relative to ttot for more consistency; we can then do f2agn/f1agn as desired.
            #fdualagn_tot[i_sim,i_mask] = t_dualagn / t_1or2agn_tot
            fdualagn_tot[i_sim,i_mask] = t_dualagn / ttot
            f1or2agn_tot[i_sim,i_mask] = t_1or2agn_tot / ttot

            oneortwoagn_histproj_dt = oneagn_histproj_dt + dualagn_histproj_dt
            oneortwoagn_histproj_dt_frac = oneortwoagn_histproj_dt / ttot_projsep_bins
            oneortwoagn_avghistproj_dt[i_sim,i_mask,:] = np.nansum(oneortwoagn_histproj_dt[i_sim,i_mask,:,:],axis=0)/(1.0*ncam)
            oneortwoagn_avghistproj_dt_frac[i_sim,i_mask,:] = np.nansum(oneortwoagn_histproj_dt_frac[i_sim,i_mask,:,:],axis=0)/(1.0*nproj_tnonzero_eachbin)
            nproj_t1or2agnnonzero_eachbin=np.array([len(oneortwoagn_histproj_dt[i_sim,i_mask,oneortwoagn_histproj_dt[i_sim,i_mask,:,k]>0,k])
                                                    for k in range(len(sep_bins)-1)])

            t_totagn_early = np.array([d.dt[(d.bhprojsep[i,:]>=latesep)&(totmask)].sum() for i in range(ncam)])            
            nproj_tagnnonzero_early = len(t_totagn_early[t_totagn_early>0])
            ftotagn_early[i_sim,i_mask] = np.nansum( t_totagn_early/t_early ) / (1.0*nproj_tnonzero_early)

            t_dualagn_early = np.array([d.dt[(d.bhprojsep[i,:]>=latesep)&(dualmask)].sum() for i in range(ncam)])
            fdualagn_early[i_sim,i_mask] = np.sum( t_dualagn_early / t_early ) / (1.0*nproj_tnonzero_early)

            t_oneagn_early = np.array([d.dt[(d.bhprojsep[i,:]>=latesep)&(onemask)].sum() for i in range(ncam)])
            fone_early = np.sum( t_oneagn_early / t_early ) / (1.0*nproj_tnonzero_early)
            t_noagn_early = np.array([d.dt[(d.bhprojsep[i,:]>=latesep)&(~dualmask)&(~onemask)].sum() for i in range(ncam)])
            fno_early = np.sum( t_noagn_early / t_early ) / (1.0*nproj_tnonzero_early)

            t_1or2agn_early = t_oneagn_early + t_dualagn_early
            nproj_t1or2agnnonzero_early = len(t_1or2agn_early[t_1or2agn_early>0])
            f1or2agn_early[i_sim,i_mask] = np.sum( t_1or2agn_early / t_early ) / (1.0*nproj_tnonzero_early)
            
            tmpfearly=fdualagn_early[i_sim,i_mask]+fone_early
            assert np.isclose(tmpfearly[tmpfearly==tmpfearly],f1or2agn_early[i_sim,i_mask]),'fdual_early+fone_early!=f1or2agn_early'
            tmpfearly=fdualagn_early[i_sim,i_mask]+fone_early+fno_early
            print "t_early=",t_early
            print "t_1or2agn_early=",t_1or2agn_early
            print "t_dualagn_early=",t_dualagn_early
            print "t_oneagn_early=",t_oneagn_early
            print "t_noagn_early=",t_noagn_early
            assert np.isclose(tmpfearly[tmpfearly==tmpfearly],1),'fdual_early+fone_early+fno_early!=1'


            t_totagn_late = np.array([d.dt[(d.bhprojsep[i,:]<latesep)&
                                           (d.bhprojsep[i,:]>0)&(totmask)].sum() for i in range(ncam)])
            nproj_tagnnonzero_late = len(t_totagn_late[t_totagn_late>0])
            ftotagn_late[i_sim,i_mask] = np.nansum( t_totagn_late/t_late ) / (1.0*nproj_tnonzero_late)
            t_totagn_obsc_late[i_sim,i_mask] = np.sum(np.array([d.dt[(totmask)&((Ngas_bh1[i,:]>1.0e23)|(Ngas_bh2[i,:]>1.0e23))&
                                                                      (d.bhprojsep[i,:]<latesep)&(d.bhprojsep[i,:]>0)].sum()
                                                                for i in range(ncam)]))/(1.0*nproj_tnonzero_late)

            t_dualagn_late = np.array([d.dt[(d.bhprojsep[i,:]<latesep)&
                                            (d.bhprojsep[i,:]>0)&(dualmask)].sum() for i in range(ncam)])
            t_oneagn_late = np.array([d.dt[(d.bhprojsep[i,:]<latesep)&
                                           (d.bhprojsep[i,:]>0)&(onemask)].sum() for i in range(ncam)])
            t_noagn_late = np.array([d.dt[(d.bhprojsep[i,:]<latesep)&(d.bhprojsep[i,:]>0)&
                                          (~dualmask)&(~onemask)].sum() for i in range(ncam)])
            t_1or2agn_late = t_oneagn_late + t_dualagn_late
            nproj_t1or2agnnonzero_late = len(t_1or2agn_late[t_1or2agn_late>0])
            f1or2agn_late[i_sim,i_mask] = np.sum( t_1or2agn_late / t_late ) / (1.0*nproj_tnonzero_late)
            fdualagn_late[i_sim,i_mask] = np.sum( t_dualagn_late / t_late ) / (1.0*nproj_tnonzero_late)
            fone_late = np.sum( t_oneagn_late / t_late ) / (1.0*nproj_tnonzero_late)
            fno_late = np.sum( t_noagn_late / t_late ) / (1.0*nproj_tnonzero_late)
            tmpflate=fdualagn_late[i_sim,i_mask]+fone_late
            assert np.isclose(tmpflate[tmpflate==tmpflate],f1or2agn_late[i_sim,i_mask]),\
                'fdual_late+fone_late!=f1or2agn_late'
            tmpflate=fdualagn_late[i_sim,i_mask]+fone_late+fno_late
            print "t_late=",t_late
            print "t_1or2agn_late=",t_1or2agn_late
            print "t_dualagn_late=",t_dualagn_late
            print "t_oneagn_late=",t_oneagn_late
            print "t_noagn_late=",t_noagn_late
            assert np.isclose(tmpflate[tmpflate==tmpflate],1),'fdual_late+fone_late+fno_late!=1'

            t_1agn_post = np.array([d.dt[(d.bhprojsep[i,:]<=0)&(onemask)].sum() for i in range(ncam)])
            print "t_1agn_post=",t_1agn_post
            print "should match this:", np.array([d.dt[(d.bhprojsep[i,:]<=0)&(totmask)].sum() for i in range(ncam)])
            t_1agn_obsc_post[i_sim,i_mask] = np.mean(np.array([d.dt[(onemask)&(Ngas_bh1[i,:]>1.0e23)&(d.bhprojsep[i,:]<=0)].sum()
                                                               for i in range(ncam)]))

            dM_agn_obsc[i_sim,i_mask] = np.nansum((mdot[(totmask)&((np.nanmedian(Ngas_bh1,axis=0)>1.0e23)|
                                                             (np.nanmedian(Ngas_bh2,axis=0)>1.0e23))] *
                                                       d.dt[(totmask)&((np.nanmedian(Ngas_bh1,axis=0)>1.0e23)|
                                                                       (np.nanmedian(Ngas_bh2,axis=0)>1.0e23))]*1.0e9))

            dM_agn_obsc_adv[i_sim,i_mask] = np.nansum((mdot[(totmask)&((np.nanmedian(Ngas_bh1,axis=0)>1.0e23)|
                                                             (np.nanmedian(Ngas_bh2,axis=0)>1.0e23))&
                                                            (d.time>=tlate0[i_sim])] * 
                                                       d.dt[(totmask)&((np.nanmedian(Ngas_bh1,axis=0)>1.0e23)|
                                                                       (np.nanmedian(Ngas_bh2,axis=0)>1.0e23))&
                                                            (d.time>=tlate0[i_sim])]*1.0e9))

            t_1or2agn_late_post_arr[i_sim,i_mask,:] = t_1or2agn_late + t_1agn_post
            t_dualagn_late_arr[i_sim,i_mask,:] = t_dualagn_late


            print "t_totagn_obsc: ",t_totagn_obsc
            print "t_totagn_obsc_late: ",t_totagn_obsc_late
            print "t_1agn_obsc_post: ",t_1agn_obsc_post
            print "t(obsc agn, adv):",t_totagn_obsc_late+t_1agn_obsc_post
            print "dM_agn_obsc_adv: ",dM_agn_obsc_adv

            nproj_t1agnnonzero_post = len(t_1agn_post[t_1agn_post>0])
            #if nproj_t1agnnonzero_post>0:
            #    f1agn_post[i_sim,i_mask] = np.nansum( t_1agn_post/t_post ) / (1.0*nproj_t1agnnonzero_post)
            f1agn_post[i_sim,i_mask] = np.sum( t_1agn_post/t_post ) / (1.0*nproj_tnonzero_post)
            print f1agn_post[i_sim,i_mask]

            ### Now for some sanity checks...
            assert_isclose(ttot,t_dualagn+t_oneagn+t_noagn,'t_dualagn+t_oneagn+t_noagn!=ttot')
            tmptotdt=dualagn_histproj_dt[i_sim,i_mask,:,:]+oneagn_histproj_dt[i_sim,i_mask,:,:]+noagn_histproj_dt[i_sim,i_mask,:,:]
            assert np.allclose((np.nansum(tmptotdt,axis=1)),ttot), 'Error: tmptotdt!=ttot'
            tmptotavgdt=dualagn_avghistproj_dt[i_sim,i_mask,:]+oneagn_avghistproj_dt+noagn_avghistproj_dt
            assert np.allclose((np.nansum(tmptotavgdt)),ttot), 'Error: tmptotavgdt!=ttot'
            tmptot=dualagn_histproj_dt_frac[i_sim,i_mask,:,:]+oneagn_histproj_dt_frac[i_sim,i_mask,:,:]+noagn_histproj_dt_frac[i_sim,i_mask,:,:]
            assert np.allclose(tmptot[tmptot==tmptot],1), 'Error: tmptot!=1'
            tmptotavg=dualagn_avghistproj_dt_frac[i_sim,i_mask,:]+oneagn_avghistproj_dt_frac+noagn_avghistproj_dt_frac
            assert np.allclose(tmptotavg[tmptotavg==tmptotavg],1), 'Error: tmptotavg!=1'
            tmptotavgb=oneortwoagn_avghistproj_dt_frac[i_sim,i_mask,:]+noagn_avghistproj_dt_frac
            assert np.allclose(tmptotavg[tmptotavgb==tmptotavgb],1), 'Error: tmptotavg!=1'
            
            assert_isclose(dualagn_hist3d_dt[i_sim,i_mask,:].sum(),np.nansum(dualagn_avghistproj_dt[i_sim,i_mask,:]),'dual agn')
            assert_isclose(oneagn_hist3d_dt[i_sim,i_mask,:].sum(),np.nansum(oneagn_avghistproj_dt),'one agn')
            assert_isclose(noagn_hist3d_dt[i_sim,i_mask,:].sum(),np.nansum(noagn_avghistproj_dt),'no agn')

            if i_mask>2:
                print "tot:"
                print t_totagn
                print t_1or2agn_tot
                print "early:"
                print t_totagn_early
                print t_1or2agn_early
                print "late:"
                print t_totagn_late
                print t_1or2agn_late
                assert t_totagn>=t_1or2agn_tot, \
                    'Error: t_totagn<t_1or2agn_tot (%g, %g)'%(t_totagn,t_1or2agn_tot)
                assert (np.float32(t_totagn_early)>=np.float32(t_1or2agn_early)).all(), \
                    'Error: t_totagn_early<t_1or2agn_early. '
                assert (np.float32(t_totagn_late)>=np.float32(t_1or2agn_late)).all(), \
                    'Error: t_totagn_late<t_1or2agn_late (%g, %g)'
                assert ftotagn_tot[i_sim,i_mask]>=f1or2agn_tot[i_sim,i_mask], \
                    'Error: ftotagn_tot<f1or2agn_tot (%g, %g)'%(ftotagn_tot[i_sim,i_mask],
                                                                f1or2agn_tot[i_sim,i_mask])
                assert (ftotagn_early[i_sim,i_mask]>=f1or2agn_early[i_sim,i_mask]).all(), \
                    'Error: ftotagn_early<f1or2agn_early (%g, %g)'%(ftotagn_early[i_sim,i_mask],
                                                                    f1or2agn_early[i_sim,i_mask])
                assert (ftotagn_late[i_sim,i_mask]>=f1or2agn_late[i_sim,i_mask]).all(), \
                    'Error: ftotagn_late<f1or2agn_late (%g, %g)'%(ftotagn_late[i_sim,i_mask],
                                                                  f1or2agn_late[i_sim,i_mask])


            ### Print some results ###
            print "\n%s:"%agn_labels[i_mask]
            print "sep [kpc]:", sep_bins[:idx_sep+1]
            print "all AGN [Myr] (total=%.4g, total(<%.2gkpc)=%.4g):"%(1000*agn_avghistproj_dt[i_sim,i_mask,:].sum(),latesep,
                                                                       1000*agn_avghistproj_dt[i_sim,i_mask,:idx_sep+1].sum())
            print 1000*agn_avghistproj_dt[i_sim,i_mask,:idx_sep+1]

            print "all dual AGN [Myr] (total=%.4g, total(<%.2gkpc)=%.4g):"%(1000*dualagn_avghistproj_dt[i_sim,i_mask,:].sum(),latesep,
                                                                          1000*dualagn_avghistproj_dt[i_sim,i_mask,:idx_sep+1].sum())
            print 1000*dualagn_avghistproj_dt[i_sim,i_mask,:idx_sep+1]

            print "1 AGN [Myr] (total=%.4g, total(<%.2gkpc)=%.4g):"%(1000*oneagn_avghistproj_dt.sum(),latesep,
                                                                   1000*oneagn_avghistproj_dt[:idx_sep+1].sum())
            print 1000*oneagn_avghistproj_dt[:idx_sep+1]

            print "1or2 AGN [Myr] (total=%.4g, total(<%.2gkpc)=%.4g):"%(1000*(oneagn_avghistproj_dt.sum()+
                                                                              dualagn_avghistproj_dt[i_sim,i_mask,:].sum()),
                                                                      latesep,1000*(oneagn_avghistproj_dt[:idx_sep+1].sum()+
                                                                                    dualagn_avghistproj_dt[i_sim,i_mask,:idx_sep+1].sum()))
            print 1000*(oneagn_avghistproj_dt[:idx_sep+1]+dualagn_avghistproj_dt[i_sim,i_mask,:idx_sep+1])
            print "should match this: "
            print "1or2 AGN [Myr] (total=%.4g, total(<%.2gkpc)=%.4g):"%(1000*(oneortwoagn_avghistproj_dt[i_sim,i_mask,:]).sum(),latesep,
                                                                        1000*(oneortwoagn_avghistproj_dt[i_sim,i_mask,:idx_sep+1].sum()))
            print 1000*(oneortwoagn_avghistproj_dt[i_sim,i_mask,:idx_sep+1])
            print "t_dualagn/(ttot): %g"%(fdualagn_tot[i_sim,i_mask])
            print "t_dualagn/(t_dualagn+t_oneagn): %g"%(t_dualagn/t_1or2agn_tot)
            print "fdualagn_tot/(f1or2agn_tot): %g"%(fdualagn_tot[i_sim,i_mask]/f1or2agn_tot[i_sim,i_mask])


        ### End of loop over AGN mask arrays ###



        if single_sim_plots and not agnx0:

            ###################################
            ### AGN HISTOGRAMS FOR EACH SIM ###
            ###################################

            print "\nplotting AGN histograms for each sim..."
            plot_dt_hist_onesim(sep_bins,dsep,histproj_dt[i_sim,:,:], 
                                dt_hists = (oneortwoagn_avghistproj_dt[i_sim,:],
                                            agn_avghistproj_dt[i_sim,:],dualagn_avghistproj_dt[i_sim,:]),
                                dtfrac_hists = (oneortwoagn_avghistproj_dt_frac[i_sim,:],
                                                agn_avghistproj_dt_frac[i_sim,:],dualagn_avghistproj_dt_frac[i_sim,:]),
                                leg_lbls=('AGN','AGN (L1+2)','Dual'),eq_bins=equal_bins,ncam=ncam,
                                plot_titles=plot_titles,path=path,subdir=subdir_arr[i_sim],
                                binstr=binstr,tmax_postmrg=tmax_postmrg,fbase='agn')
            

        print "Calculating WISE AGN timescales..."

        for i_mask in np.arange(len(agn_wmaskarr_list)):
            ### recall:
            ### agn_wmaskarr_list = fedd_wmaskarr_list + lagn_wmaskarr_list + flbol_wmaskarr_list
            ### dualagn_wmaskarr_list = dual_fedd_wmaskarr_list + dual_lagn_wmaskarr_list + dual_flbol_wmaskarr_list

            print "\n%s:\n"%wise_agn_titles[i_mask]

            tot_wmask = agn_wmaskarr_list[i_mask][0]
            tot_nowmask = agn_wmaskarr_list[i_mask][1]
            notot_wmask = agn_wmaskarr_list[i_mask][2]
            notot_nowmask = agn_wmaskarr_list[i_mask][3]
            dual_wmask = dualagn_wmaskarr_list[i_mask][0]
            dual_nowmask = dualagn_wmaskarr_list[i_mask][1]
            one_wmask = dualagn_wmaskarr_list[i_mask][2]
            one_nowmask = dualagn_wmaskarr_list[i_mask][3]
            no_wmask = dualagn_wmaskarr_list[i_mask][4]
            no_nowmask = dualagn_wmaskarr_list[i_mask][5]

            ### Time when total AGN luminosity exceeds threshold & meets WISE color criteria ###            
            t_totagn_wise = np.array([d.dt[tot_wmask[j,:]].sum() for j in range(ncam)])
            ftotagn_wise_tot[i_sim,i_mask] = np.mean(t_totagn_wise) / ttot
            wise_agn_histproj_dt[i_sim,i_mask,:,:] = np.array([np.histogram(projsepvar[j,tot_wmask[j,:]],bins=sep_bins,
                                                                            weights=d.dt[tot_wmask[j,:]])[0]
                                                               for j in range(ncam)])
            nproj_twagnnonzero_eachbin=np.array([len(wise_agn_histproj_dt[i_sim,i_mask,
                                                                          wise_agn_histproj_dt[i_sim,i_mask,:,k]>0,k]) 
                                                 for k in range(len(sep_bins)-1)])
            wise_agn_avghistproj_dt[i_sim,i_mask,:] = np.nansum(wise_agn_histproj_dt[i_sim,i_mask,:,:],axis=0) / (1.0*ncam)
            wise_agn_histproj_dt_frac[i_sim,i_mask,:,:] = wise_agn_histproj_dt[i_sim,i_mask,:,:]/ttot_projsep_bins
            wise_agn_avghistproj_dt_frac[i_sim,i_mask,:] = (np.nansum(wise_agn_histproj_dt_frac[i_sim,i_mask,:,:],axis=0)/
                                                      (1.0*nproj_tnonzero_eachbin))

            t_nototagn_wise = np.array([d.dt[notot_wmask[j,:]].sum() for j in range(ncam)])
            fnototagn_wise_tot[i_sim,i_mask] = np.mean(t_nototagn_wise) / ttot
            wise_nototagn_histproj_dt[i_sim,i_mask,:,:]=np.array([np.histogram(projsepvar[j,notot_wmask[j,:]],bins=sep_bins, 
                                                                               weights=d.dt[notot_wmask[j,:]])[0] 
                                                               for j in range(ncam)])
            print "t_nototagn_wise:" ,t_nototagn_wise
            print "fnototagn_wise_tot[i_sim,i_mask]: ",fnototagn_wise_tot[i_sim,i_mask]
            wise_nototagn_avghistproj_dt[i_sim,i_mask,:] = (np.nansum(wise_nototagn_histproj_dt[i_sim,i_mask,:,:],axis=0) / 
                                                            (1.0*ncam))
            wise_nototagn_histproj_dt_frac[i_sim,i_mask,:,:]=wise_nototagn_histproj_dt[i_sim,i_mask,:,:]/ttot_projsep_bins
            wise_nototagn_avghistproj_dt_frac[i_sim,i_mask,:] = (np.nansum(wise_nototagn_histproj_dt[i_sim,i_mask,:,:]/
                                                                           (1.0*ttot_projsep_bins),axis=0) / 
                                                                 (1.0*nproj_tnonzero_eachbin))

            t_totagn_nowise = np.array([d.dt[tot_nowmask[j,:]].sum() for j in range(ncam)])
            ftotagn_nowise_tot[i_sim,i_mask] = np.mean(t_totagn_nowise) / ttot
            t_nototagn_nowise = np.array([d.dt[notot_nowmask[j,:]].sum() for j in range(ncam)])
            agn_nowise_histproj_dt[i_sim,i_mask,:,:] = np.array([np.histogram(projsepvar[j,tot_nowmask[j,:]],bins=sep_bins,
                                                                              weights=d.dt[tot_nowmask[j,:]])[0] 
                                                                 for j in range(ncam)])
            agn_nowise_avghistproj_dt[i_sim,i_mask,:] = (np.nansum(agn_nowise_histproj_dt[i_sim,i_mask,:,:],axis=0) / 
                                                         (1.0*ncam))
            agn_nowise_histproj_dt_frac[i_sim,i_mask,:,:]=agn_nowise_histproj_dt[i_sim,i_mask,:,:]/ttot_projsep_bins
            agn_nowise_avghistproj_dt_frac[i_sim,i_mask,:] = (np.nansum(agn_nowise_histproj_dt[i_sim,i_mask,:,:]/
                                                                        (1.0*ttot_projsep_bins),axis=0) / 
                                                              (1.0*nproj_tnonzero_eachbin))


            ### Time with two AGN simultaneously active & meets WISE color criteria ###
            t_dualagn_wise = np.array([d.dt[dual_wmask[j,:]].sum() for j in range(ncam)])
            fdualagn_wise_tot[i_sim,i_mask] = np.mean(t_dualagn_wise) / ttot
            wise_dualagn_histproj_dt[i_sim,i_mask,:,:]=np.array([np.histogram(projsepvar[j,dual_wmask[j,:]],bins=sep_bins, 
                                                                              weights=d.dt[dual_wmask[j,:]])[0] 
                                                                 for j in range(ncam)])
            nproj_t2wagnnonzero_eachbin=np.array([len(wise_dualagn_histproj_dt[i_sim,i_mask,
                                                                               wise_dualagn_histproj_dt[i_sim,i_mask,:,k]>0,k]) 
                                                  for k in range(len(sep_bins)-1)])
            wise_dualagn_avghistproj_dt[i_sim,i_mask,:] = np.nansum(wise_dualagn_histproj_dt[i_sim,i_mask,:,:],axis=0) / (1.0*ncam)
            wise_dualagn_histproj_dt_frac[i_sim,i_mask,:,:] = wise_dualagn_histproj_dt[i_sim,i_mask,:,:]/ttot_projsep_bins
            wise_dualagn_avghistproj_dt_frac[i_sim,i_mask,:] = (np.nansum(wise_dualagn_histproj_dt_frac[i_sim,i_mask,:,:],axis=0)/
                                                          (1.0*nproj_tnonzero_eachbin))

            t_dualagn_nowise = np.array([d.dt[dual_nowmask[j,:]].sum() for j in range(ncam)])
            fdualagn_nowise_tot[i_sim,i_mask] = np.mean(t_dualagn_nowise) / ttot
            dualagn_nowise_histproj_dt[i_sim,i_mask,:,:]=np.array([np.histogram(projsepvar[j,dual_nowmask[j,:]],
                                                                                bins=sep_bins, 
                                                                                weights=d.dt[dual_nowmask[j,:]])[0] 
                                                                   for j in range(ncam)])
            dualagn_nowise_avghistproj_dt[i_sim,i_mask,:] = np.nansum(dualagn_nowise_histproj_dt[i_sim,i_mask,:,:],axis=0) / (1.0*ncam)
            dualagn_nowise_histproj_dt_frac[i_sim,i_mask,:,:] = dualagn_nowise_histproj_dt[i_sim,i_mask,:,:]/ttot_projsep_bins
            dualagn_nowise_avghistproj_dt_frac[i_sim,i_mask,:] = (np.nansum(dualagn_nowise_histproj_dt_frac[i_sim,i_mask,:,:],axis=0)/
                                                          (1.0*nproj_tnonzero_eachbin))


            ### Time with exactly one AGN active & meets WISE color criteria ###
            t_oneagn_wise = np.array([d.dt[one_wmask[j,:]].sum() for j in range(ncam)])
            f1agn_wise_tot[i_sim,i_mask] = np.mean(t_oneagn_wise) / ttot
            wise_oneagn_histproj_dt[i_sim,i_mask,:,:]=np.array([np.histogram(projsepvar[j,one_wmask[j,:]],bins=sep_bins, 
                                                                             weights=d.dt[one_wmask[j,:]])[0] 
                                                                for j in range(ncam)])
            nproj_t1wagnnonzero_eachbin=np.array([len(wise_oneagn_histproj_dt[i_sim,i_mask,
                                                                              wise_oneagn_histproj_dt[i_sim,i_mask,:,k]>0,k]) 
                                                  for k in range(len(sep_bins)-1)])
            wise_oneagn_avghistproj_dt[i_sim,i_mask,:] = np.nansum(wise_oneagn_histproj_dt[i_sim,i_mask,:,:],axis=0) / (1.0*ncam)
            wise_oneagn_histproj_dt_frac[i_sim,i_mask,:,:]=wise_oneagn_histproj_dt[i_sim,i_mask,:,:]/ttot_projsep_bins
            wise_oneagn_avghistproj_dt_frac[i_sim,i_mask,:] = (np.nansum(wise_oneagn_histproj_dt_frac[i_sim,i_mask,:,:],axis=0)/
                                               (1.0*nproj_tnonzero_eachbin))

            t_oneagn_nowise = np.array([d.dt[one_nowmask[j,:]].sum() for j in range(ncam)])
            f1agn_nowise_tot[i_sim,i_mask] = np.mean(t_oneagn_nowise) / ttot

            ### Time with no AGN active & meets WISE color criteria ###
            t_noagn_wise = np.array([d.dt[no_wmask[j,:]].sum() for j in range(ncam)])
            fnoagn_wise_tot[i_sim,i_mask] = np.mean(t_noagn_wise) / ttot
            print "t_noagn_wise:" ,t_noagn_wise
            print "fnoagn_wise_tot[i_sim,i_mask]: ",fnoagn_wise_tot[i_sim,i_mask]

            t_noagn_nowise = np.array([d.dt[no_nowmask[j,:]].sum() for j in range(ncam)])

            t_1or2agn_wise = t_dualagn_wise + t_oneagn_wise
            t_1or2agn_wise_arr[i_sim,i_mask,:] = t_1or2agn_wise
            f1or2agn_wise_tot[i_sim,i_mask] = np.mean(t_1or2agn_wise) / ttot

            wise_oneortwoagn_histproj_dt = wise_oneagn_histproj_dt + wise_dualagn_histproj_dt
            wise_oneortwoagn_histproj_dt_frac = wise_oneortwoagn_histproj_dt / ttot_projsep_bins
            ### Note: need to resolve issue with diferent definitions of i_mask for wise vs total agn mask arrays
            #print wise_oneortwoagn_histproj_dtagn_frac.shape
            #print wise_oneortwoagn_histproj_dt.shape
            #print oneortwoagn_histproj_dt.shape
            #wise_oneortwoagn_histproj_dtagn_frac[i_sim,i_mask,:,:] = wise_oneortwoagn_histproj_dt[i_sim,i_mask,:,:] / oneortwoagn_histproj_dt[i_sim,i_mask,:,:]
            nproj_tw1or2agnnonzero_eachbin=np.array([len(wise_oneortwoagn_histproj_dt[i_sim,i_mask,wise_oneortwoagn_histproj_dt[i_sim,i_mask,:,k]>0,k])
                                                    for k in range(len(sep_bins)-1)])
            #print wise_oneortwoagn_histproj_dtagn_frac.shape
            #print wise_oneortwoagn_avghistproj_dtagn_frac.shape
            print nproj_t1or2agnnonzero_eachbin.shape
            if np.max(nproj_tw1or2agnnonzero_eachbin)>0:
                wise_oneortwoagn_avghistproj_dt[i_sim,i_mask,:] = np.nansum(wise_oneortwoagn_histproj_dt[i_sim,i_mask,:,:],axis=0)/(1.0*ncam)
                wise_oneortwoagn_avghistproj_dt_frac[i_sim,i_mask,:] = np.nansum(wise_oneortwoagn_histproj_dt_frac[i_sim,i_mask,:,:],axis=0)/(1.0*nproj_tnonzero_eachbin)
                #wise_oneortwoagn_avghistproj_dtagn_frac[i_sim,i_mask,:] = np.nansum(wise_oneortwoagn_histproj_dtagn_frac[i_sim,i_mask,:,:],axis=0)/(1.0*nproj_t1or2agnnonzero_eachbin)

            t_1or2agn_nowise = t_dualagn_nowise + t_oneagn_nowise
            f1or2agn_nowise_tot[i_sim,i_mask] = np.mean(t_1or2agn_nowise) / ttot
            oneortwoagn_nowise_histproj_dt[i_sim,i_mask,:,:] = np.array([np.histogram(projsepvar[j,(dual_nowmask[j,:]|one_nowmask[j,:])],
                                                                                  bins=sep_bins,
                                                                                  weights=d.dt[(dual_nowmask[j,:]|one_nowmask[j,:])])[0] 
                                                                     for j in range(ncam)])
            oneortwoagn_nowise_avghistproj_dt[i_sim,i_mask,:] = np.nansum(oneortwoagn_nowise_histproj_dt[i_sim,i_mask,:,:],axis=0) / (1.0*ncam)
            oneortwoagn_nowise_histproj_dt_frac[i_sim,i_mask,:,:]=oneortwoagn_nowise_histproj_dt[i_sim,i_mask,:,:]/ttot_projsep_bins
            #oneortwoagn_nowise_histproj_dtagn_frac[i_sim,i_mask,:,:]=oneortwoagn_nowise_histproj_dt[i_sim,i_mask,:,:]/ oneortwoagn_histproj_dt[i_sim,i_mask,:]

            oneortwoagn_nowise_avghistproj_dt_frac[i_sim,i_mask,:] = (np.nansum(oneortwoagn_nowise_histproj_dt_frac[i_sim,i_mask,:,:],axis=0)/
                                                                (1.0*nproj_tnonzero_eachbin))
            #if np.max(nproj_t1or2agnnonzero_eachbin)>0:
                #oneortwoagn_nowise_avghistproj_dtagn_frac[i_sim,i_mask,:] = np.nansum(oneortwoagn_nowise_histproj_dtagn_frac[i_sim,i_mask,:,:],axis=0)/(1.0*nproj_t1or2agnnonzero_eachbin)

            t_wise = t_totagn_wise + t_nototagn_wise
            print t_wise
            print t_1or2agn_wise
            print t_noagn_wise
            assert np.allclose(t_wise,t_1or2agn_wise+t_noagn_wise),'Error: total wise tscales do not match.'
            assert np.allclose(t_wise+t_totagn_nowise+t_nototagn_nowise,ttot),'Error: total time in agn wmaskarr != ttot.'
            assert np.allclose(t_dualagn_wise+t_dualagn_nowise+t_oneagn_wise+
                               t_oneagn_nowise+t_noagn_wise+t_noagn_nowise,ttot),\
                               'Error: total time in dualagn wmaskarr != ttot.'
            
            tmptotwisedt=(wise_agn_histproj_dt[i_sim,i_mask,:,:]+wise_nototagn_histproj_dt[i_sim,i_mask,:,:])

            assert np.allclose((np.nansum(tmptotwisedt,axis=1)),t_wise), 'Error: tmptotwisedt!=t_wise'
            tmptotavgwisedt=wise_agn_avghistproj_dt[i_sim,i_mask,:]+wise_nototagn_avghistproj_dt[i_sim,i_mask,:]
            assert np.allclose((np.nansum(tmptotavgwisedt)),np.mean(t_wise)), 'Error: tmptotavgwisedt!=t_wise'



            ### Time-weighted WISE AGN column density ###
            Ngas_bh1_wiseagn_histproj[i_sim,i_mask,:,:] = np.array([np.histogram(projsepvar[j,tot_wmask[j,:]],bins=sep_bins, 
                                                                                 weights=d.dt[tot_wmask[j,:]]*Ngas_bh1[j,tot_wmask[j,:]])[0]
                                                                    for j in range(ncam)])/ttot_projsep_bins
            Ngas_bh1_wiseagn_histproj[Ngas_bh1_wiseagn_histproj==0] = np.nan

            wts = np.zeros(Ngas_bh2.shape)
            wts[Ngas_bh2==Ngas_bh2] = dt_tiled[Ngas_bh2==Ngas_bh2]*Ngas_bh2[Ngas_bh2==Ngas_bh2]
            Ngas_bh2_wiseagn_histproj[i_sim,i_mask,:,:] = np.array([np.histogram(projsepvar[j,tot_wmask[j,:]],bins=sep_bins, 
                                                                                 weights=wts[j,tot_wmask[j,:]])[0]
                                                                    for j in range(ncam)])/ttot_projsep_bins
            Ngas_bh2_wiseagn_histproj[Ngas_bh2_wiseagn_histproj==0] = np.nan

            Ngas_wiseagn_hist_postmrg[i_sim,i_mask,:,:] = np.array([np.histogram(d.time[~has2bh][tot_wmask[j,~has2bh]]-d.tmrg,bins=tpostmrg_bins, 
                                                                                 weights=d.dt[~has2bh][tot_wmask[j,~has2bh]]*Ngas_bh1[j,~has2bh][tot_wmask[j,~has2bh]])[0] 
                                                                    for j in range(ncam)])/dtbin
            Ngas_wiseagn_hist_postmrg[Ngas_wiseagn_hist_postmrg==0.0] = np.nan


            for il,llim in enumerate((1.0e44,1.0e45)):
                Ngas_bh1_wiseagn_mycut_histproj[i_sim,il,:,:] = np.array([np.histogram(projsepvar[j,(((d.lagn1>=llim)|(d.lagn2>=llim))&(mycutmask[j,:]))],bins=sep_bins,
                                                                                       weights=(d.dt[(((d.lagn1>=llim)|(d.lagn2>=llim))&(mycutmask[j,:]))]*
                                                                                                Ngas_bh1[j,(((d.lagn1>=llim)|(d.lagn2>=llim))&(mycutmask[j,:]))]))[0]
                                                                          for j in range(ncam)])/ttot_projsep_bins
                Ngas_bh1_wiseagn_mycut_histproj[Ngas_bh1_wiseagn_mycut_histproj==0] = np.nan
                
                Ngas_bh2_wiseagn_mycut_histproj[i_sim,il,:,:] = np.array([np.histogram(projsepvar[j,(((d.lagn1>=llim)|(d.lagn2>=llim))&(mycutmask[j,:]))],bins=sep_bins,
                                                                                       weights=wts[j,(((d.lagn1>=llim)|(d.lagn2>=llim))&(mycutmask[j,:]))])[0]
                                                                          for j in range(ncam)])/ttot_projsep_bins
                Ngas_bh2_wiseagn_mycut_histproj[Ngas_bh2_wiseagn_mycut_histproj==0] = np.nan
                
                Ngas_wiseagn_mycut_hist_postmrg[i_sim,il,:,:] = np.array([np.histogram(d.time[~has2bh][(((d.lagn1[~has2bh]>=llim)|(d.lagn2[~has2bh]>=llim))&
                                                                                                        (mycutmask[j,~has2bh]))]-d.tmrg, bins=tpostmrg_bins, 
                                                                                       weights=(d.dt[~has2bh][(((d.lagn1[~has2bh]>=llim)|(d.lagn2[~has2bh]>=llim))&mycutmask[j,~has2bh])]*
                                                                                                Ngas_bh1[j,~has2bh][(((d.lagn1[~has2bh]>=llim)|(d.lagn2[~has2bh]>=llim))&
                                                                                                                     mycutmask[j,~has2bh])]))[0]  for j in range(ncam)])/dtbin
                Ngas_wiseagn_mycut_hist_postmrg[Ngas_wiseagn_mycut_hist_postmrg==0.0] = np.nan



            ### Time-weighted WISE dual AGN column density ###
            Ngas_bh1_wisedualagn_histproj[i_sim,i_mask,:,:] = np.array([np.histogram(projsepvar[j,dual_wmask[j,:]],bins=sep_bins, 
                                                                                     weights=d.dt[dual_wmask[j,:]]*Ngas_bh1[j,dual_wmask[j,:]])[0]
                                                                        for j in range(ncam)])/ttot_projsep_bins
            Ngas_bh1_wisedualagn_histproj[Ngas_bh1_wisedualagn_histproj==0] = np.nan
            wts = np.zeros(Ngas_bh2.shape)
            wts[Ngas_bh2==Ngas_bh2] = dt_tiled[Ngas_bh2==Ngas_bh2]*Ngas_bh2[Ngas_bh2==Ngas_bh2]
            Ngas_bh2_wisedualagn_histproj[i_sim,i_mask,:,:] = np.array([np.histogram(projsepvar[j,dual_wmask[j,:]],bins=sep_bins, 
                                                                                     weights=wts[j,dual_wmask[j,:]])[0]
                                                                        for j in range(ncam)])/ttot_projsep_bins
            Ngas_bh2_wisedualagn_histproj[Ngas_bh2_wisedualagn_histproj==0] = np.nan

            for il,llim in enumerate((1.0e44,1.0e45)):
                Ngas_bh1_wisedualagn_mycut_histproj[i_sim,il,:,:] = np.array([np.histogram(projsepvar[j,((d.lagn1>=llim)&(d.lagn2>=llim)&(mycutmask[j,:]))],bins=sep_bins,
                                                                                           weights=(d.dt[((d.lagn1>=llim)&(d.lagn2>=llim)&(mycutmask[j,:]))]*
                                                                                                    Ngas_bh1[j,((d.lagn1>=llim)&(d.lagn2>=llim)&(mycutmask[j,:]))]))[0]
                                                                              for j in range(ncam)])/ttot_projsep_bins
                Ngas_bh1_wisedualagn_mycut_histproj[Ngas_bh1_wisedualagn_mycut_histproj==0] = np.nan
                Ngas_bh2_wisedualagn_mycut_histproj[i_sim,il,:,:] = np.array([np.histogram(projsepvar[j,((d.lagn1>=llim)&(d.lagn2>=llim)&(mycutmask[j,:]))],bins=sep_bins, 
                                                                                           weights=wts[j,((d.lagn1>=llim)&(d.lagn2>=llim)&(mycutmask[j,:]))])[0]
                                                                              for j in range(ncam)])/ttot_projsep_bins
                Ngas_bh2_wisedualagn_mycut_histproj[Ngas_bh2_wisedualagn_mycut_histproj==0] = np.nan




            ### WISE AGN by merger phase ###

            t_wise_early = np.array([d.dt[(d.bhprojsep[j,:]>=latesep)&
                                          ((tot_wmask[j,:])|(notot_wmask[j,:]))].sum() for j in range(ncam)])
            t_totagn_wise_early = np.array([d.dt[(d.bhprojsep[j,:]>=latesep)&(tot_wmask[j,:])].sum() for j in range(ncam)])
            ftotagn_wise_early[i_sim,i_mask] = np.nansum( t_totagn_wise_early / t_early ) / (1.0*nproj_tnonzero_early)

            t_totagn_nowise_early = np.array([d.dt[(d.bhprojsep[j,:]>=latesep)&
                                                   (tot_nowmask[j,:])].sum() for j in range(ncam)])
            ftotagn_nowise_early[i_sim,i_mask] = np.nansum( t_totagn_nowise_early / t_early ) / (1.0*nproj_tnonzero_early)

            t_nototagn_wise_early = np.array([d.dt[(d.bhprojsep[j,:]>=latesep)&
                                                   (notot_wmask[j,:])].sum() for j in range(ncam)])
            fnototagn_wise_early[i_sim,i_mask] = np.nansum( t_nototagn_wise_early / t_early ) / (1.0*nproj_tnonzero_early)

            t_dualagn_wise_early = np.array([d.dt[(d.bhprojsep[j,:]>=latesep)&(dual_wmask[j,:])].sum() for j in range(ncam)])
            fdualagn_wise_early[i_sim,i_mask] = np.nansum( t_dualagn_wise_early / t_early ) / (1.0*nproj_tnonzero_early)
            t_dualagn_nowise_early = np.array([d.dt[(d.bhprojsep[j,:]>=latesep)&
                                                    (dual_nowmask[j,:])].sum() for j in range(ncam)])
            fdualagn_nowise_early[i_sim,i_mask] = np.nansum( t_dualagn_nowise_early / t_early ) / (1.0*nproj_tnonzero_early)

            t_oneagn_wise_early = np.array([d.dt[(d.bhprojsep[j,:]>=latesep)&(one_wmask[j,:])].sum() for j in range(ncam)])
            f1agn_wise_early[i_sim,i_mask] = np.nansum( t_oneagn_wise_early / t_early ) / (1.0*nproj_tnonzero_early)
            t_oneagn_nowise_early = np.array([d.dt[(d.bhprojsep[j,:]>=latesep)&
                                                   (one_nowmask[j,:])].sum() for j in range(ncam)])
            f1agn_nowise_early[i_sim,i_mask] = np.nansum( t_oneagn_nowise_early / t_early ) / (1.0*nproj_tnonzero_early)

            t_noagn_wise_early = np.array([d.dt[(d.bhprojsep[j,:]>=latesep)&(no_wmask[j,:])].sum() for j in range(ncam)])
            fnoagn_wise_early[i_sim,i_mask] = np.nansum( t_noagn_wise_early / t_early ) / (1.0*nproj_tnonzero_early)

            t_1or2agn_wise_early = t_oneagn_wise_early + t_dualagn_wise_early
            f1or2agn_wise_early[i_sim,i_mask] = np.nansum( t_1or2agn_wise_early / t_early ) / (1.0*nproj_tnonzero_early)
            tmpfwearly=fdualagn_wise_early[i_sim,i_mask]+f1agn_wise_early[i_sim,i_mask]
            assert np.isclose(tmpfwearly[tmpfwearly==tmpfwearly],f1or2agn_wise_early[i_sim,i_mask]),\
                'fdual_wise_early+f1_wise_early!=f1or2agn_wise_early'
            tmpfwearly=tmpfwearly+fnoagn_wise_early[i_sim,i_mask]
            print "t_early=",t_early
            print "t_wise_early (mean=%g):"%t_wise_early.mean(),t_wise_early
            print "t_1or2agn_wise_early (mean=%g):"%t_1or2agn_wise_early.mean(),t_1or2agn_wise_early
            print "t_dualagn_wise_early=",t_dualagn_wise_early
            print "t_oneagn_wise_early=",t_oneagn_wise_early
            print "t_noagn_wise_early=",t_noagn_wise_early
            assert np.isclose(tmpfwearly[tmpfwearly==tmpfwearly],np.nansum(t_wise_early/t_early)/(1.0*nproj_tnonzero_early)),\
                'fdual_wise_early+fone_wise_early+fno_wise_early!=t_wise_early/t_early'


            t_wise_late = np.array([d.dt[(d.bhprojsep[j,:]<latesep)&(d.bhprojsep[j,:]>0)&
                                          ((tot_wmask[j,:])|(notot_wmask[j,:]))].sum() for j in range(ncam)])
            t_totagn_wise_late = np.array([d.dt[(d.bhprojsep[j,:]<latesep)&
                                                (d.bhprojsep[j,:]>0)&(tot_wmask[j,:])].sum() for j in range(ncam)])
            ftotagn_wise_late[i_sim,i_mask] = np.nansum( t_totagn_wise_late/t_late ) / (1.0*nproj_tnonzero_late)
            t_nototagn_wise_late = np.array([d.dt[(d.bhprojsep[j,:]<latesep)&
                                                  (d.bhprojsep[j,:]>0)&(notot_wmask[j,:])].sum() for j in range(ncam)])
            fnototagn_wise_late[i_sim,i_mask] = np.nansum( t_nototagn_wise_late/t_late ) / (1.0*nproj_tnonzero_late)
            t_totagn_nowise_late = np.array([d.dt[(d.bhprojsep[j,:]<latesep)&
                                                (d.bhprojsep[j,:]>0)&(tot_nowmask[j,:])].sum() for j in range(ncam)])
            ftotagn_nowise_late[i_sim,i_mask] = np.nansum( t_totagn_nowise_late/t_late ) / (1.0*nproj_tnonzero_late)

            t_dualagn_wise_late = np.array([d.dt[(d.bhprojsep[j,:]<latesep)&(d.bhprojsep[j,:]>0)&
                                                 (dual_wmask[j,:])].sum() for j in range(ncam)])
            fdualagn_wise_late[i_sim,i_mask] = np.nansum( t_dualagn_wise_late / t_late ) / (1.0*nproj_tnonzero_late)
            t_dualagn_nowise_late = np.array([d.dt[(d.bhprojsep[j,:]<latesep)&(d.bhprojsep[j,:]>0)&
                                                 (dual_nowmask[j,:])].sum() for j in range(ncam)])
            fdualagn_nowise_late[i_sim,i_mask] = np.nansum( t_dualagn_nowise_late / t_late ) / (1.0*nproj_tnonzero_late)

            t_oneagn_wise_late = np.array([d.dt[(d.bhprojsep[j,:]<latesep)&(d.bhprojsep[j,:]>0)&
                                                (one_wmask[j,:])].sum() for j in range(ncam)])
            f1agn_wise_late[i_sim,i_mask] = np.nansum( t_oneagn_wise_late / t_late ) / (1.0*nproj_tnonzero_late)
            t_oneagn_nowise_late = np.array([d.dt[(d.bhprojsep[j,:]<latesep)&(d.bhprojsep[j,:]>0)&
                                                  (one_nowmask[j,:])].sum() for j in range(ncam)])
            f1agn_nowise_late[i_sim,i_mask] = np.nansum( t_oneagn_nowise_late / t_late ) / (1.0*nproj_tnonzero_late)

            t_noagn_wise_late = np.array([d.dt[(d.bhprojsep[j,:]<latesep)&(d.bhprojsep[j,:]>0)&
                                               (no_wmask[j,:])].sum() for j in range(ncam)])
            fnoagn_wise_late[i_sim,i_mask] = np.nansum( t_noagn_wise_late / t_late ) / (1.0*nproj_tnonzero_late)

            t_1or2agn_wise_late = t_oneagn_wise_late + t_dualagn_wise_late
            f1or2agn_wise_late[i_sim,i_mask] = np.nansum( t_1or2agn_wise_late / t_late) / (1.0*nproj_tnonzero_late)
            t_1or2agn_nowise_late = t_oneagn_nowise_late + t_dualagn_nowise_late
            f1or2agn_nowise_late[i_sim,i_mask] = np.nansum( t_1or2agn_nowise_late / t_late) / (1.0*nproj_tnonzero_late)

            tmpfwlate=fdualagn_wise_late[i_sim,i_mask]+f1agn_wise_late[i_sim,i_mask]
            assert np.isclose(tmpfwlate[tmpfwlate==tmpfwlate],f1or2agn_wise_late[i_sim,i_mask]),\
                'fdual_wise_late+f1agn_wise_late!=f1or2agn_wise_late'
            tmpfwlate=tmpfwlate+fnoagn_wise_late[i_sim,i_mask]
            print "t_late=",t_late
            print "t_wise_late (mean=%g):"%t_wise_late.mean(),t_wise_late
            print "t_1or2agn_wise_late (mean=%g):"%t_1or2agn_wise_late.mean(),t_1or2agn_wise_late
            print "t_dualagn_wise_late=",t_dualagn_wise_late
            print "t_oneagn_wise_late=",t_oneagn_wise_late
            print "t_noagn_wise_late=",t_noagn_wise_late
            assert np.isclose(tmpfwlate[tmpfwlate==tmpfwlate],np.nansum(t_wise_late/t_late)/(1.0*nproj_tnonzero_late)),\
                'fdual_wise_late+f1agn_wise_late+fnoagn_wise_late!=t_wise_late/t_late'


            print "t_1agn_post (mean=%g):"%t_1agn_post.mean(),t_1agn_post
            print "should match this:", np.array([d.dt[((d.bhprojsep[j,:]<=0)&(totmask))].sum() for j in range(ncam)])

            t_wise_post = np.array([d.dt[(d.bhprojsep[j,:]<=0)&
                                         ((tot_wmask[j,:])|(notot_wmask[j,:]))].sum() for j in range(ncam)])
            print "t_wise_post (mean=%g):"%t_wise_post.mean(),t_wise_post
            print "should match this: ", np.array([d.dt[(d.bhprojsep[j,:]<=0)&
                                                        ((one_wmask[j,:])|(no_wmask[j,:]))].sum() for j in range(ncam)])

            t_1agn_wise_post = np.array([d.dt[(d.bhprojsep[j,:]<=0)&(one_wmask[j,:])].sum() for j in range(ncam)])
            f1agn_wise_post[i_sim,i_mask] = np.nansum(t_1agn_wise_post / t_post) / (1.0*nproj_tnonzero_post)

            print "t_1agn_wise_post (mean=%g):"%t_1agn_wise_post.mean(),t_wise_post

            t_1agn_nowise_post = np.array([d.dt[(d.bhprojsep[j,:]<=0)&(one_nowmask[j,:])].sum() for j in range(ncam)])
            f1agn_nowise_post[i_sim,i_mask] = np.nansum(t_1agn_nowise_post / t_post) / (1.0*nproj_tnonzero_post)

            t_nototagn_wise_post = np.array([d.dt[(d.bhprojsep[j,:]<=0)&(notot_wmask[j,:])].sum() for j in range(ncam)])
            fnototagn_wise_post[i_sim,i_mask] = np.nansum(t_nototagn_wise_post / t_post) / (1.0*nproj_tnonzero_post)

            t_noagn_wise_post = np.array([d.dt[(d.bhprojsep[j,:]<=0)&(no_wmask[j,:])].sum() for j in range(ncam)])
            fnoagn_wise_post[i_sim,i_mask] = np.nansum( t_noagn_wise_post / t_post ) / (1.0*nproj_tnonzero_post)

            print t_1agn_wise_post, t_1agn_nowise_post,t_nototagn_wise_post
            assert np.allclose(t_1agn_wise_post+t_nototagn_wise_post,t_wise_post),\
                'Error: t_1agn_wise_post+t_nototagn_wise_post!=t_wise_post'

            t_1or2agn_wise_late_post_arr[i_sim,i_mask,:] = t_1or2agn_wise_late + t_1agn_wise_post


        ########################################
        ### WISE AGN HISTOGRAMS FOR EACH SIM ###
        ########################################

        if single_sim_plots and not ngas_only:

            print "plotting WISE AGN histograms for each sim..."

            for iplot,ix in enumerate([ix_w05,ix_w08,ix_j11]):
                print '%s plots...'%wfbase[iplot]

                if not agnx0:
                    lbls=tuple([wisemask_labels[iplot]+s for s in (' AGN',' AGN (L1+2)',' dual AGN')])
                    plot_dt_hist_onesim(sep_bins,dsep,histproj_dt[i_sim,:,:], 
                                        dt_hists = (wise_oneortwoagn_avghistproj_dt[i_sim,ix,:],
                                                    wise_agn_avghistproj_dt[i_sim,ix,:],
                                                    wise_dualagn_avghistproj_dt[i_sim,ix,:]),
                                        dtfrac_hists = (wise_oneortwoagn_avghistproj_dt_frac[i_sim,ix,:],
                                                        wise_agn_avghistproj_dt_frac[i_sim,ix,:],
                                                        wise_dualagn_avghistproj_dt_frac[i_sim,ix,:]),
                                        leg_lbls=lbls,eq_bins=equal_bins,ncam=ncam,plot_titles=plot_titles,
                                        path=path,subdir=subdir_arr[i_sim],binstr=binstr,
                                        tmax_postmrg=tmax_postmrg,fbase='%s_agn'%wfbase[iplot])
                    
                    #lbls=tuple([wisemask_labels[iplot]+s for s in (' AGN',)])
                    #plot_dt_hist_onesim(sep_bins,dsep,histproj_dt[i_sim,:,:], 
                    #                    dtfrac_hists = (wise_oneortwoagn_avghistproj_dtagn_frac[i_sim,ix,:]),
                    #                    leg_lbls=lbls,eq_bins=equal_bins,ncam=ncam,plot_titles=plot_titles,
                    #                    path=path,subdir=subdir_arr[i_sim],binstr=binstr,
                    #                    tmax_postmrg=tmax_postmrg,fbase='%s_fagn'%wfbase[iplot])


                lbls=(wisemask_labels[iplot]+' AGN','AGN, not '+wisemask_labels[iplot],
                      wisemask_labels[iplot]+' no AGN')
                plot_dt_hist_onesim(sep_bins,dsep,histproj_dt[i_sim,:,:], 
                                    dt_hists = (wise_agn_avghistproj_dt[i_sim,ix,:],
                                                agn_nowise_avghistproj_dt[i_sim,ix,:],
                                                wise_nototagn_avghistproj_dt[i_sim,ix,:]),
                                    dtfrac_hists = (wise_agn_avghistproj_dt_frac[i_sim,ix,:],
                                                    agn_nowise_avghistproj_dt_frac[i_sim,ix,:],
                                                    wise_nototagn_avghistproj_dt_frac[i_sim,ix,:]),
                                    leg_lbls=lbls,eq_bins=equal_bins,ncam=ncam,plot_titles=plot_titles,
                                    path=path,subdir=subdir_arr[i_sim],binstr=binstr,
                                    tmax_postmrg=tmax_postmrg,fbase='%s_agn_noagn'%wfbase[iplot],
                                    carr=['darkblue','green','darkorange'],lsarr=['-','-',':'],lwarr=[2,1,2])
                
                if not agnx0:
                    #lbls=(wisemask_labels[iplot]+' dual AGN','Dual AGN, not '+wisemask_labels[iplot],
                    #      wisemask_labels[iplot]+' AGN, not dual')
                    #lbls=('Dual MIR AGN','Dual AGN, not MIR','MIR AGN, not dual')
                    lbls=('Dual MIR AGN','Dual, not MIR','MIR, not dual')
                    plot_dt_hist_onesim(sep_bins,dsep,histproj_dt[i_sim,:,:], 
                                        dt_hists = (wise_dualagn_avghistproj_dt[i_sim,ix,:],
                                                    dualagn_nowise_avghistproj_dt[i_sim,ix,:],
                                                    wise_oneagn_avghistproj_dt[i_sim,ix,:]),
                                        dtfrac_hists = (wise_dualagn_avghistproj_dt_frac[i_sim,ix,:],
                                                        dualagn_nowise_avghistproj_dt_frac[i_sim,ix,:],
                                                        wise_oneagn_avghistproj_dt_frac[i_sim,ix,:]),
                                        leg_lbls=lbls,eq_bins=equal_bins,ncam=ncam,plot_titles=plot_titles,
                                        path=path,subdir=subdir_arr[i_sim],binstr=binstr,
                                        title='%s (%s)'%(subdir_arr[i_sim],wisemask_labels[iplot]),
                                        tmax_postmrg=tmax_postmrg,fbase='%s_dual_nodual'%wfbase[iplot],
                                        carr=['m','green','k'],lsarr=['-','-',':'],lwarr=[2,1,2])

                    print subdir_arr[i_sim]
                    print get_sim_name(subdir_arr[i_sim])
                    plot_dt_hist_onesim(sep_bins,dsep,histproj_dt[i_sim,:,:], 
                                        dt_hists = (wise_dualagn_avghistproj_dt[i_sim,ix,:],
                                                    dualagn_nowise_avghistproj_dt[i_sim,ix,:],
                                                    wise_oneagn_avghistproj_dt[i_sim,ix,:]),
                                        dtfrac_hists = (wise_dualagn_avghistproj_dt_frac[i_sim,ix,:],
                                                        dualagn_nowise_avghistproj_dt_frac[i_sim,ix,:],
                                                        wise_oneagn_avghistproj_dt_frac[i_sim,ix,:]),
                                        leg_lbls=lbls,eq_bins=equal_bins,ncam=ncam,leg_eachplot=True,
                                        #plot_titles=[pt+wisemask_labels[iplot] for pt in plot_titles],
                                        plot_titles=[get_sim_name(subdir_arr[i_sim])],
                                        path=path,subdir=subdir_arr[i_sim],binstr=binstr,
                                        #title='%s (%s)'%(subdir_arr[i_sim],wisemask_labels[iplot]),
                                        title=' ', tmax_postmrg=tmax_postmrg,
                                        fbase='pubstyle_%s_dual_nodual'%wfbase[iplot],
                                        carr=['r','b','k'],lsarr=['-','-',':'],lwarr=[2,1,2])

        with open(maindir+'/global_sim_info.txt','a') as gfp: 
            gfp.write('\n\n%s'%path)
            gfp.write("\nMax SFR: %g"%sfr_max[i_sim])
            gfp.write("\nMax log sSFR: %g"%np.log10(ssfr_max[i_sim]))
            gfp.write("\nMax fedd: %g"%fedd_max[i_sim])
            gfp.write("\nMax log lagn: %g"%np.log10(lagn_max[i_sim]))
            gfp.write("\nMax flbol: %g"%flbol_max[i_sim])
            gfp.write("\nMax log ltot: %g"%np.log10(ltot_max[i_sim]))
            gfp.write("\nMax log NH (BH1, LOS med): %g"%np.log10(nh1_max_med[i_sim]))
            gfp.write("\nMax log NH (BH2, LOS med): %g"%np.log10(nh2_max_med[i_sim]))
            gfp.write("\nMax log NH (BH1): %g"%np.log10(nh1_max[i_sim]))
            gfp.write("\nMax log NH (BH2): %g"%np.log10(nh2_max[i_sim]))
            gfp.write("\nInit Mstar: %g"%mstar_init[i_sim])
            gfp.write("\nFinal M*: %g"%mstar_final[i_sim])
            gfp.write("\nDelta M*/Init M*: %g"%((mstar_final[i_sim]-mstar_init[i_sim])/mstar_init[i_sim]))
            gfp.write("\nDelta M*/Final M*: %g"%((mstar_final[i_sim]-mstar_init[i_sim])/mstar_final[i_sim]))
            gfp.write("\ntlate0: %g"%tlate0[i_sim])
            gfp.write("\ntmrg: %g"%tmrg[i_sim])
            gfp.write("\nInitial MBH: %g"%mbh_init[i_sim])
            gfp.write("\nMBH(tlate): %g"%mbh_tlate[i_sim])
            gfp.write("\nMBH(tmrg): %g"%mbh_tmrg[i_sim])
            gfp.write("\nFinal MBH: %g"%mbh_final[i_sim])
            gfp.write("\nDelta MBH/Init MBH: %g"%((mbh_final[i_sim]-mbh_init[i_sim])/mbh_init[i_sim]))
            gfp.write("\nDelta MBH/Final MBH: %g"%((mbh_final[i_sim]-mbh_init[i_sim])/mbh_final[i_sim]))
            gfp.write("\nMBH(tmax)-MBH(tlate)/MBH(tlate): %g"%((mbh_final[i_sim]-mbh_tlate[i_sim])/
                                                               mbh_tlate[i_sim]))
            gfp.write("\nMBH(tmax)-MBH(tlate)/MBH(tmax): %g"%((mbh_final[i_sim]-mbh_tlate[i_sim])/
                                                               mbh_final[i_sim]))
            ###fedd_lims = [0.01,0.05,0.1]
            ###lagn_lims = [1.0e43, 1.0e44, 1.0e45]
            ###flbol_lims = [0.1, 0.3, 0.5]

            gfp.write("\ndM_agn_obsc (1e44): %g"%dM_agn_obsc[i_sim,4])
            gfp.write("\ndM_agn_obsc (1e44)/Final MBH: %g"%(dM_agn_obsc[i_sim,4]/mbh_final[i_sim]))
            gfp.write("\ndM_agn_obsc (1e45): %g"%dM_agn_obsc[i_sim,5])
            gfp.write("\ndM_agn_obsc (1e45)/Final MBH: %g"%(dM_agn_obsc[i_sim,5]/mbh_final[i_sim]))
            gfp.write("\ndM_agn_obsc (fL0.1): %g"%dM_agn_obsc[i_sim,6])
            gfp.write("\ndM_agn_obsc (fL0.1)/Final MBH: %g"%(dM_agn_obsc[i_sim,6]/mbh_final[i_sim]))
            gfp.write("\ndM_agn_obsc_adv (1e44): %g"%dM_agn_obsc_adv[i_sim,4])
            gfp.write("\ndM_agn_obsc_adv (1e44)/Final MBH: %g"%(dM_agn_obsc_adv[i_sim,4]/mbh_final[i_sim]))
            gfp.write("\ndM_agn_obsc_adv (1e45): %g"%dM_agn_obsc_adv[i_sim,5])
            gfp.write("\ndM_agn_obsc_adv (1e45)/Final MBH: %g"%(dM_agn_obsc_adv[i_sim,5]/mbh_final[i_sim]))
            gfp.write("\ndM_agn_obsc_adv (fL0.1): %g"%dM_agn_obsc_adv[i_sim,6])
            gfp.write("\ndM_agn_obsc_adv (fL0.1)/Final MBH: %g"%(dM_agn_obsc_adv[i_sim,6]/mbh_final[i_sim]))
            
            gfp.write("\nt_totagn_obsc (1e44): %g"%t_totagn_obsc[i_sim,4])
            gfp.write("\nt_totagn_obsc,adv (1e44): %g"%(t_totagn_obsc_late[i_sim,4]+t_1agn_obsc_post[i_sim,4]))
            gfp.write("\nt_totagn_obsc (1e45): %g"%t_totagn_obsc[i_sim,5])
            gfp.write("\nt_totagn_obsc,adv (1e45): %g"%(t_totagn_obsc_late[i_sim,5]+t_1agn_obsc_post[i_sim,5]))
            gfp.write("\nt_totagn_obsc (fL0.1): %g"%t_totagn_obsc[i_sim,6])
            gfp.write("\nt_totagn_obsc,adv (fL0.1): %g"%(t_totagn_obsc_late[i_sim,6]+t_1agn_obsc_post[i_sim,6]))
            gfp.write("\nt_late+t_post: %g"%(t_late_arr[i_sim].mean()+t_post_arr[i_sim].mean()))



        ### Loop over hi-res wise mask arrays: ###
        agn_wmaskarr_hires_l44 = copy(l44_wmaskarr_hires)
        dualagn_wmaskarr_hires_l44 = copy(dual_l44_wmaskarr_hires)

        agn_wmaskarr_hires_l45 = copy(l45_wmaskarr_hires)
        dualagn_wmaskarr_hires_l45 = copy(dual_l45_wmaskarr_hires)


        #print len(agn_wmaskarr_hires)
        for i_wmask in range(len(agn_wmaskarr_hires_l44)):
        
            ### Lagn>1e44 ###
            dual_wmask_hires_l44 = dualagn_wmaskarr_hires_l44[i_wmask][0]
            one_wmask_hires_l44 = dualagn_wmaskarr_hires_l44[i_wmask][2]

            tot_wmask_hires_l44 = agn_wmaskarr_hires_l44[i_wmask][0]
            notot_wmask_hires_l44 = agn_wmaskarr_hires_l44[i_wmask][2]

            t_dualagn_wise_hires_l44_tot[i_sim,i_wmask,:] = np.array([d.dt[dual_wmask_hires_l44[j,:]].sum() for j in range(ncam)])
            t_oneagn_wise_hires_l44_tot = np.array([d.dt[one_wmask_hires_l44[j,:]].sum() for j in range(ncam)])
            t_1or2agn_wise_hires_l44_tot[i_sim,i_wmask,:] = (t_dualagn_wise_hires_l44_tot[i_sim,i_wmask,:] + t_oneagn_wise_hires_l44_tot)
            #f1or2agn_wise_hires_tot[i_sim,i_wmask] = np.mean(t_1or2agn_wise_hires_tot) / ttot

            t_dualagn_wise_hires_l44_late[i_sim,i_wmask,:] = np.array([d.dt[(d.bhprojsep[j,:]<latesep)&(d.bhprojsep[j,:]>0)&
                                                                            (dual_wmask_hires_l44[j,:])].sum() for j in range(ncam)])
            t_oneagn_wise_hires_l44_late = np.array([d.dt[(d.bhprojsep[j,:]<latesep)&(d.bhprojsep[j,:]>0)&
                                                          (one_wmask_hires_l44[j,:])].sum() for j in range(ncam)])
            t_1or2agn_wise_hires_l44_late = t_oneagn_wise_hires_l44_late + t_dualagn_wise_hires_l44_late[i_sim,i_wmask,:]
            
            t_1agn_wise_hires_l44_post = np.array([d.dt[(d.bhprojsep[j,:]<=0)&(one_wmask_hires_l44[j,:])].sum() for j in range(ncam)])

            t_1or2agn_wise_hires_l44_late_post[i_sim,i_wmask,:] = t_1or2agn_wise_hires_l44_late + t_1agn_wise_hires_l44_post
            

            t_wise_hires_l44_tot[i_sim,i_wmask,:] = np.array([d.dt[(tot_wmask_hires_l44[j,:])|
                                                                   (notot_wmask_hires_l44[j,:])].sum() for j in range(ncam)])
            t_wise_hires_l44_late = np.array([d.dt[(d.bhprojsep[j,:]<latesep)&(d.bhprojsep[j,:]>0)&
                                                   ((tot_wmask_hires_l44[j,:])|(notot_wmask_hires_l44[j,:]))].sum() 
                                              for j in range(ncam)])
            t_wise_hires_l44_post = np.array([d.dt[(d.bhprojsep[j,:]<=0)&((tot_wmask_hires_l44[j,:])|
                                                                          (notot_wmask_hires_l44[j,:]))].sum() 
                                              for j in range(ncam)])
            t_wise_hires_l44_late_post[i_sim,i_wmask,:] = t_wise_hires_l44_late + t_wise_hires_l44_post


            ### Lagn>1e45 ###
            dual_wmask_hires_l45 = dualagn_wmaskarr_hires_l45[i_wmask][0]
            one_wmask_hires_l45 = dualagn_wmaskarr_hires_l45[i_wmask][2]

            tot_wmask_hires_l45 = agn_wmaskarr_hires_l45[i_wmask][0]
            notot_wmask_hires_l45 = agn_wmaskarr_hires_l45[i_wmask][2]

            t_dualagn_wise_hires_l45_tot[i_sim,i_wmask,:] = np.array([d.dt[dual_wmask_hires_l45[j,:]].sum() for j in range(ncam)])
            t_oneagn_wise_hires_l45_tot = np.array([d.dt[one_wmask_hires_l45[j,:]].sum() for j in range(ncam)])
            t_1or2agn_wise_hires_l45_tot[i_sim,i_wmask,:] = (t_dualagn_wise_hires_l45_tot[i_sim,i_wmask,:] + t_oneagn_wise_hires_l45_tot)
            #f1or2agn_wise_hires_tot[i_sim,i_wmask] = np.mean(t_1or2agn_wise_hires_tot) / ttot

            t_dualagn_wise_hires_l45_late[i_sim,i_wmask,:] = np.array([d.dt[(d.bhprojsep[j,:]<latesep)&(d.bhprojsep[j,:]>0)&
                                                                            (dual_wmask_hires_l45[j,:])].sum() for j in range(ncam)])
            t_oneagn_wise_hires_l45_late = np.array([d.dt[(d.bhprojsep[j,:]<latesep)&(d.bhprojsep[j,:]>0)&
                                                (one_wmask_hires_l45[j,:])].sum() for j in range(ncam)])
            t_1or2agn_wise_hires_l45_late = t_oneagn_wise_hires_l45_late + t_dualagn_wise_hires_l45_late[i_sim,i_wmask,:]
            
            t_1agn_wise_hires_l45_post = np.array([d.dt[(d.bhprojsep[j,:]<=0)&(one_wmask_hires_l45[j,:])].sum() for j in range(ncam)])
            t_1or2agn_wise_hires_l45_late_post[i_sim,i_wmask,:] = t_1or2agn_wise_hires_l45_late + t_1agn_wise_hires_l45_post


            t_wise_hires_l45_tot[i_sim,i_wmask,:] = np.array([d.dt[(tot_wmask_hires_l45[j,:])|
                                                                   (notot_wmask_hires_l45[j,:])].sum() for j in range(ncam)])
            t_wise_hires_l45_late = np.array([d.dt[(d.bhprojsep[j,:]<latesep)&(d.bhprojsep[j,:]>0)&
                                                   ((tot_wmask_hires_l45[j,:])|(notot_wmask_hires_l45[j,:]))].sum() 
                                              for j in range(ncam)])
            t_wise_hires_l45_post = np.array([d.dt[(d.bhprojsep[j,:]<=0)&((tot_wmask_hires_l45[j,:])|
                                                                          (notot_wmask_hires_l45[j,:]))].sum() 
                                              for j in range(ncam)])
            t_wise_hires_l45_late_post[i_sim,i_wmask,:] = t_wise_hires_l45_late + t_wise_hires_l45_post

            #f1or2agn_wise_hires_late_post[i_sim,i_wmask] = (np.sum((t_1or2agn_wise_hires_late+t_1agn_wise_hires_post)/
            #                                                       (t_late+t_post))/
            #                                                (1.0*np.maximum(nproj_tnonzero_late,nproj_tnonzero_post)))

        ### End loop over hi-res wise mask arrays ###
            

        t_j11_tot[i_sim,:] = np.array([d.dt[(j11mask[j,:])].sum() for j in range(ncam)])
        t_j11_late_post[i_sim,:] = np.array([d.dt[(d.bhprojsep[j,:]<latesep)&(j11mask[j,:])].sum() 
                                             for j in range(ncam)])

        t_1or2agn_j11_l44_tot[i_sim,:] = np.array([d.dt[((d.lagn1>=1.0e44)|(d.lagn2>=1.0e44))&
                                                        (j11mask[j,:])].sum() for j in range(ncam)])
        t_1or2agn_j11_l44_late_post[i_sim,:] = np.array([d.dt[(d.bhprojsep[j,:]<latesep)&((d.lagn1>=1.0e44)|(d.lagn2>=1.0e44))&
                                                              (j11mask[j,:])].sum() for j in range(ncam)])
        t_1or2agn_j11_l45_tot[i_sim,:] = np.array([d.dt[((d.lagn1>=1.0e45)|(d.lagn2>=1.0e45))&(j11mask[j,:])].sum() for j in range(ncam)])
        t_1or2agn_j11_l45_late_post[i_sim,:] = np.array([d.dt[(d.bhprojsep[j,:]<latesep)&((d.lagn1>=1.0e45)|(d.lagn2>=1.0e45))&
                                                              (j11mask[j,:])].sum() for j in range(ncam)])

        t_dualagn_j11_l44_tot[i_sim,:] = np.array([d.dt[(d.lagn1>=1.0e44)&(d.lagn2>=1.0e44)&
                                                        (j11mask[j,:])].sum() for j in range(ncam)])
        t_dualagn_j11_l44_late[i_sim,:] = np.array([d.dt[(d.bhprojsep[j,:]<latesep)&(d.lagn1>=1.0e44)&(d.lagn2>=1.0e44)&
                                                              (j11mask[j,:])].sum() for j in range(ncam)])
        t_dualagn_j11_l45_tot[i_sim,:] = np.array([d.dt[(d.lagn1>=1.0e45)&(d.lagn2>=1.0e45)&
                                                        (j11mask[j,:])].sum() for j in range(ncam)])
        t_dualagn_j11_l45_late[i_sim,:] = np.array([d.dt[(d.bhprojsep[j,:]<latesep)&(d.lagn1>=1.0e45)&(d.lagn2>=1.0e45)&
                                                              (j11mask[j,:])].sum() for j in range(ncam)])

        t_mycut_tot[i_sim,:] = np.array([d.dt[(mycutmask[j,:])].sum() for j in range(ncam)])
        t_mycut_late_post[i_sim,:] = np.array([d.dt[(d.bhprojsep[j,:]<latesep)&(mycutmask[j,:])].sum() 
                                             for j in range(ncam)])

        t_1or2agn_mycut_l44_tot[i_sim,:] = np.array([d.dt[((d.lagn1>=1.0e44)|(d.lagn2>=1.0e44))&(mycutmask[j,:])].sum() 
                                                     for j in range(ncam)])
        t_1or2agn_mycut_l44_late_post[i_sim,:] = np.array([d.dt[(d.bhprojsep[j,:]<latesep)&((d.lagn1>=1.0e44)|(d.lagn2>=1.0e44))&
                                                              (mycutmask[j,:])].sum() for j in range(ncam)])
        t_1or2agn_mycut_l45_tot[i_sim,:] = np.array([d.dt[((d.lagn1>=1.0e45)|(d.lagn2>=1.0e45))&(mycutmask[j,:])].sum() 
                                                     for j in range(ncam)])
        t_1or2agn_mycut_l45_late_post[i_sim,:] = np.array([d.dt[(d.bhprojsep[j,:]<latesep)&((d.lagn1>=1.0e45)|(d.lagn2>=1.0e45))&
                                                              (mycutmask[j,:])].sum() for j in range(ncam)])

        t_dualagn_mycut_l44_tot[i_sim,:] = np.array([d.dt[(d.lagn1>=1.0e44)&(d.lagn2>=1.0e44)&(mycutmask[j,:])].sum() 
                                                     for j in range(ncam)])
        t_dualagn_mycut_l44_late[i_sim,:] = np.array([d.dt[(d.bhprojsep[j,:]<latesep)&(d.lagn1>=1.0e44)&(d.lagn2>=1.0e44)&
                                                              (mycutmask[j,:])].sum() for j in range(ncam)])
        t_dualagn_mycut_l45_tot[i_sim,:] = np.array([d.dt[(d.lagn1>=1.0e45)&(d.lagn2>=1.0e45)&(mycutmask[j,:])].sum() 
                                                     for j in range(ncam)])
        t_dualagn_mycut_l45_late[i_sim,:] = np.array([d.dt[(d.bhprojsep[j,:]<latesep)&(d.lagn1>=1.0e45)&(d.lagn2>=1.0e45)&
                                                           (mycutmask[j,:])].sum() for j in range(ncam)])

    ### end of loop over path_arr ###

    #print f1or2agn_wise_hires_late_post.shape
    #print f1or2agn_wise_hires_tot.shape

    ###############################
    ### GLOBAL QUANTITIES PLOTS ###
    ###############################

    if nsim > 1 and multisim_plots and not agnx0 and not ngas_only:
        print "Making plots of global quantities..."
        #xarrs = np.array([np.log10(ssfr_max),sfr_max,np.log10(lagn_max),np.log10(ltot_max),flbol_max])
        xarrs = np.array([np.log10(ssfr_max),np.log10(sfr_max),
                          np.log10(lagn_max),np.log10(ltot_max),flbol_max])

        ### new tAGN/ttot plot averaged over all sims:
            #f1or2agn_tot_hires[i_sim,i_mask_hires] = t_1or2agn_tot_hires / ttot
        
        f1or2agn_tot_hires_eachsim = np.array([ t_1or2agn_tot_hires[:,i]/ttot_arr 
                                                for i in range(t_1or2agn_tot_hires.shape[1]) ])
        f1or2agn_late_post_hires_eachsim = np.array([ np.mean(t_1or2agn_late_post_hires[:,i,:]/
                                                              (t_late_arr+t_post_arr),axis=1)
                                                      for i in range(t_1or2agn_late_post_hires.shape[1]) ])
        print f1or2agn_tot_hires_eachsim.shape
        print f1or2agn_late_post_hires_eachsim.shape

        if valtype_allsim == 'median':
            f1or2agn_tot_hires_allsim = np.array([ np.median(t_1or2agn_tot_hires[:,i]/ttot_arr) 
                                                   for i in range(t_1or2agn_tot_hires.shape[1]) ])
            f1or2agn_late_post_hires_allsim = np.array([ np.nanmedian( (t_1or2agn_late_post_hires[:,i,:])/(t_late_arr+t_post_arr) )
                                                         for i in range(t_1or2agn_late_post_hires.shape[1])])
            fdualagn_agn_tot_hires_allsim = np.array([ np.nanmedian(t_dualagn_tot_hires[:,i]/t_1or2agn_tot_hires[:,i]) 
                                                       for i in range(t_1or2agn_tot_hires.shape[1]) ])
            fdualagn_agn_late_post_hires_allsim = np.array([ np.nanmedian( (t_dualagn_late_hires[:,i,:])/
                                                                           (t_1or2agn_late_post_hires[:,i,:]) )
                                                             for i in range(t_1or2agn_late_post_hires.shape[1])])
        else:
            f1or2agn_tot_hires_allsim = np.array([ np.mean(t_1or2agn_tot_hires[:,i]/ttot_arr) 
                                                   for i in range(t_1or2agn_tot_hires.shape[1]) ])
            f1or2agn_late_post_hires_allsim = np.array([ np.sum( (t_1or2agn_late_post_hires[:,i,:])/(t_late_arr+t_post_arr) )/
                                                         (1.0*np.maximum(nproj_tnonzero_late,nproj_tnonzero_post))/nsim
                                                         for i in range(t_1or2agn_late_post_hires.shape[1])])
            fdualagn_agn_tot_hires_allsim = np.array([ np.nanmean(t_dualagn_tot_hires[:,i]/t_1or2agn_tot_hires[:,i]) 
                                                       for i in range(t_1or2agn_tot_hires.shape[1]) ])
            fdualagn_agn_late_post_hires_allsim = np.array([ np.nansum( (t_dualagn_late_hires[:,i,:])/
                                                                        (t_1or2agn_late_post_hires[:,i,:]) )/
                                                             (1.0*np.maximum(nproj_tnonzero_late,nproj_tnonzero_post))/nsim
                                                             for i in range(t_1or2agn_late_post_hires.shape[1]) ])
        yerrs_tot = errbars(np.array([ (t_1or2agn_tot_hires[:,i]/ttot_arr) 
                                       for i in range(t_1or2agn_tot_hires.shape[1]) ]).transpose(),
                            valtype=valtype_allsim,errtype=errtype)
        yerrs_late_post = errbars(np.array([ ((t_1or2agn_late_post_hires[:,i,:])/(t_late_arr+t_post_arr)).flatten()
                                             for i in range(t_1or2agn_late_post_hires.shape[1])]).transpose(), 
                                  valtype=valtype_allsim,errtype=errtype)
        yerrs_dual_tot = errbars(np.array([ (t_dualagn_tot_hires[:,i]/t_1or2agn_tot_hires[:,i]) 
                                            for i in range(t_1or2agn_tot_hires.shape[1]) ]).transpose(),
                                 valtype=valtype_allsim,errtype=errtype)
        yerrs_dual_late_post = errbars(np.array([ ((t_dualagn_late_hires[:,i,:])/(t_1or2agn_late_post_hires[:,i,:])).flatten()
                                                  for i in range(t_1or2agn_late_post_hires.shape[1])]).transpose(), 
                                       valtype=valtype_allsim,errtype=errtype)
        
        print "ttot_arr:",ttot_arr
        print "t_1or2agn_tot_hires:",t_1or2agn_tot_hires
        print f1or2agn_tot_hires_allsim.shape
        print "f1or2agn_tot_hires_allsim:"
        print f1or2agn_tot_hires_allsim
        print f1or2agn_late_post_hires_allsim.shape
        print "f1or2agn_late_post_hires_allsim:"
        print f1or2agn_late_post_hires_allsim
        print len(lagn_lims_hires)
        print "Tot:"
        print lagn_lims_hires
        print t_1or2agn_tot_hires[:,3]
        print lagn_lims
        print t_1or2agn_tot_arr[:,4]

        print "\n"
        print t_1or2agn_late_post_hires[:,0,:]/(t_late_arr+t_post_arr)
        print f1or2agn_tot_hires_allsim.shape
        print yerrs_tot[0].shape
        print yerrs_tot[1].shape
        print f1or2agn_late_post_hires_allsim.shape
        print yerrs_late_post[0].shape
        print yerrs_late_post[1].shape

        mpl.rcParams.update({'font.size': 12})

        #f1or2agn_tot_hires_allsim = np.median(f1or2agn_tot_hires,axis=0)
        #f1or2agn_late_post_hires_allsim = np.median(f1or2agn_late_post_hires,axis=0)
        plt.clf()
        plt.cla()
        plt.close()
        #fig=plt.figure(figsize=(5,5))
        fig=plt.figure(figsize=(5,3.5))
        ax=fig.add_subplot(111)
        ax.set_xscale('log')
        ax.set_xlim(1.0e42,3.0e46)
        ax.set_ylim(0,1)
        #yerrs = errbars(f1or2agn_tot_hires,valtype='mean',errtype='iqr')
        #yerrs = errbars(f1or2agn_tot_hires,valtype='mean',errtype='full')
        #yerrs = errbars(f1or2agn_tot_hires,valtype='median',errtype='mad')
        p1,=ax.plot((lagn_lims_hires),f1or2agn_late_post_hires_allsim,'*',markersize=9,mew=1.6,mec='r',color='none')
        ax.errorbar(lagn_lims_hires,f1or2agn_late_post_hires_allsim,yerr=yerrs_late_post,color='r',
                    ecolor='r',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)
        p0,=ax.plot((lagn_lims_hires),f1or2agn_tot_hires_allsim,'o',markersize=6,mew=1.2,mec='b',color='none')
        ax.errorbar(lagn_lims_hires,f1or2agn_tot_hires_allsim,yerr=yerrs_tot,color='b',ecolor='b',ls='None',capsize=4)
        ax.set_xlabel(r'min L$_{\rm AGN}$ [erg s$^{-1}$]')
        ax.set_ylabel(r't$_{\rm AGN}$ / t$_{\rm tot}$')
        ax.legend((p0,p1), ('Total','Late+Post'), fontsize=10, loc='upper right', 
                  numpoints=1, handletextpad=0.1)
        fig.subplots_adjust(bottom=0.16,left=0.14,right=0.92,top=0.9)
        fig.savefig(outdir+'/agn_fraction_allsim_%s%s_%s_%s.pdf'%(name,zstr,valtype_allsim,errtype))


        plt.clf()
        plt.cla()
        plt.close()
        #fig=plt.figure(figsize=(5,5))
        fig=plt.figure(figsize=(5,3.5))
        ax=fig.add_subplot(111)
        ax.set_xscale('log')
        ax.set_xlim(1.0e42,3.0e46)
        ax.set_ylim(0,1)
        #yerrs = errbars(f1or2agn_tot_hires,valtype='mean',errtype='iqr')
        #yerrs = errbars(f1or2agn_tot_hires,valtype='mean',errtype='full')
        #yerrs = errbars(f1or2agn_tot_hires,valtype='median',errtype='mad')
        p1,=ax.plot((lagn_lims_hires),fdualagn_agn_late_post_hires_allsim,'*',markersize=9,mew=1.6,mec='r',color='none')
        ax.errorbar(lagn_lims_hires,fdualagn_agn_late_post_hires_allsim,yerr=yerrs_dual_late_post,color='r',
                    ecolor='r',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)
        p0,=ax.plot((lagn_lims_hires),fdualagn_agn_tot_hires_allsim,'o',markersize=6,mew=1.2,mec='b',color='none')
        ax.errorbar(lagn_lims_hires,fdualagn_agn_tot_hires_allsim,yerr=yerrs_dual_tot,color='b',ecolor='b',ls='None',capsize=4)
        ax.set_xlabel(r'min L$_{\rm AGN}$ [erg s$^{-1}$]')
        ax.set_ylabel(r't$_{\rm 2AGN}$ / t$_{\rm AGN}$')
        ax.legend((p0,p1), ('Total','Late+Post'), fontsize=10, loc='upper right', 
                  numpoints=1, handletextpad=0.1)
        fig.subplots_adjust(bottom=0.16,left=0.14,right=0.92,top=0.9)
        fig.savefig(outdir+'/dualagn_fraction_allsim_%s%s_%s_%s.pdf'%(name,zstr,valtype_allsim,errtype))

        #f1or2agn_tot_hires_allsim = np.median(f1or2agn_tot_hires,axis=0)
        #f1or2agn_late_post_hires_allsim = np.median(f1or2agn_late_post_hires,axis=0)
        plt.clf()
        plt.cla()
        plt.close()
        #fig=plt.figure(figsize=(5,5))
        fig=plt.figure(figsize=(5,3.5))
        ax=fig.add_subplot(111)
        ax.set_xscale('log')
        ax.set_xlim(1.0e42,3.0e46)
        ax.set_ylim(0,1)
        #yerrs = errbars(f1or2agn_tot_hires,valtype='mean',errtype='iqr')
        #yerrs = errbars(f1or2agn_tot_hires,valtype='mean',errtype='full')
        #yerrs = errbars(f1or2agn_tot_hires,valtype='median',errtype='mad')
        for i_s in range(nsim):
            if name == 'allres_d1d0':
                lw_lp=0.75*i_s+0.5
                lw_tot=0.75*i_s+1
            else:
                lw_lp=2
                lw_tot=1
            ax.plot((lagn_lims_hires),f1or2agn_late_post_hires_eachsim[:,i_s],
                    '*-',markersize=9,mew=1.6,mec='r',lw=lw_lp,color='r')
            ax.plot((lagn_lims_hires),f1or2agn_tot_hires_eachsim[:,i_s],
                    'o-',markersize=6,mew=1.2,mec='b',lw=lw_tot,color='b')
        ax.set_xlabel(r'min L$_{\rm AGN}$ [erg s$^{-1}$]')
        ax.set_ylabel(r't$_{\rm AGN}$ / t$_{\rm tot}$')
        ax.legend((p0,p1), ('Total','Late+Post'), fontsize=10, loc='upper right', 
                  numpoints=1, handletextpad=0.1)
        fig.subplots_adjust(bottom=0.16,left=0.14,right=0.92,top=0.9)
        fig.savefig(outdir+'/agn_fraction_eachsim_%s%s.pdf'%(name,zstr))


        ###
        print t_1or2agn_tot_arr.shape
        print t_1or2agn_tot_arr[:,ixagn_l44[0]].shape
        print np.tile(t_1or2agn_tot_arr[:,ixagn_l44[0]],(ncam,1)).shape
        print np.tile(t_1or2agn_tot_arr[:,ixagn_l44[0]],(ncam,1)).transpose().shape

        n_nonzero_l44=np.max(np.array([np.nonzero(t_1or2agn_late_post_arr[:,ixagn_l44[0],c])[0].shape 
                                       for c in range(ncam)]))
        n_nonzero_l45=np.max(np.array([np.nonzero(t_1or2agn_late_post_arr[:,ixagn_l45[0],c])[0].shape 
                                       for c in range(ncam)]))

        f1or2agn_wise_hires_l44_totagn_eachsim = np.array([ np.mean(t_1or2agn_wise_hires_l44_tot[:,i,:],axis=1)/
                                                            t_1or2agn_tot_arr[:,ixagn_l44[0]] 
                                                            for i in range(t_1or2agn_wise_hires_l44_tot.shape[1]) ])
        f1or2agn_wise_hires_l44_late_post_agn_eachsim = np.array([ np.mean(t_1or2agn_wise_hires_l44_late_post[:,i,:]/
                                                                           t_1or2agn_late_post_arr[:,ixagn_l44[0],:], axis=1 )
                                                                   for i in range(t_1or2agn_wise_hires_l44_late_post.shape[1]) ])
        print f1or2agn_wise_hires_l44_totagn_eachsim.shape
        print f1or2agn_wise_hires_l44_totagn_eachsim
        print f1or2agn_wise_hires_l44_late_post_agn_eachsim.shape
        print f1or2agn_wise_hires_l44_late_post_agn_eachsim
        #return
        
        f1or2agn_wise_hires_l45_totagn_eachsim = np.array([ np.mean(t_1or2agn_wise_hires_l45_tot[:,i,:],axis=1)/
                                                            t_1or2agn_tot_arr[:,ixagn_l45[0]] 
                                                            for i in range(t_1or2agn_wise_hires_l45_tot.shape[1]) ])
        f1or2agn_wise_hires_l45_late_post_agn_eachsim = np.array([ np.mean(t_1or2agn_wise_hires_l45_late_post[:,i,:]/
                                                                           t_1or2agn_late_post_arr[:,ixagn_l45[0],:], axis=1 )
                                                                   for i in range(t_1or2agn_wise_hires_l45_late_post.shape[1]) ])

        f1or2agn_j11_l44_totagn_eachsim =  np.mean(t_1or2agn_j11_l44_tot,axis=1) / t_1or2agn_tot_arr[:,ixagn_l44[0]]
        f1or2agn_j11_l45_totagn_eachsim =  np.mean(t_1or2agn_j11_l45_tot,axis=1) / t_1or2agn_tot_arr[:,ixagn_l45[0]]

        f1or2agn_j11_l44_late_post_agn_eachsim = np.mean( (t_1or2agn_j11_l44_late_post /
                                                           t_1or2agn_late_post_arr[:,ixagn_l44[0],:]), axis=1)
        f1or2agn_j11_l45_late_post_agn_eachsim = np.mean( (t_1or2agn_j11_l45_late_post /
                                                           t_1or2agn_late_post_arr[:,ixagn_l45[0],:]), axis=1)

        f1or2agn_mycut_l44_totagn_eachsim =  np.mean(t_1or2agn_mycut_l44_tot,axis=1) / t_1or2agn_tot_arr[:,ixagn_l44[0]]
        f1or2agn_mycut_l45_totagn_eachsim =  np.mean(t_1or2agn_mycut_l45_tot,axis=1) / t_1or2agn_tot_arr[:,ixagn_l45[0]]

        f1or2agn_mycut_l44_late_post_agn_eachsim = np.mean( (t_1or2agn_mycut_l44_late_post /
                                                             t_1or2agn_late_post_arr[:,ixagn_l44[0],:]), axis=1)
        f1or2agn_mycut_l45_late_post_agn_eachsim = np.mean( (t_1or2agn_mycut_l45_late_post /
                                                             t_1or2agn_late_post_arr[:,ixagn_l45[0],:]), axis=1)

        if valtype_allsim == 'median':
            ### t_wise / ttot ###
            f1or2agn_wise_hires_l44_tot_allsim = np.array([ np.median( (t_1or2agn_wise_hires_l44_tot[:,i,:])/
                                                                    np.tile(ttot_arr,(ncam,1)).transpose() )
                                                            for i in range(t_1or2agn_wise_hires_l44_tot.shape[1])])
            f1or2agn_wise_hires_l44_late_post_allsim = np.array([ np.nanmedian( (t_1or2agn_wise_hires_l44_late_post[:,i,:])/
                                                                                (t_late_arr+t_post_arr) )
                                                                  for i in range(t_1or2agn_wise_hires_l44_late_post.shape[1])])
            ### t_wise / t_agn ###
            f1or2agn_wise_hires_l44_totagn_allsim = np.array([ np.nanmedian( (t_1or2agn_wise_hires_l44_tot[:,i,:])/
                                                                          np.tile(t_1or2agn_tot_arr[:,ixagn_l44[0]],(ncam,1)).transpose() )
                                                               for i in range(t_1or2agn_wise_hires_l44_tot.shape[1])])
            f1or2agn_wise_hires_l45_totagn_allsim = np.array([ np.nanmedian( (t_1or2agn_wise_hires_l45_tot[:,i,:])/
                                                                             np.tile(t_1or2agn_tot_arr[:,ixagn_l45[0]],(ncam,1)).transpose() )
                                                               for i in range(t_1or2agn_wise_hires_l45_tot.shape[1])])
            f1or2agn_wise_hires_l44_late_post_agn_allsim = np.array([ np.nanmedian( (t_1or2agn_wise_hires_l44_late_post[:,i,:])/
                                                                                 (t_1or2agn_late_post_arr[:,ixagn_l44[0],:]) )
                                                                      for i in range(t_1or2agn_wise_hires_l44_late_post.shape[1])])
            f1or2agn_wise_hires_l45_late_post_agn_allsim = np.array([ np.nanmedian( (t_1or2agn_wise_hires_l45_late_post[:,i,:])/
                                                                                 (t_1or2agn_late_post_arr[:,ixagn_l45[0],:]) )
                                                                      for i in range(t_1or2agn_wise_hires_l45_late_post.shape[1])])

            f1or2agn_j11_l44_totagn_allsim = np.nanmedian( t_1or2agn_j11_l44_tot/np.tile(t_1or2agn_tot_arr[:,ixagn_l44[0]],
                                                                                         (ncam,1)).transpose() )
            f1or2agn_j11_l44_late_post_agn_allsim = np.nanmedian( t_1or2agn_j11_l44_late_post/
                                                                  t_1or2agn_late_post_arr[:,ixagn_l44[0],:] )
            f1or2agn_j11_l45_totagn_allsim = np.nanmedian( t_1or2agn_j11_l45_tot/np.tile(t_1or2agn_tot_arr[:,ixagn_l45[0]],
                                                                                         (ncam,1)).transpose() )
            f1or2agn_j11_l45_late_post_agn_allsim = np.nanmedian( t_1or2agn_j11_l45_late_post/
                                                                  t_1or2agn_late_post_arr[:,ixagn_l45[0],:] )

            f1or2agn_mycut_l44_totagn_allsim = np.nanmedian( t_1or2agn_mycut_l44_tot/np.tile(t_1or2agn_tot_arr[:,ixagn_l44[0]],
                                                                                             (ncam,1)).transpose() )
            f1or2agn_mycut_l44_late_post_agn_allsim = np.nanmedian( t_1or2agn_mycut_l44_late_post/
                                                                    t_1or2agn_late_post_arr[:,ixagn_l44[0],:] )
            f1or2agn_mycut_l45_totagn_allsim = np.nanmedian( t_1or2agn_mycut_l45_tot/np.tile(t_1or2agn_tot_arr[:,ixagn_l45[0]],
                                                                   (ncam,1)).transpose() )
            f1or2agn_mycut_l45_late_post_agn_allsim = np.nanmedian( t_1or2agn_mycut_l45_late_post/
                                                                    t_1or2agn_late_post_arr[:,ixagn_l45[0],:] )

            if (f1or2agn_mycut_l44_late_post_agn_allsim.max() > 1.0 or f1or2agn_mycut_l44_totagn_allsim.max() > 1.0 or 
                f1or2agn_mycut_l45_late_post_agn_allsim.max() > 1.0 or f1or2agn_mycut_l45_totagn_allsim.max() > 1.0):

                print "\n"
                print "l44:"
                print t_1or2agn_mycut_l44_late_post.shape
                print t_1or2agn_mycut_l44_late_post
                print t_1or2agn_late_post_arr.shape
                print t_1or2agn_late_post_arr[:,ixagn_l44[0],:].shape
                print t_1or2agn_late_post_arr[:,ixagn_l44[0],:]
                print t_1or2agn_mycut_l44_late_post/t_1or2agn_late_post_arr[:,ixagn_l44[0],:]
                print f1or2agn_mycut_l44_late_post_agn_allsim
                print "\n"
                print t_1or2agn_mycut_l44_tot.shape
                print t_1or2agn_mycut_l44_tot
                print t_1or2agn_tot_arr.shape
                print t_1or2agn_tot_arr[:,ixagn_l44[0]].shape
                print t_1or2agn_tot_arr[:,ixagn_l44[0]]
                print f1or2agn_mycut_l44_totagn_allsim

                print "\n"

                print "l45:"
                print t_1or2agn_mycut_l45_late_post
                print t_1or2agn_late_post_arr[:,ixagn_l45[0],:]
                print f1or2agn_mycut_l45_late_post_agn_allsim
                print "\n"
                print t_1or2agn_mycut_l45_tot
                print t_1or2agn_tot_arr[:,ixagn_l45[0]]
                print f1or2agn_mycut_l45_totagn_allsim

                print "Error: at least one f1or2agn_mycut value > 1.0! "
                return


            ### t_dual_wise / t_dual ###
            fdualagn_wise_hires_l44_totagn_allsim = np.array([ np.nanmedian( (t_dualagn_wise_hires_l44_tot[:,i,:])/
                                                                             np.tile(t_dualagn_tot_arr[:,ixagn_l44[0]],(ncam,1)).transpose() )
                                                               for i in range(t_dualagn_wise_hires_l44_tot.shape[1])])
            fdualagn_wise_hires_l45_totagn_allsim = np.array([ np.nanmedian( (t_dualagn_wise_hires_l45_tot[:,i,:])/
                                                                             np.tile(t_dualagn_tot_arr[:,ixagn_l45[0]],(ncam,1)).transpose() )
                                                               for i in range(t_dualagn_wise_hires_l45_tot.shape[1])])
            fdualagn_wise_hires_l44_late_post_agn_allsim = np.array([ np.nanmedian( (t_dualagn_wise_hires_l44_late[:,i,:])/
                                                                                    (t_dualagn_late_arr[:,ixagn_l44[0],:]) )
                                                                      for i in range(t_dualagn_wise_hires_l44_late.shape[1])])
            fdualagn_wise_hires_l45_late_post_agn_allsim = np.array([ np.nanmedian( (t_dualagn_wise_hires_l45_late[:,i,:])/
                                                                                    (t_dualagn_late_arr[:,ixagn_l45[0],:]) )
                                                                      for i in range(t_dualagn_wise_hires_l45_late.shape[1])])

            fdualagn_j11_l44_totagn_allsim = np.nanmedian( t_dualagn_j11_l44_tot/
                                                           np.tile(t_dualagn_tot_arr[:,ixagn_l44[0]],(ncam,1)).transpose() )
            fdualagn_j11_l44_late_post_agn_allsim = np.nanmedian( t_dualagn_j11_l44_late/t_dualagn_late_arr[:,ixagn_l44[0],:] )
            fdualagn_j11_l45_totagn_allsim = np.nanmedian( t_dualagn_j11_l45_tot/
                                                           np.tile(t_dualagn_tot_arr[:,ixagn_l45[0]],(ncam,1)).transpose() )
            fdualagn_j11_l45_late_post_agn_allsim = np.nanmedian( t_dualagn_j11_l45_late/t_dualagn_late_arr[:,ixagn_l45[0],:] )

            fdualagn_mycut_l44_totagn_allsim = np.nanmedian( t_dualagn_mycut_l44_tot/
                                                             np.tile(t_dualagn_tot_arr[:,ixagn_l44[0]],(ncam,1)).transpose() )
            fdualagn_mycut_l44_late_post_agn_allsim = np.nanmedian( t_dualagn_mycut_l44_late/t_dualagn_late_arr[:,ixagn_l44[0],:] )
            fdualagn_mycut_l45_totagn_allsim = np.nanmedian( t_dualagn_mycut_l45_tot/
                                                             np.tile(t_dualagn_tot_arr[:,ixagn_l45[0]],(ncam,1)).transpose() )
            fdualagn_mycut_l45_late_post_agn_allsim = np.nanmedian( t_dualagn_mycut_l45_late/t_dualagn_late_arr[:,ixagn_l45[0],:] )

            ### t_dual_wise / t_wise ###
            fdualwagn_to_wise_hires_l44_tot_allsim = np.array([ np.nanmedian( t_dualagn_wise_hires_l44_tot[:,i,:]/
                                                                              t_wise_hires_l44_tot[:,i,:] )
                                                                for i in range(t_dualagn_wise_hires_l44_tot.shape[1]) ])
            fdualwagn_to_wise_hires_l45_tot_allsim = np.array([ np.nanmedian( t_dualagn_wise_hires_l45_tot[:,i,:]/
                                                                              t_wise_hires_l44_tot[:,i,:] )
                                                                for i in range(t_dualagn_wise_hires_l45_tot.shape[1]) ])
            fdualwagn_to_wise_hires_l44_late_post_allsim = np.array([ np.nanmedian(t_dualagn_wise_hires_l44_late[:,i,:]/
                                                                                   t_wise_hires_l44_late_post[:,i,:])
                                                                      for i in range(t_dualagn_wise_hires_l44_late.shape[1]) ])
            fdualwagn_to_wise_hires_l45_late_post_allsim = np.array([ np.nanmedian(t_dualagn_wise_hires_l45_late[:,i,:]/
                                                                                   t_wise_hires_l44_late_post[:,i,:])
                                                                      for i in range(t_dualagn_wise_hires_l45_late.shape[1]) ])

            fdualwagn_to_wise_j11_l44_tot_allsim = np.nanmedian( t_dualagn_j11_l44_tot/t_j11_tot )
            fdualwagn_to_wise_j11_l44_late_post_allsim = np.nanmedian( t_dualagn_j11_l44_late/t_j11_late_post )
            fdualwagn_to_wise_j11_l45_tot_allsim = np.nanmedian( t_dualagn_j11_l45_tot/t_j11_tot )
            fdualwagn_to_wise_j11_l45_late_post_allsim = np.nanmedian( t_dualagn_j11_l45_late/t_j11_late_post )

            fdualwagn_to_wise_mycut_l44_tot_allsim = np.nanmedian( t_dualagn_mycut_l44_tot/t_mycut_tot )
            fdualwagn_to_wise_mycut_l44_late_post_allsim = np.nanmedian( t_dualagn_mycut_l44_late/t_mycut_late_post )
            fdualwagn_to_wise_mycut_l45_tot_allsim = np.nanmedian( t_dualagn_mycut_l45_tot/t_mycut_tot )
            fdualwagn_to_wise_mycut_l45_late_post_allsim = np.nanmedian( t_dualagn_mycut_l45_late/t_mycut_late_post )

        else:
            ### t_wise / ttot ###
            f1or2agn_wise_hires_l44_tot_allsim = np.array([ np.sum( (t_1or2agn_wise_hires_l44_tot[:,i,:])/
                                                                    np.tile(ttot_arr,(ncam,1)).transpose() ) / ncam / nsim
                                                            for i in range(t_1or2agn_wise_hires_l44_tot.shape[1])])
            f1or2agn_wise_hires_l44_late_post_allsim = np.array([ np.sum( (t_1or2agn_wise_hires_l44_late_post[:,i,:])/
                                                                          (t_late_arr+t_post_arr) )/
                                                                  (1.0*np.maximum(nproj_tnonzero_late,nproj_tnonzero_post))/nsim 
                                                                  for i in range(t_1or2agn_wise_hires_l44_late_post.shape[1])])
            ### t_wise / t_agn ###
            f1or2agn_wise_hires_l44_totagn_allsim = np.array([ np.nansum( (t_1or2agn_wise_hires_l44_tot[:,i,:])/
                                                                          np.tile(t_1or2agn_tot_arr[:,ixagn_l44[0]],(ncam,1)).transpose() ) 
                                                               / ncam / np.nonzero(t_1or2agn_tot_arr[:,ixagn_l44[0]])[0].shape
                                                               for i in range(t_1or2agn_wise_hires_l44_tot.shape[1])])
            f1or2agn_wise_hires_l45_totagn_allsim = np.array([ np.nansum( (t_1or2agn_wise_hires_l45_tot[:,i,:])/
                                                                          np.tile(t_1or2agn_tot_arr[:,ixagn_l45[0]],(ncam,1)).transpose() ) 
                                                               / ncam / np.nonzero(t_1or2agn_tot_arr[:,ixagn_l45[0]])[0].shape
                                                               for i in range(t_1or2agn_wise_hires_l45_tot.shape[1])])
            f1or2agn_wise_hires_l44_late_post_agn_allsim = np.array([ np.nansum( (t_1or2agn_wise_hires_l44_late_post[:,i,:])/
                                                                                 (t_1or2agn_late_post_arr[:,ixagn_l44[0],:]) )/
                                                                      (1.0*np.maximum(nproj_tnonzero_late,nproj_tnonzero_post))/ n_nonzero_l44
                                                                      for i in range(t_1or2agn_wise_hires_l44_late_post.shape[1])])
            f1or2agn_wise_hires_l45_late_post_agn_allsim = np.array([ np.nansum( (t_1or2agn_wise_hires_l45_late_post[:,i,:])/
                                                                                 (t_1or2agn_late_post_arr[:,ixagn_l45[0],:]) )/
                                                                      (1.0*np.maximum(nproj_tnonzero_late,nproj_tnonzero_post))/ n_nonzero_l45
                                                                      for i in range(t_1or2agn_wise_hires_l45_late_post.shape[1])])

            f1or2agn_j11_l44_totagn_allsim = ( np.nansum( t_1or2agn_j11_l44_tot/
                                                          np.tile(t_1or2agn_tot_arr[:,ixagn_l44[0]],(ncam,1)).transpose() )
                                               / ncam / np.nonzero(t_1or2agn_tot_arr[:,ixagn_l44[0]])[0].shape )
            f1or2agn_j11_l45_totagn_allsim = ( np.nansum( t_1or2agn_j11_l45_tot/
                                                          np.tile(t_1or2agn_tot_arr[:,ixagn_l45[0]],(ncam,1)).transpose() )
                                               / ncam / np.nonzero(t_1or2agn_tot_arr[:,ixagn_l45[0]])[0].shape )
            f1or2agn_j11_l44_late_post_agn_allsim = ( np.nansum( t_1or2agn_j11_l44_late_post/
                                                                 t_1or2agn_late_post_arr[:,ixagn_l44[0],:] )/
                                                      (1.0*np.maximum(nproj_tnonzero_late,nproj_tnonzero_post))/ n_nonzero_l44 )
            f1or2agn_j11_l45_late_post_agn_allsim = ( np.nansum( t_1or2agn_j11_l45_late_post/
                                                                 t_1or2agn_late_post_arr[:,ixagn_l45[0],:] )/
                                                      (1.0*np.maximum(nproj_tnonzero_late,nproj_tnonzero_post))/ n_nonzero_l45 )

            f1or2agn_mycut_l44_totagn_allsim = ( np.nansum( t_1or2agn_mycut_l44_tot/
                                                          np.tile(t_1or2agn_tot_arr[:,ixagn_l44[0]],(ncam,1)).transpose() )
                                               / ncam / np.nonzero(t_1or2agn_tot_arr[:,ixagn_l44[0]])[0].shape )
            f1or2agn_mycut_l45_totagn_allsim = ( np.nansum( t_1or2agn_mycut_l45_tot/
                                                          np.tile(t_1or2agn_tot_arr[:,ixagn_l45[0]],(ncam,1)).transpose() )
                                               / ncam / np.nonzero(t_1or2agn_tot_arr[:,ixagn_l45[0]])[0].shape )
            f1or2agn_mycut_l44_late_post_agn_allsim = ( np.nansum( t_1or2agn_mycut_l44_late_post/
                                                                 t_1or2agn_late_post_arr[:,ixagn_l44[0],:] )/
                                                      (1.0*np.maximum(nproj_tnonzero_late,nproj_tnonzero_post))/ n_nonzero_l44 )
            f1or2agn_mycut_l45_late_post_agn_allsim = ( np.nansum( t_1or2agn_mycut_l45_late_post/
                                                                 t_1or2agn_late_post_arr[:,ixagn_l45[0],:] )/
                                                      (1.0*np.maximum(nproj_tnonzero_late,nproj_tnonzero_post))/ n_nonzero_l45 )

            if (f1or2agn_mycut_l44_late_post_agn_allsim.max() > 1.0 or f1or2agn_mycut_l44_totagn_allsim.max() > 1.0 or
                f1or2agn_mycut_l45_late_post_agn_allsim.max() > 1.0 or f1or2agn_mycut_l45_totagn_allsim.max() > 1.0):

                print "l44:"
                print t_1or2agn_mycut_l44_late_post
                print t_1or2agn_late_post_arr[:,ixagn_l44[0],:]
                print (1.0*np.maximum(nproj_tnonzero_late,nproj_tnonzero_post))/ n_nonzero_l44 
                print f1or2agn_mycut_l44_late_post_agn_allsim
                print "\n"
                print t_1or2agn_mycut_l44_tot
                print t_1or2agn_tot_arr[:,ixagn_l44[0]]
                np.nonzero(t_1or2agn_tot_arr[:,ixagn_l44[0]])[0].shape
                print f1or2agn_mycut_l44_tot_agn_allsim
                print "\n"

                print t_1or2agn_mycut_l45_late_post
                print t_1or2agn_late_post_arr[:,ixagn_l45[0],:]
                print (1.0*np.maximum(nproj_tnonzero_late,nproj_tnonzero_post))/ n_nonzero_l45 
                print f1or2agn_mycut_l45_late_post_agn_allsim
                print "\n"
                print t_1or2agn_mycut_l45_tot
                print t_1or2agn_tot_arr[:,ixagn_l45[0]]
                np.nonzero(t_1or2agn_tot_arr[:,ixagn_l45[0]])[0].shape
                print f1or2agn_mycut_l45_tot_agn_allsim

                print "Error: at least one f1or2agn_mycut value > 1.0! "
                return

            ### t_dual_wise / t_dual ###
            fdualagn_wise_hires_l44_totagn_allsim = np.array([ np.nansum( (t_dualagn_wise_hires_l44_tot[:,i,:])/
                                                                          np.tile(t_dualagn_tot_arr[:,ixagn_l44[0]],(ncam,1)).transpose() )
                                                               / ncam / np.nonzero(t_dualagn_tot_arr[:,ixagn_l44[0]])[0].shape
                                                               for i in range(t_dualagn_wise_hires_l44_tot.shape[1])])
            fdualagn_wise_hires_l45_totagn_allsim = np.array([ np.nansum( (t_dualagn_wise_hires_l45_tot[:,i,:])/
                                                                          np.tile(t_dualagn_tot_arr[:,ixagn_l45[0]],(ncam,1)).transpose() )
                                                               / ncam / np.nonzero(t_dualagn_tot_arr[:,ixagn_l45[0]])[0].shape
                                                               for i in range(t_dualagn_wise_hires_l45_tot.shape[1])])
            fdualagn_wise_hires_l44_late_post_agn_allsim = np.array([ np.nansum( (t_dualagn_wise_hires_l44_late[:,i,:])/
                                                                                 (t_dualagn_late_arr[:,ixagn_l44[0],:]) )/
                                                                      (1.0*nproj_tnonzero_late)/ n_nonzero_l44
                                                                      for i in range(t_dualagn_wise_hires_l44_late.shape[1])])
            fdualagn_wise_hires_l45_late_post_agn_allsim = np.array([ np.nansum( (t_dualagn_wise_hires_l45_late[:,i,:])/
                                                                                 (t_dualagn_late_arr[:,ixagn_l45[0],:]) )/
                                                                      (1.0*nproj_tnonzero_late)/ n_nonzero_l45
                                                                      for i in range(t_dualagn_wise_hires_l45_late.shape[1])])

            fdualagn_j11_l44_totagn_allsim = (np.nansum( t_dualagn_j11_l44_tot/np.tile(t_dualagn_tot_arr[:,ixagn_l44[0]],(ncam,1)).transpose() )
                                              / ncam / np.nonzero(t_dualagn_tot_arr[:,ixagn_l44[0]])[0].shape )
            fdualagn_j11_l44_late_post_agn_allsim = ( np.nansum( t_dualagn_j11_l44_late/t_dualagn_late_arr[:,ixagn_l44[0],:] )/
                                                      (1.0*nproj_tnonzero_late)/ n_nonzero_l44 )
            fdualagn_j11_l45_totagn_allsim = (np.nansum( t_dualagn_j11_l45_tot/np.tile(t_dualagn_tot_arr[:,ixagn_l45[0]],(ncam,1)).transpose() )
                                              / ncam / np.nonzero(t_dualagn_tot_arr[:,ixagn_l45[0]])[0].shape )
            fdualagn_j11_l45_late_post_agn_allsim = ( np.nansum( t_dualagn_j11_l45_late/t_dualagn_late_arr[:,ixagn_l45[0],:] )/
                                                      (1.0*nproj_tnonzero_late)/ n_nonzero_l45 )

            fdualagn_mycut_l44_totagn_allsim = (np.nansum( t_dualagn_mycut_l44_tot/np.tile(t_dualagn_tot_arr[:,ixagn_l44[0]],(ncam,1)).transpose() )
                                              / ncam / np.nonzero(t_dualagn_tot_arr[:,ixagn_l44[0]])[0].shape )
            fdualagn_mycut_l44_late_post_agn_allsim = ( np.nansum( t_dualagn_mycut_l44_late/t_dualagn_late_arr[:,ixagn_l44[0],:] )/
                                                      (1.0*nproj_tnonzero_late)/ n_nonzero_l44 )
            fdualagn_mycut_l45_totagn_allsim = (np.nansum( t_dualagn_mycut_l45_tot/np.tile(t_dualagn_tot_arr[:,ixagn_l45[0]],(ncam,1)).transpose() )
                                              / ncam / np.nonzero(t_dualagn_tot_arr[:,ixagn_l45[0]])[0].shape )
            fdualagn_mycut_l45_late_post_agn_allsim = ( np.nansum( t_dualagn_mycut_l45_late/t_dualagn_late_arr[:,ixagn_l45[0],:] )/
                                                      (1.0*nproj_tnonzero_late)/ n_nonzero_l45 )

            ### t_dual_wise / t_wise ###
            fdualwagn_to_wise_hires_l44_tot_allsim = np.array([ np.nanmean( t_dualagn_wise_hires_l44_tot[:,i,:]/
                                                                            t_wise_hires_l44_tot[:,i,:] )
                                                                for i in range(t_dualagn_wise_hires_l44_tot.shape[1]) ])
            fdualwagn_to_wise_hires_l45_tot_allsim = np.array([ np.nanmean( t_dualagn_wise_hires_l45_tot[:,i,:]/
                                                                            t_wise_hires_l44_tot[:,i,:] )
                                                                for i in range(t_dualagn_wise_hires_l45_tot.shape[1]) ])
            fdualwagn_to_wise_hires_l44_late_post_allsim = np.array([ np.nanmean(t_dualagn_wise_hires_l44_late[:,i,:]/
                                                                                 t_wise_hires_l44_late_post[:,i,:])
                                                                      for i in range(t_dualagn_wise_hires_l44_late.shape[1]) ])
            fdualwagn_to_wise_hires_l45_late_post_allsim = np.array([ np.nanmean(t_dualagn_wise_hires_l45_late[:,i,:]/
                                                                                 t_wise_hires_l44_late_post[:,i,:])
                                                                      for i in range(t_dualagn_wise_hires_l45_late.shape[1]) ])

            fdualwagn_to_wise_j11_l44_tot_allsim = np.nanmean( t_dualagn_j11_l44_tot/t_j11_tot )
            fdualwagn_to_wise_j11_l44_late_post_allsim = np.nanmean( t_dualagn_j11_l44_late/t_j11_late_post )
            fdualwagn_to_wise_j11_l45_tot_allsim = np.nanmean( t_dualagn_j11_l45_tot/t_j11_tot )
            fdualwagn_to_wise_j11_l45_late_post_allsim = np.nanmean( t_dualagn_j11_l45_late/t_j11_late_post )

            fdualwagn_to_wise_mycut_l44_tot_allsim = np.nanmean( t_dualagn_mycut_l44_tot/t_mycut_tot )
            fdualwagn_to_wise_mycut_l44_late_post_allsim = np.nanmean( t_dualagn_mycut_l44_late/t_mycut_late_post )
            fdualwagn_to_wise_mycut_l45_tot_allsim = np.nanmean( t_dualagn_mycut_l45_tot/t_mycut_tot )
            fdualwagn_to_wise_mycut_l45_late_post_allsim = np.nanmean( t_dualagn_mycut_l45_late/t_mycut_late_post )

            
        print "\nt_1or2agn_late_post_arr:"
        print t_1or2agn_late_post_arr.shape
        print t_1or2agn_late_post_arr[:,ixagn_l44[0],:].shape
        print t_1or2agn_late_post_arr[:,ixagn_l45[0],:].shape
        print "\n"

        yerrs_wise_l44_tot = errbars(np.array([ ( (t_1or2agn_wise_hires_l44_tot[:,i,:])/
                                                  np.tile(ttot_arr,(ncam,1)).transpose() ).flatten()
                                                for i in range(t_1or2agn_wise_hires_l44_tot.shape[1])]).transpose(),
                                     valtype=valtype_allsim,errtype=errtype)
        yerrs_wise_l44_late_post = errbars(np.array([ ( (t_1or2agn_wise_hires_l44_late_post[:,i,:])/
                                                        (t_late_arr+t_post_arr) ).flatten()
                                                      for i in range(t_1or2agn_wise_hires_l44_late_post.shape[1])]).transpose(),
                                           valtype=valtype_allsim,errtype=errtype)
        yerrs_wise_l44_totagn = errbars(np.array([ ( (t_1or2agn_wise_hires_l44_tot[:,i,:])/
                                                  np.tile(t_1or2agn_tot_arr[:,ixagn_l44[0]],(ncam,1)).transpose() ).flatten()
                                                for i in range(t_1or2agn_wise_hires_l44_tot.shape[1])]).transpose(),
                                       valtype=valtype_allsim,errtype=errtype)
        yerrs_wise_l45_totagn = errbars(np.array([ ( (t_1or2agn_wise_hires_l45_tot[:,i,:])/
                                                  np.tile(t_1or2agn_tot_arr[:,ixagn_l45[0]],(ncam,1)).transpose() ).flatten()
                                                for i in range(t_1or2agn_wise_hires_l45_tot.shape[1])]).transpose(),
                                       valtype=valtype_allsim,errtype=errtype)
        yerrs_wise_l44_late_post_agn = errbars(np.array([ ( (t_1or2agn_wise_hires_l44_late_post[:,i,:])/
                                                        (t_1or2agn_late_post_arr[:,ixagn_l44[0],:]) ).flatten()
                                                      for i in range(t_1or2agn_wise_hires_l44_late_post.shape[1])]).transpose(),
                                       valtype=valtype_allsim,errtype=errtype)
        yerrs_wise_l45_late_post_agn = errbars(np.array([ ( (t_1or2agn_wise_hires_l45_late_post[:,i,:])/
                                                        (t_1or2agn_late_post_arr[:,ixagn_l45[0],:]) ).flatten()
                                                      for i in range(t_1or2agn_wise_hires_l45_late_post.shape[1])]).transpose(),
                                       valtype=valtype_allsim,errtype=errtype)
        

        yerrs_dual_wise_l44_totagn = errbars(np.array([ ( (t_dualagn_wise_hires_l44_tot[:,i,:])/
                                                          np.tile(t_dualagn_tot_arr[:,ixagn_l44[0]],(ncam,1)).transpose() ).flatten()
                                                        for i in range(t_dualagn_wise_hires_l44_tot.shape[1])]).transpose(),
                                             valtype=valtype_allsim,errtype=errtype)
        yerrs_dual_wise_l45_totagn = errbars(np.array([ ( (t_dualagn_wise_hires_l45_tot[:,i,:])/
                                                          np.tile(t_dualagn_tot_arr[:,ixagn_l45[0]],(ncam,1)).transpose() ).flatten()
                                                        for i in range(t_dualagn_wise_hires_l45_tot.shape[1])]).transpose(),
                                             valtype=valtype_allsim,errtype=errtype)
        yerrs_dual_wise_l44_late_post_agn = errbars(np.array([ ( (t_dualagn_wise_hires_l44_late[:,i,:])/
                                                                 (t_dualagn_late_arr[:,ixagn_l44[0],:]) ).flatten()
                                                               for i in range(t_dualagn_wise_hires_l44_late.shape[1])]).transpose(),
                                                    valtype=valtype_allsim,errtype=errtype)
        yerrs_dual_wise_l45_late_post_agn = errbars(np.array([ ( (t_dualagn_wise_hires_l45_late[:,i,:])/
                                                                 (t_dualagn_late_arr[:,ixagn_l45[0],:]) ).flatten()
                                                               for i in range(t_dualagn_wise_hires_l45_late.shape[1])]).transpose(),
                                                    valtype=valtype_allsim,errtype=errtype)


        yerrs_dualwagn_to_wise_l44_tot = errbars(np.array([ ( t_dualagn_wise_hires_l44_tot[:,i,:]/
                                                              t_wise_hires_l44_tot[:,i,:] ).flatten()
                                                            for i in range(t_dualagn_wise_hires_l44_tot.shape[1])]).transpose(),
                                                 valtype=valtype_allsim,errtype=errtype)
        yerrs_dualwagn_to_wise_l45_tot = errbars(np.array([ ( t_dualagn_wise_hires_l45_tot[:,i,:]/
                                                              t_wise_hires_l45_tot[:,i,:] ).flatten()
                                                            for i in range(t_dualagn_wise_hires_l45_tot.shape[1])]).transpose(),
                                                 valtype=valtype_allsim,errtype=errtype)
        yerrs_dualwagn_to_wise_l44_late_post = errbars(np.array([ ( t_dualagn_wise_hires_l44_late[:,i,:]/
                                                                    t_wise_hires_l44_late_post[:,i,:] ).flatten()
                                                                  for i in range(t_dualagn_wise_hires_l44_late.shape[1])]).transpose(),
                                                       valtype=valtype_allsim,errtype=errtype)

        #tmp=(t_dualagn_wise_hires_l44_late[6:,:]/t_wise_hires_l44_late_post[6:,:]).flatten()
        #print tmp[tmp.argsort()]
        #print t_dualagn_wise_hires_l44_tot[6:,:]
        #print tmp.shape
        #print tmp[tmp==tmp].shape
        #print tmp[tmp==0].shape
        #print np.array([np.nanmedian((t_dualagn_wise_hires_l44_late[:,i,:]/t_wise_hires_l44_late_post[:,i,:] ))
        #                for i in range(t_dualagn_wise_hires_l44_late.shape[1])])
        #print np.array([np.nanmin((t_dualagn_wise_hires_l44_late[:,i,:]/t_wise_hires_l44_late_post[:,i,:] ))
        #                for i in range(t_dualagn_wise_hires_l44_late.shape[1])])
        #print np.array([np.nanmax((t_dualagn_wise_hires_l44_late[:,i,:]/t_wise_hires_l44_late_post[:,i,:] ))
        #                for i in range(t_dualagn_wise_hires_l44_late.shape[1])])
        #print fdualwagn_to_wise_hires_l44_late_post_allsim
        #print yerrs_dualwagn_to_wise_l44_late_post
        yerrs_dualwagn_to_wise_l45_late_post = errbars(np.array([ ( t_dualagn_wise_hires_l45_late[:,i,:]/
                                                                    t_wise_hires_l45_late_post[:,i,:] ).flatten()
                                                                  for i in range(t_dualagn_wise_hires_l45_late.shape[1])]).transpose(),
                                                       valtype=valtype_allsim,errtype=errtype)


        yerrs_j11_l44_totagn = errbars((t_1or2agn_j11_l44_tot/
                                        np.tile(t_1or2agn_tot_arr[:,ixagn_l44[0]],(ncam,1)).transpose()).flatten(),
                                        valtype=valtype_allsim,errtype=errtype)
        yerrs_j11_l44_late_post_agn = errbars((t_1or2agn_j11_l44_late_post/
                                               t_1or2agn_late_post_arr[:,ixagn_l44[0],:]).flatten(),
                                               valtype=valtype_allsim,errtype=errtype)        
        yerrs_j11_l45_totagn = errbars((t_1or2agn_j11_l45_tot/
                                        np.tile(t_1or2agn_tot_arr[:,ixagn_l45[0]],(ncam,1)).transpose()).flatten(),
                                        valtype=valtype_allsim,errtype=errtype)
        yerrs_j11_l45_late_post_agn = errbars((t_1or2agn_j11_l45_late_post/
                                               t_1or2agn_late_post_arr[:,ixagn_l45[0],:]).flatten(),
                                              valtype=valtype_allsim,errtype=errtype)

        yerrs_dual_j11_l44_totagn = errbars((t_dualagn_j11_l44_tot/
                                             np.tile(t_dualagn_tot_arr[:,ixagn_l44[0]],(ncam,1)).transpose()).flatten(),
                                            valtype=valtype_allsim,errtype=errtype)
        yerrs_dual_j11_l44_late_post_agn = errbars((t_dualagn_j11_l44_late/
                                                    t_dualagn_late_arr[:,ixagn_l44[0],:]).flatten(),
                                                   valtype=valtype_allsim,errtype=errtype)        
        yerrs_dual_j11_l45_totagn = errbars((t_dualagn_j11_l45_tot/
                                             np.tile(t_dualagn_tot_arr[:,ixagn_l45[0]],(ncam,1)).transpose()).flatten(),
                                            valtype=valtype_allsim,errtype=errtype)
        yerrs_dual_j11_l45_late_post_agn = errbars((t_dualagn_j11_l45_late/
                                                    t_dualagn_late_arr[:,ixagn_l45[0],:]).flatten(),
                                                   valtype=valtype_allsim,errtype=errtype)        

        yerrs_dualwagn_to_wise_j11_l44_tot = errbars((t_dualagn_j11_l44_tot/t_j11_tot).flatten(),
                                                     valtype=valtype_allsim,errtype=errtype)
        yerrs_dualwagn_to_wise_j11_l45_tot = errbars((t_dualagn_j11_l45_tot/t_j11_tot).flatten(),
                                                     valtype=valtype_allsim,errtype=errtype)
        yerrs_dualwagn_to_wise_j11_l44_late_post = errbars((t_dualagn_j11_l44_late/t_j11_late_post).flatten(),
                                                           valtype=valtype_allsim,errtype=errtype)
        yerrs_dualwagn_to_wise_j11_l45_late_post = errbars((t_dualagn_j11_l45_late/t_j11_late_post).flatten(),
                                                           valtype=valtype_allsim,errtype=errtype)



        yerrs_mycut_l44_totagn = errbars((t_1or2agn_mycut_l44_tot/
                                          np.tile(t_1or2agn_tot_arr[:,ixagn_l44[0]],(ncam,1)).transpose()).flatten(),
                                         valtype=valtype_allsim,errtype=errtype)
        yerrs_mycut_l44_late_post_agn = errbars((t_1or2agn_mycut_l44_late_post/
                                                 t_1or2agn_late_post_arr[:,ixagn_l44[0],:]).flatten(),
                                                valtype=valtype_allsim,errtype=errtype)        
        yerrs_mycut_l45_totagn = errbars((t_1or2agn_mycut_l45_tot/
                                          np.tile(t_1or2agn_tot_arr[:,ixagn_l45[0]],(ncam,1)).transpose()).flatten(),
                                         valtype=valtype_allsim,errtype=errtype)
        yerrs_mycut_l45_late_post_agn = errbars((t_1or2agn_mycut_l45_late_post/
                                                 t_1or2agn_late_post_arr[:,ixagn_l45[0],:]).flatten(),
                                                valtype=valtype_allsim,errtype=errtype)

        yerrs_dual_mycut_l44_totagn = errbars((t_dualagn_mycut_l44_tot/
                                               np.tile(t_dualagn_tot_arr[:,ixagn_l44[0]],(ncam,1)).transpose()).flatten(),
                                              valtype=valtype_allsim,errtype=errtype)
        yerrs_dual_mycut_l44_late_post_agn = errbars((t_dualagn_mycut_l44_late/
                                                      t_dualagn_late_arr[:,ixagn_l44[0],:]).flatten(),
                                                   valtype=valtype_allsim,errtype=errtype)        
        yerrs_dual_mycut_l45_totagn = errbars((t_dualagn_mycut_l45_tot/
                                               np.tile(t_dualagn_tot_arr[:,ixagn_l45[0]],(ncam,1)).transpose()).flatten(),
                                              valtype=valtype_allsim,errtype=errtype)
        yerrs_dual_mycut_l45_late_post_agn = errbars((t_dualagn_mycut_l45_late/
                                                      t_dualagn_late_arr[:,ixagn_l45[0],:]).flatten(),
                                                     valtype=valtype_allsim,errtype=errtype)        

        yerrs_dualwagn_to_wise_mycut_l44_tot = errbars((t_dualagn_mycut_l44_tot/t_mycut_tot).flatten(),
                                                     valtype=valtype_allsim,errtype=errtype)
        yerrs_dualwagn_to_wise_mycut_l45_tot = errbars((t_dualagn_mycut_l45_tot/t_mycut_tot).flatten(),
                                                     valtype=valtype_allsim,errtype=errtype)
        yerrs_dualwagn_to_wise_mycut_l44_late_post = errbars((t_dualagn_mycut_l44_late/t_mycut_late_post).flatten(),
                                                           valtype=valtype_allsim,errtype=errtype)
        yerrs_dualwagn_to_wise_mycut_l45_late_post = errbars((t_dualagn_mycut_l45_late/t_mycut_late_post).flatten(),
                                                           valtype=valtype_allsim,errtype=errtype)

        #f1or2agn_late_post = f1or2agn_late+f1agn_post

        xlim=(0.28,1.06)

        plt.clf()
        plt.cla()
        plt.close()
        #fig=plt.figure(figsize=(5,5))
        fig=plt.figure(figsize=(5,3.5))
        ax=fig.add_subplot(111)
        ax.set_xlim(xlim)
        ax.set_ylim(0,1)
        ax.set_xlabel(r'min $W1 - W2$')
        ax.set_ylabel(r't$_{\rm WISEAGN}$ / t$_{\rm tot}$')
        #ax.plot(w12_lims_hires,np.mean(f1or2agn_wise_hires_late_post/f1or2agn_late_post[:,ixagn_l44],axis=0),'*',
        #        markersize=9,mew=1.3,mec='r',color='None')
        p1,=ax.plot(w12_lims_hires,f1or2agn_wise_hires_l44_late_post_allsim,'*',
                    markersize=9,mew=1.6,mec='r',color='None')
        ax.errorbar(w12_lims_hires,f1or2agn_wise_hires_l44_late_post_allsim,yerr=yerrs_wise_l44_late_post,color='r',
                    ecolor='r',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)
        p0,=ax.plot(w12_lims_hires,f1or2agn_wise_hires_l44_tot_allsim,'o',
                markersize=6,mew=1.2,mec='b',color='None')
        ax.errorbar(w12_lims_hires,f1or2agn_wise_hires_l44_tot_allsim,yerr=yerrs_wise_l44_tot,color='b',
                    ecolor='b',ls='None',capsize=4)
        ax.legend((p0,p1), ('Total','Late+Post'), fontsize=10, loc='upper right', 
                  numpoints=1, handletextpad=0.1)
        fig.subplots_adjust(bottom=0.16,left=0.14,right=0.92,top=0.9)
        fig.savefig(outdir+'/wiseagn_totfraction_allsim_l44_%s%s_%s_%s.pdf'%(name,zstr,valtype_allsim,errtype))

        #print yerrs_wise_l44_late_post_agn[0]
        #print yerrs_wise_l44_late_post_agn[1]
        #print yerrs_wise_l44_late_post_agn[0][0]
        #print yerrs_wise_l44_late_post_agn[1][3]
        #print f1or2agn_wise_hires_l44_totagn_allsim.shape
        #print f1or2agn_wise_hires_l44_totagn_allsim
        #print "entering loop"

        with open('%s/completeness_%sz%g_%s.txt'%(outdir,wflim_str,z,name),'w') as fc:
            for jj,wlimh in enumerate(w12_lims_hires):

                #print wlimh
                #print f1or2agn_wise_hires_l44_late_post_agn_allsim[jj]
                #print yerrs_wise_l44_late_post_agn[0][jj]
                #print yerrs_wise_l44_late_post_agn[1][jj]
                #print f1or2agn_wise_hires_l44_totagn_allsim[jj]
                #print yerrs_wise_l44_totagn[0][jj]
                #print yerrs_wise_l44_totagn[1][jj]


                fc.write(((13*"%g ")+"\n")%(wlimh,f1or2agn_wise_hires_l44_late_post_agn_allsim[jj],
                                            yerrs_wise_l44_late_post_agn[0][jj],yerrs_wise_l44_late_post_agn[1][jj],
                                            f1or2agn_wise_hires_l44_totagn_allsim[jj],
                                            yerrs_wise_l44_totagn[0][jj],yerrs_wise_l44_totagn[1][jj],
                                            f1or2agn_wise_hires_l45_late_post_agn_allsim[jj],
                                            yerrs_wise_l45_late_post_agn[0][jj],yerrs_wise_l45_late_post_agn[1][jj],
                                            f1or2agn_wise_hires_l45_totagn_allsim[jj],
                                            yerrs_wise_l45_totagn[0][jj],yerrs_wise_l45_totagn[1][jj]))

            fc.write(((13*"%g ")+"\n")%(-11,f1or2agn_j11_l44_late_post_agn_allsim,
                                         yerrs_j11_l44_late_post_agn[0],yerrs_j11_l44_late_post_agn[1],
                                         f1or2agn_j11_l44_totagn_allsim,
                                         yerrs_j11_l44_totagn[0],yerrs_j11_l44_totagn[1],
                                         f1or2agn_j11_l45_late_post_agn_allsim,
                                         yerrs_j11_l45_late_post_agn[0],yerrs_j11_l45_late_post_agn[1],
                                         f1or2agn_j11_l45_totagn_allsim,
                                         yerrs_j11_l45_totagn[0],yerrs_j11_l45_totagn[1]))
            fc.write(((13*"%g ")+"\n")%(-12,f1or2agn_mycut_l44_late_post_agn_allsim,
                                         yerrs_mycut_l44_late_post_agn[0],yerrs_mycut_l44_late_post_agn[1],
                                         f1or2agn_mycut_l44_totagn_allsim,
                                         yerrs_mycut_l44_totagn[0],yerrs_mycut_l44_totagn[1],
                                         f1or2agn_mycut_l45_late_post_agn_allsim,
                                         yerrs_mycut_l45_late_post_agn[0],yerrs_mycut_l45_late_post_agn[1],
                                         f1or2agn_mycut_l45_totagn_allsim,
                                         yerrs_mycut_l45_totagn[0],yerrs_mycut_l45_totagn[1]))
            fc.close()
            #f1or2agn_wise_hires_l45_late_post_agn_allsim
            #yerrs_wise_l45_late_post_agn
            #f1or2agn_wise_hires_l45_totagn_allsim
            #yerrs_wise_l45_totagn

            #f1or2agn_j11_l44_late_post_agn_allsim
            #yerrs_j11_l44_late_post_agn
            #f1or2agn_mycut_l44_late_post_agn_allsim
            #yerrs_mycut_l44_late_post_agn

            #f1or2agn_j11_l45_late_post_agn_allsim
            #yerrs_j11_l45_late_post_agn
            #f1or2agn_mycut_l45_late_post_agn_allsim
            #yerrs_mycut_l45_late_post_agn


        #return

        if plot_wedge_allsim:


            ################ LAGN > 1e44, ALLSIM ################
            plt.clf()
            plt.cla()
            plt.close()
            fig=plt.figure(figsize=(5,3.5))
            ax=plt.subplot2grid((1,5),(0,0),colspan=4)
            ax.set_xlim(xlim)
            ax.set_ylim(0,1)
            ax.set_xlabel(r'min $W1 - W2$')
            ax.set_ylabel(r't$_{\rm WISEAGN}$ / t$_{\rm AGN}$')
            p1,=ax.plot(w12_lims_hires,f1or2agn_wise_hires_l44_late_post_agn_allsim,'*',
                        markersize=9,mew=1.6,mec='r',color='None')
            ax.errorbar(w12_lims_hires,f1or2agn_wise_hires_l44_late_post_agn_allsim,yerr=yerrs_wise_l44_late_post_agn,
                        ecolor='r',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)
            
            p0,=ax.plot(w12_lims_hires,f1or2agn_wise_hires_l44_totagn_allsim,'o',
                        markersize=6,mew=1.2,mec='b',color='None')
            ax.errorbar(w12_lims_hires,f1or2agn_wise_hires_l44_totagn_allsim,yerr=yerrs_wise_l44_totagn,color='b',
                        ecolor='b',ls='None',capsize=4)
            ax.legend((p0,p1), ('Total','Late+Post'), fontsize=10, loc='lower left', 
                      numpoints=1, handletextpad=0.1)
            
            x1=0.25
            x2=0.75
            ax2=plt.subplot2grid((1,5),(0,4))
            ax2.set_xlim(0,1)
            ax2.set_ylim(0,1)
            ax2.set_xticks(np.array([x1,x2]))
            ax2.set_xticklabels(np.array(['J11','This\nwork']),fontsize=10)
            ax2.set_yticklabels(np.repeat('',10))
            
            ax2.plot(x1,f1or2agn_j11_l44_late_post_agn_allsim,"*",
                     markersize=9,mew=1.6,mec='darkorange',color='None')
            ax2.errorbar(x1,f1or2agn_j11_l44_late_post_agn_allsim,yerr=yerrs_j11_l44_late_post_agn,
                         ecolor='darkorange',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)
            ax2.plot(x2,f1or2agn_mycut_l44_late_post_agn_allsim,"*",
                     markersize=9,mew=1.6,mec='darkorange',color='None')
            ax2.errorbar(x2,f1or2agn_mycut_l44_late_post_agn_allsim,yerr=yerrs_mycut_l44_late_post_agn,
                         ecolor='darkorange',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)
            
            ax2.plot(x1,f1or2agn_j11_l44_totagn_allsim,"o",ms=6,mew=1.2,mec='c',color='None')
            ax2.errorbar(x1,f1or2agn_j11_l44_totagn_allsim,yerr=yerrs_j11_l44_totagn,color='c',
                         ecolor='c',ls='None',capsize=4)
            ax2.plot(x2,f1or2agn_mycut_l44_totagn_allsim,"o",ms=6,mew=1.2,mec='c',color='None')
            ax2.errorbar(x2,f1or2agn_mycut_l44_totagn_allsim,yerr=yerrs_mycut_l44_totagn,color='c',
                         ecolor='c',ls='None',capsize=4)
            
            fig.suptitle(r'min L$_{\rm AGN} = 10^{44}$ erg s$^{-1}$')
            fig.subplots_adjust(bottom=0.16,left=0.14,right=0.92,top=0.9,wspace=0,hspace=0)
            wstr='_alt' if alt_wedge else ''
            fig.savefig(outdir+'/wiseagn_fraction_allsim_l44_%s%s_%s_%s%s.pdf'%(name,zstr,valtype_allsim,
                                                                                errtype,wstr))

            ################ DUAL LAGN > 1e44, ALLSIM ################
            plt.clf()
            plt.cla()
            plt.close()
            fig=plt.figure(figsize=(5,3.5))
            ax=plt.subplot2grid((1,5),(0,0),colspan=4)
            ax.set_xlim(xlim)
            ax.set_ylim(0,1)
            ax.set_xlabel(r'min $W1 - W2$')
            ax.set_ylabel(r't$_{\rm WISE,Dual}$ / t$_{\rm Dual}$')
            p1,=ax.plot(w12_lims_hires,fdualagn_wise_hires_l44_late_post_agn_allsim,'*',
                        markersize=9,mew=1.6,mec='r',color='None')
            ax.errorbar(w12_lims_hires,fdualagn_wise_hires_l44_late_post_agn_allsim,yerr=yerrs_dual_wise_l44_late_post_agn,
                        ecolor='r',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)
            
            p0,=ax.plot(w12_lims_hires,fdualagn_wise_hires_l44_totagn_allsim,'o',
                        markersize=6,mew=1.2,mec='b',color='None')
            ax.errorbar(w12_lims_hires,fdualagn_wise_hires_l44_totagn_allsim,yerr=yerrs_dual_wise_l44_totagn,color='b',
                        ecolor='b',ls='None',capsize=4)
            ax.legend((p0,p1), ('Total','Late'), fontsize=10, loc='lower left', 
                      numpoints=1, handletextpad=0.1)
            
            x1=0.25
            x2=0.75
            ax2=plt.subplot2grid((1,5),(0,4))
            ax2.set_xlim(0,1)
            ax2.set_ylim(0,1)
            ax2.set_xticks(np.array([x1,x2]))
            ax2.set_xticklabels(np.array(['J11','This\nwork']),fontsize=10)
            ax2.set_yticklabels(np.repeat('',10))
            
            ax2.plot(x1,fdualagn_j11_l44_late_post_agn_allsim,"*",
                     markersize=9,mew=1.6,mec='darkorange',color='None')
            ax2.errorbar(x1,fdualagn_j11_l44_late_post_agn_allsim,yerr=yerrs_dual_j11_l44_late_post_agn,
                         ecolor='darkorange',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)
            ax2.plot(x2,fdualagn_mycut_l44_late_post_agn_allsim,"*",
                     markersize=9,mew=1.6,mec='darkorange',color='None')
            ax2.errorbar(x2,fdualagn_mycut_l44_late_post_agn_allsim,yerr=yerrs_dual_mycut_l44_late_post_agn,
                         ecolor='darkorange',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)
            
            ax2.plot(x1,fdualagn_j11_l44_totagn_allsim,"o",ms=6,mew=1.2,mec='c',color='None')
            ax2.errorbar(x1,fdualagn_j11_l44_totagn_allsim,yerr=yerrs_dual_j11_l44_totagn,color='c',
                         ecolor='c',ls='None',capsize=4)
            ax2.plot(x2,fdualagn_mycut_l44_totagn_allsim,"o",ms=6,mew=1.2,mec='c',color='None')
            ax2.errorbar(x2,fdualagn_mycut_l44_totagn_allsim,yerr=yerrs_dual_mycut_l44_totagn,color='c',
                         ecolor='c',ls='None',capsize=4)
            
            fig.suptitle(r'min L$_{\rm AGN} = 10^{44}$ erg s$^{-1}$')
            fig.subplots_adjust(bottom=0.16,left=0.14,right=0.92,top=0.9,wspace=0,hspace=0)
            wstr='_alt' if alt_wedge else ''
            fig.savefig(outdir+'/wise_dual_agn_fraction_allsim_l44_%s%s_%s_%s%s.pdf'%(name,zstr,valtype_allsim,
                                                                                      errtype,wstr))

            ################ DUAL / WISE  LAGN > 1e44, ALLSIM ################
            plt.clf()
            plt.cla()
            plt.close()
            fig=plt.figure(figsize=(5,3.5))
            ax=plt.subplot2grid((1,5),(0,0),colspan=4)
            ax.set_xlim(xlim)
            ax.set_ylim(0,1)
            ax.set_xlabel(r'min $W1 - W2$')
            ax.set_ylabel(r't$_{\rm WISE,Dual}$ / t$_{\rm WISE}$')
            p1,=ax.plot(w12_lims_hires,fdualwagn_to_wise_hires_l44_late_post_allsim,'*',
                        markersize=9,mew=1.6,mec='r',color='None')
            ax.errorbar(w12_lims_hires,fdualwagn_to_wise_hires_l44_late_post_allsim,yerr=yerrs_dualwagn_to_wise_l44_late_post,
                        ecolor='r',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)
            
            p0,=ax.plot(w12_lims_hires,fdualwagn_to_wise_hires_l44_tot_allsim,'o',
                        markersize=6,mew=1.2,mec='b',color='None')
            ax.errorbar(w12_lims_hires,fdualwagn_to_wise_hires_l44_tot_allsim,yerr=yerrs_dualwagn_to_wise_l44_tot,color='b',
                        ecolor='b',ls='None',capsize=4)
            ax.legend((p0,p1), ('Total','Late+Post'), fontsize=10, loc='lower left', 
                      numpoints=1, handletextpad=0.1)
            
            x1=0.25
            x2=0.75
            ax2=plt.subplot2grid((1,5),(0,4))
            ax2.set_xlim(0,1)
            ax2.set_ylim(0,1)
            ax2.set_xticks(np.array([x1,x2]))
            ax2.set_xticklabels(np.array(['J11','This\nwork']),fontsize=10)
            ax2.set_yticklabels(np.repeat('',10))
            
            ax2.plot(x1,fdualwagn_to_wise_j11_l44_late_post_allsim,"*",
                     markersize=9,mew=1.6,mec='darkorange',color='None')
            ax2.errorbar(x1,fdualwagn_to_wise_j11_l44_late_post_allsim,yerr=yerrs_dualwagn_to_wise_j11_l44_late_post,
                         ecolor='darkorange',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)
            ax2.plot(x2,fdualwagn_to_wise_mycut_l44_late_post_allsim,"*",
                     markersize=9,mew=1.6,mec='darkorange',color='None')
            ax2.errorbar(x2,fdualwagn_to_wise_mycut_l44_late_post_allsim,yerr=yerrs_dualwagn_to_wise_mycut_l44_late_post,
                         ecolor='darkorange',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)
            
            ax2.plot(x1,fdualwagn_to_wise_j11_l44_tot_allsim,"o",ms=6,mew=1.2,mec='c',color='None')
            ax2.errorbar(x1,fdualwagn_to_wise_j11_l44_tot_allsim,yerr=yerrs_dualwagn_to_wise_j11_l44_tot,color='c',
                         ecolor='c',ls='None',capsize=4)
            ax2.plot(x2,fdualwagn_to_wise_mycut_l44_tot_allsim,"o",ms=6,mew=1.2,mec='c',color='None')
            ax2.errorbar(x2,fdualwagn_to_wise_mycut_l44_tot_allsim,yerr=yerrs_dualwagn_to_wise_mycut_l44_tot,color='c',
                         ecolor='c',ls='None',capsize=4)
            
            fig.suptitle(r'min L$_{\rm AGN} = 10^{44}$ erg s$^{-1}$')
            fig.subplots_adjust(bottom=0.16,left=0.14,right=0.92,top=0.9,wspace=0,hspace=0)
            wstr='_alt' if alt_wedge else ''
            fig.savefig(outdir+'/wise_dual_agn_to_wise_fraction_allsim_l44_%s%s_%s_%s%s.pdf'%(name,zstr,valtype_allsim,
                                                                                              errtype,wstr))

            ################ LAGN > 1e45, ALLSIM ################
            plt.clf()
            plt.cla()
            plt.close()
            fig=plt.figure(figsize=(5,3.5))
            ax=plt.subplot2grid((1,5),(0,0),colspan=4)
            ax.set_xlim(xlim)
            ax.set_ylim(0,1)
            ax.set_xlabel(r'min $W1 - W2$')
            ax.set_ylabel(r't$_{\rm WISEAGN}$ / t$_{\rm AGN}$')
            p1,=ax.plot(w12_lims_hires,f1or2agn_wise_hires_l45_late_post_agn_allsim,'*',
                        markersize=9,mew=1.6,mec='r',color='None')
            ax.errorbar(w12_lims_hires,f1or2agn_wise_hires_l45_late_post_agn_allsim,yerr=yerrs_wise_l45_late_post_agn,
                        ecolor='r',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)            
            p0,=ax.plot(w12_lims_hires,f1or2agn_wise_hires_l45_totagn_allsim,'o',
                        markersize=6,mew=1.2,mec='b',color='None')
            ax.errorbar(w12_lims_hires,f1or2agn_wise_hires_l45_totagn_allsim,yerr=yerrs_wise_l45_totagn,color='b',
                        ecolor='b',ls='None',capsize=4)
            ax.legend((p0,p1), ('Total','Late+Post'), fontsize=10, loc='lower left', 
                      numpoints=1, handletextpad=0.1)
            x1=0.25
            x2=0.75
            ax2=plt.subplot2grid((1,5),(0,4))
            ax2.set_xlim(0,1)
            ax2.set_ylim(0,1)
            ax2.set_xticks(np.array([x1,x2]))
            ax2.set_xticklabels(np.array(['J11','This\nwork']),fontsize=10)
            ax2.set_yticklabels(np.repeat('',10))
            
            ax2.plot(x1,f1or2agn_j11_l45_late_post_agn_allsim,"*",
                     markersize=9,mew=1.6,mec='darkorange',color='None')
            ax2.errorbar(x1,f1or2agn_j11_l45_late_post_agn_allsim,yerr=yerrs_j11_l45_late_post_agn,
                         ecolor='darkorange',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)
            ax2.plot(x2,f1or2agn_mycut_l45_late_post_agn_allsim,"*",
                     markersize=9,mew=1.6,mec='darkorange',color='None')
            ax2.errorbar(x2,f1or2agn_mycut_l45_late_post_agn_allsim,yerr=yerrs_mycut_l45_late_post_agn,
                         ecolor='darkorange',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)
            
            ax2.plot(x1,f1or2agn_j11_l45_totagn_allsim,"o",ms=6,mew=1.2,mec='c',color='None')
            ax2.errorbar(x1,f1or2agn_j11_l45_totagn_allsim,yerr=yerrs_j11_l45_totagn,color='c',
                         ecolor='c',ls='None',capsize=4)
            ax2.plot(x2,f1or2agn_mycut_l45_totagn_allsim,"o",ms=6,mew=1.2,mec='c',color='None')
            ax2.errorbar(x2,f1or2agn_mycut_l45_totagn_allsim,yerr=yerrs_mycut_l45_totagn,color='c',
                         ecolor='c',ls='None',capsize=4)
            
            fig.suptitle(r'min L$_{\rm AGN} = 10^{45}$ erg s$^{-1}$')
            fig.subplots_adjust(bottom=0.16,left=0.14,right=0.92,top=0.9,wspace=0,hspace=0)
            wstr='_alt' if alt_wedge else ''
            fig.savefig(outdir+'/wiseagn_fraction_allsim_l45_%s%s_%s_%s%s.pdf'%(name,zstr,valtype_allsim,
                                                                                errtype,wstr))


            ################ DUAL LAGN > 1e45, ALLSIM ################
            plt.clf()
            plt.cla()
            plt.close()
            fig=plt.figure(figsize=(5,3.5))
            ax=plt.subplot2grid((1,5),(0,0),colspan=4)
            ax.set_xlim(xlim)
            ax.set_ylim(0,1)
            ax.set_xlabel(r'min $W1 - W2$')
            ax.set_ylabel(r't$_{\rm WISE,Dual}$ / t$_{\rm Dual}$')
            p1,=ax.plot(w12_lims_hires,fdualagn_wise_hires_l45_late_post_agn_allsim,'*',
                        markersize=9,mew=1.6,mec='r',color='None')
            ax.errorbar(w12_lims_hires,fdualagn_wise_hires_l45_late_post_agn_allsim,yerr=yerrs_dual_wise_l45_late_post_agn,
                        ecolor='r',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)            
            p0,=ax.plot(w12_lims_hires,fdualagn_wise_hires_l45_totagn_allsim,'o',
                        markersize=6,mew=1.2,mec='b',color='None')
            ax.errorbar(w12_lims_hires,fdualagn_wise_hires_l45_totagn_allsim,yerr=yerrs_dual_wise_l45_totagn,color='b',
                        ecolor='b',ls='None',capsize=4)
            ax.legend((p0,p1), ('Total','Late'), fontsize=10, loc='lower left', 
                      numpoints=1, handletextpad=0.1)
            x1=0.25
            x2=0.75
            ax2=plt.subplot2grid((1,5),(0,4))
            ax2.set_xlim(0,1)
            ax2.set_ylim(0,1)
            ax2.set_xticks(np.array([x1,x2]))
            ax2.set_xticklabels(np.array(['J11','This\nwork']),fontsize=10)
            ax2.set_yticklabels(np.repeat('',10))
            
            ax2.plot(x1,fdualagn_j11_l45_late_post_agn_allsim,"*",
                     markersize=9,mew=1.6,mec='darkorange',color='None')
            ax2.errorbar(x1,fdualagn_j11_l45_late_post_agn_allsim,yerr=yerrs_dual_j11_l45_late_post_agn,
                         ecolor='darkorange',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)
            ax2.plot(x2,fdualagn_mycut_l45_late_post_agn_allsim,"*",
                     markersize=9,mew=1.6,mec='darkorange',color='None')
            ax2.errorbar(x2,fdualagn_mycut_l45_late_post_agn_allsim,yerr=yerrs_dual_mycut_l45_late_post_agn,
                         ecolor='darkorange',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)
            
            ax2.plot(x1,fdualagn_j11_l45_totagn_allsim,"o",ms=6,mew=1.2,mec='c',color='None')
            ax2.errorbar(x1,fdualagn_j11_l45_totagn_allsim,yerr=yerrs_dual_j11_l45_totagn,color='c',
                         ecolor='c',ls='None',capsize=4)
            ax2.plot(x2,fdualagn_mycut_l45_totagn_allsim,"o",ms=6,mew=1.2,mec='c',color='None')
            ax2.errorbar(x2,fdualagn_mycut_l45_totagn_allsim,yerr=yerrs_dual_mycut_l45_totagn,color='c',
                         ecolor='c',ls='None',capsize=4)
            
            fig.suptitle(r'min L$_{\rm AGN} = 10^{45}$ erg s$^{-1}$')
            fig.subplots_adjust(bottom=0.16,left=0.14,right=0.92,top=0.9,wspace=0,hspace=0)
            wstr='_alt' if alt_wedge else ''
            fig.savefig(outdir+'/wise_dual_agn_fraction_allsim_l45_%s%s_%s_%s%s.pdf'%(name,zstr,valtype_allsim,
                                                                                      errtype,wstr))


            ################ DUAL / WISE  LAGN > 1e45, ALLSIM ################
            plt.clf()
            plt.cla()
            plt.close()
            fig=plt.figure(figsize=(5,3.5))
            ax=plt.subplot2grid((1,5),(0,0),colspan=4)
            ax.set_xlim(xlim)
            ax.set_ylim(0,1)
            ax.set_xlabel(r'min $W1 - W2$')
            ax.set_ylabel(r't$_{\rm WISE,Dual}$ / t$_{\rm WISE}$')
            p1,=ax.plot(w12_lims_hires,fdualwagn_to_wise_hires_l45_late_post_allsim,'*',
                        markersize=9,mew=1.6,mec='r',color='None')
            ax.errorbar(w12_lims_hires,fdualwagn_to_wise_hires_l45_late_post_allsim,yerr=yerrs_dualwagn_to_wise_l45_late_post,
                        ecolor='r',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)
            
            p0,=ax.plot(w12_lims_hires,fdualwagn_to_wise_hires_l45_tot_allsim,'o',
                        markersize=6,mew=1.2,mec='b',color='None')
            ax.errorbar(w12_lims_hires,fdualwagn_to_wise_hires_l45_tot_allsim,yerr=yerrs_dualwagn_to_wise_l45_tot,color='b',
                        ecolor='b',ls='None',capsize=4)
            ax.legend((p0,p1), ('Total','Late+Post'), fontsize=10, loc='lower left', 
                      numpoints=1, handletextpad=0.1)
            
            x1=0.25
            x2=0.75
            ax2=plt.subplot2grid((1,5),(0,4))
            ax2.set_xlim(0,1)
            ax2.set_ylim(0,1)
            ax2.set_xticks(np.array([x1,x2]))
            ax2.set_xticklabels(np.array(['J11','This\nwork']),fontsize=10)
            ax2.set_yticklabels(np.repeat('',10))
            
            ax2.plot(x1,fdualwagn_to_wise_j11_l45_late_post_allsim,"*",
                     markersize=9,mew=1.6,mec='darkorange',color='None')
            ax2.errorbar(x1,fdualwagn_to_wise_j11_l45_late_post_allsim,yerr=yerrs_dualwagn_to_wise_j11_l45_late_post,
                         ecolor='darkorange',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)
            ax2.plot(x2,fdualwagn_to_wise_mycut_l45_late_post_allsim,"*",
                     markersize=9,mew=1.6,mec='darkorange',color='None')
            ax2.errorbar(x2,fdualwagn_to_wise_mycut_l45_late_post_allsim,yerr=yerrs_dualwagn_to_wise_mycut_l45_late_post,
                         ecolor='darkorange',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)
            
            ax2.plot(x1,fdualwagn_to_wise_j11_l45_tot_allsim,"o",ms=6,mew=1.2,mec='c',color='None')
            ax2.errorbar(x1,fdualwagn_to_wise_j11_l45_tot_allsim,yerr=yerrs_dualwagn_to_wise_j11_l45_tot,color='c',
                         ecolor='c',ls='None',capsize=4)
            ax2.plot(x2,fdualwagn_to_wise_mycut_l45_tot_allsim,"o",ms=6,mew=1.2,mec='c',color='None')
            ax2.errorbar(x2,fdualwagn_to_wise_mycut_l45_tot_allsim,yerr=yerrs_dualwagn_to_wise_mycut_l45_tot,color='c',
                         ecolor='c',ls='None',capsize=4)
            
            fig.suptitle(r'min L$_{\rm AGN} = 10^{45}$ erg s$^{-1}$')
            fig.subplots_adjust(bottom=0.16,left=0.14,right=0.92,top=0.9,wspace=0,hspace=0)
            wstr='_alt' if alt_wedge else ''
            fig.savefig(outdir+'/wise_dual_agn_to_wise_fraction_allsim_l45_%s%s_%s_%s%s.pdf'%(name,zstr,valtype_allsim,
                                                                                              errtype,wstr))


            #return


            ################ LAGN > 1e44, EACHSIM ################
            plt.clf()
            plt.cla()
            plt.close()
            fig=plt.figure(figsize=(5,3.5))
            ax=plt.subplot2grid((1,5),(0,0),colspan=4)
            ax.set_xlim(xlim)
            ax.set_ylim(0,1)
            ax.set_xlabel(r'min $W1 - W2$')
            ax.set_ylabel(r't$_{\rm WISEAGN}$ / t$_{\rm AGN}$')
            for i_s in range(nsim):
                if name == 'allres_d1d0':
                    lw_lp=0.75*i_s+0.5
                    lw_tot=0.75*i_s+1
                else:
                    lw_lp=2
                    lw_tot=1
                ax.plot((w12_lims_hires),f1or2agn_wise_hires_l44_late_post_agn_eachsim[:,i_s],
                        '*-',markersize=9,mew=1.6,mec='r',lw=lw_lp,color='r')
                ax.plot((w12_lims_hires),f1or2agn_wise_hires_l44_totagn_eachsim[:,i_s],
                        'o-',markersize=6,mew=1.2,mec='b',lw=lw_tot,color='b')
            ax.legend((p0,p1), ('Total','Late+Post'), fontsize=10, loc='lower left', 
                      numpoints=1, handletextpad=0.1)
            
            x1=0.25
            x2=0.75
            ax2=plt.subplot2grid((1,5),(0,4))
            ax2.set_xlim(0,1)
            ax2.set_ylim(0,1)
            ax2.set_xticks(np.array([x1,x2]))
            ax2.set_xticklabels(np.array(['J11','This\nwork']),fontsize=10)
            ax2.set_yticklabels(np.repeat('',10))
                        
            for i_s in range(nsim):
                if name == 'allres_d1d0':
                    mark_lp=i_s+9
                    mark_tot=i_s+6
                else:
                    mark_lp=9
                    mark_tot=6
                ax2.plot(x1,f1or2agn_j11_l44_late_post_agn_eachsim[i_s],"*",
                         ms=mark_lp,mew=1.6,mec='darkorange',color='None')
                ax2.plot(x2,f1or2agn_mycut_l44_late_post_agn_eachsim[i_s],"*",
                         ms=mark_lp,mew=1.6,mec='darkorange',color='None')
                ax2.plot(x1,f1or2agn_j11_l44_totagn_eachsim[i_s],"o",ms=mark_tot,mew=1.2,mec='c',color='None')
                ax2.plot(x2,f1or2agn_mycut_l44_totagn_eachsim[i_s],"o",ms=mark_tot,mew=1.2,mec='c',color='None')
            
            fig.suptitle(r'min L$_{\rm AGN} = 10^{44}$ erg s$^{-1}$')
            fig.subplots_adjust(bottom=0.16,left=0.14,right=0.92,top=0.9,wspace=0,hspace=0)
            wstr='_alt' if alt_wedge else ''
            fig.savefig(outdir+'/wiseagn_fraction_eachsim_l44_%s%s%s.pdf'%(name,zstr,wstr))

            print "f1or2agn_j11_l44_late_post_agn_eachsim:", f1or2agn_j11_l44_late_post_agn_eachsim
            print "f1or2agn_j11_l44_totagn_eachsim:", f1or2agn_j11_l44_totagn_eachsim
            print "f1or2agn_mycut_l44_late_post_agn_eachsim:", f1or2agn_mycut_l44_late_post_agn_eachsim
            print "f1or2agn_mycut_l44_totagn_eachsim:", f1or2agn_mycut_l44_totagn_eachsim

            ###########


        else:

            plt.clf()
            plt.cla()
            plt.close()
            fig=plt.figure(figsize=(5,3.5))
            ax=fig.add_subplot(111)
            ax.set_xlim(xlim)
            ax.set_ylim(0,1)
            ax.set_xlabel(r'min $W1 - W2$')
            ax.set_ylabel(r't$_{\rm WISEAGN}$ / t$_{\rm AGN}$')
            
            p1,=ax.plot(w12_lims_hires,f1or2agn_wise_hires_l44_late_post_agn_allsim,'*',
                        markersize=9,mew=1.6,mec='r',color='None')
            ax.errorbar(w12_lims_hires,f1or2agn_wise_hires_l44_late_post_agn_allsim,yerr=yerrs_wise_l44_late_post_agn,
                        ecolor='r',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)
            
            p0,=ax.plot(w12_lims_hires,f1or2agn_wise_hires_l44_totagn_allsim,'o',
                        markersize=6,mew=1.2,mec='b',color='None')
            ax.errorbar(w12_lims_hires,f1or2agn_wise_hires_l44_totagn_allsim,yerr=yerrs_wise_l44_totagn,color='b',
                        ecolor='b',ls='None',capsize=4)

            ax.legend((p0,p1), ('Total','Late+Post'), fontsize=10, loc='lower left', 
                      numpoints=1, handletextpad=0.1)
            fig.suptitle(r'min L$_{\rm AGN} = 10^{44}$ erg s$^{-1}$')
            fig.subplots_adjust(bottom=0.16,left=0.14,right=0.92,top=0.9)
            fig.savefig(outdir+'/wiseagn_fraction_allsim_nowedge_l44_%s%s_%s_%s.pdf'%(name,zstr,valtype_allsim,errtype))
            

            plt.clf()
            plt.cla()
            plt.close()
            fig=plt.figure(figsize=(5,3.5))
            ax=fig.add_subplot(111)
            ax.set_xlim(xlim)
            ax.set_ylim(0,1)
            ax.set_xlabel(r'min $W1 - W2$')
            ax.set_ylabel(r't$_{\rm WISEAGN}$ / t$_{\rm AGN}$')
            p1,=ax.plot(w12_lims_hires,f1or2agn_wise_hires_l45_late_post_agn_allsim,'*',
                        markersize=9,mew=1.6,mec='r',color='None')
            ax.errorbar(w12_lims_hires,f1or2agn_wise_hires_l45_late_post_agn_allsim,yerr=yerrs_wise_l45_late_post_agn,
                        ecolor='r',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)                
            p0,=ax.plot(w12_lims_hires,f1or2agn_wise_hires_l45_totagn_allsim,'o',
                        markersize=6,mew=1.2,mec='b',color='None')
            ax.errorbar(w12_lims_hires,f1or2agn_wise_hires_l45_totagn_allsim,yerr=yerrs_wise_l45_totagn,color='b',
                        ecolor='b',ls='None',capsize=4)
            ax.legend((p0,p1), ('Total','Late+Post'), fontsize=10, loc='lower left', 
                      numpoints=1, handletextpad=0.1)
            fig.suptitle(r'min L$_{\rm AGN} = 10^{45}$ erg s$^{-1}$')
            fig.subplots_adjust(bottom=0.16,left=0.14,right=0.92,top=0.9)
            fig.savefig(outdir+'/wiseagn_fraction_allsim_nowedge_l45_%s%s_%s_%s.pdf'%(name,zstr,valtype_allsim,errtype))


        print "f1or2agn_j11_l44_totagn_allsim",f1or2agn_j11_l44_totagn_allsim
        print "f1or2agn_j11_l45_totagn_allsim",f1or2agn_j11_l45_totagn_allsim
        print "t_1or2agn_j11_l45_tot",t_1or2agn_j11_l45_tot
        print "f1or2agn_j11_l44_late_post_agn_allsim",f1or2agn_j11_l44_late_post_agn_allsim
        print "f1or2agn_j11_l45_late_post_agn_allsim",f1or2agn_j11_l45_late_post_agn_allsim

        if allsim_only: return


        plot_multisim_totals((fdualagn_tot/(1.0*f1or2agn_tot),fdualagn_early/(1.0*f1or2agn_early),
                              fdualagn_late/(1.0*f1or2agn_late)),
                             ylbl = r'<t$_{\rm 2AGN}$/t$_{\rm AGN}$>',outdir=outdir,
                             extra=extra,tmax_postmrg=tmax_postmrg,latesep=latesep,xlbl_eachline=True,
                             leg_eachplot=True,leg_loc='lower left',
                             fbase='dual_agn_fraction',plot_titles=plot_titles,title='Dual AGN fraction')

        plot_multisim_totals((fdualagn_early/(1.0*f1or2agn_early),fdualagn_late/(1.0*f1or2agn_late)),
                             ylbl = r'<t$_{\rm 2AGN}$/t$_{\rm AGN}$>',outdir=outdir,
                             extra=extra,tmax_postmrg=tmax_postmrg,latesep=latesep,
                             leg_lbls=('Early','Late'),carr=['g','m'],msarr=['^','s'],
                             mewarr=[1.5,1.5],mszarr=[6,6],
                             fbase='dual_agn_fraction_notot',plot_titles=plot_titles,
                             title='Dual AGN fraction',leg_eachline=True)

        plot_multisim_totals((f1or2agn_tot,ftotagn_tot,f1or2agn_early,ftotagn_early,
                              f1or2agn_late,ftotagn_late,f1agn_post),
                             ylbl = r'<t$_{\rm AGN}$/t$_{\rm tot}$>',outdir=outdir,
                             extra=extra,tmax_postmrg=tmax_postmrg,latesep=latesep,
                             fbase='agn_fraction',plot_titles=plot_titles,
                             leg_lbls=('Total','Total (L1+2)','Early','Early (L1+2)',
                                       'Late','Late (L1+2)','Post'),
                             carr=2*['k']+2*['g']+2*['m']+['orange'],
                             msarr=2*['o']+2*['^']+2*['s']+['v'],
                             mewarr=2*[1.5]+5*[1.0],mszarr=[6,4,6,4,6,4,6],title='AGN fraction')

        plot_multisim_totals((f1or2agn_early,f1or2agn_late,f1agn_post),
                             ylbl = r'<t$_{\rm AGN}$/t$_{\rm tot}$>',outdir=outdir,
                             extra=extra,tmax_postmrg=tmax_postmrg,latesep=latesep,
                             fbase='agn_fraction_notot',plot_titles=plot_titles,xlbl_eachline=True,
                             leg_lbls=('Early','Late','Post'),carr=['g','m','orange'],msarr=['^','s','v'],
                             mewarr=3*[1.5],mszarr=[6,6,6],title='AGN fraction',leg_eachline=True)


        for iv,v in enumerate(('lgssfr','lgsfr','lglagn','lgltot','flbol')):
            plot_multisim_totals((f1or2agn_tot,ftotagn_tot,f1or2agn_early,ftotagn_early,
                                  f1or2agn_late,ftotagn_late,f1agn_post),xvar=v,xarr=xarrs[iv],
                                 ylbl = r'<t$_{\rm AGN}$/t$_{\rm tot}$>',outdir=outdir,
                                 extra=extra,tmax_postmrg=tmax_postmrg,latesep=latesep,
                                 fbase='agn_fraction_vs_max%s'%v,plot_titles=plot_titles,
                                 leg_lbls=('Total','Total (L1+2)','Early','Early (L1+2)',
                                       'Late','Late (L1+2)','Post'),
                                 carr=2*['k']+2*['g']+2*['m']+['orange'],
                                 msarr=2*['o']+2*['^']+2*['s']+['v'],
                                 mewarr=2*[1.5]+5*[1.0],mszarr=[6,4,6,4,6,4,6],title='AGN fraction')

    if nsim > 1 and multisim_plots and not ngas_only:
        for iplot,ix in enumerate([ix_w05,ix_w08,ix_j11]):
            print '%s plots...'%wfbase[iplot]

            #if agnx0:
                #plot_multisim_totals((fnoagn_wise_tot[:,ix],fnoagn_wise_early[:,ix],
                #                      fnoagn_wise_late[:,ix],fnoagn_wise_post[:,ix]),
                #                     ylbl = r' <t$_{\rm WISE}$/t$_{\rm tot}$>',
                #                     outdir=outdir,extra=extra,tmax_postmrg=tmax_postmrg,latesep=latesep,
                #                     fbase='wise_noagn_fraction_%s'%wfbase[iplot],plot_titles=plot_titles,
                #                     leg_lbls=('Total','Early','Late','Post'),
                #                     carr=['k','g','m','orange'],msarr=['o']+['^']+['s']+['v'],
                #                     mewarr=[1.5]+3*[1.0],mszarr=4*[6],
                #                     title='WISE fraction (%s) - No AGN'%wisemask_labels[iplot])

            if not agnx0:
                plot_multisim_totals((f1or2agn_wise_tot[:,ix]/f1or2agn_tot,
                                      f1or2agn_wise_early[:,ix]/f1or2agn_early,
                                      f1or2agn_wise_late[:,ix]/f1or2agn_late,
                                      f1agn_wise_post[:,ix]/f1agn_post),
                                     ylbl = r' <t$_{\rm WISE AGN}$/t$_{\rm AGN}$>',
                                     outdir=outdir,extra=extra,tmax_postmrg=tmax_postmrg,latesep=latesep,
                                     fbase='wise_agn_fraction_%s'%wfbase[iplot],plot_titles=plot_titles,
                                     leg_lbls=('Total','Early','Late','Post'),xlbl_eachline=True,
                                     carr=['k','g','m','orange'],msarr=['o']+['^']+['s']+['v'],
                                     mewarr=[1.5]+3*[1.0],mszarr=4*[6],
                                     title='WISE AGN fraction (%s)'%wisemask_labels[iplot])
                
                plot_multisim_totals((f1or2agn_wise_early[:,ix]/f1or2agn_early,
                                      f1or2agn_wise_late[:,ix]/f1or2agn_late,
                                      f1agn_wise_post[:,ix]/f1agn_post),
                                     ylbl = r' <t$_{\rm WISE AGN}$/t$_{\rm AGN}$>',
                                     outdir=outdir,extra=extra,tmax_postmrg=tmax_postmrg,latesep=latesep,
                                     fbase='wise_agn_fraction_notot_%s'%wfbase[iplot],plot_titles=plot_titles,
                                     leg_lbls=('Early','Late','Post'),
                                     carr=['g','m','orange'],msarr=['^']+['s']+['v'],
                                     mewarr=3*[1.5],mszarr=3*[6],
                                     title='WISE AGN fraction (%s)'%wisemask_labels[iplot],leg_eachline=True)

                ix_w44 = np.array([k for k in ix if '44' in wise_agn_titles[k]])
                ix_44 = np.array([k for k in range(nagnmask) if '44' in agn_labels[k]])
                plot_multisim_totals((f1or2agn_wise_early[:,ix_w44:ix_w44+2]/f1or2agn_early[:,ix_44:ix_44+2],
                                      f1or2agn_wise_late[:,ix_w44:ix_w44+2]/f1or2agn_late[:,ix_44:ix_44+2],
                                      f1agn_wise_post[:,ix_w44:ix_w44+2]/f1agn_post[:,ix_44:ix_44+2]),
                                     nplots=2, ylbl = r' <t$_{\rm WISE AGN}$/t$_{\rm AGN}$>',
                                     outdir=outdir,extra=extra,tmax_postmrg=tmax_postmrg,latesep=latesep,
                                     fbase='wise_agn_fraction_notot_2panel_%s'%wfbase[iplot],
                                     plot_titles=plot_titles[ix_44:ix_44+2],leg_lbls=('Early','Late','Post'),
                                     xlbl_eachline=True,carr=['g','m','orange'],msarr=['^']+['s']+['v'],
                                     mewarr=3*[1.5],mszarr=3*[6], title='%s'%wisemask_labels[iplot],
                                     leg_eachline=True,leg_loc='lower right')

                
                plot_multisim_totals((fdualagn_wise_tot[:,ix]/fdualagn_tot,
                                      fdualagn_wise_early[:,ix]/fdualagn_early,
                                      fdualagn_wise_late[:,ix]/fdualagn_late),
                                     ylbl = r' <t$_{\rm WISE dual}$/t$_{\rm dual}$>',
                                     outdir=outdir,extra=extra,tmax_postmrg=-1,latesep=latesep,
                                     fbase='wise_dual_agn_fraction_%s'%wfbase[iplot],plot_titles=plot_titles,
                                     title='WISE dual AGN fraction (%s)'%wisemask_labels[iplot])
                
                plot_multisim_totals((fdualagn_wise_early[:,ix]/fdualagn_early,
                                      fdualagn_wise_late[:,ix]/fdualagn_late),
                                     ylbl = r' <t$_{\rm WISE dual}$/t$_{\rm dual}$>',
                                     outdir=outdir,extra=extra,tmax_postmrg=-1,latesep=latesep,
                                     fbase='wise_dual_agn_fraction_notot_%s'%wfbase[iplot],
                                     plot_titles=plot_titles,leg_lbls=('Early','Late'),xlbl_eachline=True,
                                     carr=['g','m'],msarr=['^']+['s'],mewarr=2*[1.5],mszarr=2*[6],
                                     title='WISE dual AGN fraction (%s)'%wisemask_labels[iplot],
                                     leg_loc='lower right',leg_eachplot=True)                

                plot_multisim_totals((fdualagn_wise_early[:,ix_w44]/fdualagn_early[:,ix_44],
                                      fdualagn_wise_late[:,ix_w44]/fdualagn_late[:,ix_44]),
                                     nplots=1, ylbl = r' <t$_{\rm WISE dual}$/t$_{\rm dual}$>',
                                     outdir=outdir,extra=extra,tmax_postmrg=-1,latesep=latesep,
                                     fbase='wise_dual_agn_fraction_notot_1panel_lagn44_%s'%wfbase[iplot],
                                     plot_titles=[wisemask_labels[iplot]],leg_lbls=('Early','Late'),xlbl_eachline=True,
                                     carr=['g','m'],msarr=['^']+['s'],mewarr=2*[1.5],mszarr=2*[6],
                                     title=wisemask_labels[iplot],leg_loc='lower right',leg_eachplot=True) 
                
                plot_multisim_totals((fdualagn_wise_tot[:,ix]/(f1or2agn_wise_tot[:,ix]+fnoagn_wise_tot[:,ix]),
                                      fdualagn_wise_early[:,ix]/(f1or2agn_wise_early[:,ix]+fnoagn_wise_early[:,ix]),
                                      fdualagn_wise_late[:,ix]/(f1or2agn_wise_late[:,ix]+fnoagn_wise_late[:,ix])),
                                     ylbl = r' <t$_{\rm WISE dual}$/t$_{\rm WISE}$>',
                                     outdir=outdir,extra=extra,tmax_postmrg=-1,latesep=latesep,
                                     fbase='wise_dual_agn_to_wise_fraction_%s'%wfbase[iplot],
                                     plot_titles=plot_titles,leg_lbls=('Total','Early','Late'),
                                     leg_eachplot=True,leg_loc='lower left',
                                     carr=['k','g','m'],msarr=['o','^','s'],mewarr=3*[1.5],mszarr=3*[6],
                                     title='WISE dual AGN fraction (%s)'%wisemask_labels[iplot])

                plot_multisim_totals((fdualagn_wise_tot[:,ix_w44]/(f1or2agn_wise_tot[:,ix_w44]+
                                                                   fnoagn_wise_tot[:,ix_w44]),
                                      fdualagn_wise_early[:,ix_w44]/(f1or2agn_wise_early[:,ix_w44]+
                                                                     fnoagn_wise_early[:,ix_w44]),
                                      fdualagn_wise_late[:,ix_w44]/(f1or2agn_wise_late[:,ix_w44]+
                                                                    fnoagn_wise_late[:,ix_w44])),
                                     nplots=1, ylbl = r' <t$_{\rm WISE dual}$/t$_{\rm WISE}$>',
                                     outdir=outdir,extra=extra,tmax_postmrg=-1,latesep=latesep,
                                     fbase='wise_dual_agn_to_wise_fraction_1panel_lagn44_%s'%wfbase[iplot],
                                     leg_lbls=('Total','Early','Late'),leg_eachplot=True,leg_loc='lower left',
                                     carr=['k','g','m'],msarr=['o','^','s'],mewarr=3*[1.5],mszarr=3*[6],
                                     plot_titles=[wisemask_labels[iplot]],title=wisemask_labels[iplot])

                for iv,v in enumerate(('lgssfr','lgsfr','lglagn','lgltot','flbol')):
                    plot_multisim_totals((f1or2agn_wise_tot[:,ix]/f1or2agn_tot,
                                          f1or2agn_wise_early[:,ix]/f1or2agn_early,
                                          f1or2agn_wise_late[:,ix]/f1or2agn_late,
                                          f1agn_wise_post[:,ix]/f1agn_post),
                                         ylbl = r' <t$_{\rm WISE AGN}$/t$_{\rm AGN}$>',xvar=v,xarr=xarrs[iv],
                                         outdir=outdir,extra=extra,tmax_postmrg=tmax_postmrg,latesep=latesep,
                                         fbase='wise_agn_fraction_%s_vs_max%s'%(wfbase[iplot],v),plot_titles=plot_titles,
                                         leg_lbls=('Total','Early','Late','Post'),
                                         carr=['k','g','m','orange'],msarr=['o']+['^']+['s']+['v'],
                                         mewarr=[1.5]+3*[1.0],mszarr=4*[6],
                                         title='WISE AGN fraction (%s)'%wisemask_labels[iplot])

                    plot_multisim_totals((fdualagn_wise_tot[:,ix]/fdualagn_tot,
                                          fdualagn_wise_early[:,ix]/fdualagn_early,
                                          fdualagn_wise_late[:,ix]/fdualagn_late),
                                         ylbl = r' <t$_{\rm WISE dual}$/t$_{\rm dual}$>',xvar=v,xarr=xarrs[iv],
                                         outdir=outdir,extra=extra,tmax_postmrg=-1,latesep=latesep,
                                         fbase='wise_dual_agn_fraction_%s_vs_max%s'%(wfbase[iplot],v),plot_titles=plot_titles,
                                         title='WISE dual AGN fraction (%s)'%wisemask_labels[iplot])
                    
                    plot_multisim_totals((fdualagn_wise_tot[:,ix]/(f1or2agn_wise_tot[:,ix]+fnoagn_wise_tot[:,ix]),
                                          fdualagn_wise_early[:,ix]/(f1or2agn_wise_early[:,ix]+fnoagn_wise_early[:,ix]),
                                          fdualagn_wise_late[:,ix]/(f1or2agn_wise_late[:,ix]+fnoagn_wise_late[:,ix])),
                                         ylbl = r' <t$_{\rm WISE dual}$/t$_{\rm WISE}$>',xvar=v,xarr=xarrs[iv],
                                         outdir=outdir,extra=extra,tmax_postmrg=-1,latesep=latesep,
                                         fbase='wise_dual_agn_to_wise_fraction_%s_vs_max%s'%(wfbase[iplot],v),
                                         plot_titles=plot_titles,leg_lbls=('Total','Early','Late'),
                                         carr=['k','g','m'],msarr=['o','^','s'],mewarr=3*[1.5],mszarr=3*[6],
                                         title='WISE dual AGN fraction (%s)'%wisemask_labels[iplot])


            #plot_multisim_totals((fdualagn_wise_tot[:,ix],fdualagn_nowise_tot[:,ix],f1agn_wise_tot[:,ix],
            #                      fdualagn_wise_early[:,ix],fdualagn_nowise_early[:,ix],f1agn_wise_early[:,ix],
            #                      fdualagn_wise_late[:,ix],fdualagn_nowise_late[:,ix],f1agn_wise_late[:,ix],
            #                      f1agn_nowise_post[:,ix],f1agn_wise_post[:,ix]),
            #                     ylbl = wisemask_labels[iplot]+r' <t/t$_{\rm tot}$>',
            #                     outdir=outdir,extra=extra,tmax_postmrg=tmax_postmrg,latesep=latesep,
            #                     fbase='wise_nowise_dual_agn_fraction_%s'%wfbase[iplot],plot_titles=plot_titles,
            #                     leg_lbls=('WISE dual (Total)', 'Dual, no WISE (Total)', 'WISE 1AGN (Total)',
            #                               'WISE dual (Early)', 'Dual, no WISE (Early)', 'WISE 1AGN (Early)',
            #                               'WISE dual (Late)', 'Dual, no WISE (Late)', 'WISE 1AGN (Late)',
            #                               'AGN, no WISE (Post)', 'WISE AGN (Post)'),
            #                     carr=3*['darkblue','g','c']+['g','c'],msarr=3*['o']+3*['^']+3*['s']+2*['v'],
            #                     mewarr=3*[1.5]+8*[1.0],mszarr=11*[6])

            ##plot_multisim_totals((f1or2agn_wise_tot[:,ix],f1or2agn_wise_early[:,ix],f1or2agn_wise_late[:,ix],
            ##                      f1agn_wise_post[:,ix],f1or2agn_nowise_tot[:,ix],f1or2agn_nowise_early[:,ix],
            ##                      f1or2agn_nowise_late[:,ix],f1agn_nowise_post[:,ix],fnoagn_wise_tot[:,ix],
            ##                      fnoagn_wise_early[:,ix],fnoagn_wise_late[:,ix],fnoagn_wise_post[:,ix]),
            ##                     ylbl = wisemask_labels[iplot]+r' <t/t$_{\rm tot}$>',
            ##                     outdir=outdir,extra=extra,tmax_postmrg=tmax_postmrg,
            ##                     fbase='wise_nowise_agn_fraction_%s'%wfbase[iplot],plot_titles=plot_titles,
            ##                     leg_lbls=(),
            ##                     carr=['g']+['m']+['orange']+['k'],
            ##                     msarr=['^']+['s']+['v']+['o'],
            ##                     mewarr=3*[1.0]+[1.5],mszarr=4*[6])

            #plot_multisim_totals((ftotagn_wise_tot[:,ix],ftotagn_wise_early[:,ix],ftotagn_wise_late[:,ix],
            #                      f1agn_wise_post[:,ix],ftotagn_nowise_tot[:,ix],ftotagn_nowise_early[:,ix],
            #                      ftotagn_nowise_late[:,ix],f1agn_nowise_post[:,ix],fnototagn_wise_tot[:,ix],
            #                      fnototagn_wise_early[:,ix],fnototagn_wise_late[:,ix],fnoagn_wise_post[:,ix]),
            #                     ylbl = wisemask_labels[iplot]+r' <t/t$_{\rm tot}$>',
            #                     outdir=outdir,extra=extra,tmax_postmrg=tmax_postmrg,
            #                     fbase='wise_nowise_agn_fraction_%s'%wfbase[iplot],plot_titles=plot_titles,
            #                     leg_lbls=(),
            #                     carr=['g']+['m']+['orange']+['k'],
            #                     msarr=['^']+['s']+['v']+['o'],
            #                     mewarr=3*[1.0]+[1.5],mszarr=4*[6])



    ###############################
    ### MERGER AGN PHASES PLOTS ###
    ###############################

    if nsim>1 and not multisim_plots: return

    if not agnx0 and not ngas_only:
        print "Making merger AGN phase plots..."

        bw=0.3

        if equal_bins:
            xlim = (-2*bw,len(sep_bins[sep_bins<=latesep])-1+2*bw)
            xticklbl = ["%g"%(sep_bins[j]) for j in np.arange(len(sep_bins[sep_bins<=latesep]))]
        else:
            xticklbl = ["%g-%g"%(sep_bins[j],sep_bins[j+1]) for j in np.arange(len(sep_bins)-1)]
            if not use_logbins: xticklbl[0] = 'post-mrg'
        if plot_3dsep:
            xlabel = '3D sep [%skpc]'%('log ' if use_logbins else '')
        else:
            xlabel = 'proj. sep. [%skpc]'%('log ' if use_logbins else '')

        plt.close('all')
        plt.ioff()        
        for ptype in ['dt','dtfrac']:

            ylim = (0.0,1.41) if ptype=='dt' else (0.0,1.05)
            ylabel = r'$\Delta$t [Gyr]' if ptype=='dt' else r'$\Delta$t/t$_{\rm bin}$'
            fig = plt.figure(figsize=(9,6))
            for i in range(6):
                ax = fig.add_subplot(231+i)
                ax.set_ylim(ylim)
                ax.set_ylabel(ylabel) if  i in [0,3] else ax.set_ylabel('')
                if equal_bins:
                    ax.set_xlim(xlim)
                    ax.set_xticks(np.arange(len(sep_bins[sep_bins<latesep])+1))
                else:
                    ax.set_xticks(np.arange(len(sep_bins)-1))
                ax.set_xticklabels(xticklbl,fontsize=9)
                if ptype == 'dt':
                    var = hist3d_dt if plot_3dsep else histproj_dt
                    tot=cmph.makebar(ax,var,nbins=len(sep_bins)-1,color='c',width=bw,xoffset=-bw,
                                     val='median',errtype=errtype)
                    agn_var = agn_hist3d[:,i,:] if plot_3dsep else agn_histproj_dt[:,i,:,:]
                    agn=cmph.makebar(ax,agn_var,nbins=len(sep_bins)-1,color='m',width=bw,
                                     xoffset=0,val='median',errtype=errtype)
                else:
                    agn_var = agn_hist3d_dt_frac[:,i,:] if plot_3dsep else agn_histproj_dt_frac[:,i,:,:]
                    agn=cmph.makebar(ax,agn_var,nbins=len(sep_bins)-1,color='m',
                                     width=2*bw,xoffset=-bw,val='median',errtype=errtype)
                ax.set_xlabel(xlabel) if i in [3,4,5] else ax.set_xlabel('')
                if ptype=='dt': ax.legend((tot,agn),('Total','AGN'),fontsize=9,loc='upper left')
                ax.set_title(plot_titles[i],fontsize=10)
            fig.subplots_adjust(wspace=0.2,hspace=0.2,left=0.08,right=0.94,bottom=0.1,top=0.92)

            fbase_extra = '' if ptype=='dt' else '_'+ptype
            plotname='%s/merger_agn_phases%s%s%s%s_tpost%g.pdf'%(outdir,fbase_extra,sepstring,binstr,extra,tmax_postmrg)
            fig.savefig(plotname,format='pdf')
            plt.cla()
            plt.clf()

    ###############################
    ###   NGAS VS BHSEP PLOTS   ###
    ###############################

    if equal_bins: 
        print "equal_bins=True; skipping Ngas plots."
        return

    print "Making Ngas vs BH sep plots..."
    ngastype='_grid' if sfrhist_w_grid else '_aux'
    if not agnx0:

        bw=1.0
        xlabel = r'projected sep. [kpc]'
        xlabel_post = r't(post-merger) [Myr]'
        xlim_post = (-(len(tpostmrg_bins)-2)-0.5*bw-0.2,0+0.5*bw+0.2)
        xticklbl_post = [r"%g"%(t*1000) for t in tpostmrg_bins[1:]][::-1]
        xticks_post = np.arange(-(len(tpostmrg_bins)-2),1)
        if equal_bins:
        #xlim_pre = (-2*bw,len(sep_bins[sep_bins<=latesep])-1+2*bw)
            xlim_pre = (1-0.5*bw-0.2,(len(sep_bins[sep_bins<=latesep])-1)+0.5*bw+0.2)
            xticklbl_pre = ["%g"%(sep_bins[j]) for j in np.arange(1,len(sep_bins[sep_bins<=latesep]))]
            xticks_pre = np.arange(1,len(sep_bins[sep_bins<=latesep]))
        else:
            xlim_pre = (1-0.5*bw-0.2,(len(sep_bins)-3)+0.5*bw+0.2)
            xticklbl_pre = ["%g-%g"%(sep_bins[j],sep_bins[j+1]) for j in np.arange(1,len(sep_bins)-2)]
            xticks_pre = np.arange(1,len(sep_bins)-2)
        if (len(xticklbl_post)!=len(xticks_post)) or (len(xticklbl_pre)!=len(xticks_pre)):
            print 'Error: length of xticklbl doesnt match xticks:'
            print xticklbl_post
            print xticks_post
            print xticklbl_pre
            print xticks_pre
            sys.exit()

        for j_def,agndef in enumerate(['fedd','lagn','flbol']):


            fig = plt.figure(figsize=(6.5,8))
            for i in range(3):

                axa = fig.add_subplot(321+i*2)
                axa.set_yscale('log')
                axa.set_ylim(3.0e20,3.0e24)
                #axa.set_ylabel(r'N$_{\rm gas}$ [cm$^{-2}$]',fontsize=12)
                axa.set_ylabel(r'N$_{\rm H}$ [cm$^{-2}$]',fontsize=12)
                axa.set_xlim(xlim_post)
                axa.set_xticks(xticks_post)
                axa.set_xticklabels(xticklbl_post,fontsize=10)
                axa.set_xlabel(xlabel_post,fontsize=11)
                axa.set_title(plot_titles[3*j_def+i],fontsize=11)
                ngpost=cmph.makebar(axa,Ngas_agn_hist_postmrg[:,3*j_def+i,:,:],nbins=len(tpostmrg_bins)-1,ecolor='k',
                                    color='k',width=0.5*bw,xoffset=0,val='median',errtype=errtype,
                                    verbose=verbose,binpoints_only=True,reverse=True)

                axb = fig.add_subplot(321+i*2+1)
                axb.set_yscale('log')
                axb.set_ylim(3.0e20,3.0e24)
                axb.set_ylabel('')
                axb.set_xlim(xlim_pre)
                axb.set_xticks(xticks_pre)
                axb.set_xticklabels(xticklbl_pre,fontsize=10)
                axb.set_yticklabels('')
                axb.set_xlabel(xlabel,fontsize=11)
                axb.set_title(plot_titles[3*j_def+i],fontsize=11)
                ng1=cmph.makebar(axb,Ngas_bh1_agn_histproj[:,3*j_def+i,:,1:-1],xmin=1,
                                 nbins=(len(sep_bins)-1)-2,ecolor='b',
                                 color='b',width=0.5*bw,xoffset=-0.05*bw,val='median',errtype=errtype,
                                 verbose=verbose,binpoints_only=True)
                ng2=cmph.makebar(axb,Ngas_bh2_agn_histproj[:,3*j_def+i,:,1:-1],xmin=1,
                                 nbins=(len(sep_bins)-1)-2,ecolor='m',
                                 color='m',width=0.5*bw,xoffset=0.05*bw,val='median',errtype=errtype,
                                 verbose=verbose, binpoints_only=True)
                axb.legend((ng1,ng2),('BH1','BH2'),fontsize=9,loc='upper right',numpoints=1)
            
            fig.subplots_adjust(hspace=0.42,wspace=0.08,left=0.14,right=0.96,bottom=0.08,top=0.94)
            plotname='%s/Ngas%s%s_%s_merger_phases%s%s_tpost%g_%s%s.pdf'%(outdir,ngastype,rstr,agndef,sepstring,binstr,
                                                                            tmax_postmrg,errtype,extra)
            fig.savefig(plotname,format='pdf')
            plt.cla()
            plt.clf()
            plt.close()
            
            ### START HERE FOR WISE AGN NGAS PLOTS ###
            #ix_fedd = np.array([k for k in range(nwagnmask) if 'fEdd' in wise_agn_titles[k]])
            #ix_lagn = np.array([k for k in range(nwagnmask) if 'Lagn' in wise_agn_titles[k]])
            #ix_flbol = np.array([k for k in range(nwagnmask) if 'fLbol' in wise_agn_titles[k]])
            #ix_w05 = np.array([k for k in range(nwagnmask) if 'W1W2>0.5' in wise_agn_titles[k]])
            #ix_w08 = np.array([k for k in range(nwagnmask) if 'W1W2>0.8' in wise_agn_titles[k]])
            #ix_j11 = np.array([k for k in range(nwagnmask) if 'J11' in wise_agn_titles[k]])

            if agndef=='lagn':
                ix_w05_lagn44 = np.intersect1d(ix_w05,ix_lagn)[1]
                ix_w08_lagn44 = np.intersect1d(ix_w08,ix_lagn)[1]
                ix_j11_lagn44 = np.intersect1d(ix_j11,ix_lagn)[1]
                ix_w05_lagn45 = np.intersect1d(ix_w05,ix_lagn)[2]
                ix_w08_lagn45 = np.intersect1d(ix_w08,ix_lagn)[2]
                ix_j11_lagn45 = np.intersect1d(ix_j11,ix_lagn)[2]

                fig = plt.figure(figsize=(6.5,8))

                for i,windex in enumerate((ix_w05_lagn44,ix_w08_lagn44,ix_j11_lagn44)):
                    titles=[r"%s; L$_{\rm AGN} > 10^{44}$ erg s$^{-1}$"%s for s in ('W1-W2>0.5','W1-W2>0.8','J11')]

                    print i,windex
                    axa = fig.add_subplot(321+i*2)
                    axa.set_yscale('log')
                    axa.set_ylim(3.0e20,3.0e24)
                #axa.set_ylabel(r'N$_{\rm gas}$ [cm$^{-2}$]',fontsize=12)
                    axa.set_ylabel(r'N$_{\rm H}$ [cm$^{-2}$]',fontsize=12)
                    axa.set_xlim(xlim_post)
                    axa.set_xticks(xticks_post)
                    axa.set_xticklabels(xticklbl_post,fontsize=10)
                    axa.set_xlabel(xlabel_post,fontsize=11)
                    axa.set_title(titles[i],fontsize=11)
                    ngpost=cmph.makebar(axa,Ngas_agn_hist_postmrg[:,3*2+1,:,:],nbins=len(tpostmrg_bins)-1,ecolor='c',
                                    color='c',width=0.5*bw,xoffset=-0.05*bw,val='median',errtype=errtype,
                                    verbose=verbose,binpoints_only=True,reverse=True,elinewidth=1.2,capthick=1.2)
                    ngpost=cmph.makebar(axa,Ngas_wiseagn_hist_postmrg[:,windex,:,:],nbins=len(tpostmrg_bins)-1,ecolor='k',
                                        color='k',width=0.5*bw,xoffset=0.05*bw,val='median',errtype=errtype,
                                        verbose=verbose,binpoints_only=True,reverse=True)
                    
                    axb = fig.add_subplot(321+i*2+1)
                    axb.set_yscale('log')
                    axb.set_ylim(3.0e20,3.0e24)
                    axb.set_ylabel('')
                    axb.set_xlim(xlim_pre)
                    axb.set_xticks(xticks_pre)
                    axb.set_xticklabels(xticklbl_pre,fontsize=10)
                    axb.set_yticklabels('')
                    axb.set_xlabel(xlabel,fontsize=11)
                    axb.set_title(titles[i],fontsize=11)
                    ng1=cmph.makebar(axb,np.mean([Ngas_bh1_agn_histproj[:,3*2+1,:,1:-1],
                                                  Ngas_bh2_agn_histproj[:,3*2+1,:,1:-1]],axis=0),
                                     xmin=1,nbins=(len(sep_bins)-1)-2,ecolor='c',
                                     color='c',width=0.5*bw,xoffset=-0.1*bw,val='median',errtype=errtype,
                                     verbose=verbose,binpoints_only=True,elinewidth=1.2,capthick=1.2)
                    ng2=cmph.makebar(axb,np.mean([Ngas_bh1_wiseagn_histproj[:,windex,:,1:-1],
                                                  Ngas_bh2_wiseagn_histproj[:,windex,:,1:-1]],axis=0),
                                     xmin=1,nbins=(len(sep_bins)-1)-2,ecolor='k',
                                     color='k',width=0.5*bw,xoffset=0,val='median',errtype=errtype,
                                     verbose=verbose,binpoints_only=True)
                    #ng2=cmph.makebar(axb,Ngas_bh2_wiseagn_histproj[:,windex,:,1:-1],xmin=1,
                    ##                 nbins=(len(sep_bins)-1)-2,ecolor='m',
                    #                 color='m',width=0.5*bw,xoffset=0.05*bw,val='median',errtype=errtype,
                    #                 verbose=verbose, binpoints_only=True)
                    ng3=cmph.makebar(axb,np.mean([Ngas_bh1_wisedualagn_histproj[:,windex,:,1:-1],
                                                  Ngas_bh2_wisedualagn_histproj[:,windex,:,1:-1]],axis=0),
                                     xmin=1, nbins=(len(sep_bins)-1)-2,ecolor='r',
                                     color='r',width=0.5*bw,xoffset=0.1*bw,val='median',errtype=errtype,
                                     verbose=verbose, binpoints_only=True,elinewidth=2,capthick=2)
                    #ng4=cmph.makebar(axb,Ngas_bh2_wisedualagn_histproj[:,windex,:,1:-1],xmin=1,
                    #                 nbins=(len(sep_bins)-1)-2,ecolor='m',
                    #                 color='m',width=0.5*bw,xoffset=0.05*bw,val='median',errtype=errtype,
                    #                 verbose=verbose, binpoints_only=True)
                    axb.legend((ng1,ng2,ng3),('AGN','MIR AGN','Dual MIR AGN'),fontsize=9,loc='upper right',
                               numpoints=1,handletextpad=0.08)
                    
                fig.subplots_adjust(hspace=0.42,wspace=0.08,left=0.14,right=0.96,bottom=0.08,top=0.94)
                plotname='%s/Ngas%s%s_wise_lagn44_merger_phases%s%s_tpost%g_%s%s.pdf'%(outdir,ngastype,rstr,sepstring,binstr,
                                                                                       tmax_postmrg,errtype,extra)
                fig.savefig(plotname,format='pdf')
                plt.cla()
                plt.clf()
                plt.close()


                for i,windex in enumerate((ix_w05_lagn45,ix_w08_lagn45,ix_j11_lagn45)):
                    titles=[r"%s; L$_{\rm AGN} > 10^{45}$ erg s$^{-1}$"%s for s in ('W1-W2>0.5','W1-W2>0.8','J11')]

                    print i,windex
                    axa = fig.add_subplot(321+i*2)
                    axa.set_yscale('log')
                    axa.set_ylim(3.0e20,3.0e24)
                #axa.set_ylabel(r'N$_{\rm gas}$ [cm$^{-2}$]',fontsize=12)
                    axa.set_ylabel(r'N$_{\rm H}$ [cm$^{-2}$]',fontsize=12)
                    axa.set_xlim(xlim_post)
                    axa.set_xticks(xticks_post)
                    axa.set_xticklabels(xticklbl_post,fontsize=10)
                    axa.set_xlabel(xlabel_post,fontsize=11)
                    axa.set_title(titles[i],fontsize=11)
                    ngpost=cmph.makebar(axa,Ngas_agn_hist_postmrg[:,3*2+2,:,:],nbins=len(tpostmrg_bins)-1,ecolor='c',
                                    color='c',width=0.5*bw,xoffset=-0.05*bw,val='median',errtype=errtype,
                                    verbose=verbose,binpoints_only=True,reverse=True,elinewidth=1.2,capthick=1.2)
                    ngpost=cmph.makebar(axa,Ngas_wiseagn_hist_postmrg[:,windex,:,:],nbins=len(tpostmrg_bins)-1,ecolor='k',
                                        color='k',width=0.5*bw,xoffset=0.05*bw,val='median',errtype=errtype,
                                        verbose=verbose,binpoints_only=True,reverse=True)
                    
                    axb = fig.add_subplot(321+i*2+1)
                    axb.set_yscale('log')
                    axb.set_ylim(3.0e20,3.0e24)
                    axb.set_ylabel('')
                    axb.set_xlim(xlim_pre)
                    axb.set_xticks(xticks_pre)
                    axb.set_xticklabels(xticklbl_pre,fontsize=10)
                    axb.set_yticklabels('')
                    axb.set_xlabel(xlabel,fontsize=11)
                    axb.set_title(titles[i],fontsize=11)
                    ng1=cmph.makebar(axb,np.mean([Ngas_bh1_agn_histproj[:,3*2+2,:,1:-1],
                                                  Ngas_bh2_agn_histproj[:,3*2+2,:,1:-1]],axis=0),
                                     xmin=1,nbins=(len(sep_bins)-1)-2,ecolor='c',
                                     color='c',width=0.5*bw,xoffset=-0.1*bw,val='median',errtype=errtype,
                                     verbose=verbose,binpoints_only=True,elinewidth=1.2,capthick=1.2)
                    ng2=cmph.makebar(axb,np.mean([Ngas_bh1_wiseagn_histproj[:,windex,:,1:-1],
                                                  Ngas_bh2_wiseagn_histproj[:,windex,:,1:-1]],axis=0),
                                     xmin=1,nbins=(len(sep_bins)-1)-2,ecolor='k',
                                     color='k',width=0.5*bw,xoffset=0,val='median',errtype=errtype,
                                     verbose=verbose,binpoints_only=True)
                    #ng2=cmph.makebar(axb,Ngas_bh2_wiseagn_histproj[:,windex,:,1:-1],xmin=1,
                    ##                 nbins=(len(sep_bins)-1)-2,ecolor='m',
                    #                 color='m',width=0.5*bw,xoffset=0.05*bw,val='median',errtype=errtype,
                    #                 verbose=verbose, binpoints_only=True)
                    ng3=cmph.makebar(axb,np.mean([Ngas_bh1_wisedualagn_histproj[:,windex,:,1:-1],
                                                  Ngas_bh2_wisedualagn_histproj[:,windex,:,1:-1]],axis=0),
                                     xmin=1, nbins=(len(sep_bins)-1)-2,ecolor='r',
                                     color='r',width=0.5*bw,xoffset=0.1*bw,val='median',errtype=errtype,
                                     verbose=verbose, binpoints_only=True,elinewidth=2,capthick=2)
                    #ng4=cmph.makebar(axb,Ngas_bh2_wisedualagn_histproj[:,windex,:,1:-1],xmin=1,
                    #                 nbins=(len(sep_bins)-1)-2,ecolor='m',
                    #                 color='m',width=0.5*bw,xoffset=0.05*bw,val='median',errtype=errtype,
                    #                 verbose=verbose, binpoints_only=True)
                    axb.legend((ng1,ng2,ng3),('AGN','MIR AGN','Dual MIR AGN'),fontsize=9,loc='upper right',
                               numpoints=1,handletextpad=0.08)
                    
                fig.subplots_adjust(hspace=0.42,wspace=0.08,left=0.14,right=0.96,bottom=0.08,top=0.94)
                plotname='%s/Ngas%s%s_wise_lagn45_merger_phases%s%s_tpost%g_%s%s.pdf'%(outdir,ngastype,rstr,sepstring,binstr,
                                                                                       tmax_postmrg,errtype,extra)
                fig.savefig(plotname,format='pdf')
                plt.cla()
                plt.clf()
                plt.close()



                for i,llim in enumerate((44,45)):
                    ##suptit='This work'
                    #suptit=''
                    if name=='hires':
                        suptit='High-resolution mergers'
                    elif name=='fid_compare_to_hires':
                        suptit='Fiducial mergers (with high-res equivalent)'
                    else: suptit=''

                    if suptit:
                        fig = plt.figure(figsize=(6.5,3.1))
                    else:
                        fig = plt.figure(figsize=(6.5,3))
                    title=r"L$_{\rm AGN} > 10^{%d}$ erg s$^{-1}$"%llim

                    axa = fig.add_subplot(121)
                    axa.set_yscale('log')
                    axa.set_ylim(1.0e21,3.0e24)
                    #axa.set_ylim(3.0e20,3.0e24)
                    axa.set_ylabel(r'N$_{\rm H}$ [cm$^{-2}$]',fontsize=12)
                    axa.set_xlim(xlim_post)
                    axa.set_xticks(xticks_post)
                    axa.set_xticklabels(xticklbl_post,fontsize=10)
                    axa.set_xlabel(xlabel_post,fontsize=11)
                    axa.set_title(title,fontsize=11)
                    ngpost=cmph.makebar(axa,Ngas_agn_hist_postmrg[:,3*2+i+1,:,:],nbins=len(tpostmrg_bins)-1,ecolor='c',
                                    color='c',width=0.5*bw,xoffset=-0.05*bw,val='median',errtype=errtype,
                                    verbose=verbose,binpoints_only=True,reverse=True,elinewidth=1.2,capthick=1.2)
                    ngpost=cmph.makebar(axa,Ngas_wiseagn_mycut_hist_postmrg[:,i,:,:],nbins=len(tpostmrg_bins)-1,ecolor='k',
                                        color='k',width=0.5*bw,xoffset=0.05*bw,val='median',errtype=errtype,
                                        verbose=verbose,binpoints_only=True,reverse=True)
                    
                    axb = fig.add_subplot(122)
                    axb.set_yscale('log')
                    axb.set_ylim(1.0e21,3.0e24)
                    #axb.set_ylim(3.0e20,3.0e24)
                    axb.set_ylabel('')
                    axb.set_xlim(xlim_pre)
                    axb.set_xticks(xticks_pre)
                    axb.set_xticklabels(xticklbl_pre,fontsize=10)
                    axb.set_yticklabels('')
                    axb.set_xlabel(xlabel,fontsize=11)
                    axb.set_title(title,fontsize=11)
                    ng1=cmph.makebar(axb,np.mean([Ngas_bh1_agn_histproj[:,3*2+i+1,:,1:-1],
                                                  Ngas_bh2_agn_histproj[:,3*2+i+1,:,1:-1]],axis=0),
                                     xmin=1,nbins=(len(sep_bins)-1)-2,ecolor='c',
                                     color='c',width=0.5*bw,xoffset=-0.1*bw,val='median',errtype=errtype,
                                     verbose=verbose,binpoints_only=True,elinewidth=1.2,capthick=1.2)
                    ng2=cmph.makebar(axb,np.mean([Ngas_bh1_wiseagn_mycut_histproj[:,i,:,1:-1],
                                                  Ngas_bh2_wiseagn_mycut_histproj[:,i,:,1:-1]],axis=0),
                                     xmin=1,nbins=(len(sep_bins)-1)-2,ecolor='k',
                                     color='k',width=0.5*bw,xoffset=0,val='median',errtype=errtype,
                                     verbose=verbose,binpoints_only=True)
                    ng3=cmph.makebar(axb,np.mean([Ngas_bh1_wisedualagn_mycut_histproj[:,i,:,1:-1],
                                                  Ngas_bh2_wisedualagn_mycut_histproj[:,i,:,1:-1]],axis=0),
                                     xmin=1, nbins=(len(sep_bins)-1)-2,ecolor='r',
                                     color='r',width=0.5*bw,xoffset=0.1*bw,val='median',errtype=errtype,
                                     verbose=verbose, binpoints_only=True,elinewidth=2,capthick=2)
                    axb.legend((ng1,ng2,ng3),('AGN','MIR AGN','Dual MIR AGN'),fontsize=9,loc='upper right',
                               numpoints=1,handletextpad=0.1)
                    #          numpoints=1,handletextpad=0.08,borderpad=0.2)
                    
                    if suptit:
                        fig.suptitle(suptit)
                        fig.subplots_adjust(hspace=0.42,wspace=0.06,left=0.12,right=0.98,bottom=0.18,top=0.83)
                    else:
                        fig.subplots_adjust(hspace=0.42,wspace=0.06,left=0.12,right=0.98,bottom=0.18,top=0.88)
                    #fig.subplots_adjust(hspace=0.42,wspace=0.08,left=0.14,right=0.96,bottom=0.16,top=0.84)
                    ##fig.subplots_adjust(hspace=0.42,wspace=0.08,left=0.14,right=0.96,bottom=0.08,top=0.94)
                    plotname='%s/Ngas%s%s_wise_mycut_lagn%d_merger_phases%s%s_tpost%g_%s%s.pdf'%(outdir,ngastype,rstr,llim,sepstring,binstr,
                                                                                                 tmax_postmrg,errtype,extra)
                    fig.savefig(plotname,format='pdf')
                    plt.cla()
                    plt.clf()
                    plt.close()




            #Ngas_bh1_wiseagn_histproj[i_sim,i_mask,:,:] = np.array([np.histogram(projsepvar[j,tot_wmask[j,:]],bins=sep_bins, 
            #                                                                     weights=d.dt[tot_wmask[j,:]]*Ngas_bh1[j,tot_wmask[j,:]])[0]
            #                                                        for j in range(ncam)])/ttot_projsep_bins
            #Ngas_bh2_wiseagn_histproj[i_sim,i_mask,:,:] = np.array([np.histogram(projsepvar[j,tot_wmask[j,:]],bins=sep_bins, 
            #                                                                     weights=wts[j,tot_wmask[j,:]])[0]
            #                                                        for j in range(ncam)])/ttot_projsep_bins
            #Ngas_wiseagn_hist_postmrg[i_sim,i_mask,:,:] = np.array([np.histogram(d.time[~has2bh][tot_wmask[j,~has2bh]]-d.tmrg,bins=tpostmrg_bins, 
            #                                                                     weights=d.dt[~has2bh][tot_wmask[j,~has2bh]]*Ngas_bh1[j,~has2bh][tot_wmask[j,~has2bh]])[0] 
            #                                                        for j in range(ncam)])/dtbin


            ### Time-weighted WISE dual AGN column density ###
            #Ngas_bh1_wisedualagn_histproj[i_sim,i_mask,:,:] = np.array([np.histogram(projsepvar[j,dual_wmask[j,:]],bins=sep_bins, 
            #                                                                         weights=d.dt[dual_wmask[j,:]]*Ngas_bh1[j,dual_wmask[j,:]])[0]
            #                                                            for j in range(ncam)])/ttot_projsep_bins
            #Ngas_bh2_wisedualagn_histproj[i_sim,i_mask,:,:] = np.array([np.histogram(projsepvar[j,dual_wmask[j,:]],bins=sep_bins, 
            #                                                                     weights=wts[j,dual_wmask[j,:]])[0]
            #                                                        for j in range(ncam)])/ttot_projsep_bins




    ### Ngas plot with no AGN weighting
    fig = plt.figure(figsize=(6.5,2.7))
    axa = fig.add_subplot(121)
    axa.set_yscale('log')
    axa.set_ylim(3.0e20,3.0e24)
    #axa.set_ylabel(r'N$_{\rm gas}$ [cm$^{-2}$]',fontsize=12)
    axa.set_ylabel(r'N$_{\rm H}$ [cm$^{-2}$]',fontsize=12)
    axa.set_xlim(xlim_post)
    axa.set_xticks(xticks_post)
    axa.set_xticklabels(xticklbl_post,fontsize=10)
    axa.set_xlabel(xlabel_post,fontsize=11)
    ngpost=cmph.makebar(axa,Ngas_hist_postmrg,nbins=len(tpostmrg_bins)-1,ecolor='k',
                        color='k',width=0.5*bw,xoffset=0,val='median',errtype=errtype,
                        verbose=verbose,binpoints_only=True,reverse=True)
    
    axb = fig.add_subplot(122)
    axb.set_yscale('log')
    axb.set_ylim(3.0e20,3.0e24)
    axb.set_ylabel('')
    axb.set_xlim(xlim_pre)
    axb.set_xticks(xticks_pre)
    axb.set_xticklabels(xticklbl_pre,fontsize=10)
    axb.set_yticklabels('')
    axb.set_xlabel(xlabel,fontsize=11)
    ng1=cmph.makebar(axb,Ngas_bh1_histproj[:,:,1:-1],xmin=1,
                     nbins=(len(sep_bins)-1)-2,ecolor='b',
                     color='b',width=0.5*bw,xoffset=-0.05*bw,val='median',errtype=errtype,
                     verbose=verbose,binpoints_only=True)
    ng2=cmph.makebar(axb,Ngas_bh2_histproj[:,:,1:-1],xmin=1,
                     nbins=(len(sep_bins)-1)-2,ecolor='m',
                     color='m',width=0.5*bw,xoffset=0.05*bw,val='median',errtype=errtype,
                     verbose=verbose, binpoints_only=True)
    axb.legend((ng1,ng2),('BH1','BH2'),fontsize=9,loc='upper right',numpoints=1)
    fig.suptitle(ngas_title,fontsize=11)

    fig.subplots_adjust(hspace=0.42,wspace=0.08,left=0.14,right=0.96,bottom=0.16,top=0.84)
    plotname='%s/Ngas%s%s_merger_phases%s%s_tpost%g_%s%s.pdf'%(outdir,ngastype,rstr,sepstring,binstr,
                                                             tmax_postmrg,errtype,extra)
    fig.savefig(plotname,format='pdf')
    plt.cla()
    plt.clf()

   
    
def merger_phases_agnx0compare(maindir='/oasis/projects/nsf/hvd115/lblecha/',subdir_arr='test', 
                               z=0,tmax_postmrg=0.1, latesep=10, xlim=(), ylim=(), 
                               skip_snap0=True, equal_bins=True,alt_wedge=True,
                               dsep=10,cam=range(7),dust_to_metal_ratio=0.4, plot_3dsep=False,
                               use_logbins=False, ylog=False,extra='',errtype='mad',verbose=False,
                               multisim_plots=True,sfrhist_w_grid=False,grid_nh_res=0,valtype='median',
                               plot_wedge_allsim=False,wise_fluxlims=[0.0, 0.0, 0.0, 0.0],wise_fluxlim_type='',
                               omega_m=0.308, h=0.678,contam_file_only=False):


    wflim_str = ''
    if z>0:
        omega_l = 1.0-omega_m
        cosmo = pc.Cosmology(omega_m, omega_l, 0.0, -1.0, h)
        dL = cosmo.lum_dist(z) * 1000.0*ac.KPC ## Mpc to cm

        if wise_fluxlim_type in ('SN3','SN10'):
            print 'Imposing cuts for WISE sensitivity limits for W1 & W2, assuming %s'%wise_fluxlim_type
            if wise_fluxlim_type=='SN3':
                ### total band fluxes (haven't coded for redshifted SEDs yet)
                #W1fluxlim = 2.157e-15 ## erg/s/cm^2
                #W2fluxlim = 3.605e-15 ## erg/s/cm^2
                ### monochromatic fluxes F_lambda 
                W1fluxlim = 3.256e-16 ## W/m/cm^2
                W2fluxlim = 3.45873e-16 ## W/m/cm^2
            elif wise_fluxlim_type=='SN10':
                ### total band fluxes (haven't coded for redshifted SEDs yet)
                #W1fluxlim = 7.690e-15 ## erg/s/cm^2
                #W2fluxlim = 1.285e-14 ## erg/s/cm^2
                ### monochromatic fluxes F_lambda 
                W1fluxlim = 1.1606e-15 ## W/m/cm^2
                W2fluxlim = 1.23287e-15 ## W/m/cm^2
            w1lum_lim=W1fluxlim * 4*np.pi*dL*dL ##W/m (sunrise units)
            w2lum_lim=W2fluxlim * 4*np.pi*dL*dL 
            wflim_str = wise_fluxlim_type+'_'
            
        elif not wise_fluxlim_type and np.max(wise_fluxlims)>0.0:
            print 'Imposing cuts based on input minimum WISE band fluxes:'
            print wise_fluxlims
            w1lum_lim = wise_fluxlims[0] * 4*np.pi*dL*dL ## W/m (sunrise units)
            w2lum_lim = wise_fluxlims[1] * 4*np.pi*dL*dL
            wflim_str = 'man_wflims_'

        else:
            print "No cuts imposed on minimum WISE band fluxes."
            w1lum_lim=0.0
            w2lum_lim=0.0

        print "z = %g"%z
        print "dL = %g cm"%dL
        print "w1lum_lim = %g"%w1lum_lim
        print "w2lum_lim = %g"%w2lum_lim
        #return
    else:
        if wise_fluxlim_type or np.max(wise_fluxlims)>0.0:
            print "Error: cannot set WISE fluxlims for rest-frame (z=0) calculations."
            return
        w1lum_lim=0.0
        w2lum_lim=0.0


    ncam=len(cam)
    print "plot_3dsep?", plot_3dsep

    if z==0:
        path_arr = define_sunrise_path_arr(basepath=maindir,name=subdir_arr,sfrhist_w_grid=sfrhist_w_grid)
        path_0_arr = define_sunrise_path_arr(basepath=maindir,name=subdir_arr+'_agnx0',sfrhist_w_grid=sfrhist_w_grid)
        name=subdir_arr[0] if isinstance(subdir_arr,list) else subdir_arr
    else:
        path_arr = define_sunrise_path_arr_z(basepath=maindir,name=subdir_arr,z=z)
        path_0_arr = define_sunrise_path_arr_z(basepath=maindir,name=subdir_arr,z=z,agnx0=True)
        name="%s_z%.1f"%(subdir_arr[0],z) if isinstance(subdir_arr,list) else "%s_z%.1f"%(subdir_arr,z)

    print path_arr
    print path_0_arr
    #if len(path_arr)>1 or len(path_0_arr)>1:
    #    assert len(path_arr)==len(path_0_arr),'Error! mismatch in agnx1 and agnx0 dir list.'
    nsim = len(path_arr)
    #extra = '_%dsims'%nsim if nsim>1 else ''
    extra=''
    #outdir = maindir if nsim>1 else subdir_arr[0]
    outdir=maindir

    ### set limits for agn masks ###
    fedd_lims = [0.01,0.05,0.1]
    flbol_lims = [0.1, 0.3, 0.5]
    lagn_lims = [1.0e43, 1.0e44, 1.0e45]
    nagnmask = len(fedd_lims)+len(lagn_lims)+(len(flbol_lims))
    nwagnmask = nagnmask * 3

    w12_lims_hires = [0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
    w23_lims_hires = [2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5]
    lagn_lims_hires = [10**x for x in [42.5,43.0,43.5,44.0,44.5,45.0,45.5,46.0]]
    
    sepstr = '_bhsep' if plot_3dsep else '_bhprojsep'
    lgsepstr = '_lgbhsep' if plot_3dsep else '_lgbhprojsep'
    sepstring = lgsepstr if use_logbins else sepstr
    binstr = '_eqbins%d'%dsep if equal_bins else ''


    fedd_plot_titles = [r'f$_{\rm Edd}$>%g'%(lim) for lim in fedd_lims]
    lagn_plot_titles = [r'L$_{\rm AGN}$>10$^{%g}$ [erg s$^{-1}$]'%(np.log10(lim))  for lim in lagn_lims]
    flbol_plot_titles = [r'L$_{\rm AGN}$/L$_{\rm tot}$>%g'%(lim) for lim in flbol_lims]
    agn_plot_titles = plot_titles = fedd_plot_titles + lagn_plot_titles + flbol_plot_titles
    fedd_labels = ['fEdd>%g'%(lim) for lim in fedd_lims]
    lagn_labels = ['Lagn>%g'%(np.log10(lim))  for lim in lagn_lims]
    flbol_labels = ['fLbol>%g'%(lim) for lim in flbol_lims]
    agn_labels = fedd_labels + lagn_labels + flbol_labels
    
    wisemask_labels = ('W1W2>0.5','W1W2>0.8','J11 wedge')
    fedd_wise_titles = ['%s, fEdd>%g'%(lbl,lim) for lbl in wisemask_labels for lim in fedd_lims]
    lagn_wise_titles = ['%s, Lagn>%g'%(lbl,lim) for lbl in wisemask_labels for lim in lagn_lims]
    flbol_wise_titles = ['%s, fLbol>%g'%(lbl,lim) for lbl in wisemask_labels for lim in flbol_lims]
    wise_agn_titles = fedd_wise_titles + lagn_wise_titles + flbol_wise_titles
    wise_selection_labels = ['MIR AGN', 'MIR false neg', 'MIR false pos','']
    wise_duals_labels = ['MIR dual AGN', 'dual AGN, not MIR', 'MIR AGN, not dual','']
    print "wise_agn_titles:",wise_agn_titles
    
    ix_fedd = np.array([k for k in range(nwagnmask) if 'fEdd' in wise_agn_titles[k]])
    ix_lagn = np.array([k for k in range(nwagnmask) if 'Lagn' in wise_agn_titles[k]])
    ix_flbol = np.array([k for k in range(nwagnmask) if 'fLbol' in wise_agn_titles[k]])
    ix_w05 = np.array([k for k in range(nwagnmask) if 'W1W2>0.5' in wise_agn_titles[k]])
    ix_w08 = np.array([k for k in range(nwagnmask) if 'W1W2>0.8' in wise_agn_titles[k]])
    ix_j11 = np.array([k for k in range(nwagnmask) if 'J11' in wise_agn_titles[k]])
    wfbase = ('w05','w08','j11')
    

    #wagn_global0 = np.zeros((nsim,nwagnmask))
    wagn_global0 = np.zeros((nsim,nwagnmask))

    f_wise_tot = np.copy(wagn_global0)
    f_wise_early = np.copy(wagn_global0)
    f_wise_late = np.copy(wagn_global0)
    f_wise_post = np.copy(wagn_global0)

    f_wise_tot_0 = np.copy(wagn_global0)
    f_wise_early_0 = np.copy(wagn_global0)
    f_wise_late_0 = np.copy(wagn_global0)
    f_wise_post_0 = np.copy(wagn_global0)

    f_wise_agnx0x1_tot = np.copy(wagn_global0)
    f_wise_agnx0x1_early = np.copy(wagn_global0)
    f_wise_agnx0x1_late = np.copy(wagn_global0)
    f_wise_agnx0x1_post = np.copy(wagn_global0)

    f_wise_agnx0x1_noagn_tot = np.copy(wagn_global0)
    f_wise_agnx0x1_noagn_early = np.copy(wagn_global0)
    f_wise_agnx0x1_noagn_late = np.copy(wagn_global0)
    f_wise_agnx0x1_noagn_post = np.copy(wagn_global0)

    tw_hires = np.zeros((nsim,len(w12_lims_hires),ncam))
    tw_hires_late_post = np.zeros((nsim,len(w12_lims_hires),ncam))
    tw_hires_0 = np.zeros((nsim,len(w12_lims_hires),ncam))
    tw_hires_late_post_0 = np.zeros((nsim,len(w12_lims_hires),ncam))
    tw_hires_noagn_0 = np.zeros((nsim,len(w12_lims_hires),ncam))

    tj11 = np.zeros((nsim,ncam))
    tj11_late_post = np.zeros((nsim,ncam))
    tj11_0 = np.zeros((nsim,ncam))
    tj11_late_post_0 = np.zeros((nsim,ncam))
    tj11_noagn_l44_0 = np.zeros((nsim,ncam))

    tmycut = np.zeros((nsim,ncam))
    tmycut_late_post = np.zeros((nsim,ncam))
    tmycut_0 = np.zeros((nsim,ncam))
    tmycut_late_post_0 = np.zeros((nsim,ncam))
    tmycut_noagn_l44_0 = np.zeros((nsim,ncam))

    sfr_max = np.zeros((nsim))
    ssfr_max = np.zeros((nsim))
    fedd_max = np.zeros((nsim))
    lagn_max = np.zeros((nsim))
    ltot_max = np.zeros((nsim))
    flbol_max = np.zeros((nsim))
    mstar_final = np.zeros((nsim))

    xtypes = np.array(['lgssfr','lgsfr','lglagn','lgltot','flbol'])   

    ### load & process agn metadata ###
    for i_sim,path in enumerate(path_arr):

        print "\nloading agn metadata for sim %d (%s)..."%(i_sim,path)
        d = agn_metadata(path, ncam=ncam, skip_snap0=skip_snap0,
                         tmax_postmrg = tmax_postmrg, agnx0=False, skip_nh=True)
        print d.dt.shape, d.posbh1.shape, d.lbol_grid.shape, d.bhprojsep.shape, d.w2w3.shape, d.nsnaps

        print "\nloading agn metadata for sim %d (%s)..."%(i_sim,path_0_arr[i_sim])
        #d0 = agn_metadata(path_0_arr[i_sim], ncam=ncam, skip_snap0=skip_snap0,
        #                  tmax_postmrg = tmax_postmrg, agnx0=True, grid_nh_res=grid_nh_res)
        d0 = agn_metadata(path_0_arr[i_sim], ncam=ncam, skip_snap0=skip_snap0,
                          tmax_postmrg = tmax_postmrg, agnx0=True, skip_nh=True)
        print d0.dt.shape, d0.posbh1.shape, d0.lbol_grid.shape, d0.bhprojsep.shape, d0.w2w3.shape, d0.nsnaps

        mask = np.in1d(d.snap,d0.snap)
        snap_nomatch=d.snap[~mask]

        mask0 = np.in1d(d0.snap,d.snap)
        snap0_nomatch=d0.snap[~mask0]

        if len(snap_nomatch)>0: 
            print "agnx1 nsnaps (raw): %d"%d.snap.size
            print "snap_nomatch:",snap_nomatch
            mask[d.snap==snap_nomatch]=False
            d.mask_snaps(mask,skip_nh=True)
            print "agnx1 nsnaps (retained): %d"%d.snap.size

        if len(snap0_nomatch)>0: 
            print "agnx0 nsnaps (raw): %d"%d0.snap.size
            print "snap0_nomatch:",snap0_nomatch
            mask0[d0.snap==snap0_nomatch]=False
            d0.mask_snaps(mask0,skip_nh=True)
            print "agnx0 nsnaps (retained): %d"%d0.snap.size
        

        ### AGNx1 ###
    

        has2bh = (d.nbh==2) 
        ttot = d.dt.sum()

        sfr_max[i_sim] = np.max(d.sfr)
        ssfr_max[i_sim] = np.max(d.sfr/d.mstar)
        fedd_max[i_sim] = np.max(d.fedd)
        lagn_max[i_sim] = np.max(d.lagn)
        ltot_max[i_sim] = np.max(d.lbol_grid)
        flbol_max[i_sim] = np.max(d.fagn_lbol)
        mstar_final[i_sim] = d.mstar[-1]

        meanw12 = np.mean(d.w1w2,axis=0)
        meanw23 = np.mean(d.w2w3,axis=0)

        j11mask = j11_wedge_cut(d.w1w2,d.w2w3,
                                w1lum_lim=w1lum_lim,w2lum_lim=w2lum_lim,w1lum=d.w1lum,w2lum=d.w2lum)
        j11mask_mean = j11_wedge_cut(meanw12,meanw23,
                                     w1lum_lim=w1lum_lim,w2lum_lim=w2lum_lim,w1lum=d.w1lum,w2lum=d.w2lum)

        mycutmask = my_wedge_cut(d.w1w2,d.w2w3,alt=alt_wedge,
                                 w1lum_lim=w1lum_lim,w2lum_lim=w2lum_lim,w1lum=d.w1lum,w2lum=d.w2lum)
        mycutmask_mean = my_wedge_cut(meanw12,meanw23,alt=alt_wedge,
                                      w1lum_lim=w1lum_lim,w2lum_lim=w2lum_lim,w1lum=d.w1lum,w2lum=d.w2lum)

        w05mask = (d.w1w2>0.5)
        w05mask_mean = (meanw12>0.5)
        w08mask = (d.w1w2>0.8)
        w08mask_mean = (meanw12>0.8)

        fedd_maskarr,dual_fedd_maskarr = set_agn_maskarr_list(d.fedd,d.fedd1,d.fedd2,fedd_lims)
        lagn_maskarr,dual_lagn_maskarr = set_agn_maskarr_list(d.lagn,d.lagn1,d.lagn2,lagn_lims)
        flbol_maskarr,dual_flbol_maskarr = set_agn_maskarr_list(d.fagn_lbol,d.fagn1_lbol,
                                                                d.fagn2_lbol,flbol_lims)
        agn_maskarr_list = fedd_maskarr + lagn_maskarr + flbol_maskarr
        dualagn_maskarr_list = dual_fedd_maskarr + dual_lagn_maskarr + dual_flbol_maskarr

        avg_wisemask_list = [w05mask_mean, w08mask_mean, j11mask_mean]
        wisemask_list = [w05mask, w08mask, j11mask]
        assert len(wisemask_list)==3, 'Error: wise agn arrays need wisemask_list of len 3, not %d'%len(wisemask_list)

        if w1lum_lim>0.0 or w2lum_lim>0.0:
            fedd_wmaskarr_list,dual_fedd_wmaskarr_list = set_agn_wmaskarr_list(d.fedd,d.fedd1,d.fedd2,fedd_lims,wisemask_list,
                                                                               w1lum=d.w1lum, w2lum=d.w2lum, w1lum_lim=w1lum_lim,
                                                                               w2lum_lim=w2lum_lim )
            flbol_wmaskarr_list,dual_flbol_wmaskarr_list = set_agn_wmaskarr_list(d.fagn_lbol,d.fagn1_lbol,d.fagn2_lbol,
                                                                                 flbol_lims,wisemask_list,
                                                                                 w1lum=d.w1lum, w2lum=d.w2lum, w1lum_lim=w1lum_lim,
                                                                                 w2lum_lim=w2lum_lim )
            lagn_wmaskarr_list,dual_lagn_wmaskarr_list = set_agn_wmaskarr_list(d.lagn,d.lagn1,d.lagn2,lagn_lims,wisemask_list,
                                                                               w1lum=d.w1lum, w2lum=d.w2lum, w1lum_lim=w1lum_lim,
                                                                               w2lum_lim=w2lum_lim )

        else:
            fedd_wmaskarr_list,dual_fedd_wmaskarr_list = set_agn_wmaskarr_list(d.fedd,d.fedd1,d.fedd2,fedd_lims,wisemask_list)
            flbol_wmaskarr_list,dual_flbol_wmaskarr_list = set_agn_wmaskarr_list(d.fagn_lbol,d.fagn1_lbol,d.fagn2_lbol,
                                                                                 flbol_lims,wisemask_list)
            lagn_wmaskarr_list,dual_lagn_wmaskarr_list = set_agn_wmaskarr_list(d.lagn,d.lagn1,d.lagn2,lagn_lims,wisemask_list)

        agn_wmaskarr_list = fedd_wmaskarr_list + lagn_wmaskarr_list + flbol_wmaskarr_list
        dualagn_wmaskarr_list = dual_fedd_wmaskarr_list + dual_lagn_wmaskarr_list + dual_flbol_wmaskarr_list

        wisemask_hires_list = [(d.w1w2>wlim) for wlim in w12_lims_hires]

        if w1lum_lim>0.0 or w2lum_lim>0.0:
            l44_wmaskarr_hires,dual_l44_wmaskarr_hires = set_agn_wmaskarr_list(d.lagn,d.lagn1,d.lagn2,[1.0e44],wisemask_hires_list,
                                                                               w1lum=d.w1lum, w2lum=d.w2lum, w1lum_lim=w1lum_lim,
                                                                               w2lum_lim=w2lum_lim )
        else:
            l44_wmaskarr_hires,dual_l44_wmaskarr_hires = set_agn_wmaskarr_list(d.lagn,d.lagn1,d.lagn2,[1.0e44],wisemask_hires_list)

        ttot = d.dt.sum()
        t_early = np.array([d.dt[d.bhprojsep[i,:]>=latesep].sum() for i in range(ncam)])
        t_late = np.array([d.dt[(d.bhprojsep[i,:]<latesep)&
                                (d.bhprojsep[i,:]>0)].sum() for i in range(ncam)])
        t_post = np.array([d.dt[d.bhprojsep[i,:]<=0].sum() for i in range(ncam)])
        print "\nt_early:",t_early
        print "\nt_late:",t_late
        print "\nt_post:",t_post
        print "\nttot:",ttot
        print "max proj. sep for each cam:",[d.bhprojsep[i,:].max() for i in range(ncam)]
        assert t_early.all()>0 and t_late.all()>0 and t_post.all()>0,'Error: at least one projection has a merger phase with t=0. Consider a different value of latesep (%g) or tmax_postmrg (%g).'%(latesep,tmax_postmrg)

        assert np.allclose(t_early+t_late+t_post,ttot),'Error: time in merger phases !=ttot'
        nproj_tnonzero_early = t_early[t_early>0].size
        nproj_tnonzero_late = t_late[t_late>0].size
        nproj_tnonzero_post = t_post[t_post>0].size
        print "nproj_tnonzero_[early,late,post]: %d %d %d"%(nproj_tnonzero_early,nproj_tnonzero_late,nproj_tnonzero_post)

        ### AGNx0 ###

        has2bh_0 = (d0.nbh==2) 
        ttot_0 = d0.dt.sum()

        meanw12_0 = np.mean(d0.w1w2,axis=0)
        meanw23_0 = np.mean(d0.w2w3,axis=0)

        j11mask_0 = j11_wedge_cut(d0.w1w2,d0.w2w3,
                                  w1lum_lim=w1lum_lim,w2lum_lim=w2lum_lim,w1lum=d0.w1lum,w2lum=d0.w2lum)
        j11mask_mean_0 = j11_wedge_cut(meanw12_0,meanw23_0,
                                       w1lum_lim=w1lum_lim,w2lum_lim=w2lum_lim,w1lum=d0.w1lum,w2lum=d0.w2lum)

        mycutmask_0 = my_wedge_cut(d0.w1w2,d0.w2w3,alt=alt_wedge,
                                   w1lum_lim=w1lum_lim,w2lum_lim=w2lum_lim,w1lum=d0.w1lum,w2lum=d0.w2lum)
        mycutmask_mean_0 = my_wedge_cut(meanw12_0,meanw23_0,alt=alt_wedge,
                                        w1lum_lim=w1lum_lim,w2lum_lim=w2lum_lim,w1lum=d0.w1lum,w2lum=d0.w2lum)


        w05mask_0 = (d0.w1w2>0.5)
        w05mask_mean_0 = (meanw12_0>0.5)
        w08mask_0 = (d0.w1w2>0.8)
        w08mask_mean_0 = (meanw12_0>0.8)

        avg_wisemask_list_0 = [w05mask_mean_0, w08mask_mean_0, j11mask_mean_0]
        wisemask_list_0 = [w05mask_0, w08mask_0, j11mask_0]
        assert len(wisemask_list_0)==3, 'Error: wise agn arrays need wisemask_list of len 3, not %d'%len(wisemask_list_0)

        if w1lum_lim>0.0 or w2lum_lim>0.0:
            fedd_wmaskarr_list_0,dual_fedd_wmaskarr_list_0 = set_agn_wmaskarr_list(d0.fedd,d0.fedd1,d0.fedd2,fedd_lims,wisemask_list_0,
                                                                               w1lum=d.w1lum, w2lum=d.w2lum, w1lum_lim=w1lum_lim,
                                                                               w2lum_lim=w2lum_lim )
            flbol_wmaskarr_list_0,dual_flbol_wmaskarr_list_0 = set_agn_wmaskarr_list(d0.fagn_lbol,d0.fagn1_lbol,d0.fagn2_lbol,
                                                                                     flbol_lims,wisemask_list_0,
                                                                                     w1lum=d.w1lum, w2lum=d.w2lum, w1lum_lim=w1lum_lim,
                                                                                     w2lum_lim=w2lum_lim )
            lagn_wmaskarr_list_0,dual_lagn_wmaskarr_list_0 = set_agn_wmaskarr_list(d0.lagn,d0.lagn1,d0.lagn2,lagn_lims,wisemask_list_0,
                                                                                   w1lum=d.w1lum, w2lum=d.w2lum, w1lum_lim=w1lum_lim,
                                                                                   w2lum_lim=w2lum_lim )
        else:
            fedd_wmaskarr_list_0,dual_fedd_wmaskarr_list_0 = set_agn_wmaskarr_list(d0.fedd,d0.fedd1,d0.fedd2,fedd_lims,wisemask_list_0)
            flbol_wmaskarr_list_0,dual_flbol_wmaskarr_list_0 = set_agn_wmaskarr_list(d0.fagn_lbol,d0.fagn1_lbol,d0.fagn2_lbol,
                                                                                     flbol_lims,wisemask_list_0)
            lagn_wmaskarr_list_0,dual_lagn_wmaskarr_list_0 = set_agn_wmaskarr_list(d0.lagn,d0.lagn1,d0.lagn2,lagn_lims,wisemask_list_0)

        agn_wmaskarr_list_0 = fedd_wmaskarr_list_0 + lagn_wmaskarr_list_0 + flbol_wmaskarr_list_0
        dualagn_wmaskarr_list_0 = dual_fedd_wmaskarr_list_0 + dual_lagn_wmaskarr_list_0 + dual_flbol_wmaskarr_list_0

        wisemask_hires_list_0 = [(d0.w1w2>wlim) for wlim in w12_lims_hires]
        if w1lum_lim>0.0 or w2lum_lim>0.0:
            l44_wmaskarr_hires_0,dual_l44_wmaskarr_hires_0 = set_agn_wmaskarr_list(d0.lagn,d0.lagn1,d0.lagn2,[1.0e44],
                                                                                   wisemask_hires_list_0,
                                                                                   w1lum=d0.w1lum, w2lum=d0.w2lum, 
                                                                                   w1lum_lim=w1lum_lim,w2lum_lim=w2lum_lim )
        else:
            l44_wmaskarr_hires_0,dual_l44_wmaskarr_hires_0 = set_agn_wmaskarr_list(d0.lagn,d0.lagn1,d0.lagn2,[1.0e44],
                                                                                   wisemask_hires_list_0)

        ttot_0 = d0.dt.sum()
        t_early_0 = np.array([d0.dt[d0.bhprojsep[i,:]>=latesep].sum() for i in range(ncam)])
        t_late_0 = np.array([d0.dt[(d0.bhprojsep[i,:]<latesep)&
                                (d0.bhprojsep[i,:]>0)].sum() for i in range(ncam)])
        t_post_0 = np.array([d0.dt[d0.bhprojsep[i,:]<=0].sum() for i in range(ncam)])
        print "\nt_early_0:",t_early_0
        print "\nt_late_0:",t_late_0
        print "\nt_post_0:",t_post_0
        print "\nttot_0:",ttot_0
        print "max proj. sep for each cam:",[d.bhprojsep[i,:].max() for i in range(ncam)]
        assert t_early.all()>0 and t_late.all()>0 and t_post.all()>0,'Error: at least one projection has a merger phase with t=0. Consider a different value of latesep (%g) or tmax_postmrg (%g).'%(latesep,tmax_postmrg)

        assert np.allclose(t_early_0+t_late_0+t_post_0,ttot_0),'Error: time in merger phases !=ttot_0'
        nproj_tnonzero_early_0 = t_early_0[t_early_0>0].size
        nproj_tnonzero_late_0 = t_late_0[t_late_0>0].size
        nproj_tnonzero_post_0 = t_post_0[t_post_0>0].size
        print "nproj_tnonzero_[early,late,post]_0: %d %d %d"%(nproj_tnonzero_early_0,nproj_tnonzero_late_0,nproj_tnonzero_post_0)

        for i_mask in np.arange(len(agn_wmaskarr_list)):
            ### recall:
            ### agn_wmaskarr_list = fedd_wmaskarr_list + lagn_wmaskarr_list + flbol_wmaskarr_list
            ### dualagn_wmaskarr_list = dual_fedd_wmaskarr_list + dual_lagn_wmaskarr_list + dual_flbol_wmaskarr_list

            #print "\n%s:\n"%wise_agn_titles[i_mask]

            tot_wmask = agn_wmaskarr_list[i_mask][0]
            tot_nowmask = agn_wmaskarr_list[i_mask][1]
            notot_wmask = agn_wmaskarr_list[i_mask][2]
            notot_nowmask = agn_wmaskarr_list[i_mask][3]
            dual_wmask = dualagn_wmaskarr_list[i_mask][0]
            dual_nowmask = dualagn_wmaskarr_list[i_mask][1]
            one_wmask = dualagn_wmaskarr_list[i_mask][2]
            one_nowmask = dualagn_wmaskarr_list[i_mask][3]
            no_wmask = dualagn_wmaskarr_list[i_mask][4]
            no_nowmask = dualagn_wmaskarr_list[i_mask][5]

            notot_wmask_0 = agn_wmaskarr_list_0[i_mask][2]
            notot_nowmask_0 = agn_wmaskarr_list_0[i_mask][3]
            no_wmask_0 = dualagn_wmaskarr_list_0[i_mask][4]
            no_nowmask_0 = dualagn_wmaskarr_list_0[i_mask][5]

            tw = np.array([d.dt[((tot_wmask[j,:])|(notot_wmask[j,:]))].sum() for j in range(ncam)])
            f_wise_tot[i_sim,i_mask] = np.sum( tw / ttot ) / (1.0*ncam)

            tw_early = np.array([d.dt[(d.bhprojsep[j,:]>=latesep)&
                                      ((tot_wmask[j,:])|(notot_wmask[j,:]))].sum() for j in range(ncam)])
            f_wise_early[i_sim,i_mask] = np.nansum( tw_early / t_early ) / (1.0*nproj_tnonzero_early)

            tw_late = np.array([d.dt[(d.bhprojsep[j,:]<latesep)&(d.bhprojsep[j,:]>0)&
                                     ((tot_wmask[j,:])|(notot_wmask[j,:]))].sum() for j in range(ncam)])
            f_wise_late[i_sim,i_mask] = np.nansum( tw_late / t_late ) / (1.0*nproj_tnonzero_late)

            tw_post = np.array([d.dt[(d.bhprojsep[j,:]<=0)&((tot_wmask[j,:])|(notot_wmask[j,:]))].sum() 
                                for j in range(ncam)])
            f_wise_post[i_sim,i_mask] = np.nansum( tw_post / t_post ) / (1.0*nproj_tnonzero_post)

            tw_0 = np.array([d0.dt[(notot_wmask_0[j,:])].sum() for j in range(ncam)])
            f_wise_tot_0[i_sim,i_mask] = np.sum( tw_0 / ttot_0 ) / (1.0*ncam)

            tw_early_0 = np.array([d0.dt[(d0.bhprojsep[j,:]>=latesep)&(notot_wmask_0[j,:])].sum() 
                                   for j in range(ncam)])
            f_wise_early_0[i_sim,i_mask] = np.nansum( tw_early_0 / t_early_0 ) / (1.0*nproj_tnonzero_early_0)

            tw_late_0 = np.array([d0.dt[(d0.bhprojsep[j,:]<latesep)&(d0.bhprojsep[j,:]>0)&
                                        (notot_wmask_0[j,:])].sum() for j in range(ncam)])
            f_wise_late_0[i_sim,i_mask] = np.nansum( tw_late_0 / t_late_0 ) / (1.0*nproj_tnonzero_late_0)

            tw_post_0 = np.array([d0.dt[(d0.bhprojsep[j,:]<=0)&(notot_wmask_0[j,:])].sum() 
                                  for j in range(ncam)])
            f_wise_post_0[i_sim,i_mask] = np.nansum( tw_post_0 / t_post_0 ) / (1.0*nproj_tnonzero_post_0)

            f_wise_agnx0x1_tot[i_sim,i_mask] = np.nanmean( tw_0/tw ) 
            f_wise_agnx0x1_early[i_sim,i_mask] = np.nanmean( tw_early_0/tw_early )
            f_wise_agnx0x1_late[i_sim,i_mask] = np.nanmean( tw_late_0/tw_late ) 
            f_wise_agnx0x1_post[i_sim,i_mask] = np.nanmean( tw_post_0/tw_post ) 

            
            tw_noagn_0 = np.array([d0.dt[(notot_wmask_0[j,:])&
                                         ((notot_wmask[j,:])|(notot_nowmask[j,:]))].sum() for j in range(ncam)])
            f_wise_agnx0x1_noagn_tot[i_sim,i_mask] = np.nanmean( tw_noagn_0/tw ) 
            

            tw_noagn_early_0 = np.array([d0.dt[(d0.bhprojsep[j,:]>=latesep)&(notot_wmask_0[j,:])&
                                               ((notot_wmask[j,:])|(notot_nowmask[j,:]))].sum() for j in range(ncam)])
            f_wise_agnx0x1_noagn_early[i_sim,i_mask] = np.nanmean( tw_noagn_early_0/tw_early ) 

            tw_noagn_late_0 = np.array([d0.dt[(d0.bhprojsep[j,:]<latesep)&(d0.bhprojsep[j,:]>0)&
                                              (notot_wmask_0[j,:])&
                                              ((notot_wmask[j,:])|(notot_nowmask[j,:]))].sum() for j in range(ncam)])
            f_wise_agnx0x1_noagn_late[i_sim,i_mask] = np.nanmean( tw_noagn_late_0/tw_late ) 

            tw_noagn_post_0 = np.array([d0.dt[(d0.bhprojsep[j,:]<=0)&(notot_wmask_0[j,:])&
                                              ((notot_wmask[j,:])|(notot_nowmask[j,:]))].sum() 
                                  for j in range(ncam)])
            f_wise_agnx0x1_noagn_post[i_sim,i_mask] = np.nanmean( tw_noagn_post_0/tw_post ) 
            if tw_late_0.max()>0:
                print "t_late: ",t_late
                print "t_late_0: ",t_late_0
                print "tw_late: ",tw_late
                print "tw_late_0: ",tw_late_0
                print "tw_noagn_late_0: ",tw_noagn_late_0
                print "f_wise_agnx0x1_late: ",f_wise_agnx0x1_late[i_sim,i_mask]
                print "f_wise_agnx0x1_noagn_late: ",f_wise_agnx0x1_noagn_late[i_sim,i_mask]


        ### end loop over wmaskarr ###

        if tw_noagn_0.max() > 0:
            print "WARNING!! tw_noagn_0.max() > 0!"
            print "tot, late, post:"
            print tw_noagn_0.max(),tw_noagn_late_0.max(),tw_noagn_post_0.max()
            #return

        ### Loop over hi-res wise mask arrays: ###
        agn_wmaskarr_hires = copy(l44_wmaskarr_hires)
        dualagn_wmaskarr_hires = copy(dual_l44_wmaskarr_hires)
        agn_wmaskarr_hires_0 = copy(l44_wmaskarr_hires_0)
        dualagn_wmaskarr_hires_0 = copy(dual_l44_wmaskarr_hires_0)

        print len(agn_wmaskarr_hires)
        for i_wmask in range(len(agn_wmaskarr_hires)):
        
            tot_wmask_hires = agn_wmaskarr_hires[i_wmask][0]
            tot_nowmask_hires = agn_wmaskarr_hires[i_wmask][1]
            notot_wmask_hires = agn_wmaskarr_hires[i_wmask][2]
            notot_nowmask_hires = agn_wmaskarr_hires[i_wmask][3]
            dual_wmask_hires = dualagn_wmaskarr_hires[i_wmask][0]
            dual_nowmask_hires = dualagn_wmaskarr_hires[i_wmask][1]
            one_wmask_hires = dualagn_wmaskarr_hires[i_wmask][2]
            one_nowmask_hires = dualagn_wmaskarr_hires[i_wmask][3]
            no_wmask_hires = dualagn_wmaskarr_hires[i_wmask][4]
            no_nowmask_hires = dualagn_wmaskarr_hires[i_wmask][5]

            notot_wmask_hires_0 = agn_wmaskarr_hires_0[i_wmask][2]
            notot_nowmask_hires_0 = agn_wmaskarr_hires_0[i_wmask][3]
            no_wmask_hires_0 = dualagn_wmaskarr_hires_0[i_wmask][4]
            no_nowmask_hires_0 = dualagn_wmaskarr_hires_0[i_wmask][5]

            tw_hires[i_sim,i_wmask,:] = np.array([d.dt[((tot_wmask_hires[j,:])|(notot_wmask_hires[j,:]))].sum() 
                                                  for j in range(ncam)])

            tw_hires_late = np.array([d.dt[(d.bhprojsep[j,:]<latesep)&(d.bhprojsep[j,:]>0)&
                                     ((tot_wmask_hires[j,:])|(notot_wmask_hires[j,:]))].sum() for j in range(ncam)])
            #f_wise_late[i_sim,i_mask] = np.nansum( tw_late / t_late ) / (1.0*nproj_tnonzero_late)

            tw_hires_post = np.array([d.dt[(d.bhprojsep[j,:]<=0)&((tot_wmask_hires[j,:])|(notot_wmask_hires[j,:]))].sum() 
                                for j in range(ncam)])
            #f_wise_post[i_sim,i_mask] = np.nansum( tw_post / t_post ) / (1.0*nproj_tnonzero_post)
            tw_hires_late_post[i_sim,i_wmask,:] = tw_hires_late + tw_hires_post
            

            tw_hires_0[i_sim,i_wmask,:] = np.array([d0.dt[(notot_wmask_hires_0[j,:])].sum() for j in range(ncam)])
            #f_wise_agnx0x1_tot[i_sim,i_mask] = np.nanmean( tw_0/tw ) 

            #f_wise_tot[i_sim,i_mask] = np.sum( tw / ttot ) / (1.0*ncam)

            #f1or2agn_wise_hires_tot[i_sim,i_wmask] = np.mean(t_1or2agn_wise_hires_tot) / ttot

            #f_wise_agnx0x1_late[i_sim,i_mask] = np.nanmean( tw_late_0/tw_late ) 
            #f_wise_agnx0x1_post[i_sim,i_mask] = np.nanmean( tw_post_0/tw_post ) 

            ### TESTING ###
            #tw_hires_noagn_0[i_sim,i_wmask,:] = np.array([d0.dt[(notot_wmask_hires_0[j,:])&(notot_wmask_hires[j,:])].sum() 
            #                                              for j in range(ncam)])
            tw_hires_noagn_0[i_sim,i_wmask,:] = np.array([d0.dt[(notot_wmask_hires_0[j,:])&((notot_wmask_hires[j,:])|
                                                                                            (notot_nowmask_hires[j,:]))].sum() 
                                                          for j in range(ncam)])
            #f_wise_agnx0x1_noagn_tot[i_sim,i_mask] = np.nanmean( tw_noagn_0/tw ) 


            tw_hires_late_0 = np.array([d0.dt[(d0.bhprojsep[j,:]<latesep)&(d0.bhprojsep[j,:]>0)&
                                        (notot_wmask_hires_0[j,:])].sum() for j in range(ncam)])
            #f_wise_late_0[i_sim,i_mask] = np.nansum( tw_late_0 / t_late_0 ) / (1.0*nproj_tnonzero_late_0)

            tw_hires_post_0 = np.array([d0.dt[(d0.bhprojsep[j,:]<=0)&(notot_wmask_hires_0[j,:])].sum() 
                                  for j in range(ncam)])
            #f_wise_post_0[i_sim,i_mask] = np.nansum( tw_post_0 / t_post_0 ) / (1.0*nproj_tnonzero_post_0)

            tw_hires_late_post_0[i_sim,i_wmask,:] = tw_hires_late_0 + tw_hires_post_0


        tj11[i_sim,:] = np.array([d.dt[(j11mask[j,:])].sum() for j in range(ncam)])
        tj11_late_post[i_sim,:] = np.array([d.dt[(d.bhprojsep[j,:]<latesep)&
                                                 (j11mask[j,:])].sum() for j in range(ncam)])
        tj11_0[i_sim,:] = np.array([d0.dt[(j11mask_0[j,:])].sum() for j in range(ncam)])
        tj11_late_post_0[i_sim,:] = np.array([d0.dt[(d0.bhprojsep[j,:]<latesep)&
                                                    (j11mask_0[j,:])].sum() for j in range(ncam)])
        tj11_noagn_l44_0[i_sim,:] = np.array([d0.dt[(d.lagn<1.0e44)&(j11mask_0[j,:])].sum() 
                                              for j in range(ncam)])
        
        tmycut[i_sim,:] = np.array([d.dt[(mycutmask[j,:])].sum() for j in range(ncam)])
        tmycut_late_post[i_sim,:] = np.array([d.dt[(d.bhprojsep[j,:]<latesep)&
                                                   (mycutmask[j,:])].sum() for j in range(ncam)])
        tmycut_0[i_sim,:] = np.array([d0.dt[(mycutmask_0[j,:])].sum() for j in range(ncam)])
        tmycut_late_post_0[i_sim,:] = np.array([d0.dt[(d0.bhprojsep[j,:]<latesep)&
                                                      (mycutmask_0[j,:])].sum() for j in range(ncam)])
        tmycut_noagn_l44_0[i_sim,:] = np.array([d0.dt[(d.lagn<1.0e44)&(mycutmask_0[j,:])].sum() 
                                                for j in range(ncam)])
        
    
    ### end loop over path_arr ###

    if (tw_hires_noagn_0.max()>0) or (tj11_noagn_l44_0.max()>0) or (tmycut_noagn_l44_0.max()>0):
            print "WARNING!! nonzero noagn_0 found!"
            print "tw_hires_noagn_0.max():",tw_hires_noagn_0.max() 
            print "tj11_noagn_l44_0.max():",tj11_noagn_l44_0.max()
            print "tmycut_noagn_l44_0.max():",tmycut_noagn_l44_0.max()
            #return


    if nsim>0 and multisim_plots:
    #if nsim>1 and multisim_plots:

        if valtype=='median':
            f_wise_hires_agnx0x1_tot_allsim = np.array([ np.nanmedian( (tw_hires_0[:,i,:]/tw_hires[:,i,:]) )
                                                         for i in range(tw_hires.shape[1]) ])
            f_wise_hires_agnx0x1_late_post_allsim = np.array([ np.nanmedian( (tw_hires_late_post_0[:,i,:]/
                                                                              tw_hires_late_post[:,i,:]) )
                                                               for i in range(tw_hires_late_post.shape[1]) ])
            f_wise_hires_agnx0x1_noagn_tot_allsim = np.array([ np.nanmedian( (tw_hires_noagn_0[:,i,:]/tw_hires[:,i,:]) )
                                                               for i in range(tw_hires.shape[1]) ])
            f_j11_agnx0x1_tot_allsim = np.nanmedian(tj11_0/tj11)
            f_j11_agnx0x1_late_post_allsim = np.nanmedian(tj11_late_post_0/tj11_late_post)
            f_j11_agnx0x1_noagn_tot_allsim = np.nanmedian(tj11_noagn_l44_0/tj11)
            f_mycut_agnx0x1_tot_allsim = np.nanmedian(tmycut_0/tmycut)
            f_mycut_agnx0x1_late_post_allsim = np.nanmedian(tmycut_late_post_0/tmycut_late_post)
            f_mycut_agnx0x1_noagn_tot_allsim = np.nanmedian(tmycut_noagn_l44_0/tmycut)
        else:
            f_wise_hires_agnx0x1_tot_allsim = np.array([ np.nanmean( (tw_hires_0[:,i,:]/tw_hires[:,i,:]) )
                                                         for i in range(tw_hires.shape[1]) ])
            f_wise_hires_agnx0x1_late_post_allsim = np.array([ np.nanmean( (tw_hires_late_post_0[:,i,:]/
                                                                            tw_hires_late_post[:,i,:]) )
                                                               for i in range(tw_hires_late_post.shape[1]) ])
            f_wise_hires_agnx0x1_noagn_tot_allsim = np.array([ np.nanmean( (tw_hires_noagn_0[:,i,:]/tw_hires[:,i,:]) )
                                                               for i in range(tw_hires.shape[1]) ])
            f_j11_agnx0x1_tot_allsim = np.nanmean(tj11_0/tj11)
            f_j11_agnx0x1_late_post_allsim = np.nanmean(tj11_late_post_0/tj11_late_post)
            f_j11_agnx0x1_noagn_tot_allsim = np.nanmean(tj11_noagn_l44_0/tj11)
            f_mycut_agnx0x1_tot_allsim = np.nanmean(tmycut_0/tmycut)
            f_mycut_agnx0x1_late_post_allsim = np.nanmean(tmycut_late_post_0/tmycut_late_post)
            f_mycut_agnx0x1_noagn_tot_allsim = np.nanmean(tmycut_noagn_l44_0/tmycut)


        yerrs_wise_agnx0x1_tot = errbars(np.array([ ( (tw_hires_0[:,i,:]/tw_hires[:,i,:]) ).flatten()
                                                    for i in range(tw_hires.shape[1]) ]).transpose(), 
                                         valtype=valtype,errtype=errtype)
        yerrs_wise_agnx0x1_late_post = errbars( np.array([ ( (tw_hires_late_post_0[:,i,:]/
                                                              tw_hires_late_post[:,i,:]) ).flatten()
                                                           for i in range(tw_hires_late_post.shape[1]) ]).transpose(), 
                                                valtype=valtype,errtype=errtype)
        yerrs_wise_agnx0x1_noagn_tot = errbars(np.array([ ( (tw_hires_noagn_0[:,i,:]/tw_hires[:,i,:]) ).flatten()
                                                          for i in range(tw_hires.shape[1]) ]).transpose(), 
                                               valtype=valtype,errtype=errtype)


        yerrs_j11_agnx0x1_tot = errbars(np.array([ (tj11_0/tj11).flatten() ]).transpose(), 
                                        valtype=valtype,errtype=errtype)
        yerrs_j11_agnx0x1_late_post = errbars(np.array([ (tj11_late_post_0/tj11_late_post).flatten() ]).transpose(),
                                              valtype=valtype,errtype=errtype)
        yerrs_j11_agnx0x1_noagn_tot = errbars(np.array([ (tj11_noagn_l44_0/tj11).flatten() ]).transpose(), 
                                              valtype=valtype,errtype=errtype)
            
        yerrs_mycut_agnx0x1_tot = errbars(np.array([ (tmycut_0/tmycut).flatten() ]).transpose(),
                                          valtype=valtype,errtype=errtype)
        yerrs_mycut_agnx0x1_late_post = errbars(np.array([ (tmycut_late_post_0/
                                                            tmycut_late_post).flatten() ]).transpose(),
                                                valtype=valtype,errtype=errtype)
        yerrs_mycut_agnx0x1_noagn_tot = errbars(np.array([ (tmycut_noagn_l44_0/tmycut).flatten() ]).transpose(),
                                          valtype=valtype,errtype=errtype)
            
        print "tw_hires:"
        print tw_hires
        print "tw_hires_0:"
        print tw_hires_0
        print f_wise_hires_agnx0x1_tot_allsim
        print f_wise_hires_agnx0x1_late_post_allsim

        print "\nf_wise_hires_agnx0x1_noagn_tot_allsim:"
        print f_wise_hires_agnx0x1_noagn_tot_allsim


        with open('%s/contamination_%sz%g_%s.txt'%(outdir,wflim_str,z,name),'w') as fc:
            for jj,wlimh in enumerate(w12_lims_hires):

                ### twise,agnx0/twise
                #w12_lims_hires
                #f_wise_hires_agnx0x1_late_post_allsim
                #yerrs_wise_agnx0x1_late_post
                #f_wise_hires_agnx0x1_tot_allsim
                #yerrs_wise_agnx0x1_tot

                ## twise,sf/twise
                #f_wise_hires_agnx0x1_noagn_tot_allsim


                fc.write(((10*"%g ")+"\n")%(wlimh,f_wise_hires_agnx0x1_late_post_allsim[jj],
                                            yerrs_wise_agnx0x1_late_post[0][jj],yerrs_wise_agnx0x1_late_post[1][jj],
                                            f_wise_hires_agnx0x1_tot_allsim[jj],
                                            yerrs_wise_agnx0x1_tot[0][jj], yerrs_wise_agnx0x1_tot[1][jj],
                                            f_wise_hires_agnx0x1_noagn_tot_allsim[jj],
                                            yerrs_wise_agnx0x1_noagn_tot[0][jj], yerrs_wise_agnx0x1_noagn_tot[1][jj]))
            ## twise,agnx0/twise
            #f_j11_agnx0x1_late_post_allsim
            #yerrs_j11_agnx0x1_late_post
            #f_j11_agnx0x1_tot_allsim
            #yerr=yerrs_j11_agnx0x1_tot

            #f_mycut_agnx0x1_late_post_allsim
            #yerrs_mycut_agnx0x1_late_post
            #f_mycut_agnx0x1_tot_allsim
            #yerrs_mycut_agnx0x1_tot

            ## twise,sf/twise
            #f_j11_agnx0x1_noagn_tot_allsim
            #f_mycut_agnx0x1_noagn_tot_allsim

            fc.write(((10*"%g ")+"\n")%(-11,f_j11_agnx0x1_late_post_allsim,
                                        yerrs_j11_agnx0x1_late_post[0],yerrs_j11_agnx0x1_late_post[1],
                                        f_j11_agnx0x1_tot_allsim,
                                        yerrs_j11_agnx0x1_tot[0], yerrs_j11_agnx0x1_tot[1],
                                        f_j11_agnx0x1_noagn_tot_allsim,
                                        yerrs_j11_agnx0x1_noagn_tot[0], yerrs_j11_agnx0x1_noagn_tot[1]))

            fc.write(((10*"%g ")+"\n")%(-12,f_mycut_agnx0x1_late_post_allsim,
                                        yerrs_mycut_agnx0x1_late_post[0],yerrs_mycut_agnx0x1_late_post[1],
                                        f_mycut_agnx0x1_tot_allsim,
                                        yerrs_mycut_agnx0x1_tot[0], yerrs_mycut_agnx0x1_tot[1],
                                        f_mycut_agnx0x1_noagn_tot_allsim,
                                        yerrs_mycut_agnx0x1_noagn_tot[0], yerrs_mycut_agnx0x1_noagn_tot[1]))

            fc.close()


        if contam_file_only: return
        
        ### make plots ###
        mpl.rcParams.update({'font.size': 12})

        if plot_wedge_allsim:
            ymax=np.max(np.append(f_wise_hires_agnx0x1_tot_allsim+yerrs_wise_agnx0x1_tot,
                                  f_wise_hires_agnx0x1_late_post_allsim+yerrs_wise_agnx0x1_late_post))
            if not ylim: ylim=(0,ymax+0.12)
            xlim= (0.28,1.06)
            plt.clf()
            plt.cla()
            plt.close()
            #fig=plt.figure(figsize=(5,3.5))
            fig=plt.figure(figsize=(5,5))
            ax=plt.subplot2grid((3,5),(0,0),colspan=4,rowspan=2)
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            ax.set_xlabel(r'min $W1 - W2$')
            ax.set_ylabel(r't$_{\rm WISE,AGNx0}$ / t$_{\rm WISE}$',fontsize=13)
            #ax.set_ylabel('t(WISE,AGNx0)/\nt(WISE)',fontsize=11)
            p1,=ax.plot(w12_lims_hires,f_wise_hires_agnx0x1_late_post_allsim,'*',
                        markersize=9,mew=1.6,mec='r',color='None')
            ax.errorbar(w12_lims_hires,f_wise_hires_agnx0x1_late_post_allsim,yerr=yerrs_wise_agnx0x1_late_post,
                        color='r',ecolor='r',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)
                
            p0,=ax.plot(w12_lims_hires,f_wise_hires_agnx0x1_tot_allsim,'o',
                        markersize=6,mew=1.2,mec='b',color='None')
            ax.errorbar(w12_lims_hires,f_wise_hires_agnx0x1_tot_allsim,yerr=yerrs_wise_agnx0x1_tot,
                        color='b',ecolor='b',ls='None',capsize=4)
            ax.legend((p0,p1), ('Total','Late+Post'), fontsize=10, loc='upper left', 
                      numpoints=1, handletextpad=0.1)

            ax2=plt.subplot2grid((3,5),(0,4),rowspan=2)
            x1=0.25
            x2=0.75
            ax2.set_xlim(0,1)
            ax2.set_ylim(ylim)
            ax2.set_xticks(np.array([x1,x2]))
            ax2.set_xticklabels(np.array(['J11','This\nwork']),fontsize=10)
            ax2.set_yticklabels(np.repeat('',10))
            ax2.plot(x1,f_j11_agnx0x1_late_post_allsim,"*",markersize=9,mew=1.6,mec='darkorange',color='None')
            ax2.errorbar(x1,f_j11_agnx0x1_late_post_allsim,yerr=yerrs_j11_agnx0x1_late_post,
                        color='darkorange',ecolor='darkorange',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)
            ax2.plot(x2,f_mycut_agnx0x1_late_post_allsim,"*",markersize=9,mew=1.6,mec='darkorange',color='None')
            ax2.errorbar(x2,f_mycut_agnx0x1_late_post_allsim,yerr=yerrs_mycut_agnx0x1_late_post,
                        color='darkorange',ecolor='darkorange',ls='None',
                        capsize=4,elinewidth=2,markeredgewidth=2)
            ax2.plot(x1,f_j11_agnx0x1_tot_allsim,"o",markersize=6,mew=1.2,mec='c',color='None')
            ax2.errorbar(x1,f_j11_agnx0x1_tot_allsim,yerr=yerrs_j11_agnx0x1_tot,
                        color='c',ecolor='c',ls='None',capsize=4)
            ax2.plot(x2,f_mycut_agnx0x1_tot_allsim,"o",markersize=6,mew=1.2,mec='c',color='None')
            ax2.errorbar(x2,f_mycut_agnx0x1_tot_allsim,yerr=yerrs_mycut_agnx0x1_tot,
                        color='c',ecolor='c',ls='None',capsize=4)

            
            ax3=plt.subplot2grid((3,5),(2,0),colspan=4)
            ax3.set_xlabel(r'min $W1 - W2$')
            ax3.set_ylabel(r't$_{\rm WISE,SF}$ / t$_{\rm WISE}$',fontsize=13)
            #ax3.set_ylabel('t(WISE,SFonly)/\nt(WISE)',fontsize=11)
            ax3.set_xlim(xlim)
            ax3.set_ylim(0,0.18)
            ax3.set_yticks(np.array([0,0.1]))
            p2,=ax3.plot(w12_lims_hires,f_wise_hires_agnx0x1_noagn_tot_allsim,'H',markersize=9,color='k')
            ax3.legend((p2,), ('Total (SF only)',), fontsize=10, loc='upper left', 
                       numpoints=1, handletextpad=0.1)

            ax4=plt.subplot2grid((3,5),(2,4))
            ax4.set_xlim(0,1)
            ax4.set_ylim(0,0.18)
            ax4.set_xticks(np.array([x1,x2]))
            ax4.set_xticklabels(np.array(['J11','This\nwork']),fontsize=10)
            ax4.set_yticks(np.array([0,0.1]))
            ax4.set_yticklabels(np.repeat('',2))
            ax4.plot(x1,f_j11_agnx0x1_noagn_tot_allsim,'H',ms=9,color='m',mec='m')
            ax4.plot(x2,f_mycut_agnx0x1_noagn_tot_allsim,'H',ms=9,color='m',mec='m')

            fig.subplots_adjust(bottom=0.1,left=0.14,right=0.94,top=0.96,wspace=0,hspace=0.55)
            wstr='_alt' if alt_wedge else ''
            fig.savefig(outdir+'/wise_fraction_agnx0compare_allsim_%s_%s_%s%s.pdf'%(name,valtype,
                                                                                    errtype,wstr))

        else:
            plt.clf()
            plt.cla()
            plt.close()
        #fig=plt.figure(figsize=(5,5))
            fig=plt.figure(figsize=(5,3.5))
            ax=fig.add_subplot(111)
        #ax.set_xlim(0.28,1.02)
            xlim=(0.28,1.22) if plot_wedge_allsim else (0.28,1.02)
            ax.set_xlim(xlim)
        #ax.set_ylim(0,1)
        #ylim=(0,0.28)
            ymax=np.max(np.append(f_wise_hires_agnx0x1_tot_allsim+yerrs_wise_agnx0x1_tot,
                                  f_wise_hires_agnx0x1_late_post_allsim+yerrs_wise_agnx0x1_late_post))
            ylim=(0,ymax+0.12)
            ax.set_ylim(ylim)
            ax.set_xlabel(r'min $W1 - W2$')
            ax.set_ylabel(r't$_{\rm WISE,AGNx0}$ / t$_{\rm WISE}$')
            p1,=ax.plot(w12_lims_hires,f_wise_hires_agnx0x1_late_post_allsim,'*',
                        markersize=9,mew=1.6,mec='r',color='None')
            ax.errorbar(w12_lims_hires,f_wise_hires_agnx0x1_late_post_allsim,yerr=yerrs_wise_agnx0x1_late_post,
                        color='r',ecolor='r',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)
            if plot_wedge_allsim:
                ax.plot(1.1,f_j11_agnx0x1_late_post_allsim,"*",markersize=9,mew=1.6,mec='m',color='None')
                ax.errorbar(1.1,f_j11_agnx0x1_late_post_allsim,yerr=yerrs_j11_agnx0x1_late_post,
                            color='m',ecolor='m',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)
                ax.plot(1.2,f_mycut_agnx0x1_late_post_allsim,"*",markersize=9,mew=1.6,mec='darkorange',color='None')
                ax.errorbar(1.2,f_mycut_agnx0x1_late_post_allsim,yerr=yerrs_mycut_agnx0x1_late_post,
                            color='darkorange',ecolor='darkorange',ls='None',
                            capsize=4,elinewidth=2,markeredgewidth=2)
                
            p0,=ax.plot(w12_lims_hires,f_wise_hires_agnx0x1_tot_allsim,'o',
                markersize=6,mew=1.2,mec='b',color='None')
            ax.errorbar(w12_lims_hires,f_wise_hires_agnx0x1_tot_allsim,yerr=yerrs_wise_agnx0x1_tot,
                        color='b',ecolor='b',ls='None',capsize=4)
            if plot_wedge_allsim:
                ax.plot(1.1,f_j11_agnx0x1_tot_allsim,"o",markersize=6,mew=1.2,mec='g',color='None')
                ax.errorbar(1.1,f_j11_agnx0x1_tot_allsim,yerr=yerrs_j11_agnx0x1_tot,
                            color='g',ecolor='g',ls='None',capsize=4)
                ax.plot(1.2,f_mycut_agnx0x1_tot_allsim,"o",markersize=6,mew=1.2,mec='c',color='None')
                ax.errorbar(1.2,f_mycut_agnx0x1_tot_allsim,yerr=yerrs_mycut_agnx0x1_tot,
                            color='c',ecolor='c',ls='None',capsize=4)
            ax.legend((p0,p1), ('Total','Late+Post'), fontsize=10, loc='upper left', 
                      numpoints=1, handletextpad=0.1)

            ax2=ax.twinx()
            p2,=ax2.plot(w12_lims_hires,f_wise_hires_agnx0x1_noagn_tot_allsim,'H',markersize=9,color='k')
            if plot_wedge_allsim:
                ax.plot(1.1,f_j11_agnx0x1_noagn_tot_allsim,'H',ms=9,color='g',mec='g')
                ax.plot(1.2,f_mycut_agnx0x1_noagn_tot_allsim,'H',ms=9,color='c',mec='c')
                
            ax2.set_ylabel(r't$_{\rm WISE,SFonly}$ / t$_{\rm WISE}$')
            ax2.set_xlim(xlim)
            ax2.set_ylim(ylim)
            ax2.legend((p2,), ('Total, SF only',), fontsize=10, loc='upper right', 
                       numpoints=1, handletextpad=0.1)
            fig.subplots_adjust(bottom=0.14,left=0.14,right=0.86,top=0.88)
            fig.savefig(outdir+'/wise_fraction_agnx0compare_allsim_%s_%s_%s.pdf'%(name,valtype,errtype))
            
        print f_j11_agnx0x1_tot_allsim
        print f_j11_agnx0x1_late_post_allsim
        print f_j11_agnx0x1_noagn_tot_allsim
        #return

        #xarrs = np.array([np.log10(ssfr_max),sfr_max,np.log10(lagn_max),np.log10(ltot_max),flbol_max])
        #xarrs=np.append(np.arange(nsim).reshape(1,2),xarrs,axis=0)
        #xarrs = np.array([np.arange(nsim),np.log10(ssfr_max),sfr_max,np.log10(lagn_max),
        #                  np.log10(ltot_max),flbol_max])
        xarrs = np.array([np.arange(nsim),np.log10(ssfr_max),np.log10(sfr_max),np.log10(lagn_max),
                          np.log10(ltot_max),flbol_max])

        for iplot,ix in enumerate([ix_w05[0],ix_w08[0],ix_j11[0]]):
            print '%s plots...'%wfbase[iplot]

            ymax=0.35 if latesep>=10 else 0.6
            plot_multisim_totals((f_wise_tot_0[:,ix],f_wise_early_0[:,ix],
                                  f_wise_late_0[:,ix],f_wise_post_0[:,ix]),
                                 ylbl = r' <t$_{\rm WISE}$/t$_{\rm tot}$>', nplots=3,
                                 xvar = np.append('',xtypes[:2]),xarr=xarrs[:3,:], ylim=(0.0,ymax),
                                 outdir=outdir,extra=extra,tmax_postmrg=tmax_postmrg,latesep=latesep,
                                 fbase='wise_fraction_agnx0_%s_%s'%(wfbase[iplot],name),
                                 plot_titles=np.repeat('',6),leg_lbls=('Total','Early','Late','Post'),
                                 carr=['k','g','m','orange'],msarr=['o']+['^']+['s']+['v'],
                                 mewarr=[1.5]+3*[1.0],mszarr=4*[6],leg_loc='upper left',
                                 title='WISE fraction (%s) - No AGN'%wisemask_labels[iplot])
            
            ymax=0.41 if latesep>=10 else 0.8
            plot_multisim_totals((f_wise_agnx0x1_tot[:,ix],f_wise_agnx0x1_early[:,ix],
                                  f_wise_agnx0x1_late[:,ix],f_wise_agnx0x1_post[:,ix]),
                                 ylbl = r' <t$_{\rm WISE,AGNX0}$/t$_{\rm WISE}$>', nplots=3,
                                 xvar = np.append('',xtypes[:2]),xarr=xarrs[:3,:],outdir=outdir,extra=extra,
                                 tmax_postmrg=tmax_postmrg,latesep=latesep,ylim=(0.0,ymax),
                                 fbase='wise_fraction_agnx0compare_%s_%s'%(wfbase[iplot],name),
                                 plot_titles=np.repeat('',6),leg_lbls=('Total','Early','Late','Post'),
                                 carr=['k','g','m','orange'],msarr=['o']+['^']+['s']+['v'],
                                 mewarr=[1.5]+3*[1.0],mszarr=4*[6],leg_loc='upper left',
                                 title='WISE lifetime (AGNx0 / AGNx1) (%s)'%wisemask_labels[iplot])

        for iplot,ix in enumerate([ix_w05,ix_w08,ix_j11]):
            print '%s plots...'%wfbase[iplot]

            ymax=1.0
            plot_multisim_totals((f_wise_agnx0x1_noagn_tot[:,ix],f_wise_agnx0x1_noagn_early[:,ix],
                                  f_wise_agnx0x1_noagn_late[:,ix],f_wise_agnx0x1_noagn_post[:,ix]),
                                 ylbl = r' <t$_{\rm WISE,SFonly}$/t$_{\rm WISE}$>',
                                 tmax_postmrg=tmax_postmrg,latesep=latesep,ylim=(0.0,ymax),
                                 fbase='wise_fraction_agnx0compare_noagn_%s_%s'%(wfbase[iplot],name),
                                 plot_titles=plot_titles,leg_lbls=('Total','Early','Late','Post'),
                                 carr=['k','g','m','orange'],msarr=['o']+['^']+['s']+['v'],
                                 mewarr=[1.5]+3*[1.0],mszarr=4*[6],
                                 title='WISE lifetime (AGNX0,SF only / AGN) (%s)'%wisemask_labels[iplot])
            


def completeness_vs_z(maindir='/oasis/projects/nsf/hvd115/lblecha/',subdir_arr='fid',wise_fluxlim_type=''):

    
    wflim_str=wise_fluxlim_type+'_' if wise_fluxlim_type else wise_fluxlim_type

    compfiles = np.array(glob.glob(maindir+'/completeness_%sz*_%s.txt'%(wflim_str,subdir_arr)))
    z_arr = np.array([ np.float((s.split('z')[1]).split('_')[0]) for s in compfiles ])
    ixsort=z_arr.argsort(kind='mergesort')
    z_arr=z_arr[ixsort]
    compfiles=compfiles[ixsort]
    print z_arr

    fl44_late_w05 = np.nan*np.ones((z_arr.size))
    fl44_tot_w05 = np.nan*np.ones((z_arr.size))
    yerr_fl44_late_w05 = np.nan*np.ones((2,z_arr.size))
    yerr_fl44_tot_w05 = np.nan*np.ones((2,z_arr.size))

    fl45_late_w05 = np.nan*np.ones((z_arr.size))
    fl45_tot_w05 = np.nan*np.ones((z_arr.size))
    yerr_fl45_late_w05 = np.nan*np.ones((2,z_arr.size))
    yerr_fl45_tot_w05 = np.nan*np.ones((2,z_arr.size))

    fl44_late_w08 = np.nan*np.ones((z_arr.size))
    fl44_tot_w08 = np.nan*np.ones((z_arr.size))
    yerr_fl44_late_w08 = np.nan*np.ones((2,z_arr.size))
    yerr_fl44_tot_w08 = np.nan*np.ones((2,z_arr.size))

    fl45_late_w08 = np.nan*np.ones((z_arr.size))
    fl45_tot_w08 = np.nan*np.ones((z_arr.size))
    yerr_fl45_late_w08 = np.nan*np.ones((2,z_arr.size))
    yerr_fl45_tot_w08 = np.nan*np.ones((2,z_arr.size))
    
    fl44_late_j11 = np.nan*np.ones((z_arr.size))
    fl44_tot_j11 = np.nan*np.ones((z_arr.size))
    yerr_fl44_late_j11 = np.nan*np.ones((2,z_arr.size))
    yerr_fl44_tot_j11 = np.nan*np.ones((2,z_arr.size))

    fl45_late_j11 = np.nan*np.ones((z_arr.size))
    fl45_tot_j11 = np.nan*np.ones((z_arr.size))
    yerr_fl45_late_j11 = np.nan*np.ones((2,z_arr.size))
    yerr_fl45_tot_j11 = np.nan*np.ones((2,z_arr.size))
    
    fl44_late_mycut = np.nan*np.ones((z_arr.size))
    fl44_tot_mycut = np.nan*np.ones((z_arr.size))
    yerr_fl44_late_mycut = np.nan*np.ones((2,z_arr.size))
    yerr_fl44_tot_mycut = np.nan*np.ones((2,z_arr.size))

    fl45_late_mycut = np.nan*np.ones((z_arr.size))
    fl45_tot_mycut = np.nan*np.ones((z_arr.size))
    yerr_fl45_late_mycut = np.nan*np.ones((2,z_arr.size))
    yerr_fl45_tot_mycut = np.nan*np.ones((2,z_arr.size))

    for cf in compfiles:
        print "reading file %s"%cf
        with open(cf,'r') as f:
            z = np.float((cf.split('z')[1]).split('_')[0])
            tmpdata = np.loadtxt(f,unpack=True,dtype='float64')
            wlim,fl44_late,yerr0_fl44_late,yerr1_fl44_late,fl44_tot,yerr0_fl44_tot,yerr1_fl44_tot,fl45_late,yerr0_fl45_late,yerr1_fl45_late,fl45_tot,yerr0_fl45_tot,yerr1_fl45_tot = tmpdata

            ### W1W2 > 0.5 ###
            fl44_late_w05[z_arr==z] = fl44_late[wlim==0.5]
            fl44_tot_w05[z_arr==z] = fl44_tot[wlim==0.5]
            yerr_fl44_late_w05[:,z_arr==z] = np.array([yerr0_fl44_late[wlim==0.5],yerr1_fl44_late[wlim==0.5]])
            yerr_fl44_tot_w05[:,z_arr==z] = np.array([yerr0_fl44_tot[wlim==0.5],yerr1_fl44_tot[wlim==0.5]])

            fl45_late_w05[z_arr==z] = fl45_late[wlim==0.5]
            fl45_tot_w05[z_arr==z] = fl45_tot[wlim==0.5]
            yerr_fl45_late_w05[:,z_arr==z] = np.array([yerr0_fl45_late[wlim==0.5],yerr1_fl45_late[wlim==0.5]])
            yerr_fl45_tot_w05[:,z_arr==z] = np.array([yerr0_fl45_tot[wlim==0.5],yerr1_fl45_tot[wlim==0.5]])

            ### W1W2 > 0.8 ###
            fl44_late_w08[z_arr==z] = fl44_late[wlim==0.8]
            fl44_tot_w08[z_arr==z] = fl44_tot[wlim==0.8]
            yerr_fl44_late_w08[:,z_arr==z] = np.array([yerr0_fl44_late[wlim==0.8],yerr1_fl44_late[wlim==0.8]])
            yerr_fl44_tot_w08[:,z_arr==z] = np.array([yerr0_fl44_tot[wlim==0.8],yerr1_fl44_tot[wlim==0.8]])

            fl45_late_w08[z_arr==z] = fl45_late[wlim==0.8]
            fl45_tot_w08[z_arr==z] = fl45_tot[wlim==0.8]
            yerr_fl45_late_w08[:,z_arr==z] = np.array([yerr0_fl45_late[wlim==0.8],yerr1_fl45_late[wlim==0.8]])
            yerr_fl45_tot_w08[:,z_arr==z] = np.array([yerr0_fl45_tot[wlim==0.8],yerr1_fl45_tot[wlim==0.8]])
            
            ### J11 ###
            fl44_late_j11[z_arr==z] = fl44_late[wlim==-11]
            fl44_tot_j11[z_arr==z] = fl44_tot[wlim==-11]
            yerr_fl44_late_j11[:,z_arr==z] = np.array([yerr0_fl44_late[wlim==-11],yerr1_fl44_late[wlim==-11]])
            yerr_fl44_tot_j11[:,z_arr==z] = np.array([yerr0_fl44_tot[wlim==-11],yerr1_fl44_tot[wlim==-11]])

            fl45_late_j11[z_arr==z] = fl45_late[wlim==-11]
            fl45_tot_j11[z_arr==z] = fl45_tot[wlim==-11]
            yerr_fl45_late_j11[:,z_arr==z] = np.array([yerr0_fl45_late[wlim==-11],yerr1_fl45_late[wlim==-11]])
            yerr_fl45_tot_j11[:,z_arr==z] = np.array([yerr0_fl45_tot[wlim==-11],yerr1_fl45_tot[wlim==-11]])

            ### my cut ###
            fl44_late_mycut[z_arr==z] = fl44_late[wlim==-12]
            fl44_tot_mycut[z_arr==z] = fl44_tot[wlim==-12]
            yerr_fl44_late_mycut[:,z_arr==z] = np.array([yerr0_fl44_late[wlim==-12],yerr1_fl44_late[wlim==-12]])
            yerr_fl44_tot_mycut[:,z_arr==z] = np.array([yerr0_fl44_tot[wlim==-12],yerr1_fl44_tot[wlim==-12]])

            fl45_late_mycut[z_arr==z] = fl45_late[wlim==-12]
            fl45_tot_mycut[z_arr==z] = fl45_tot[wlim==-12]
            yerr_fl45_late_mycut[:,z_arr==z] = np.array([yerr0_fl45_late[wlim==-12],yerr1_fl45_late[wlim==-12]])
            yerr_fl45_tot_mycut[:,z_arr==z] = np.array([yerr0_fl45_tot[wlim==-12],yerr1_fl45_tot[wlim==-12]])

    

    print "making l44 plots"
    plt.clf()
    plt.cla()
    plt.close()
    mpl.rcParams.update({'font.size': 12})
    #fig=plt.figure(figsize=(5,3.5))
    fig=plt.figure(figsize=(8,6))
    #ylim=(0,1.02)
    ylim=(0,1.0)

    ax1=fig.add_subplot(221)
    ax1.set_xlim(0-0.05,z_arr.max()+0.05)
    ax1.set_ylim(ylim)
    ax1.set_xlabel('z')
    ax1.set_ylabel(r't$_{\rm WISEAGN}$ / t$_{\rm AGN}$')
    print "plotting mycut"
    p11,=ax1.plot(z_arr,fl44_late_mycut,'*-',color='darkorange',lw=2.5,ms=9)
    ax1.errorbar(z_arr,fl44_late_mycut,yerr=yerr_fl44_late_mycut,
                ecolor='darkorange',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)            
    p01,=ax1.plot(z_arr,fl44_tot_mycut,'c',marker='o',ls='-',lw=2.5,ms=5)
    ax1.errorbar(z_arr,fl44_tot_mycut,yerr=yerr_fl44_tot_mycut,
                ecolor='c',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)            
    #ax1.set_title(r'2-color (this work); L$_{\rm AGN}>10^{44}$ erg s$^{-1}$',fontsize=11)
    ax1.set_title(r'2-color (this work)',fontsize=11)
    ax1.legend((p01,p11), ('Total','Late+Post'), fontsize=10, loc='lower left', 
              numpoints=1, handletextpad=0.1)


    ax2=fig.add_subplot(222)
    ax2.set_xlim(0-0.05,z_arr.max()+0.05)
    ax2.set_ylim(ylim)
    ax2.set_xlabel('z')
    ax2.set_ylabel(r't$_{\rm WISEAGN}$ / t$_{\rm AGN}$')
    print "plotting w08"
    p12,=ax2.plot(z_arr,fl44_late_w08,'r*-',lw=2,ms=9)
    ax2.errorbar(z_arr,fl44_late_w08,yerr=yerr_fl44_late_w08,
                ecolor='r',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)            
    p02,=ax2.plot(z_arr,fl44_tot_w08,'bo-',lw=2,ms=5)
    ax2.errorbar(z_arr,fl44_tot_w08,yerr=yerr_fl44_tot_w08,
                ecolor='b',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)            
    #ax2.set_title(r'W1-W2>0.8; L$_{\rm AGN}>10^{44}$ erg s$^{-1}$',fontsize=11)
    ax2.set_title(r'W1-W2>0.8',fontsize=11)
    ax2.legend((p02,p12), ('Total','Late+Post'), fontsize=10, loc='lower left', 
              numpoints=1, handletextpad=0.1)


    ax3=fig.add_subplot(223)
    ax3.set_xlim(0-0.05,z_arr.max()+0.05)
    ax3.set_ylim(ylim)
    ax3.set_xlabel('z')
    ax3.set_ylabel(r't$_{\rm WISEAGN}$ / t$_{\rm AGN}$')
    print "plotting j11"
    p13,=ax3.plot(z_arr,fl44_late_j11,'k',marker='*',ls='-',lw=2.5,ms=9)
    ax3.errorbar(z_arr,fl44_late_j11,yerr=yerr_fl44_late_j11,
                ecolor='k',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)            
    p03,=ax3.plot(z_arr,fl44_tot_j11,'g',marker='o',ls='-',lw=2.5,ms=5)
    ax3.errorbar(z_arr,fl44_tot_j11,yerr=yerr_fl44_tot_j11,
                ecolor='g',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)            
    #ax3.set_title(r'2-color (J11); L$_{\rm AGN}>10^{44}$ erg s$^{-1}$',fontsize=11)
    ax3.set_title(r'2-color (J11)',fontsize=11)
    ax3.legend((p03,p13), ('Total','Late+Post'), fontsize=10, loc='lower left', 
              numpoints=1, handletextpad=0.1)

    ax4=fig.add_subplot(224)
    ax4.set_xlim(0-0.05,z_arr.max()+0.05)
    ax4.set_ylim(ylim)
    ax4.set_xlabel('z')
    ax4.set_ylabel(r't$_{\rm WISEAGN}$ / t$_{\rm AGN}$')
    print "plotting w05"
    p14,=ax4.plot(z_arr,fl44_late_w05,'m*-',ms=9)
    ax4.errorbar(z_arr,fl44_late_w05,yerr=yerr_fl44_late_w05,
                ecolor='m',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)            
    p04,=ax4.plot(z_arr,fl44_tot_w05,'o-',color='darkblue',ms=5)
    ax4.errorbar(z_arr,fl44_tot_w05,yerr=yerr_fl44_tot_w05,
                 ecolor='darkblue',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)            
    ax4.legend((p04,p14), ('Total','Late+Post'), fontsize=10, loc='lower left', 
              numpoints=1, handletextpad=0.1)
    #ax4.set_title(r'W1-W2>0.5; L$_{\rm AGN}>10^{44}$ erg s$^{-1}$',fontsize=11)
    ax4.set_title(r'W1-W2>0.5',fontsize=11)


            
    fig.suptitle(r'L$_{\rm AGN}>10^{44}$ erg s$^{-1}$',fontsize=12)
    fig.subplots_adjust(bottom=0.08,left=0.1,right=0.96,top=0.9,wspace=0.25,hspace=0.35)
    fig.savefig(maindir+'/completeness_vs_z_l44_%s%s.pdf'%(wflim_str,subdir_arr))
    #plt.clf()
    #plt.cla()
    #plt.close()



    print "making l45 plots"
    plt.clf()
    plt.cla()
    plt.close()
    mpl.rcParams.update({'font.size': 12})
    #fig=plt.figure(figsize=(5,3.5))
    fig=plt.figure(figsize=(8,6))
    #ylim=(0,1.02)
    ylim=(0,1.0)

    ax1=fig.add_subplot(221)
    ax1.set_xlim(0-0.05,z_arr.max()+0.05)
    ax1.set_ylim(ylim)
    ax1.set_xlabel('z')
    ax1.set_ylabel(r't$_{\rm WISEAGN}$ / t$_{\rm AGN}$')
    print "plotting mycut"
    p11,=ax1.plot(z_arr,fl45_late_mycut,'*-',color='darkorange',lw=2.5,ms=9)
    ax1.errorbar(z_arr,fl45_late_mycut,yerr=yerr_fl45_late_mycut,
                ecolor='darkorange',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)            
    p01,=ax1.plot(z_arr,fl45_tot_mycut,'c',marker='o',ls='-',lw=2.5,ms=5)
    ax1.errorbar(z_arr,fl45_tot_mycut,yerr=yerr_fl45_tot_mycut,
                ecolor='c',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)            
    #ax1.set_title(r'2-color (this work); L$_{\rm AGN}>10^{45}$ erg s$^{-1}$',fontsize=11)
    ax1.set_title(r'2-color (this work)',fontsize=11)
    ax1.legend((p01,p11), ('Total','Late+Post'), fontsize=10, loc='lower left', 
              numpoints=1, handletextpad=0.1)


    ax2=fig.add_subplot(222)
    ax2.set_xlim(0-0.05,z_arr.max()+0.05)
    ax2.set_ylim(ylim)
    ax2.set_xlabel('z')
    ax2.set_ylabel(r't$_{\rm WISEAGN}$ / t$_{\rm AGN}$')
    print "plotting w08"
    p12,=ax2.plot(z_arr,fl45_late_w08,'r*-',lw=2,ms=9)
    ax2.errorbar(z_arr,fl45_late_w08,yerr=yerr_fl45_late_w08,
                ecolor='r',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)            
    p02,=ax2.plot(z_arr,fl45_tot_w08,'bo-',lw=2,ms=5)
    ax2.errorbar(z_arr,fl45_tot_w08,yerr=yerr_fl45_tot_w08,
                ecolor='b',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)            
    #ax2.set_title(r'W1-W2>0.8; L$_{\rm AGN}>10^{45}$ erg s$^{-1}$',fontsize=11)
    ax2.set_title(r'W1-W2>0.8',fontsize=11)
    ax2.legend((p02,p12), ('Total','Late+Post'), fontsize=10, loc='lower left', 
              numpoints=1, handletextpad=0.1)


    ax3=fig.add_subplot(223)
    ax3.set_xlim(0-0.05,z_arr.max()+0.05)
    ax3.set_ylim(ylim)
    ax3.set_xlabel('z')
    ax3.set_ylabel(r't$_{\rm WISEAGN}$ / t$_{\rm AGN}$')
    print "plotting j11"
    p13,=ax3.plot(z_arr,fl45_late_j11,'k',marker='*',ls='-',lw=2.5,ms=9)
    ax3.errorbar(z_arr,fl45_late_j11,yerr=yerr_fl45_late_j11,
                ecolor='k',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)            
    p03,=ax3.plot(z_arr,fl45_tot_j11,'g',marker='o',ls='-',lw=2.5,ms=5)
    ax3.errorbar(z_arr,fl45_tot_j11,yerr=yerr_fl45_tot_j11,
                ecolor='g',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)            
    #ax3.set_title(r'2-color (J11); L$_{\rm AGN}>10^{45}$ erg s$^{-1}$',fontsize=11)
    ax3.set_title(r'2-color (J11)',fontsize=11)
    ax3.legend((p03,p13), ('Total','Late+Post'), fontsize=10, loc='lower left', 
              numpoints=1, handletextpad=0.1)

    ax4=fig.add_subplot(224)
    ax4.set_xlim(0-0.05,z_arr.max()+0.05)
    ax4.set_ylim(ylim)
    ax4.set_xlabel('z')
    ax4.set_ylabel(r't$_{\rm WISEAGN}$ / t$_{\rm AGN}$')
    print "plotting w05"
    p14,=ax4.plot(z_arr,fl45_late_w05,'m*-',ms=9)
    ax4.errorbar(z_arr,fl45_late_w05,yerr=yerr_fl45_late_w05,
                ecolor='m',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)            
    p04,=ax4.plot(z_arr,fl45_tot_w05,'o-',color='darkblue',ms=5)
    ax4.errorbar(z_arr,fl45_tot_w05,yerr=yerr_fl45_tot_w05,
                 ecolor='darkblue',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)            
    ax4.legend((p04,p14), ('Total','Late+Post'), fontsize=10, loc='lower left', 
              numpoints=1, handletextpad=0.1)
    #ax4.set_title(r'W1-W2>0.5; L$_{\rm AGN}>10^{45}$ erg s$^{-1}$',fontsize=11)
    ax4.set_title(r'W1-W2>0.5',fontsize=11)


            
    fig.suptitle(r'L$_{\rm AGN}>10^{45}$ erg s$^{-1}$',fontsize=12)
    fig.subplots_adjust(bottom=0.08,left=0.1,right=0.96,top=0.9,wspace=0.25,hspace=0.35)
    fig.savefig(maindir+'/completeness_vs_z_l45_%s%s.pdf'%(wflim_str,subdir_arr))
    #plt.clf()
    #plt.cla()
    #plt.close()


def contamination_vs_z(maindir='/oasis/projects/nsf/hvd115/lblecha/',subdir_arr='fid',wise_fluxlim_type=''):

    
    wflim_str=wise_fluxlim_type+'_' if wise_fluxlim_type else wise_fluxlim_type

    confiles = np.array(glob.glob(maindir+'/contamination_%sz*_%s*.txt'%(wflim_str,subdir_arr)))
    z_arr = np.array([ np.float((s.split('z')[1]).split('_')[0]) for s in confiles ])
    ixsort=z_arr.argsort(kind='mergesort')
    z_arr=z_arr[ixsort]
    confiles=confiles[ixsort]

    fagnx0x1_late_w05 = np.nan*np.ones((z_arr.size))
    fagnx0x1_tot_w05 = np.nan*np.ones((z_arr.size))
    fagnx0x1_noagn_tot_w05 = np.nan*np.ones((z_arr.size))
    yerr_fagnx0x1_late_w05 = np.nan*np.ones((2,z_arr.size))
    yerr_fagnx0x1_tot_w05 = np.nan*np.ones((2,z_arr.size))
    yerr_fagnx0x1_noagn_tot_w05 = np.nan*np.ones((2,z_arr.size))

    fagnx0x1_late_w08 = np.nan*np.ones((z_arr.size))
    fagnx0x1_tot_w08 = np.nan*np.ones((z_arr.size))
    fagnx0x1_noagn_tot_w08 = np.nan*np.ones((z_arr.size))
    yerr_fagnx0x1_late_w08 = np.nan*np.ones((2,z_arr.size))
    yerr_fagnx0x1_tot_w08 = np.nan*np.ones((2,z_arr.size))
    yerr_fagnx0x1_noagn_tot_w08 = np.nan*np.ones((2,z_arr.size))

    fagnx0x1_late_j11 = np.nan*np.ones((z_arr.size))
    fagnx0x1_tot_j11 = np.nan*np.ones((z_arr.size))
    fagnx0x1_noagn_tot_j11 = np.nan*np.ones((z_arr.size))
    yerr_fagnx0x1_late_j11 = np.nan*np.ones((2,z_arr.size))
    yerr_fagnx0x1_tot_j11 = np.nan*np.ones((2,z_arr.size))
    yerr_fagnx0x1_noagn_tot_j11 = np.nan*np.ones((2,z_arr.size))

    fagnx0x1_late_mycut = np.nan*np.ones((z_arr.size))
    fagnx0x1_tot_mycut = np.nan*np.ones((z_arr.size))
    fagnx0x1_noagn_tot_mycut = np.nan*np.ones((z_arr.size))
    yerr_fagnx0x1_late_mycut = np.nan*np.ones((2,z_arr.size))
    yerr_fagnx0x1_tot_mycut = np.nan*np.ones((2,z_arr.size))
    yerr_fagnx0x1_noagn_tot_mycut = np.nan*np.ones((2,z_arr.size))


    for cf in confiles:
        print "reading file %s"%cf
        with open(cf,'r') as f:
            z = np.float((cf.split('z')[1]).split('_')[0])
            tmpdata = np.loadtxt(f,unpack=True,dtype='float64')
            wlim,fagnx0x1_late,yerr0_fagnx0x1_late,yerr1_fagnx0x1_late,fagnx0x1_tot,yerr0_fagnx0x1_tot,yerr1_fagnx0x1_tot,fagnx0x1_noagn_tot,yerr0_fagnx0x1_noagn_tot,yerr1_fagnx0x1_noagn_tot = tmpdata

            ### W1W2 > 0.5 ###

            fagnx0x1_late_w05[z_arr==z] = fagnx0x1_late[wlim==0.5]
            fagnx0x1_tot_w05[z_arr==z] = fagnx0x1_tot[wlim==0.5]
            fagnx0x1_noagn_tot_w05[z_arr==z] = fagnx0x1_noagn_tot[wlim==0.5]
            yerr_fagnx0x1_late_w05[:,z_arr==z] = np.array([yerr0_fagnx0x1_late[wlim==0.5],yerr1_fagnx0x1_late[wlim==0.5]])
            yerr_fagnx0x1_tot_w05[:,z_arr==z] = np.array([yerr0_fagnx0x1_tot[wlim==0.5],yerr1_fagnx0x1_tot[wlim==0.5]])
            yerr_fagnx0x1_noagn_tot_w05[:,z_arr==z] = np.array([yerr0_fagnx0x1_noagn_tot[wlim==0.5],
                                                                yerr1_fagnx0x1_noagn_tot[wlim==0.5]])

            fagnx0x1_late_w08[z_arr==z] = fagnx0x1_late[wlim==0.8]
            fagnx0x1_tot_w08[z_arr==z] = fagnx0x1_tot[wlim==0.8]
            fagnx0x1_noagn_tot_w08[z_arr==z] = fagnx0x1_noagn_tot[wlim==0.8]
            yerr_fagnx0x1_late_w08[:,z_arr==z] = np.array([yerr0_fagnx0x1_late[wlim==0.8],yerr1_fagnx0x1_late[wlim==0.8]])
            yerr_fagnx0x1_tot_w08[:,z_arr==z] = np.array([yerr0_fagnx0x1_tot[wlim==0.8],yerr1_fagnx0x1_tot[wlim==0.8]])
            yerr_fagnx0x1_noagn_tot_w08[:,z_arr==z] = np.array([yerr0_fagnx0x1_noagn_tot[wlim==0.8],
                                                                yerr1_fagnx0x1_noagn_tot[wlim==0.8]])

            fagnx0x1_late_j11[z_arr==z] = fagnx0x1_late[wlim==-11]
            fagnx0x1_tot_j11[z_arr==z] = fagnx0x1_tot[wlim==-11]
            fagnx0x1_noagn_tot_j11[z_arr==z] = fagnx0x1_noagn_tot[wlim==-11]
            yerr_fagnx0x1_late_j11[:,z_arr==z] = np.array([yerr0_fagnx0x1_late[wlim==-11],yerr1_fagnx0x1_late[wlim==-11]])
            yerr_fagnx0x1_tot_j11[:,z_arr==z] = np.array([yerr0_fagnx0x1_tot[wlim==-11],yerr1_fagnx0x1_tot[wlim==-11]])
            yerr_fagnx0x1_noagn_tot_j11[:,z_arr==z] = np.array([yerr0_fagnx0x1_noagn_tot[wlim==-11],
                                                                yerr1_fagnx0x1_noagn_tot[wlim==-11]])

            fagnx0x1_late_mycut[z_arr==z] = fagnx0x1_late[wlim==-12]
            fagnx0x1_tot_mycut[z_arr==z] = fagnx0x1_tot[wlim==-12]
            fagnx0x1_noagn_tot_mycut[z_arr==z] = fagnx0x1_noagn_tot[wlim==-12]
            yerr_fagnx0x1_late_mycut[:,z_arr==z] = np.array([yerr0_fagnx0x1_late[wlim==-12],yerr1_fagnx0x1_late[wlim==-12]])
            yerr_fagnx0x1_tot_mycut[:,z_arr==z] = np.array([yerr0_fagnx0x1_tot[wlim==-12],yerr1_fagnx0x1_tot[wlim==-12]])
            yerr_fagnx0x1_noagn_tot_mycut[:,z_arr==z] = np.array([yerr0_fagnx0x1_noagn_tot[wlim==-12],
                                                                  yerr1_fagnx0x1_noagn_tot[wlim==-12]])



    print fagnx0x1_noagn_tot_w05.max()
    print fagnx0x1_noagn_tot_w08.max()
    print fagnx0x1_noagn_tot_j11.max()
    print fagnx0x1_noagn_tot_mycut.max()

    print "making contam plots"
    plt.clf()
    plt.cla()
    plt.close()
    mpl.rcParams.update({'font.size': 12})
    #fig=plt.figure(figsize=(5,3.5))
    fig=plt.figure(figsize=(8,6))
    #ylim=(0,1.02)
    ylim=(0,1.0)

    ax1=fig.add_subplot(221)
    ax1.set_xlim(0-0.05,z_arr.max()+0.05)
    ax1.set_ylim(ylim)
    ax1.set_xlabel('z')
    ax1.set_ylabel(r't$_{\rm WISE,AGNx0}$ / t$_{\rm WISE}$')
    print "plotting mycut"
    p11,=ax1.plot(z_arr,fagnx0x1_late_mycut,'*-',color='darkorange',lw=2.5,ms=9)
    ax1.errorbar(z_arr,fagnx0x1_late_mycut,yerr=yerr_fagnx0x1_late_mycut,
                ecolor='darkorange',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)            
    p01,=ax1.plot(z_arr,fagnx0x1_tot_mycut,'c',marker='o',ls='-',lw=2.5,ms=5)
    ax1.errorbar(z_arr,fagnx0x1_tot_mycut,yerr=yerr_fagnx0x1_tot_mycut,
                ecolor='c',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)            
    ax1.set_title(r'2-color (this work)',fontsize=11)
    ax1.legend((p01,p11), ('Total','Late+Post'), fontsize=10, loc='upper left', 
              numpoints=1, handletextpad=0.1)


    ax2=fig.add_subplot(222)
    ax2.set_xlim(0-0.05,z_arr.max()+0.05)
    ax2.set_ylim(ylim)
    ax2.set_xlabel('z')
    ax2.set_ylabel(r't$_{\rm WISE,AGNx0}$ / t$_{\rm WISE}$')
    print "plotting w08"
    p12,=ax2.plot(z_arr,fagnx0x1_late_w08,'r*-',lw=2,ms=9)
    ax2.errorbar(z_arr,fagnx0x1_late_w08,yerr=yerr_fagnx0x1_late_w08,
                ecolor='r',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)            
    p02,=ax2.plot(z_arr,fagnx0x1_tot_w08,'bo-',lw=2,ms=5)
    ax2.errorbar(z_arr,fagnx0x1_tot_w08,yerr=yerr_fagnx0x1_tot_w08,
                ecolor='b',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)            
    ax2.set_title(r'W1-W2>0.8',fontsize=11)
    ax2.legend((p02,p12), ('Total','Late+Post'), fontsize=10, loc='upper left', 
              numpoints=1, handletextpad=0.1)


    ax3=fig.add_subplot(223)
    ax3.set_xlim(0-0.05,z_arr.max()+0.05)
    ax3.set_ylim(ylim)
    ax3.set_xlabel('z')
    ax3.set_ylabel(r't$_{\rm WISE,AGNx0}$ / t$_{\rm WISE}$')
    print "plotting j11"
    p13,=ax3.plot(z_arr,fagnx0x1_late_j11,'k',marker='*',ls='-',lw=2.5,ms=9)
    ax3.errorbar(z_arr,fagnx0x1_late_j11,yerr=yerr_fagnx0x1_late_j11,
                ecolor='k',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)            
    p03,=ax3.plot(z_arr,fagnx0x1_tot_j11,'g',marker='o',ls='-',lw=2.5,ms=5)
    ax3.errorbar(z_arr,fagnx0x1_tot_j11,yerr=yerr_fagnx0x1_tot_j11,
                ecolor='g',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)            
    ax3.set_title(r'2-color (J11)',fontsize=11)
    ax3.legend((p03,p13), ('Total','Late+Post'), fontsize=10, loc='upper left', 
              numpoints=1, handletextpad=0.1)

    ax4=fig.add_subplot(224)
    ax4.set_xlim(0-0.05,z_arr.max()+0.05)
    ax4.set_ylim(ylim)
    ax4.set_xlabel('z')
    ax4.set_ylabel(r't$_{\rm WISE,AGNx0}$ / t$_{\rm WISE}$')
    print "plotting w05"
    p14,=ax4.plot(z_arr,fagnx0x1_late_w05,'m*-',ms=9)
    ax4.errorbar(z_arr,fagnx0x1_late_w05,yerr=yerr_fagnx0x1_late_w05,
                ecolor='m',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)            
    p04,=ax4.plot(z_arr,fagnx0x1_tot_w05,'o-',color='darkblue',ms=5)
    ax4.errorbar(z_arr,fagnx0x1_tot_w05,yerr=yerr_fagnx0x1_tot_w05,
                 ecolor='darkblue',ls='None',capsize=4,elinewidth=2,markeredgewidth=2)            
    ax4.legend((p04,p14), ('Total','Late+Post'), fontsize=10, loc='upper left', 
              numpoints=1, handletextpad=0.1)
    ax4.set_title(r'W1-W2>0.5',fontsize=11)


            
    fig.suptitle(r'',fontsize=12)
    fig.subplots_adjust(bottom=0.08,left=0.1,right=0.96,top=0.9,wspace=0.25,hspace=0.35)
    
    if len(subdir_arr)==1:
        fig.savefig(maindir+'/contamination_vs_z_%s%s.pdf'%(wflim_str,subdir_arr[0]))
    else:
        fig.savefig(maindir+'/contamination_vs_z_%s%s.pdf'%(wflim_str,subdir_arr))
    #plt.clf()
    #plt.cla()
    #plt.close()






def koss_comparison_plots(maindir='/oasis/projects/nsf/hvd115/lblecha/',
                          subdir_arr=None, d=None, cam=range(7), plot_3dsep=False, 
                          use_logbins=False, ylog=False, skip_snap0=True, tmax_postmrg=0.1, 
                          agnx0=False, errtype='mad',verbose=False):

    if subdir_arr:
        if not isinstance(subdir_arr,list): subdir_arr = [subdir_arr]
        #extra = '_'+'_'.join([s for s in subdir_arr]) 
        extra = '_%dsims' if len(subdir_arr)>1 else ''
        print "Using subdir list (%d sims):"%(len(subdir_arr))
        print subdir_arr

        #data = agn_metadata(maindir+'/'+subdir, ncam=ncam, skip_snap0=skip_snap0,
        #                    tmax_postmrg = tmax_postmrg, agnx0=agnx0)
        #plottitle = subdir
        #outdir = maindir+'/'+subdir
        #nsim = 1
    else:
        test_subdir_arr = ['q1_fg0.3_sunruns/all','q0.5_fg0.1_sunruns/all']
        subdir_arr = test_subdir_arr

    ncam = len(cam)
    nsim = len(subdir_arr)
    outdir = maindir if nsim>1 else maindir+'/'+subdir_arr[0]
    path_arr = ['%s/%s'%(maindir,sub) for sub in subdir_arr]

    dsep=10
    sep_bins = np.array([-3.0, 0, 3, 10, 30, 100, 300])
    #sep_bins = np.array([-3.0, 0, 1, 3, 10, 30, 100, 300])
    dlgsep = 0.5    
    lgsep_bins = np.array([-4.0, -2, -1, 0, 0.5, 1, 1.5, 2, 2.5]) 
    bw=0.3

    ### set limits for agn masks ###
    hilo_lagn_lims = [10**43.2, 10**44.2]

    hist3d_dt = np.zeros((nsim,len(sep_bins)-1))
    histproj_dt = np.zeros((nsim,ncam,(len(sep_bins)-1)))

    obsc_agn_hist3d_dt = np.zeros((nsim,4,len(sep_bins)-1))
    obsc_agn_hist3d_dt_frac = np.zeros((nsim,4,len(sep_bins)-1))
    obsc_agn_histproj_dt = np.zeros((nsim,4,ncam,(len(sep_bins)-1)))
    obsc_agn_histproj_dt_frac = np.zeros((nsim,4,ncam,(len(sep_bins)-1)))    
    

    ### load & process agn metadata ###
    for i_sim,path in enumerate(path_arr):

        print "\nloading agn metadata for sim %d (%s)..."%(i_sim,path)
        if not d:
            d = agn_metadata(path, ncam=ncam, skip_snap0=skip_snap0,
                             tmax_postmrg = tmax_postmrg, agnx0=agnx0)
    
        print d.dt.shape, d.posbh1.shape, d.lbol_grid.shape, d.bhprojsep.shape, d.w2w3.shape

        has2bh = (d.nbh==2) 
        ttot = d.dt.sum()

        meanw12 = np.mean(d.w1w2,axis=0)
        meanw23 = np.mean(d.w2w3,axis=0)

        j11mask = j11_wedge_cut(d.w1w2,d.w2w3,
                                w1lum_lim=w1lum_lim,w2lum_lim=w2lum_lim,w1lum=w1lum,w2lum=w2lum)
        j11mask_mean = j11_wedge_cut(meanw12,meanw23,
                                     w1lum_lim=w1lum_lim,w2lum_lim=w2lum_lim,w1lum=w1lum,w2lum=w2lum)


        w05mask = (d.w1w2>0.5)
        w05mask_mean = (meanw12>0.5)
        w08mask = (d.w1w2>0.8)
        w08mask_mean = (meanw12>0.8)
        
        ### calc histograms for each simulation ###
        print "calculating histograms for sim %d..."%i_sim

        dt_tiled = d.dt if plot_3dsep else np.tile(d.dt,(ncam,1))
        ttot_avg = d.dt.sum()/(1.0*nsim)
        lagn_tiled = d.lagn if plot_3dsep else np.tile(d.lagn,(ncam,1))
        Ngas_bh1 = np.nanmedian(d.Ngas_bh1_grid,axis=0) if plot_3dsep else d.Ngas_bh1_grid
        Ngas_bh2 = np.nanmedian(d.Ngas_bh2_grid,axis=0) if plot_3dsep else d.Ngas_bh2_grid
        Ngas_avg = np.nanmean(np.array([Ngas_bh1,Ngas_bh2]),axis=0)

        sepvar = np.log10(d.bhsep) if use_logbins else d.bhsep
        projsepvar = np.log10(d.bhprojsep) if use_logbins else d.bhprojsep
        
        if len(sep_bins[sep_bins<0])>1 or sep_bins.min()>d.bhprojsep.min():
            print "Error: invalid definition of sep_bins:",sep_bins
            sys.exit()                                           
        idx_sep_bins = np.digitize(d.bhsep,bins=sep_bins,right=True)
        print "idx_sep_bins:",idx_sep_bins.shape
        #print sep_bins.shape, d.bhsep.shape                
        idx_projsep_bins = np.array([ np.digitize(d.bhprojsep[i,:],bins=sep_bins,
                                                  right=True) for i in range(ncam) ])
        print "idx_projsep_bins:",idx_projsep_bins.shape
        #print sep_bins.shape, d.bhprojsep.shape
        #print sep_bins[idx_projsep_bins].shape
        ttot_sep_bins = np.histogram(d.bhsep,weights=d.dt,bins=sep_bins)[0]
        ttot_projsep_bins = np.array([ np.histogram(d.bhprojsep[i,:],weights=d.dt,
                                                    bins=sep_bins)[0] for i in range(ncam) ])
        #print sep_bins
        #print ttot_projsep_bins
        print "ttot = ",ttot.sum()
        print "ttot_sep_bins.sum() = ",ttot_sep_bins.sum()
        print "ttot_projsep_bins.sum(axis=1) =",ttot_projsep_bins.sum(axis=1)
        #print ttot_projsep_bins.shape
        #print idx_projsep_bins.shape

        idx_lgsep_bins = np.digitize(np.log10(d.bhsep),bins=lgsep_bins,right=True)
        idx_lgsep_bins[idx_lgsep_bins==len(lgsep_bins)] = 1
        idx_lgprojsep_bins = np.array([ np.digitize(np.log10(d.bhprojsep[i,:]),bins=lgsep_bins,
                                                    right=True) for i in range(ncam) ])
        idx_lgprojsep_bins[idx_lgprojsep_bins==len(lgsep_bins)] = 1

        ttot_lgsep_bins = np.histogram(np.log10(d.bhsep),weights=d.dt,bins=lgsep_bins)[0]
        ttot_lgprojsep_bins = np.array([ np.histogram(np.log10(d.bhprojsep[i,:]),weights=d.dt,
                                                      bins=lgsep_bins)[0] for i in range(ncam) ])

        dt_frac_sep_bins = d.dt/ttot_sep_bins[idx_sep_bins-1]
        dt_frac_projsep_bins = dt_tiled/np.array([ttot_projsep_bins[i,idx_projsep_bins[i,:]-1]
                                                  for i in range(ncam)])
        #print dt_frac_projsep_bins.shape
        dt_frac_lgsep_bins = np.zeros(d.dt.shape)
        dt_frac_lgsep_bins[has2bh] = d.dt[has2bh]/ttot_lgsep_bins[idx_lgsep_bins[has2bh]-1]
        dt_frac_lgprojsep_bins = np.zeros(dt_tiled.shape)
        dt_frac_lgprojsep_bins[:,has2bh] = dt_tiled[:,has2bh]/np.array([ttot_lgprojsep_bins[i,idx_lgprojsep_bins[i,has2bh]-1] for i in range(ncam)])

        hilo_lagn_maskarr = [(d.lagn<hilo_lagn_lims[0]),
                             ((d.lagn>=hilo_lagn_lims[0])&(d.lagn<hilo_lagn_lims[1])),
                             (d.lagn>=hilo_lagn_lims[1])]

        obsc_hilo_lagn_maskarr = [(lagn_tiled<hilo_lagn_lims[0]),
                                  ((lagn_tiled>=hilo_lagn_lims[0])&(lagn_tiled<hilo_lagn_lims[1])),
                                  ((lagn_tiled>=hilo_lagn_lims[1])&
                                   ((Ngas_bh1<1.0e22)|(Ngas_bh2<1.0e22))),
                                  ((lagn_tiled>=hilo_lagn_lims[1])&(Ngas_bh1>=1.0e22)&
                                   ((Ngas_bh2>=1.0e22)|(Ngas_bh2!=Ngas_bh2)))]

        hist3d_dt[i_sim,:] = np.histogram(sepvar, bins=sep_bins, weights=d.dt)[0]
        histproj_dt[i_sim,:,:] = np.array([np.histogram(projsepvar[j,:],
                                                        bins=sep_bins, weights=d.dt)[0]
                                           for j in range(ncam)])
        histproj_dt_frac[i_sim,:,:] = histproj_dt[i_sim,:,:]/ttot_projsep_bins

        for i_mask,mask in enumerate(obsc_hilo_lagn_maskarr):
            if plot_3dsep:
                obsc_agn_hist3d_dt[i_sim,i_mask,:] = np.histogram(sepvar[mask],bins=sep_bins,
                                                                  weights=d.dt[mask])[0]
                obsc_agn_hist3d_dt_frac[i_sim,i_mask,:] = obsc_agn_hist3d_dt[i_sim,i_mask,:]/ttot_sep_bins
            else:
                #print len(mask), len(mask[0]), len(mask[1])
                #print projsepvar.shape
                #print dt_tiled.shape
                #print [len(mask[j]) for j in range(ncam)]
                obsc_agn_histproj_dt[i_sim,i_mask,:,:]=np.array([np.histogram(projsepvar[j,mask[j]],bins=sep_bins, 
                                                                              weights=dt_tiled[j,mask[j]])[0] 
                                                                 for j in range(ncam)])
                obsc_agn_histproj_dt_frac[i_sim,i_mask,:,:] = obsc_agn_histproj_dt[i_sim,i_mask,:,:]/ttot_projsep_bins   

    #### end of loop over path_arr ####

         
    if plot_3dsep:
        tmp_dttot_ratio_3d = hist3d_dt[:,2]/hist3d_dt[:,1]
        tmp_ratio_3d = obsc_agn_hist3d_dt[:,:,2]/obsc_agn_hist3d_dt[:,:,1]

        print "\n\ntotal dt[3-10kpc]/dt[0-3kpc]:"
        print_arr_info(tmp_dttot_ratio_3d,nanvals=True)
        print tmp_dttot_ratio_3d
        print "\nobscured,lum agn: dt[3-10kpc]/dt[0-3kpc]:"
        print_arr_info(tmp_ratio_3d[:,3],nanvals=True)
        print tmp_ratio_3d[:,3]
        #print "unobscured,log(lagn)>=%g: dt[3-10kpc]/dt[0-3kpc]:"%hilo_lagn_lims[1]
        #print_arr_info(tmp_ratio_3d[:,2],nanvals=True)
        print "low-lum agn: dt[3-10kpc]/dt[0-3kpc]:"
        print_arr_info(tmp_ratio_3d[:,1],nanvals=True)
        print "inactive dt[3-10kpc]/dt[0-3kpc]"
        print_arr_info(tmp_ratio_3d[:,0],nanvals=True)
    
        tmp_ratio_3d = obsc_agn_hist3d_dt_frac[:,:,2]/obsc_agn_hist3d_dt_frac[:,:,1]

        print "\nobscured,lum agn: dt_frac[3-10kpc]/dt_frac[0-3kpc]:"
        print_arr_info(tmp_ratio_3d[:,3],nanvals=True)
        print tmp_ratio_3d[:,3]
        #print "unobscured,log(lagn)>=%g: dt_frac[3-10kpc]/dt_frac[0-3kpc]:"%hilo_lagn_lims[1]
        #print_arr_info(tmp_ratio_3d[:,2],nanvals=True)
        print "low-lum agn: dt_frac[3-10kpc]/dt_frac[0-3kpc]:"
        print_arr_info(tmp_ratio_3d[:,1],nanvals=True)
        print "inactive: dt_frac[3-10kpc]/dt_frac[0-3kpc]:"
        print_arr_info(tmp_ratio_3d[:,0],nanvals=True)

    else:

        tmp_dttot_ratio_proj = histproj_dt[:,:,2]/histproj_dt[:,:,1]
        tmp_ratio_proj = obsc_agn_histproj_dt[:,:,:,2]/obsc_agn_histproj_dt[:,:,:,1]

        print "\n\ntotal dt[3-10kpc]/dt[0-3kpc]:"
        print_arr_info(tmp_dttot_ratio_proj,nanvals=True)
        print tmp_dttot_ratio_proj
        print "\nobscured,lum agn: dt[3-10kpc]/dt[0-3kpc]:"
        print_arr_info(tmp_ratio_proj[:,3,:],nanvals=True)
        print tmp_ratio_proj[:,3,:]
        #print "unobscured,log(lagn)>=%g: dt[3-10kpc]/dt[0-3kpc]:"%hilo_lagn_lims[1]
        #print_arr_info(tmp_ratio_proj[:,2,:],nanvals=True)
        print "low-lum agn: dt[3-10kpc]/dt[0-3kpc]:"
        print_arr_info(tmp_ratio_proj[:,1,:],nanvals=True)
        print "inactive: dt[3-10kpc]/dt[0-3kpc]:"
        print_arr_info(tmp_ratio_proj[:,0,:],nanvals=True)

        tmp_ratio_proj = obsc_agn_histproj_dt_frac[:,:,:,2]/obsc_agn_histproj_dt_frac[:,:,:,1]

        print "\nobscured,lum agn: dt_frac[3-10kpc]/dt_frac[0-3kpc]:"
        print_arr_info(tmp_ratio_proj[:,3,:],nanvals=True)
        print tmp_ratio_proj[:,3,:]
        #print "unobscured,log(lagn)>=%g: dt_frac[3-10kpc]/dt_frac[0-3kpc]:"%hilo_lagn_lims[1]
        #print_arr_info(tmp_ratio_proj[:,2,:],nanvals=True)
        #print "%g<=log(lagn)<%g: dt_frac[3-10kpc]/dt_frac[0-3kpc]:"%(hilo_lagn_lims[0],hilo_lagn_lims[1])
        print "low-lum agn: dt_frac[3-10kpc]/dt_frac[0-3kpc]:"
        print_arr_info(tmp_ratio_proj[:,1,:],nanvals=True)
        print "inactive: dt_frac[3-10kpc]/dt_frac[0-3kpc]:"
        print_arr_info(tmp_ratio_proj[:,0,:],nanvals=True)



    ### separate high- and low-lum AGN phases:
    bw = 0.18
    plotbins = sep_bins[sep_bins<=latesep]
    nplotbins = len(plotbins)-1
    xticklbl_gap = ["post-merger"," "]+["%g-%g"%(sep_bins[j],sep_bins[j+1]) for j in np.arange(1,len(sep_bins)-1)]
    xticklbl = ["%g-%g"%(sep_bins[j],sep_bins[j+1]) for j in np.arange(len(sep_bins)-1)]
    #if not use_logbins: xticklbl[0] = '0'
    if not use_logbins: xticklbl[0] = 'post-merger'
    if plot_3dsep:
        xlabel = '3D sep [%skpc]'%('log ' if use_logbins else '')
    else:
        xlabel = 'proj. sep. [%skpc]'%('log ' if use_logbins else '')
    sepstr = '_bhsep' if plot_3dsep else '_bhprojsep'
    lgsepstr = '_lgbhsep' if plot_3dsep else '_lgbhprojsep'
    sepstring = lgsepstr if use_logbins else sepstr

    fig = plt.figure(figsize=(6,6))

    ax1 = fig.add_subplot(211)
        #ax1.set_ylim(0.0,1.35)
    ax1.set_xlim(-0.75,4.75)
    ax1.set_ylim(0.0,0.42)
    ax1.set_yticks(np.arange(5)*0.1)
    ax1.set_ylabel(r'$\Delta$t [Gyr]')
    ax1.set_xticks(np.append(0,np.arange(2,nplotbins+1)))
    ax1.set_xticklabels(xticklbl)

        #var = hist3d_dt[:nplotbins] if plot_3dsep else histproj_dt[:nplotbins]
    hiobsc_var = obsc_agn_hist3d_dt[:,3,:nplotbins] if plot_3dsep else obsc_agn_histproj_dt[:,3,:,:nplotbins]
    hiagn_var = obsc_agn_hist3d_dt[:,2,:nplotbins] if plot_3dsep else obsc_agn_histproj_dt[:,2,:,:nplotbins]
    loagn_var = obsc_agn_hist3d_dt[:,1,:nplotbins] if plot_3dsep else obsc_agn_histproj_dt[:,1,:,:nplotbins]
    noagn_var = obsc_agn_hist3d_dt[:,0,:nplotbins] if plot_3dsep else obsc_agn_histproj_dt[:,0,:,:nplotbins]
        #tot=cmph.makebar(ax1,var,nbins=len(plotbins),color='c',width=bw,xoffset=-3*bw,
        #                 val='median',errtype=errtype)
    hiobscagn=cmph.makebar(ax1,hiobsc_var,nbins=nplotbins,color='b',width=bw, xoffset=-1.5*bw,
                           xgap_index=0,val='median',errtype=errtype,nsim=nsim)
        #hiagn=cmph.makebar(ax1,hiagn_var,nbins=nplotbins,
        #                   color='gray', width=bw,xoffset=-bw,val='median',errtype=errtype)
    loagn=cmph.makebar(ax1,loagn_var,nbins=nplotbins,color='r',width=bw, xoffset=-0.5*bw,
                       xgap_index=0,val='median',errtype=errtype,nsim=nsim)
    noagn=cmph.makebar(ax1,noagn_var,nbins=nplotbins,color='g',width=bw,xoffset=0.5*bw,
                       xgap_index=0,val='median',errtype=errtype,nsim=nsim)

    ax1.set_xlabel(xlabel)
    ax1.legend((hiobscagn,loagn,noagn),
               ('Obsc. luminous AGN','Low-lum AGN','Inactive'),
               fontsize=9,loc='upper left')
    
    ax2 = fig.add_subplot(212)
    ax2.set_ylim(0.0,1.2)
    ax2.set_xlim(-0.75,4.75)
    ax2.set_ylabel(r'$\Delta$t/t$_{\rm bin}$',fontsize=11)
        #ax2.set_xticks(np.arange(nplotbins))
    ax2.set_xticks(np.append(0,np.arange(2,nplotbins+1)))
    ax2.set_xticklabels(xticklbl)
    hiobsc_var=obsc_agn_hist3d_dt_frac[:,3,:nplotbins] if plot_3dsep else \
        obsc_agn_histproj_dt_frac[:,3,:,:nplotbins]
    hiagn_var=obsc_agn_hist3d_dt_frac[:,2,:nplotbins] if plot_3dsep else \
        obsc_agn_histproj_dt_frac[:,2,:,:nplotbins]
    loagn_var=obsc_agn_hist3d_dt_frac[:,1,:nplotbins] if plot_3dsep else \
        obsc_agn_histproj_dt_frac[:,1,:,:nplotbins]
    noagn_var=obsc_agn_hist3d_dt_frac[:,0,:nplotbins] if plot_3dsep else \
        obsc_agn_histproj_dt_frac[:,0,:,:nplotbins]
    print "agn_hist3d_dt_frac(obscured,log(lagn)>=%g):"%hilo_lagn_lims[1]
    hiobscagn=cmph.makebar(ax2,hiobsc_var,nbins=nplotbins,color='b',width=bw,xoffset=-1.5*bw,
                           xgap_index=0,val='median',errtype=errtype,verbose=verbose,nsim=nsim)
        #print "agn_hist3d_dt_frac(log(lagn)>=%g):"%hilo_lagn_lims[1]
    hiagn=cmph.makebar(ax2,hiagn_var,nbins=nplotbins,color='gray',width=bw,xoffset=-bw,
                       val='median', errtype=errtype,verbose=verbose)
    print "agn_hist3d_dt_frac(%g<=log(lagn)<%g):"%(hilo_lagn_lims[0],hilo_lagn_lims[1])
    loagn=cmph.makebar(ax2,loagn_var,nbins=nplotbins,color='r',width=bw,xoffset=-0.5*bw,
                       xgap_index=0,val='median',errtype=errtype,verbose=verbose,nsim=nsim)
    print "agn_hist3d_dt_frac(log(lagn)<%g):"%hilo_lagn_lims[0]
    noagn=cmph.makebar(ax2,noagn_var,nbins=nplotbins,color='g',width=bw,xoffset=0.5*bw,
                       xgap_index=0,val='median',errtype=errtype,verbose=verbose,nsim=nsim)

    ax2.set_xlabel(xlabel)
    ax2.legend((hiobscagn,loagn,noagn),
               ('Obsc. luminous AGN','Low-lum AGN','Inactive'),
               fontsize=9,loc='upper right')
    figtitle=" ".join(subdir_arr) if nsim==1 else extra
    fig.suptitle(figtitle,fontsize=9) 
    
    fig.subplots_adjust(wspace=0.2,hspace=0.25,left=0.1,right=0.96,bottom=0.1,top=0.94)
    if nsim==1:
        plotname='%s/%s/merger_obsc_lum_agn_phases%s_%s_tpost%g.pdf'%(maindir,subdir_arr[0],sepstring,errtype,tmax_postmrg)
    else:
        plotname='%s/merger_obsc_lum_agn_phases%s_%s%s_tpost%g.pdf'%(maindir,sepstring,errtype,extra,tmax_postmrg)
    fig.savefig(plotname,format='pdf')
    plt.cla()
    plt.clf()
        

def make_all_barplots(maindir='/oasis/projects/nsf/hvd115/lblecha',name='fid',
                      tmax_postmrg=0.1,agnx0=False,xlbl_eachline=True,
                      sfrhist_w_grid=True,grid_nh_res=0.136,pubstyle=False):

    path_arr = define_sunrise_path_arr(basepath=maindir,name=name,
                                       sfrhist_w_grid=sfrhist_w_grid)

    for i,path in enumerate(path_arr):

        print "\n\n Making WISE AGN barplots for sunrise sim %s (%d of %d)..."%(path,i+1,len(path_arr))
        wise_agn_barplots(path,tmax_postmrg=tmax_postmrg,agnx0=agnx0,xlbl_eachline=xlbl_eachline,
                          sfrhist_w_grid=sfrhist_w_grid,grid_nh_res=grid_nh_res,pubstyle=pubstyle)


def wise_agn_barplots(path, d=None, cam=range(7), plot_3dsep=False, use_logbins=False,
                      ylog=False, skip_snap0=True, tmax_postmrg=None, agnx0=False,
                      xlbl_eachline=True,sfrhist_w_grid=True,grid_nh_res=0.136,pubstyle=False):


    ###
    # required input values: path (to sunrise sim directory) 
    # optional input: d (agn_metadata class instance, if not given will read from data in path)
    # optional input: cam (indices of sunrise camera array)
    # other optional input (needed to instantiate agn_metadata): skip_snap0, tmax_postmrg, agnx0
    # makes merger-phase bar plots for WISE AGN phases, outputs some relevant quantities
    ###

    ncam=len(cam)
    if not d:
        d = agn_metadata(path, ncam=ncam, skip_snap0=skip_snap0,
                         tmax_postmrg=tmax_postmrg, agnx0=agnx0)
    
    sepstr = '_bhsep' if plot_3dsep else '_bhprojsep'
    lgsepstr = '_lgbhsep' if plot_3dsep else '_lgbhprojsep'
    sepstring = lgsepstr if use_logbins else sepstr
    extra='_agnx0' if agnx0 else ''
    plottitle=path

    dsep=10
    sep_bins = np.arange(-dsep,200,dsep)
    dlgsep = 0.5    
    lgsep_bins = np.arange(-4,3,dlgsep)

    dt_tiled = d.dt if plot_3dsep else np.tile(d.dt,(ncam,1))

    idx_sep_bins = np.digitize(d.bhsep,bins=sep_bins,right=True)
    ttot_sep_bins = np.histogram(d.bhsep,weights=d.dt,bins=sep_bins)[0]
    idx_projsep_bins = np.array([ np.digitize(d.bhprojsep[i,:],bins=sep_bins,
                                                   right=True) for i in range(ncam) ])
    ttot_projsep_bins = np.array([ np.histogram(d.bhprojsep[i,:],weights=d.dt,
                                                     bins=sep_bins)[0] for i in range(ncam) ])

    idx_lgsep_bins = np.digitize(np.log10(d.bhsep),bins=lgsep_bins,right=True)
    idx_lgsep_bins[idx_lgsep_bins==len(lgsep_bins)] = 1
    ttot_lgsep_bins = np.histogram(np.log10(d.bhsep),weights=d.dt,
                                        bins=lgsep_bins)[0]
    idx_lgprojsep_bins = np.array([ np.digitize(np.log10(d.bhprojsep[i,:]),
                                                     bins=lgsep_bins,
                                                     right=True) for i in range(ncam) ])
    idx_lgprojsep_bins[idx_lgprojsep_bins==len(lgsep_bins)] = 1
    ttot_lgprojsep_bins = np.array([ np.histogram(np.log10(d.bhprojsep[i,:]),
                                                       weights=d.dt,bins=lgsep_bins)[0] 
                                          for i in range(ncam) ])

    dt_frac_sep_bins = d.dt/ttot_sep_bins[idx_sep_bins-1]
    dt_frac_projsep_bins = dt_tiled/np.array([ttot_projsep_bins[i,idx_projsep_bins[i,:]-1]
                                                   for i in range(ncam)])
    dt_frac_lgsep_bins = np.zeros(d.dt.shape)
    dt_frac_lgsep_bins[d.nbh==2] = d.dt[d.nbh==2]/ttot_lgsep_bins[idx_lgsep_bins[d.nbh==2]-1]
    dt_frac_lgprojsep_bins = np.zeros(dt_tiled.shape)
    dt_frac_lgprojsep_bins[:,d.nbh==2] = dt_tiled[:,d.nbh==2]/np.array([ttot_lgprojsep_bins[i,idx_lgprojsep_bins[i,d.nbh==2]-1] for i in range(ncam)])


    meanw12 = np.mean(d.w1w2,axis=0)
    meanw23 = np.mean(d.w2w3,axis=0)
    j11mask = j11_wedge_cut(d.w1w2,d.w2w3,
                            w1lum_lim=w1lum_lim,w2lum_lim=w2lum_lim,w1lum=w1lum,w2lum=w2lum)
    j11mask_mean = j11_wedge_cut(meanw12,meanw23,
                                 w1lum_lim=w1lum_lim,w2lum_lim=w2lum_lim,w1lum=w1lum,w2lum=w2lum)
    w05mask = (d.w1w2>0.5)
    w05mask_mean = (meanw12>0.5)
    w08mask = (d.w1w2>0.8)
    w08mask_mean = (meanw12>0.8)

    fedd_lims = [0.01,0.05,0.1]
    flbol_lims = [0.1, 0.3, 0.5]
    lagn_lims = [1.0e43, 1.0e44, 1.0e45]

    print "Making WISE AGN plots..."

    #print [ttot_projsep_bins[i,ttot_projsep_bins[i,:]>0].shape for i in range(ncam)]
    #for i in range(ncam):
    #    print len(sep_bins[:-1][ttot_projsep_bins[i,:]>0])
    #for k in range(len(sep_bins)-1):
    #    print sep_bins[k]
    #    print len(ttot_projsep_bins[:,k])
    #    print len(ttot_projsep_bins[ttot_projsep_bins[:,k]>0,k])

    ##wrong!
    #nproj_tnonzero_eachbin = np.array([len(sep_bins[:-1][ttot_projsep_bins[:,k]>0]) for k in range(len(sep_bins)-1)])
    nproj_tnonzero_eachbin = np.array([len(ttot_projsep_bins[ttot_projsep_bins[:,k]>0,k]) for k in range(len(sep_bins)-1)])
    #print nproj_tnonzero_eachbin

    nproj_tnonzero_eachlgbin = np.array([len(ttot_lgprojsep_bins[ttot_lgprojsep_bins[:,k]>0,k]) 
                                         for k in range(len(lgsep_bins)-1)])
    #print nproj_tnonzero_eachlgbin

    fedd_wmaskarr_list,dual_fedd_wmaskarr_list = (),()
    flbol_wmaskarr_list,dual_flbol_wmaskarr_list = (),()
    lagn_wmaskarr_list,dual_lagn_wmaskarr_list = (),()
    avg_wisemask_list = [w05mask_mean, w08mask_mean, j11mask_mean]
    wisemask_list = [w05mask, w08mask, j11mask]
        
    fedd_wmaskarr_list,dual_fedd_wmaskarr_list = set_agn_wmaskarr_list(d.fedd,d.fedd1,d.fedd2,
                                                                       fedd_lims,wisemask_list,old_fmt=True)
    flbol_wmaskarr_list,dual_flbol_wmaskarr_list = set_agn_wmaskarr_list(d.fagn_lbol,d.fagn1_lbol,d.fagn2_lbol,
                                                                         flbol_lims,wisemask_list,old_fmt=True)
    lagn_wmaskarr_list,dual_lagn_wmaskarr_list = set_agn_wmaskarr_list(d.lagn,d.lagn1,d.lagn2,
                                                                       lagn_lims,wisemask_list,old_fmt=True)


    mpl.rcParams.update({'font.size':11})
    wisemask_labels = ('W1W2>0.5','W1W2>0.8','J11 wedge')
    fedd_wise_titles = ['%s, fEdd>%g'%(lbl,lim) for lbl in wisemask_labels for lim in fedd_lims]
    #flbol_wise_titles = ['%s, fLbol>%g'%(lbl,lim) for lbl in wisemask_labels for lim in flbol_lims]
    flbol_wise_titles = [r'%s, f$_{\rm AGN}$>%g'%(lbl,lim) for lbl in wisemask_labels for lim in flbol_lims]
    lagn_wise_titles = ['%s, Lagn>%g'%(lbl,lim) for lbl in wisemask_labels for lim in lagn_lims]
    wise_agn_titles = (fedd_wise_titles, flbol_wise_titles, lagn_wise_titles)

    #wise_selection_labels = ['MIR AGN', 'MIR false neg', 'MIR false pos','']
    #wise_selection_labels = ['MIR AGN', 'AGN, not MIR', 'MIR low-lum AGN','']
    wise_selection_labels = ['MIR AGN', 'AGN, not MIR', '','']
    wise_duals_labels = ['MIR dual AGN', 'dual AGN, not MIR', 'MIR AGN, not dual','']
    tothist_kwargs = dict(range=(sep_bins.min(),sep_bins.max()),bins=len(sep_bins)-1,
                          color='k', histtype='step')
    #stackhist_kwargs = dict(range=(sep_bins.min(),sep_bins.max()), 
    #                        bins=len(sep_bins)-1, histtype='bar', stacked=True,
    #                        color=['b','m','orange','None'], lw=0.5)
    stackhist_kwargs = dict(range=(sep_bins.min(),sep_bins.max()), 
                            bins=len(sep_bins)-1, histtype='bar', stacked=True,
                            color=['b','m','None','None'], lw=1)
    stackhist_frac_kwargs = dict(range=(sep_bins.min(),sep_bins.max()), 
                                 bins=len(sep_bins)-1, histtype='bar', stacked=True,
                                 color=['darkblue','cyan','None'], lw=1)
    totlghist_kwargs = dict(range=(lgsep_bins.min(),lgsep_bins.max()),
                            bins=len(lgsep_bins)-1, color='k', histtype='step')
    stacklghist_kwargs = dict(range=(lgsep_bins.min(),lgsep_bins.max()), 
                              bins=len(lgsep_bins)-1, histtype='bar', stacked=True,
                              color=['b','m','orange','None'], lw=0.5)
    stacklghist_frac_kwargs = dict(range=(lgsep_bins.min(),lgsep_bins.max()), 
                                   bins=len(lgsep_bins)-1, histtype='bar', stacked=True,
                                   color=['b','m','orange',], lw=0.5)
    stack_kwargs = stacklghist_kwargs if use_logbins else stackhist_kwargs
    stack_frac_kwargs = stacklghist_frac_kwargs if use_logbins else stackhist_frac_kwargs

    if agnx0: extra = extra+'_agnx0'
    if ylog: extra = extra+'_ylog'

    xlab = '3D BH sep [kpc]' if plot_3dsep else 'proj. BH sep [kpc]'
    lgxlab = '3D BH sep [log kpc]' if plot_3dsep else 'proj. BH sep [log kpc]'
    xlabel = lgxlab if use_logbins else xlab

    plotnames = ['wise_agn_vs%s_%s%s_tpost%g'%(sepstring,s,extra,tmax_postmrg) 
                 for s in ['fedd','flbol','lagn']]
    frac_plotnames = ['wise_agn_vs%s_dtfrac_%s%s_tpost%g'%(sepstring,s,extra,tmax_postmrg) 
                      for s in ['fedd','flbol','lagn']]

    #xlim = (-1,2) if use_logbins else (-1.2*dsep,100) 
    xlim = (-1,2) if use_logbins else (-1.2*dsep,50) 
    #xlim = (-1,2) if use_logbins else (-1.2*dsep,150) 
    ylim = (0.5*d.dt[d.dt>0].min(),3.0) if ylog else (0.0,0.39)
    #ylim_frac = (-3,0.1) if ylog else (0.0,1.15)
    ylim_frac = (-3,0.1) if ylog else (0.0,1.05)

    for ilist,this_maskarr_list in enumerate((fedd_wmaskarr_list, flbol_wmaskarr_list, 
                                              lagn_wmaskarr_list)):

        #plt.close()
        plt.cla()
        plt.clf()
        fig = plt.figure(figsize=(8,8))
        for i in range(9):
            
            ax = fig.add_subplot(331+i)
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            if (i in [6,7,8]) or (xlbl_eachline):
                ax.set_xlabel(xlabel)
            else:  ax.set_xlabel('')
            #ax.set_ylabel(r'$\Delta$t [Gyr]') if  i in [0,3,6] else ax.set_ylabel('')
            ax.set_ylabel(r'$\Delta$t [Gyr]')
            if plot_3dsep:
                histdata = [np.log10(d.bhsep[mask]) if use_logbins else d.bhsep[mask] 
                            for mask in this_maskarr_list[i]]
                #wts = [d.dt[mask]/(1.0*nsim) for mask in this_maskarr_list[i]]
                wts = [d.dt[mask] for mask in this_maskarr_list[i]]
            else: 
                #tmphist = np.zeros((len(this_maskarr_list[i]),len(sep_bins)-1))
                #for j,mask in enumerate(this_maskarr_list[i]):
                #    wts=dt_frac_projsep_bins[mask]/(1.0*nproj_tnonzero_eachbin*nsim)
                #    tmphist[j,:],tmp=np.histogram(d.bhprojsep[mask],bins=sep_bins,weights=wts)
                if use_logbins:
                    histdata = [np.log10(d.bhprojsep[mask]) for mask in this_maskarr_list[i]]
                    #wts = [dt_tiled[mask]/(1.0*nproj_tnonzero_eachlgbin[idx_lgprojsep_bins[mask]-1]*nsim) 
                    wts = [dt_tiled[mask]/(1.0*nproj_tnonzero_eachlgbin[idx_lgprojsep_bins[mask]-1]) 
                           for mask in this_maskarr_list[i]]
                else:
                    histdata = [d.bhprojsep[mask] for mask in this_maskarr_list[i]]
                    wts = [dt_tiled[mask]/(1.0*nproj_tnonzero_eachbin[idx_projsep_bins[mask]-1]) 
                           for mask in this_maskarr_list[i]]
            if np.max([d.bhprojsep[mask].shape for mask in this_maskarr_list[i]][:-1]) > 0:
                n,bins,patches=ax.hist(histdata, weights=wts, log=ylog, 
                                       label=wise_selection_labels, **stack_kwargs)
            else:
                print "no data for mask %d. generating empty plot and continuing..."%i
            ax.legend(fontsize=9,loc='upper left')
            ax.set_title((wise_agn_titles[ilist])[i],fontsize=10)    
            if not agnx0 and 'fLbol>0.1' in wise_agn_titles[ilist][i]:
                print "\n%s"%wise_agn_titles[ilist][i]
                print "bins:",bins[:6]
                print "WISE AGN [Myr] (total=%.4g; total(<50kpc)=%.4g):"%(n[0].sum()*1000,n[0][:6].sum()*1000)
                print "   ",n[0][:6]*1000
                print "all AGN [Myr] (total=%.4g; total(<50kpc)=%.4g):"%(n[1].sum()*1000,n[1][:6].sum()*1000)
                print "   ",n[1][:6]*1000
                print "MIR AGN false neg [Myr] (total=%.4g, total(<50kpc)=%.4g):"%((n[1]-n[0]).sum()*1000,(n[1][:6]-n[0][:6]).sum()*1000)
                print "  ",(n[1][:6]-n[0][:6])*1000
                print "MIR AGN false pos [Myr] (total=%.4g, total(<50kpc)=%.4g):"%((n[2]-n[1]).sum()*1000,(n[2][:6]-n[1][:6]).sum()*1000)
                print "  ",(n[2][:6]-n[1][:6])*1000

        fig.suptitle(plottitle)
        fig.subplots_adjust(wspace=0.34,hspace=0.3,left=0.08,
                            right=0.97,bottom=0.06,top=0.92)
        #fig.subplots_adjust(wspace=0.3,hspace=0.3,left=0.08,
        #                    right=0.94,bottom=0.06,top=0.92)
        fig.savefig("%s/%s.pdf"%(path,plotnames[ilist]),format='pdf')

    
        #print np.array(n).shape
        #print np.array(n)[-1,:]
        #print ttot_projsep_bins
        #return

        plt.close()
        plt.cla()
        plt.clf()
        fig = plt.figure(figsize=(8,8))
        for i in range(9):
            ax = fig.add_subplot(331+i)
            ax.set_xlim(xlim)
            ax.set_ylim(ylim_frac)
            if (i in [6,7,8]) or (xlbl_eachline):
                ax.set_xlabel(xlabel)
            else: ax.set_xlabel('')
            #ax.set_ylabel(r'$\Delta$t/t$_{\rm bin}$') if  i in [0,3,6] else ax.set_ylabel('')
            ax.set_ylabel(r'$\Delta$t/t$_{\rm bin}$')
            if plot_3dsep:
                histdata = [np.log10(d.bhsep[mask]) if use_logbins else d.bhsep[mask] 
                            for mask in this_maskarr_list[i]]
                wts = [dt_frac_sep_bins[mask] if use_logbins else dt_frac_lgsep_bins[mask] 
                       for mask in this_maskarr_list[i]]
                ### TO DO: doesn't work yet for log bins: ###
                ### TO DO: doesn't work for agnx0 (empty arrays, this 'fix' doesn't work, 
                ###        array has wrong shape)
                #if len(histdata)==0 and len(wts)==0: 
                #    histdata = [np.zeros(len(this_maskarr_list[i][:-1]))]
                #    wts = [np.zeros(len(this_maskarr_list[i][:-1]))]
            else: 
                if use_logbins:
                    histdata = [np.log10(d.bhprojsep[mask]) for mask in this_maskarr_list[i]][:-1]
                    wts = [dt_frac_lgprojsep_bins[mask] / \
                               (1.0*nproj_tnonzero_eachlgbin[idx_lgprojsep_bins[mask]-1]) 
                           for mask in this_maskarr_list[i]][:-1]
                else:
                    histdata = [d.bhprojsep[mask] for mask in this_maskarr_list[i]][:-1]
                    wts = [dt_frac_projsep_bins[mask] / \
                               (1.0*nproj_tnonzero_eachbin[idx_projsep_bins[mask]-1]) 
                           for mask in this_maskarr_list[i]][:-1]
            
            if np.max([d.bhprojsep[mask].shape for mask in this_maskarr_list[i]][:-1]) > 0:
                n,bins,patches=ax.hist(histdata, weights=wts, log=ylog, 
                                       label=wise_selection_labels, **stack_frac_kwargs)
                for k in range(len(patches[0])): patches[0][k].set_hatch('xxxx')
                for k in range(len(patches[1])): patches[1][k].set_hatch('////')
                #for k in range(len(patches[2])): patches[2][k].set_hatch('\\\\')
            else:
                print "no data for mask %d. generating empty plot and continuing..."%i
            print len(histdata), len(wts), len(patches), len(patches[0])
            print patches[0][0]
            ax.legend(fontsize=9,loc='upper right')
            if pubstyle:
                ax.set_title("%s (%s)"%(get_sim_name(path.split('lblecha/')[1]),
                                        (wise_agn_titles[ilist][i]).split(',')[0]),fontsize=10)
            else:
                ax.set_title((wise_agn_titles[ilist])[i],fontsize=10)    
            if not agnx0 and 'fLbol>0.1' in wise_agn_titles[ilist][i]:
                print "\n%s"%wise_agn_titles[ilist][i]
                print "bins:",bins[:6]
                print "WISE AGN:",n[0][:6]
                print "total AGN:", n[1][:6]
                print "MIR AGN false neg:", (n[1][:6]-n[0][:6])
                print "MIR AGN false pos:", (n[2][:6]-n[1][:6])

        if not pubstyle: fig.suptitle(plottitle)
        fig.subplots_adjust(wspace=0.34,hspace=0.36,left=0.08,
                            right=0.97,bottom=0.06,top=0.94)
        fig.savefig("%s/%s.pdf"%(path,frac_plotnames[ilist]),format='pdf')

    stack_kwargs.update(color=['darkblue','green','cyan','None']) 
    stack_frac_kwargs.update(color=['darkblue','green','cyan']) 

    if agnx0: return

    ### dual AGN plots ###
    plotnames = ['wise_duals_vs%s_%s%s_tpost%g'%(sepstring,s,extra,tmax_postmrg) 
                 for s in ['fedd','flbol','lagn']]
    frac_plotnames = ['wise_duals_vs%s_dtfrac_%s%s_tpost%g'%(sepstring,s,extra,tmax_postmrg) 
                      for s in ['fedd','flbol','lagn']]
    for ilist,this_maskarr_list in enumerate((dual_fedd_wmaskarr_list, dual_flbol_wmaskarr_list, 
                                              dual_lagn_wmaskarr_list)):
        #plt.close()
        plt.cla()
        plt.clf()
        fig = plt.figure(figsize=(8,8))
        for i in range(9):
            ax = fig.add_subplot(331+i)
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            ax.set_xlabel(xlabel) if i in [6,7,8] else ax.set_xlabel('')
            #ax.set_ylabel(r'$\Delta$t [Gyr]') if  i in [0,3,6] else ax.set_ylabel('')
            ax.set_ylabel(r'$\Delta$t [Gyr]') 
            if plot_3dsep:
                histdata = [np.log10(d.bhsep[mask]) if use_logbins else d.bhsep[mask] 
                            for mask in this_maskarr_list[i]]
                #wts = [d.dt[mask]/(1.0*nsim) for mask in this_maskarr_list[i]]
                wts = [d.dt[mask] for mask in this_maskarr_list[i]]
            else:
                if use_logbins:
                    histdata = [np.log10(d.bhprojsep[mask]) for mask in this_maskarr_list[i]]
                    #wts = [dt_tiled[mask]/(1.0*nproj_tnonzero_eachlgbin[idx_lgprojsep_bins[mask]-1]*nsim) 
                    wts = [dt_tiled[mask]/(1.0*nproj_tnonzero_eachlgbin[idx_lgprojsep_bins[mask]-1]) 
                           for mask in this_maskarr_list[i]]
                else:
                    histdata = [d.bhprojsep[mask] for mask in this_maskarr_list[i]]
                    #wts = [dt_tiled[mask]/(1.0*nproj_tnonzero_eachbin[idx_projsep_bins[mask]-1]*nsim) 
                    wts = [dt_tiled[mask]/(1.0*nproj_tnonzero_eachbin[idx_projsep_bins[mask]-1]) 
                           for mask in this_maskarr_list[i]]
            if np.max([d.bhprojsep[mask].shape for mask in this_maskarr_list[i]][:-1]) > 0:
                n,bins,patches=ax.hist(histdata, weights=wts, log=ylog, 
                                       label=wise_duals_labels, **stack_kwargs)
            else:
                print "no data for mask %d. generating empty plot and continuing..."%i
            ax.legend(fontsize=9,loc='upper left')
            ax.set_title((wise_agn_titles[ilist])[i],fontsize=10)    
            if 'fLbol>0.1' in wise_agn_titles[ilist][i]:
                print "\n%s"%wise_agn_titles[ilist][i]
                print "bins:",bins[:6]
                print "WISE dual AGN [Myr] (total=%.4g; total(<50kpc)=%.4g):"%(n[0].sum()*1000,n[0][:6].sum()*1000)
                print "   ",n[0][:6]*1000
                print "all dual AGN [Myr] (total=%.4g; total(<50kpc)=%.4g):"%(n[1].sum()*1000,n[1][:6].sum()*1000)
                print "   ",n[1][:6]*1000
                print "MIR AGN, not dual [Myr] (total=%.4g, total(<50kpc)=%.4g):"%((n[2]-n[1]).sum()*1000,(n[2][:6]-n[1][:6]).sum()*1000)
                print "  ",(n[2][:6]-n[1][:6])*1000

        fig.suptitle(plottitle)
        #fig.subplots_adjust(wspace=0.3,hspace=0.3,left=0.08,
        #                    right=0.94,bottom=0.06,top=0.92)
        fig.subplots_adjust(wspace=0.34,hspace=0.3,left=0.08,
                            right=0.97,bottom=0.06,top=0.92)
        fig.savefig("%s/%s.pdf"%(path,plotnames[ilist]),format='pdf')


        #plt.close()
        plt.cla()
        plt.clf()
        fig = plt.figure(figsize=(8,8))
        for i in range(9):
            ax = fig.add_subplot(331+i)
            ax.set_xlim(xlim)
            ax.set_ylim(ylim_frac)
            #ax.set_ylim(0.0,1.6)
            ax.set_xlabel(xlabel) if i in [6,7,8] else ax.set_xlabel('')
            #ax.set_ylabel(r'$\Delta$t/t$_{\rm bin}$') if  i in [0,3,6] else ax.set_ylabel('')
            ax.set_ylabel(r'$\Delta$t/t$_{\rm bin}$')
            if plot_3dsep:
                histdata = [np.log10(d.bhsep[mask]) if use_logbins else d.bhsep[mask] 
                            for mask in this_maskarr_list[i][:-1]]
                wts = [dt_frac_sep_bins[mask] for mask in this_maskarr_list[i][:-1]]
            else:
                if use_logbins:
                    histdata = [np.log10(d.bhprojsep[mask]) for mask in this_maskarr_list[i]][:-1]
                    wts = [dt_frac_lgprojsep_bins[mask] / \
                               (1.0*nproj_tnonzero_eachlgbin[idx_lgprojsep_bins[mask]-1]) 
                           for mask in this_maskarr_list[i]][:-1]
                else:
                    histdata = [d.bhprojsep[mask] for mask in this_maskarr_list[i]][:-1]
                    wts = [dt_frac_projsep_bins[mask] / \
                               (1.0*nproj_tnonzero_eachbin[idx_projsep_bins[mask]-1]) 
                           for mask in this_maskarr_list[i]][:-1]
                    
            if np.max([d.bhprojsep[mask].shape for mask in this_maskarr_list[i]][:-1]) > 0:
                n,bins,patches=ax.hist(histdata, weights=wts, log=ylog, 
                                       label=wise_duals_labels, **stack_frac_kwargs)
            else:
                print "no data for mask %d. generating empty plot and continuing..."%i

            ax.legend(fontsize=9,loc='upper right')
            ax.set_title((wise_agn_titles[ilist])[i],fontsize=10)    
        fig.suptitle(plottitle)
        #fig.subplots_adjust(wspace=0.3,hspace=0.3,left=0.08,
        #                    right=0.94,bottom=0.06,top=0.92)
        fig.subplots_adjust(wspace=0.34,hspace=0.3,left=0.08,
                            right=0.97,bottom=0.06,top=0.92)
        fig.savefig("%s/%s.pdf"%(path,frac_plotnames[ilist]),format='pdf')


def wise_agn_multisim(path, d=None, cam=range(7), plot_3dsep=False, use_logbins=False,
                      ylog=False, skip_snap0=True, tmax_postmrg=None, agnx0=False,
                      xlbl_eachline=True,sfrhist_w_grid=True,grid_nh_res=0.136,pubstyle=False):


    ###
    # required input values: path (to sunrise sim directory) 
    # optional input: d (agn_metadata class instance, if not given will read from data in path)
    # optional input: cam (indices of sunrise camera array)
    # other optional input (needed to instantiate agn_metadata): skip_snap0, tmax_postmrg, agnx0
    # makes merger-phase bar plots for WISE AGN phases, outputs some relevant quantities
    ###

    ncam=len(cam)
    if not d:
        d = agn_metadata(path, ncam=ncam, skip_snap0=skip_snap0,
                         tmax_postmrg=tmax_postmrg, agnx0=agnx0)


    twise_
        


def contour_plot(maindir='/oasis/projects/nsf/hvd115/lblecha/',
                 subdir='q1_fg0.3_allrx10_sunruns/all',
                 agnx0=False, tmax_postmrg=0.1, skip_snap0=True, 
                 cam=range(7),use_logbins=False, ylog=False,
                 extra='',verbose=False):

    ncam = len(cam)
    path = "%s/%s"%(maindir,subdir)

    d = agn_metadata(path, ncam=ncam, skip_snap0=skip_snap0,
                     tmax_postmrg = tmax_postmrg, agnx0=agnx0)
    
    
    print d.dt.shape, d.posbh1.shape, d.lbol_grid.shape, d.bhprojsep.shape, d.w2w3.shape
    
    has2bh = (d.nbh==2)
    bhprojsep_mean = d.bhprojsep.mean(axis=0)
    print bhprojsep_mean.shape
    print bhprojsep_mean
    ttot = d.dt.sum()
    dt_n0 = copy(d.dt)
    dt_n0 = np.array([d.dt[i] if d.dt[i]>0 else d.dt[i-1] for i in range(len(d.dt))])
    #dt_n0_tiled = np.tile(dt_n0,(ncam,1))
    dt_n0_tiled = np.tile(dt_n0,ncam)
    dt_tiled = np.tile(d.dt,ncam)
    ttot_tiled = dt_tiled.sum()/(1.0*ncam)

    meanw12 = np.mean(d.w1w2,axis=0)
    meanw23 = np.mean(d.w2w3,axis=0)
    
    j11mask = j11_wedge_cut(d.w1w2,d.w2w3,
                            w1lum_lim=w1lum_lim,w2lum_lim=w2lum_lim,w1lum=w1lum,w2lum=w2lum)
    j11mask_mean = j11_wedge_cut(meanw12,meanw23,
                                 w1lum_lim=w1lum_lim,w2lum_lim=w2lum_lim,w1lum=w1lum,w2lum=w2lum)
    
    w05mask = (d.w1w2>0.5)
    w05mask_mean = (meanw12>0.5)
    w08mask = (d.w1w2>0.8)
    w08mask_mean = (meanw12>0.8)
    
    ### calc histograms for each simulation ###
    
    #dt_tiled = d.dt if plot_3dsep else np.tile(d.dt,(ncam,1))
    #ttot_avg = d.dt.sum()/(1.0*nsim)
    #lagn_tiled = np.tile(d.lagn,(ncam,1))
    lagn_tiled = np.tile(d.lagn,ncam)
    lagn_ratio_tiled = np.tile(lagn_ratio,ncam)
    
    sepvar = np.log10(d.bhsep) if use_logbins else d.bhsep
    projsepvar = np.log10(d.bhprojsep) if use_logbins else d.bhprojsep

    my_cmap=copy(mpl.cm.get_cmap('hot'))
    #my_cmap=copy(mpl.cm.get_cmap('jet'))
    print my_cmap
    #my_cmap=copy(mpl.cm.get_cmap('YlOrRd'))                                                     
    my_cmap.set_bad('g')
    my_cmap.set_over('b')
    my_cmap.set_under('w')

    print 'lagn:',np.log10(d.lagn)
    #print np.log10(lagn_tiled)
    #print np.log10(lagn_tiled).reshape(ncam*len(d.lagn))
    #print d.bhprojsep.reshape(ncam*len(d.lagn))
    print 'bhprojsep_mean:',bhprojsep_mean
    print 'bhprojsep_mean[lagn>45]:',bhprojsep_mean[np.log10(d.lagn)>45]
    #print dt_n0_tiled
    bhprojsep_flat = copy(d.bhprojsep).flatten()

    print np.tile(d.time,ncam)[np.log10(lagn_tiled)>45]
    print dt_tiled[np.log10(lagn_tiled)>45]
    print d.bhprojsep.flatten()[np.log10(lagn_tiled)>45]
    print lagn_tiled[np.log10(lagn_tiled)>45]
    rlim=(-10,130)
    llim=(41,47)
    lratlim=(-2.2,2.2)
    dtmax=180
    print dt_tiled.min(),dt_tiled.max()
    print d.dt.min(), d.dt.max()
    print (1000*dt_tiled/(1.0*ncam)).min()
    print (1000*dt_tiled/(1.0*ncam)).max()
    print (1000*d.dt).min()
    print (1000*d.dt).max()
    nbins=20
    lbol_hist2d,ledg,redg = np.histogram2d(np.log10(lagn_tiled),bhprojsep_flat,
                                           bins=nbins,weights=(1000*dt_tiled/(1.0*ncam)),range=(llim,rlim))
    #lbol_hist2d,ledg,redg = np.histogram2d(np.log10(d.lagn),d.bhsep,
    #                                       bins=nbins,weights=d.dt,range=(llim,rlim))
    #lbol_hist2d,ledg,redg = np.histogram2d(np.log10(d.lagn),bhprojsep_mean,
    #                                       bins=nbins,weights=d.dt,range=(llim,rlim))
    print d.lagn.shape, lagn_tiled.shape, bhprojsep_mean.shape
    print lbol_hist2d.shape, ledg.shape, redg.shape
    print lbol_hist2d.max()
    print lbol_hist2d

    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(221)
    im = ax.imshow(lbol_hist2d,cmap=my_cmap,interpolation='nearest',origin='lower',
                   aspect=(rlim[1]-rlim[0])/(llim[1]-llim[0]),
                   extent=[redg[0],redg[-1],ledg[0],ledg[-1]],
                   vmin=0,vmax=dtmax)
                   #vmin=0,vmax=1.1*lbol_hist2d.max())
    cb = fig.colorbar(im, orientation='vertical')                                   

    lbol_hist2d,ledg,redg = np.histogram2d(np.log10(d.lagn),bhprojsep_mean,
                                           bins=nbins,weights=(1000*d.dt),range=(llim,rlim))
    print lbol_hist2d.max()
    
    ax = fig.add_subplot(222)
    im = ax.imshow(lbol_hist2d,cmap=my_cmap,interpolation='nearest',origin='lower',
                   aspect=(rlim[1]-rlim[0])/(llim[1]-llim[0]),
                   extent=[redg[0],redg[-1],ledg[0],ledg[-1]],
                   vmin=0,vmax=dtmax)
                   #vmin=0,vmax=1.1*lbol_hist2d.max())
    #ax.set_xlim(-10,130)
    #ax.set_ylim(41,47)
    #plt.plot(bhprojsep_flat,np.log10(lagn_tiled),'.')
    ##cbar_ax = fig.add_axes([0.12, 0.12, 0.78, 0.04])
    #cbar_ax = fig.add_axes([0.9, 0.15, 0.02, 0.7])                                             
    ##cb = fig.colorbar(im, cax=cbar_ax,orientation='horizontal')                                  
    cb = fig.colorbar(im, orientation='vertical')                                   

    #lbol_ratio_hist2d,lratedg,redg = np.histogram2d(np.log10(lagn_ratio),d.bhsep,
    #                                             bins=nbins,weights=d.dt,range=(lratlim,rlim))
    lbol_ratio_hist2d,lratedg,redg = np.histogram2d(np.log10(lagn_ratio_tiled),bhprojsep_flat,
                                                 bins=nbins,weights=(1000*dt_tiled/(1.0*ncam)),range=(lratlim,rlim))
    print lbol_ratio_hist2d.max()
    
    ax = fig.add_subplot(223)
    im = ax.imshow(lbol_ratio_hist2d,cmap=my_cmap,interpolation='nearest',origin='lower',
                   aspect=(rlim[1]-rlim[0])/(lratlim[1]-lratlim[0]),
                   extent=[redg[0],redg[-1],lratedg[0],lratedg[-1]],
                   vmin=0,vmax=dtmax)
                   #vmin=0,vmax=1.1*lbol_ratio_hist2d.max())
    cb = fig.colorbar(im, orientation='vertical')                                   

    lbol_ratio_hist2d,lratedg,redg = np.histogram2d(np.log10(lagn_ratio),bhprojsep_mean,
                                                 bins=nbins,weights=(1000*d.dt),range=(lratlim,rlim))
    print lbol_ratio_hist2d.max()
    
    ax = fig.add_subplot(224)
    im = ax.imshow(lbol_ratio_hist2d,cmap=my_cmap,interpolation='nearest',origin='lower',
                   aspect=(rlim[1]-rlim[0])/(lratlim[1]-lratlim[0]),
                   extent=[redg[0],redg[-1],lratedg[0],lratedg[-1]],
                   vmin=0,vmax=dtmax)
                   #vmin=0,vmax=1.1*lbol_ratio_hist2d.max())
    cb = fig.colorbar(im, orientation='vertical')                                   

    fig.savefig(path+"/hist_lbol_bhsep.pdf")
    #im.show()
    plt.clf()
    plt.cla()
    plt.close()




def wise_tscale_stacked_hist(bhsep, dt, agn_maskarr_list, nproj_tnonzero_eachbin=0, idx_projsep_bins=-1, 
                             use_logbins=False, plot_3dsep=False, xlim=(), ylim=()):

    #### WARNING: just a skeleton of a function thus far. ####

    ###
    xlabel = xxx
    ###

    plt.cla()
    plt.clf()
    fig = plt.figure(figsize=(8,8))
    for i in range(9):
        ax = fig.add_subplot(331+i)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim_frac)
        ax.set_xlabel(xlabel) if i in [6,7,8] else ax.set_xlabel('')
        ax.set_ylabel(r'$\Delta$t/t$_{\rm bin}$') if  i in [0,3,6] else ax.set_ylabel('')
        if plot_3dsep:
            histdata = [np.log10(d.bhsep[mask]) if use_logbins else d.bhsep[mask] 
                        for mask in this_maskarr_list[i][:-1]]
            wts = [dt_frac_sep_bins[mask] for mask in this_maskarr_list[i][:-1]]
        else:
            if use_logbins:
                histdata = [np.log10(d.bhprojsep[mask]) for mask in this_maskarr_list[i]][:-1]
                wts = [dt_frac_lgprojsep_bins[mask] / \
                           (1.0*nproj_tnonzero_eachlgbin[idx_lgprojsep_bins[mask]-1]) 
                       for mask in this_maskarr_list[i]][:-1]
            else:
                histdata = [d.bhprojsep[mask] for mask in this_maskarr_list[i]][:-1]
                wts = [dt_frac_projsep_bins[mask] / \
                           (1.0*nproj_tnonzero_eachbin[idx_projsep_bins[mask]-1]) 
                       for mask in this_maskarr_list[i]][:-1]
                    
        if np.max([d.bhprojsep[mask].shape for mask in this_maskarr_list[i]][:-1]) > 0:
            n,bins,patches=ax.hist(histdata, weights=wts, log=ylog, 
                                   label=wise_duals_labels, **stack_frac_kwargs)
        else:
            print "no data for mask %d. generating empty plot and continuing..."%i

        ax.legend(fontsize=9,loc='upper right')
        ax.set_title((wise_agn_titles[ilist])[i],fontsize=10)    

    fig.suptitle(plottitle)
    fig.subplots_adjust(wspace=0.3,hspace=0.3,left=0.08,
                        right=0.94,bottom=0.06,top=0.92)
    fig.savefig("%s/%s.pdf"%(path,frac_plotnames[ilist]),format='pdf')




def plot_wise_agn_cuts(ax,xlim=(-5,10),ylim=(-0.3,2.0),c='gray',new_c='b',
                       include_new=True,alt_wedge=True):

    ## Jarrett et al. 2011 criteria:
    ax.plot([2.2,4.2],[1.7,1.7],color=c,ls='--',lw=1.5)
    ax.plot([2.2,2.2],[0.6,1.7],color=c,ls='--',lw=1.5)
    ax.plot([4.2,4.2],[0.8,1.7],color=c,ls='--',lw=1.5)
    ax.plot([2.2,4.2],[0.6,0.8],color=c,ls='--',lw=1.5)
    ## Stern et al. criteria:
    ax.plot(xlim,[0.5,0.5],color=c,ls=':',lw=1.5)
    #ax.plot(xlim,[0.8,0.8],color=c)
    ax.plot(xlim,[0.8,0.8],color=c,ls='-.',lw=1.5)

    if alt_wedge:
        ## My alternate proposed criteria:
        ###ax.plot([4.61,5.39],[0.5,1.9],color=new_c)
        ###ax.plot([4.9,5.6],[0.5,1.9],color=new_c)
        ##ax.plot([xlim[0],4.8],[0.5,0.5],color=new_c)
        ##ax.plot([4.8,5.5],[0.5,1.9],color=new_c)
        #ax.plot([xlim[0],4.7],[0.5,0.5],color=new_c)

        ax.plot([2.2,2.2],[0.5,ylim[1]],color=new_c)
        ax.plot([2.2,4.7],[0.5,0.5],color=new_c)
        ax.plot([4.7,5.4],[0.5,1.9],color=new_c)
    else:
        ## My proposed criteria:
        ax.plot([xlim[0],4.9],[0.5,0.5],color=new_c)
        ax.plot([4.9,4.9],[0.5,3],color=new_c)

def j11_wedge_cut(w1w2,w2w3,w1lum_lim=0.0,w2lum_lim=0.0,
                 w1lum=np.array([]),w2lum=np.array([])):
    
    assert w1w2.shape == w2w3.shape, \
        "array size mismatch in j11_wedge_cut: %g %g"%(w1w2.shape, w2w3.shape)

    if w1lum_lim>0.0 or w2lum_lim>0.0:
        mask = ((w2w3>=2.2)&(w2w3<=4.2)&(w1w2<=1.7)&(w1w2>=0.1*w2w3+0.38)&(w1lum>=w1lum_lim)&(w2lum>=w2lum_lim))
    else:
        mask = ((w2w3>=2.2)&(w2w3<=4.2)&(w1w2<=1.7)&(w1w2>=0.1*w2w3+0.38))

    return mask

def my_wedge_cut(w1w2,w2w3,alt=True,w1lum_lim=0.0,w2lum_lim=0.0,
                 w1lum=np.array([]),w2lum=np.array([])):

    assert w1w2.shape == w2w3.shape, \
        "array size mismatch in my_wedge_cut: %g %g"%(w1w2.shape, w2w3.shape)
    if alt:
        ##mask = ((w1w2>=0.5)&(w1w2>=1.8*w2w3-7.8))
        ##mask = ((w1w2>=0.5)&(w1w2>=2.0*w2w3-9.1))
        #mask = ((w1w2>=0.5)&(w1w2>=2.0*w2w3-8.9))
        if w1lum_lim>0.0 or w2lum_lim>0.0:
            mask = ((w1w2>=0.5)&(w1w2>=2.0*w2w3-8.9)&(w2w3>=2.2)&(w1lum>=w1lum_lim)&(w2lum>=w2lum_lim))
        else:
            mask = ((w1w2>=0.5)&(w1w2>=2.0*w2w3-8.9)&(w2w3>=2.2))
    else:
        #mask = ((w2w3<=4.8)&(w1w2>=0.5))
        if w1lum_lim>0.0 or w2lum_lim>0.0:
            mask = ((w2w3<=4.9)&(w1w2>=0.5)&(w1lum>=w1lum_lim)&(w2lum>=w2lum_lim))
        else:
            mask = ((w2w3<=4.9)&(w1w2>=0.5))

    return mask


def check_axis_range(axrange,data,fatal=False):

    assert len(axrange)==2 and len(data)>0, "invalid input to check_axis_range."
    if axrange[0]>data.min() or axrange[1]<data.max():
        print "WARNING: Data outside axis range will be clipped."
        print "range: ",axrange
        print "data max/min: ",data.max(),data.min()
        if fatal: raise ValueError()


def w1w2_scatter_plot(ax,w1w2,val,altvals=(),oplot_val=(),oplot_w1w2=(),ncam=7,
                      xlim=(),ylim=(-0.1,2.0),xlabel='',ylabel='W1-W2',mec='None',
                      mark='o',c0='k',ms0=2, altc=('r',),altms=(3,),
                      plot_val_last=False,return_leg_handle=False):

    if ncam>0:
        assert w1w2.shape[0] == ncam, \
            "Expected %d camera angles for w1w2, found %d."%(ncam,w1w2.shape[0])
        
        if len(val.shape)==2:
            assert val.shape[0] == ncam, \
                "Expected %d camera angles for val, found %d."%(ncam,val.shape[0])
        else: 
            val = np.tile(val,(ncam,1))
        assert val.shape[1] == w1w2.shape[1], \
            "Data size mismatch for w1w2 & val: %d %d"%(w1w2.shape[1], val.shape[1])
    else:
        assert (w1w2.size == val.size) & len(w1w2.shape)==1, \
            "For ncam=0, expected 1D arrays for w1w2 and val."

    if len(xlim)!=2: xlim = (0.95*val.min(),1.05*val.max())
    if len(ylim)!=2: ylim = (0.95*w1w2.min(),1.05*w1w2.max())
    check_axis_range(xlim,val)
    check_axis_range(ylim,w1w2)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    if len(oplot_val)>0: 
        for v in oplot_val: ax.plot([v,v],ylim,'k',ls='--')

    #kwargs = dict(marker=mark,ls='None',color=c0,ms=ms0,mec='None')
    kwargs = dict(marker=mark,ls='None',color=c0,ms=ms0,mec=mec)
    if not plot_val_last:
        #for icam in range(ncam): ax.plot(val[icam,:],w1w2[icam,:],**kwargs)
        #ax.plot(val.flatten(),w1w2.flatten(),**kwargs)        
        if return_leg_handle:
            p0,=ax.plot(val.flatten(),w1w2.flatten(),**kwargs)
        else:
            ax.plot(val,w1w2,**kwargs)


    if len(altvals)>0:
        assert isinstance(altvals,tuple), "altvals must be a tuple."
        if len(altms)!=len(altvals): altms = altms+(altms[-1],)*(len(altvals)-len(altms))
        if len(altc)!=len(altvals): altc = altc+(altc[-1],)*(len(altvals)-len(altc))
        for j,altval in enumerate(altvals):
            if ncam>0:
                if len(altval.shape)==2:
                    assert altval.shape[0] == ncam, \
                        "Expected %d camera angles for altval, found %d."%(ncam,altval.shape[0])
                else:
                    altval = np.tile(altval,(ncam,1))    
                assert altval.shape[1] == w1w2.shape[1], \
                    "Data size mismatch for w1w2 & altval: %d %d"%(w1w2.shape[1], val.shape[1])
            else:
                assert (w1w2.size == altval.size) & len(altval.shape)==1, \
                    "For ncam=0, expected 1D arrays for w1w2 and altval."

            check_axis_range(xlim,altval)

            kwargs.update(color=altc[j],ms=altms[j],mec=altc[j])
            #for icam in range(ncam): ax.plot(altval[icam,:],w1w2[icam,:],**kwargs)
            ax.plot(altval.flatten(),w1w2.flatten(),**kwargs)

    if plot_val_last:
        kwargs.update(color=c0,ms=ms0,mec=c0)
        #for icam in range(ncam): ax.plot(val[icam,:],w1w2[icam,:],**kwargs)
        p0,=ax.plot(val.flatten(),w1w2.flatten(),**kwargs)

    if return_leg_handle:
        return p0

def w1w2_scatter_multiplot(path='/oasis/projects/nsf/hvd115/lblecha/',
                           multidir='combined_sim_plots_hires',multifbase='hires',
                           skip_snap0=True, cam=range(7), dust_to_metal_ratio=0.4, 
                           latesep=10, squish=False, grid_nh_res=0,old_nh=False):

    mpl.rcParams.update({'font.size': 13})
    carr = ['m','darkblue','b','darkcyan','cyan','darkgreen','green']
    ncam=len(cam)
    path='%s/%s/'%(path,multidir)

    ### load txt files:
    tmpdata = load_infofile(path,fbase='bh_info',fext=multifbase,skip_snap0=skip_snap0)
    snap,time,bh1x,bh1y,bh1z,bh1_mdot,bh1_lbol,bh1_lbol_sunsed,bh1_mass,bh2x,bh2y,bh2z,bh2_mdot,bh2_lbol,bh2_lbol_sunsed,bh2_mass = tmpdata
    nsnaps=len(snap)
    ix_2bh = np.where(bh2_lbol==bh2_lbol)[0]
    nsnaps_2bh=len(snap[ix_2bh])


    lagn = bh1_lbol
    lagn[ix_2bh] = bh1_lbol[ix_2bh] + bh2_lbol[ix_2bh]
    lagn_sun = bh1_lbol_sunsed
    lagn_sun[ix_2bh] = bh1_lbol_sunsed[ix_2bh] + bh2_lbol_sunsed[ix_2bh]
    lagn_sun_tiled = np.tile(lagn_sun,(ncam,1))
    bh_mass = bh1_mass
    bh_mass[ix_2bh] = bh1_mass[ix_2bh] + bh2_mass[ix_2bh]
    ledd = calc_lbol_edd(bh_mass)
    fedd = np.zeros(lagn.size)
    fedd[ledd>0] = lagn[ledd>0]/ledd[ledd>0]
    
    sepsnap,septime,bh_3d_sep,bh_proj_sep = load_infofile(path,fbase='bh_sep',
                                                          fext=multifbase,skip_snap0=skip_snap0)
    sep_ix_2bh = np.where(bh_3d_sep==bh_3d_sep)[0]
    assert nsnaps_2bh==0 or len(sep_ix_2bh) == nsnaps_2bh, \
        "bh array size mismatch: ncam=%d, nsnaps_2bh=%d (=%d in bh_sep)"%(ncam,nsnaps_2bh,len(sep_ix_2bh))

    tmpdata = load_infofile(path,fbase='star_gas_info',fext=multifbase,skip_snap0=skip_snap0)
    sgsnap,sgtime,mstar,mgas,mmetals,sfr = tmpdata
    fgas = 1.0*mgas / (mgas + mstar)
    fdust = 1.0*mmetals*dust_to_metal_ratio / mgas
    ssfr = 1.0*sfr/mstar
    assert len(snap) == len(sgsnap), \
        "array size mismatch in bh_info.txt, star_gas_info.txt: %d, %d"%(snap.size,sgsnap.size)

    tmpdata = load_infofile(path,fbase='lum_info',fext=multifbase,skip_snap0=skip_snap0)
    lsnap,ltime,lbol_grid,lbol_absorbed,ltot_out,lir = tmpdata
    assert lsnap.size == snap.size, \
        "array size mismatch in bh_info.txt, lum_info.txt: %d, %d"%(snap.size,lsnap.size)


    time_wise, snapnum_wise, w1w2, w2w3 = load_wise(maindir=path,subdir='',fext=multifbase,ncam=ncam,skip_snap0=skip_snap0)
    wmask = np.in1d(snapnum_wise,snap)
    print "retaining %d of %d elements of snapnum_wise also in snap."%(snapnum_wise[wmask].size,
                                                                       snapnum_wise.size)
    nowise = np.in1d(snap,snapnum_wise,invert=True)
    print "%d of %d elements of snap are not in snapnum_wise."%(snap[nowise].size,snap.size)    
    print snap[nowise]

    snapnum_wise = snapnum_wise[wmask]
    time_wise = time_wise[wmask]
    w1w2 = w1w2[:,wmask]
    w2w3 = w2w3[:,wmask]

    meanw12 = np.mean(w1w2,axis=0)
    meanw23 = np.mean(w2w3,axis=0)

    initarr=np.zeros((ncam,nsnaps))-1
    bh1_xpix = copy(initarr)
    bh1_ypix = copy(initarr)
    bh1_Ngas_tot_aux = copy(initarr)
    bh1_Ngas_tot_grid = copy(initarr)
    bh1_Ngas_grid = copy(initarr)
    bh2_xpix = copy(initarr)
    bh2_ypix = copy(initarr)
    bh2_Ngas_tot_aux = copy(initarr)
    bh2_Ngas_tot_grid = copy(initarr)
    bh2_Ngas_grid = copy(initarr)

    for c in cam:
        rstr='_grid_res_%.3d'%(1000*grid_nh_res) if grid_nh_res>0 else ''
        with open(path+'/Ngas_bh_cam%d%s_%s.txt'%(c,rstr,multifbase)) as f:
            tmpdata = np.loadtxt(f,unpack=True,dtype='float64')
            if skip_snap0:
                sn = tmpdata[0,:]
                snapmask = (sn>0)
                tmpdata = tmpdata[:,snapmask]
            if old_nh:
                sn,bh1_xpix[c,:],bh1_ypix[c,:],bh1_Ngas_tot_aux[c,:],bh2_xpix[c,:],bh2_ypix[c,:],bh2_Ngas_tot_aux[c,:]=tmpdata
            else:
                sn,bh1_xpix[c,:],bh1_ypix[c,:],bh1_Ngas_tot_aux[c,:],bh1_Ngas_tot_grid[c,:],bh1_Ngas_grid[c,:],bh2_xpix[c,:],bh2_ypix[c,:],bh2_Ngas_tot_aux[c,:],bh2_Ngas_tot_grid[c,:],bh2_Ngas_grid[c,:]=tmpdata

    bh1_Ngas_aux = bh1_Ngas_tot_aux/2.0
    ## don't include BHs that are outside the FOV:
    bh1_Ngas_aux[(bh1_xpix<0)|(bh1_ypix<0)]=np.nan
    med_bh1_Ngas_aux = np.nanmedian(bh1_Ngas_aux, axis=0)
    min_bh1_Ngas_aux = np.nanmin(bh1_Ngas_aux, axis=0)
    max_bh1_Ngas_aux = np.nanmax(bh1_Ngas_aux, axis=0)
    yerr_bh1_aux=[-min_bh1_Ngas_aux+med_bh1_Ngas_aux,max_bh1_Ngas_aux-med_bh1_Ngas_aux]

    med_bh1_Ngas_grid = np.nanmedian(bh1_Ngas_grid, axis=0)
    min_bh1_Ngas_grid = np.nanmin(bh1_Ngas_grid, axis=0)
    max_bh1_Ngas_grid = np.nanmax(bh1_Ngas_grid, axis=0)
    yerr_bh1_grid=[-min_bh1_Ngas_grid+med_bh1_Ngas_grid,max_bh1_Ngas_grid-med_bh1_Ngas_grid]

    if nsnaps_2bh > 0:
        bh2_Ngas_aux = bh2_Ngas_tot_aux/2.0
        bh2_Ngas_aux[(bh2_xpix<0)|(bh2_ypix<0)]=np.nan
        med_bh2_Ngas_aux = np.nanmedian(bh2_Ngas_aux, axis=0)
        min_bh2_Ngas_aux = np.nanmin(bh2_Ngas_aux, axis=0)
        max_bh2_Ngas_aux = np.nanmax(bh2_Ngas_aux, axis=0)
        yerr_bh2_aux=[-min_bh2_Ngas_aux+med_bh2_Ngas_aux,max_bh2_Ngas_aux-med_bh2_Ngas_aux]

        med_bh2_Ngas_grid = np.nanmedian(bh2_Ngas_grid, axis=0)
        min_bh2_Ngas_grid = np.nanmin(bh2_Ngas_grid, axis=0)
        max_bh2_Ngas_grid = np.nanmax(bh2_Ngas_grid, axis=0)
        yerr_bh2_grid=[-min_bh2_Ngas_grid+med_bh2_Ngas_grid,max_bh2_Ngas_grid-med_bh2_Ngas_grid]

    flbol_grid = lagn_sun/lbol_grid
    flbol_grid_tiled = np.tile(flbol_grid,(ncam,1))
    flbol_out = lagn_sun_tiled/ltot_out
    print flbol_grid.shape,flbol_out.shape
    print "\nMin fLbol with W1W2>0.8: %g"%flbol_out[w1w2>0.8].min()
    print "Min fLbol with W1W2>0.5: %g"%flbol_out[w1w2>0.5].min()
    print "Min fLbol_grid with mean W1W2>0.8: %g"%flbol_grid[meanw12>0.8].min()
    print "Min fLbol_grid with mean W1W2>0.5: %g"%flbol_grid[meanw12>0.5].min()
    print "Max flbol_out each cam:",flbol_out.max(axis=1)
    mean_flbol_out = flbol_out.mean(axis=0)
    print "\n(Note: for small number of cams, ok if mean flbol_out is slightly > 1.)"
    print "Max flbol_out averaged over cams:",mean_flbol_out.max()
    print "snaps with mean flbol_out > 1:",snap[mean_flbol_out>1]
    print "Max flbol_grid: %g\n"%flbol_grid.max()
    print "Lum range for flbol>0.3: ",lagn[flbol_grid>0.3].min(),lagn[flbol_grid>0.3].max()
    print "Lum range for 0.1<flbol<=0.3: ",lagn[(flbol_grid>0.1)&(flbol_grid<=0.3)].min(),lagn[(flbol_grid>0.1)&(flbol_grid<0.3)].max()
    print "Lum range for flbol<=0.1: ",lagn[flbol_grid<=0.1].min(),lagn[flbol_grid<=0.1].max()

    bh_proj_sep_negvals = copy(bh_proj_sep)
    bh_proj_sep_negvals[bh_proj_sep!=bh_proj_sep] = -10
    print bh_proj_sep_negvals.shape
    bh_proj_sep_mean = bh_proj_sep_negvals.mean(axis=0)
    print bh_proj_sep_mean.shape

    plt.ioff()
    fig = plt.figure(figsize=(5,5))
    plt.xlim(-12,120)
    plt.ylim(21.0,24.5)
    plt.plot(bh_proj_sep_negvals,np.log10(bh1_Ngas_grid),'bo',ms=2,mec='b')
    plt.plot(bh_proj_sep_negvals,np.log10(bh2_Ngas_grid),'mo',ms=2,mec='m')

    fig.savefig(path+'/bhsep_vs_Ngas_quickplot_%s.pdf'%multifbase,format='pdf')
    #plt.close('all')
    plt.cla()
    plt.clf()
    plt.close()
    

    fs=(8,6) if squish else (8,8)
    fig = plt.figure(figsize=(8,8))

    w1w2_scatter_plot(fig.add_subplot(321),w1w2,np.log10(lagn), xlim=(41.5,47),
                      xlabel=r'L$_{AGN}$ [log erg s$^{-1}$]',ms0=1.5,altms=(2.5,))
    #w1w2_scatter_plot(fig.add_subplot(322),w1w2,np.log10(flbol_out), 
    #                  altvals=(np.log10(flbol_grid),),
    #                  xlim=(-5,0.5), xlabel=r'log (L$_{\rm AGN}$/L$_{\rm tot}$)',ms0=1.5,altms=(2.5,))
    w1w2_scatter_plot(fig.add_subplot(322),w1w2,np.log10(flbol_grid_tiled), 
                      xlim=(-5,0.5), xlabel=r'log (L$_{\rm AGN}$/L$_{\rm tot}$)',ms0=1.5,altms=(2.5,))
    w1w2_scatter_plot(fig.add_subplot(323),w1w2,np.log10(lir), xlim=(41.5,47),
                      xlabel=r'L$_{IR}$ [log erg s$^{-1}$]',ms0=1.5,altms=(2.5,))
    w1w2_scatter_plot(fig.add_subplot(324),w1w2,np.log10(lir/ltot_out), xlim=(-5,0.5),
                      xlabel=r'log (L$_{IR}$/L$_{tot}$)',ms0=1.5,altms=(2.5,))
    w1w2_scatter_plot(fig.add_subplot(325),w1w2,np.log10(ltot_out), altvals=(np.log10(lbol_grid),),
                      xlim=(41.5,47), xlabel=r'L$_{tot}$ [log erg s$^{-1}$]',ms0=1.5,altms=(2.5,))
    #w1w2_scatter_plot(fig.add_subplot(326),w1w2,np.log10(lagn/lbol_grid), 
    #                  xlim=(-3.9,0.1), ylim=(-0.1,1.8),
    #                  xlabel=r'log (L$_{\rm AGN}$/L$_{\rm tot}$)',ms0=1.8)
    w1w2_scatter_plot(fig.add_subplot(326),w1w2,np.log10(fedd),
                      xlim=(-5,0.05), xlabel=r'f$_{\rm Edd}$',ms0=1.5)


    fig.suptitle(multifbase)
    fig.subplots_adjust(bottom=0.1,left=0.12,right=0.94,hspace=0.35,wspace=0.3,top=0.92)
    fig.savefig(path+'/w1w2_lum_scatter_plots_%s.pdf'%multifbase,format='pdf')
    plt.clf()
    plt.cla()
    plt.close()


    print "\nmaking w1w2_flbol plot\n"
    #fig = plt.figure(figsize=(5,3.5))
    fig = plt.figure(figsize=(5,3.0))
    ##w1w2_scatter_plot(fig.add_subplot(111),w1w2,np.log10(flbol_grid_tiled), 
    ##                  xlim=(-5,0.01), xlabel=r'log (L$_{\rm AGN}$/L$_{\rm tot}$)',ms0=1.5)
    #w1w2_scatter_plot(fig.add_subplot(111),w1w2[bh_proj_sep_negvals>latesep],
    #                  np.log10(flbol_grid_tiled[bh_proj_sep_negvals>latesep]), 
    #                  xlim=(-5,0.01), xlabel=r'log (L$_{\rm AGN}$/L$_{\rm tot}$)',ms0=3,c0='None',mec='b',ncam=0,mark='^')
    #w1w2_scatter_plot(fig.add_subplot(111),w1w2[(bh_proj_sep_negvals<=latesep)&(bh_proj_sep_negvals>0)],
    #                  np.log10(flbol_grid_tiled[(bh_proj_sep_negvals<=latesep)&(bh_proj_sep_negvals>0)]), 
    #                  xlim=(-5,0.01), xlabel=r'log (L$_{\rm AGN}$/L$_{\rm tot}$)',ms0=3,c0='None',mec='r',ncam=0)
    #w1w2_scatter_plot(fig.add_subplot(111),w1w2[(bh_proj_sep_negvals<=0)],
    #                  np.log10(flbol_grid_tiled[(bh_proj_sep_negvals<=0)]), 
    #                  xlim=(-5,0.01), xlabel=r'log (L$_{\rm AGN}$/L$_{\rm tot}$)',ms0=3,c0='None',mec='g',ncam=0,mark='s')
    ### note: for hires, limits xlim=(-4.15,0.01), ylim=(-0.1,1.7) work well
    ###       for lowres, xlim=(-3.9,0.01), ylim=(-0.1,1.55) works well
    p1=w1w2_scatter_plot(fig.add_subplot(111),meanw12[bh_proj_sep_mean>latesep],
                         np.log10(flbol_grid[bh_proj_sep_mean>latesep]), ylim=(-0.1,1.55),
                         xlim=(-3.9,0.01), xlabel=r'log (L$_{\rm AGN}$/L$_{\rm tot}$)',ms0=5,
                         c0='None',mec='g',ncam=0,mark='^',return_leg_handle=True)
    p2=w1w2_scatter_plot(fig.add_subplot(111),meanw12[(bh_proj_sep_mean<=latesep)&(bh_proj_sep_mean>0)],
                          np.log10(flbol_grid[(bh_proj_sep_mean<=latesep)&(bh_proj_sep_mean>0)]), 
                          xlim=(-3.9,0.01),  ylim=(-0.1,1.55), xlabel=r'log (L$_{\rm AGN}$/L$_{\rm tot}$)',
                         ms0=5,c0='None',mec='m',ncam=0,mark='s',return_leg_handle=True)
    p3=w1w2_scatter_plot(fig.add_subplot(111),meanw12[(bh_proj_sep_mean<=0)],
                         np.log10(flbol_grid[(bh_proj_sep_mean<=0)]),  ylim=(-0.1,1.55), 
                         xlim=(-3.9,0.01), xlabel=r'log (L$_{\rm AGN}$/L$_{\rm tot}$)',ms0=5,
                         c0='None',mec='darkorange',ncam=0,mark='v',return_leg_handle=True)
    
    plt.plot([-3.9,0.01],[0.8,0.8],color='gray',linestyle='-.',lw=2)
    plt.plot([-3.9,0.01],[0.5,0.5],color='gray',linestyle=':',lw=2)
    lh=(p1,p2,p3)
    print len(lh)
    print lh
    print p1, p2, p3
    ll=('Early','Late','Post')
    plt.legend(lh, ll, fontsize=11, loc='upper left',numpoints=1, handletextpad=0.1)
    #plt.legend(lh, ll, fontsize=9, loc='upper left', numpoints=1, handletextpad=0.08)
    #fig.suptitle(multifbase)
    #fig.subplots_adjust(bottom=0.15,left=0.12,right=0.94,hspace=0.35,wspace=0.3,top=0.92)
    fig.subplots_adjust(bottom=0.18,left=0.12,right=0.94,top=0.95)
    fig.savefig(path+'/w1w2_flbol_scatter_plot_%s.pdf'%multifbase,format='pdf')
    plt.clf()
    plt.cla()
    plt.close()
    #return


    print "\nmaking w1w2_lagn plot\n"
    lagn_xlim=(41.5,46.5)
    fig = plt.figure(figsize=(5,3.0))
    p1=w1w2_scatter_plot(fig.add_subplot(111),meanw12[bh_proj_sep_mean>latesep],
                         np.log10(lagn[bh_proj_sep_mean>latesep]), ylim=(-0.1,1.55),
                         xlim=lagn_xlim, xlabel=r'log L$_{\rm AGN}$ [erg s$^{-1}$]',ms0=5,
                         c0='None',mec='g',ncam=0,mark='^',return_leg_handle=True)
    p2=w1w2_scatter_plot(fig.add_subplot(111),meanw12[(bh_proj_sep_mean<=latesep)&(bh_proj_sep_mean>0)],
                          np.log10(lagn[(bh_proj_sep_mean<=latesep)&(bh_proj_sep_mean>0)]), 
                          xlim=lagn_xlim,  ylim=(-0.1,1.55), xlabel=r'log L$_{\rm AGN}$ [erg s$^{-1}$]',
                         ms0=5,c0='None',mec='m',ncam=0,mark='s',return_leg_handle=True)
    p3=w1w2_scatter_plot(fig.add_subplot(111),meanw12[(bh_proj_sep_mean<=0)],
                         np.log10(lagn[(bh_proj_sep_mean<=0)]),  ylim=(-0.1,1.55), 
                         xlim=lagn_xlim, xlabel=r'log L$_{\rm AGN}$ [erg s$^{-1}$]',ms0=5,
                         c0='None',mec='darkorange',ncam=0,mark='v',return_leg_handle=True)
    
    plt.plot([lagn_xlim[0],lagn_xlim[1]],[0.8,0.8],color='gray',linestyle='-.',lw=2)
    plt.plot([lagn_xlim[0],lagn_xlim[1]],[0.5,0.5],color='gray',linestyle=':',lw=2)
    lh=(p1,p2,p3)
    print len(lh)
    print lh
    print p1, p2, p3
    ll=('Early','Late','Post')
    plt.legend(lh, ll, fontsize=11, loc='upper left',numpoints=1, handletextpad=0.1)
    #plt.legend(lh, ll, fontsize=9, loc='upper left', numpoints=1, handletextpad=0.08)
    #fig.suptitle(multifbase)
    #fig.subplots_adjust(bottom=0.15,left=0.12,right=0.94,top=0.92)
    fig.subplots_adjust(bottom=0.18,left=0.12,right=0.94,top=0.95)
    fig.savefig(path+'/w1w2_lagn_scatter_plot_%s.pdf'%multifbase,format='pdf')
    plt.clf()
    plt.cla()
    plt.close()
    return


    fig = plt.figure(figsize=(8,8))

    w1w2_scatter_plot(fig.add_subplot(321),w1w2,np.log10(sfr), xlim=(-1,2.8),
                      xlabel=r'SFR [log M$_{\odot}$ yr$^{-1}$]',ms0=1.5,altms=(2.5,))
    w1w2_scatter_plot(fig.add_subplot(322),w1w2,np.log10(ssfr), xlim=(-11.9,-8),
                      xlabel=r'sSFR [log yr$^{-1}$]',ms0=1.5,altms=(2.5,))
    w1w2_scatter_plot(fig.add_subplot(323),w1w2,fgas, xlim=(0,0.3),
                      xlabel=r'Mgas/(Mgas+M*)',ms0=1.5,altms=(2.5,))
    #w1w2_scatter_plot(fig.add_subplot(324),w1w2,fdust, xlim=(0.004,0.015),
    #                  xlabel=r'Mdust/Mgas',ms0=1.5,altms=(2.5,))
    w1w2_scatter_plot(fig.add_subplot(324),w1w2,np.log10(mmetals*dust_to_metal_ratio), 
                      xlim=(7,9), xlabel=r'Mdust [log M$_{\odot}$]',ms0=1.5,altms=(2.5,))
    if nsnaps_2bh <= 0:
        w1w2_scatter_plot(fig.add_subplot(325),w1w2,np.log10(med_bh1_Ngas_grid),
                          xlim=(21.1,24.5), xlabel=r'$N_{\rm H}$ [log cm$^{-2}$]',ms0=1.5)
    else:
        w1w2_scatter_plot(fig.add_subplot(325),w1w2,np.log10(med_bh1_Ngas_grid),
                          altvals=(np.log10(med_bh2_Ngas_grid),), xlim=(21,24.5), 
                          xlabel=r'$N_{\rm H} [log cm^{-2}$]',ms0=1.5,altms=(2.5,))
        bh_proj_sep[:,ix_2bh]
        #w1w2_scatter_plot(fig.add_subplot(326),w1w2[:,ix_2bh],np.log10(min_psep),
        w1w2_scatter_plot(fig.add_subplot(326),w1w2[:,ix_2bh],np.log10(bh_proj_sep[:,ix_2bh]),
                          altvals=(np.log10(bh_3d_sep[ix_2bh]),), xlim=(-2,2.2),
                          xlabel=r'BH sep [log kpc]',plot_val_last=True,ms0=1.5,altms=(2.5,))

    fig.suptitle(multifbase)
    fig.subplots_adjust(bottom=0.1,left=0.12,right=0.94,hspace=0.35,wspace=0.3,top=0.9)
    fig.savefig(path+'/w1w2_other_scatter_plots_%s.pdf'%multifbase,format='pdf')
    plt.clf()
    plt.cla()
    plt.close()


    print "Making SFR vs L_AGN plot..."
    fig = plt.figure(figsize=(6,6))
    ax=fig.add_subplot(221)

    #xlim=(1.5,6.1)
    #ylim=(0,2.0)
    #plt.xlim(xlim)
    #plt.ylim(ylim)
    ax.set_xlim(41.5,46.5)
    ax.set_ylim(-1,2.8)
    ax.set_xlabel(r'log($L_{\rm AGN}$)')
    ax.set_ylabel(r'log(SFR)')

    ax.plot(np.log10(lagn),np.log10(sfr),'ko',ms=2)
    ax.plot(np.log10(lagn[meanw12>0.5]),np.log10(sfr[meanw12>0.5]),'g^',ms=4,mec='g')
    ax.plot(np.log10(lagn[meanw12>0.8]),np.log10(sfr[meanw12>0.8]),'rs',ms=4,mec='r')
    
    ax=fig.add_subplot(222)

    ax.set_xlim(-3.5,0.1)
    ax.set_xticks(np.arange(-3.5,0.1,dtype='int'))
    ax.set_ylim(-1,2.8)
    ax.set_xlabel(r'log($L_{\rm AGN}/L_{\rm tot}$)')
    #ax.set_ylabel(r'log(SFR)')

    ax.plot(np.log10(flbol_grid),np.log10(sfr),'ko',ms=2)
    ax.plot(np.log10(flbol_grid[meanw12>0.5]),
            np.log10(sfr[meanw12>0.5]),'g^',ms=4,mec='g')
    ax.plot(np.log10(flbol_grid[meanw12>0.8]),
            np.log10(sfr[meanw12>0.8]),'rs',ms=4,mec='r')   

    j11mask = j11_wedge_cut(meanw12,meanw23,
                            w1lum_lim=w1lum_lim,w2lum_lim=w2lum_lim,w1lum=w1lum,w2lum=w2lum)
    

    ax=fig.add_subplot(223)

    #ax.set_xlim(40.5,46.5)
    ax.set_xlim(41.5,46.5)
    ax.set_ylim(-1,2.8)
    ax.set_xlabel(r'log($L_{\rm AGN}$)')
    ax.set_ylabel(r'log(SFR)')

    ax.plot(np.log10(lagn),np.log10(sfr),'ko',ms=2)
    ax.plot(np.log10(lagn[j11mask]),np.log10(sfr[j11mask]),'bD',ms=4,mec='b')
    
    ax=fig.add_subplot(224)

    ax.set_xlim(-3.5,0.1)
    ax.set_xticks(np.arange(-3.5,0.1,dtype='int'))
    ax.set_ylim(-1,2.8)
    ax.set_xlabel(r'log($L_{\rm AGN}/L_{\rm tot}$)')
    #ax.set_ylabel(r'log(SFR)')

    ax.plot(np.log10(flbol_grid),np.log10(sfr),'ko',ms=2)
    ax.plot(np.log10(flbol_grid[j11mask]),
            np.log10(sfr[j11mask]),'bD',ms=4,mec='b')

    fig.suptitle(multifbase)
    fig.subplots_adjust(bottom=0.1,left=0.12,right=0.96,wspace=0.25,hspace=0.25,top=0.94)
    fig.savefig(path+'/sfr_vs_lagn_%s.pdf'%multifbase,format='pdf')
    plt.clf()
    plt.cla()
    plt.close()





if __name__ == "__main__":

    helpstr = 'column_density.py -f <func> -g <gadgetdir> -s <subdir> -S <startsnap> -t <tmax_postmrg> -e <extra> -E <errtype> -r <grid_nh_res> -a (append) -o (overwrite_all) -b (write_bhfile_only) -N (write_ngasfiles_only) -d (plot_3dsep) -n (agnx0) -l (use_logbins) -L (ylog)'
    

    func='plot'
    gadgetdir='p3new_merger_sims/q1_fg0.3'
    subdir='q1_fg0.3_sunruns/all'
    startsnap=-1
    tmax_postmrg=0.25
    append, overwrite_all = False, False
    write_bhfile_only = False
    write_ngasfiles_only = False
    plot_3dsep=False
    agnx0=False
    use_logbins=False
    ylog=False
    extra=''
    errtype='mad'
    grid_nh_res=0

    try:
        opts,args = getopt.getopt(sys.argv[1:],'f:g:s:x:S:t:e:E:r:aobNmdnlLh')
    except getopt.GetoptError:
        print helpstr
        sys.exit()
    print "len(opts)=",len(opts),"len(args)=",len(args)
    print "opts=",opts
    print "args=",args
    if len(opts)==0:
        sys.exit()
    for opt,arg in opts:
        print "opt,arg=",opt,arg
        if opt == '-h':
            print helpstr
            sys.exit()
        elif opt == '-f':
            func = arg
            if func not in ('plot','merger_phases','process_snaps','load_wise'):
                print "Error: invalid function name: %s"%func
                sys.exit()
        elif opt == '-g':
            gadgetdir = arg
        elif opt == '-s':
            subdir = arg
        elif opt == '-S':
            startsnap = np.int(arg)
        elif opt == '-t':
            tmax_postmrg = np.float(arg)
        elif opt == '-e':
            extra = arg
        elif opt == '-E':
            errtype = arg
        elif opt == '-r':
            grid_nh_res = np.float(arg)
        elif opt == '-a':
            append = True
        elif opt == '-o':
            overwrite_all = True
        elif opt == '-b':
            write_bhfile_only=True
        elif opt == '-N':
            write_ngasfiles_only=True
        elif opt == '-d':
            plot_3dsep = True
        elif opt == '-n':
            agnx0 = True
        elif opt == '-l':
            use_logbins = True
        elif opt == '-L':
            ylog = True
        else:
            print "invalid function name: ",func
            sys.exit()

    if func=='plot':
        print "Making plots for subdir %s..."%subdir
        plot(subdir=subdir)
    elif func=='merger_phases':
        print "Making merger phase plots for subdir %s..."%subdir
        print "plot_3dsep: ",plot_3dsep
        print "agnx0: ",agnx0
        print "use_logbins: ",use_logbins
        print "ylog: ",ylog
        print "errtype:", errtype
        merger_phases(subdir_arr=subdir, plot_3dsep=plot_3dsep, errtype=errtype,
                      agnx0=agnx0, use_logbins=use_logbins, ylog=ylog, extra=extra,
                      grid_nh_res=grid_nh_res)
    elif func=='process_snaps':
        print "processing snaps for gadgetdir %s, subdir %s..."%(gadgetdir, subdir)
        print "startsnap=%d"%startsnap
        print "append: ",append
        print "overwrite_all: ",overwrite_all
        print "write_bhfile_only: ",write_bhfile_only
        process_snaps(gadgetdir=gadgetdir,subdir=subdir,startsnap=startsnap,
                      append=append,overwrite_all=overwrite_all,grid_nh_res=grid_nh_res,
                      write_bhfile_only=write_bhfile_only,write_ngasfiles_only=write_ngasfiles_only)
    elif func=='load_wise':
        tmp=load_wise(subdir=subdir)
