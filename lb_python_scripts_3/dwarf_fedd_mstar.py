import numpy as np
import matplotlib as mpl
#mpl.use('agg')
import matplotlib.pyplot as plt
from copy import copy
import sys
import PyCos as pc
import astro_constants as ac
import lb_utils as lb

omega_m = 0.2726
omega_l = 0.7274
h = 0.704
cosmo = pc.Cosmology(omega_m, omega_l, 0.0, -1.0, h)


def plot(path="/n/home00/lblecha/dwarf_agn_data/"):


    with open("%s/dwarf_sample_all.dat"%path,"r") as fp:
        dtypes = np.dtype({'names':(['f{}'.format(i) for i in range(27)]),
                           'formats':['|S30','|S30','|S30',np.float,np.float,
                                      np.float,np.float,np.float,np.float,
                                      np.float,'|S30',np.float,np.float,
                                      np.float,np.float,np.float,
                                      np.float,np.float,np.float,np.float,
                                      np.float,np.float,
                                      np.float,np.float,np.float,
                                      '|S10','|S10']})
        tmpdata = np.loadtxt(fp,unpack=True,dtype=dtypes,skiprows=1)
        catID, name, altname, ra, dec, Lx, Lbol, Bmag, AGNtype, z, optspec_str, bhmass_low, bhmass_hi, mstarK11, mstarother_low, mstarother_hi, vdisp_low, vdisp_low_err, vdisp_hi, vdisp_hi_err, pctAGNrband, dist, gal_gmag, gal_rmag, gminusr, hasHST, hasChandra = tmpdata
        print "Loaded data for %d sources."%len(name)
        ix=np.where(Lx>0.0)[0]
        if len(ix) < len(Lx):
            print "cutting %d sources with negative Lx. names:"%(len(Lx)-len(ix)),name[Lx<=0]
            tmpdata = lb.apply_cuts(ix, tmpdata)
            catID, name, altname, ra, dec, Lx, Lbol, Bmag, AGNtype, z, optspec_str, bhmass_low, bhmass_hi, mstarK11, mstarother_low, mstarother_hi, vdisp_low, vdisp_low_err, vdisp_hi, vdisp_hi_err, pctAGNrband, dist, gal_gmag, gal_rmag, gminusr, hasHST, hasChandra = tmpdata
            print "Retained data for %d sources."%len(name)


    plt.clf()
    plt.cla()
    plt.close('all')
    fig = plt.figure(figsize=(6,6))
    
    ax1 = fig.add_subplot(2,2,1)
    ## stellar mass histogram
    plt.xlim(8,10.1)
    plt.ylim(0,25)
    plt.xlabel(r'log (M$_*$/M$_{\odot}$)')
    plt.ylabel('N')
    ax1.hist(mstarK11[(mstarK11>0)&(mstarK11<1.0e10)],
             histtype='step',color='k',bins=7,range=(8.25,10))

    ax2 = fig.add_subplot(2,2,2)
    ## redshift histogram
    plt.xlim(0,0.061)
    plt.ylim(0,18)
    plt.xlabel('z')
    plt.ylabel('N')
    ax2.hist(z[(mstarK11<1.0e10)],histtype='step',
             color='k',bins=6,range=(0,0.06))
    ax2.set_xticks(ax2.get_xticks()[::2])
    ax2.set_yticks(ax2.get_yticks()[::2])

    ax3 = fig.add_subplot(2,2,3)
    ## fedd histogram
    
    ax4 = fig.add_subplot(2,2,4)
    ## mstar vs Lx
    plt.xlim(8,10.1)
    plt.ylim(40.2,44.7)
    plt.xlabel(r'log (M$_*$/M$_{\odot}$)')
    plt.ylabel(r'log L$_{\rm X}$ [erg s$^{-1}$]')
    ax4.plot(mstarK11[mstarK11>0],Lx[mstarK11>0],'ko',markersize=3)
    ix = np.where(mstarK11<=0)[0]
    if len(ix)>0:
        arr=np.zeros((len(ix)))
        for xtup,ytup in zip(zip(arr+8.0,arr+8.25),zip(Lx[ix],Lx[ix])):
            ax4.plot(xtup,ytup,'r-',linewidth=0.5)

    fig.subplots_adjust(wspace=0.3,hspace=0.3)
    fig.savefig('%s/dwarf_agn_properties.eps'%path)
    plt.clf()
    plt.cla()
    plt.close('all')
    
