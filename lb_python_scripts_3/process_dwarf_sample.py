import numpy as np
import matplotlib as mpl
#mpl.use('agg')
import matplotlib.pyplot as plt
from copy import copy
import sys, os
import PyCos as pc
from itertools import islice
from astropy import units as u
from astropy.coordinates import SkyCoord
import astro_constants as ac
import lb_utils as lb
import match_catalog as matchcat

omega_m = 0.2726
omega_l = 0.7274
h = 0.704
cosmo = pc.Cosmology(omega_m, omega_l, 0.0, -1.0, h)

def read_data(path="/n/home00/lblecha/dwarf_agn_data/",
              fname='dwarf_sample_with_b13_extras.dat',
              with_lx_only=True,has_b13_extras=True,
              has_lx_a12=False):

    with open("%s/%s"%(path,fname),"r") as fp:

        ncols = 27 
        if has_b13_extras: ncols=ncols+3
        if has_lx_a12: ncols=ncols+2
        fmt_list = ['|S30','|S30','|S30',np.float64,np.float64,
                    np.float,np.float,np.float,np.float,np.float,
                    '|S30',np.float,np.float,np.float,np.float,np.float,
                    np.float,np.float,np.float,np.float,np.float,
                    np.float,np.float,np.float,np.float,'|S10','|S10']
        if has_b13_extras: fmt_list = fmt_list+[np.float,np.float,np.float]
        if has_lx_a12: fmt_list = fmt_list+[np.float,np.float]
        dtypes = np.dtype({'names':(['f{}'.format(i) for i in range(ncols)]),
                           'formats':fmt_list})
        tmpdata = np.loadtxt(fp,unpack=True,dtype=dtypes,skiprows=1)
        name=tmpdata[1]
        Lx=tmpdata[5]
        #catID, name, altname, ra, dec, Lx, Lbol, Bmag, AGNtype, z, optspec_str, bhmass_low, bhmass_hi, mstarK11, mstarother_low, mstarother_hi, vdisp_low, vdisp_low_err, vdisp_hi, vdisp_hi_err, pctAGNrband, dist, gal_gmag, gal_rmag, gminusr, hasHST, hasChandra = tmpdata[:27]
        print "Loaded data for %d dwarf sources."%len(name)
        if with_lx_only:
            ix=np.where(Lx>0.0)[0]
            if len(ix) < len(Lx):
                print "cutting %d sources with negative Lx. names:"%(len(Lx)-len(ix)),name[Lx<=0]
                tmpdata = lb.apply_cuts(ix, tmpdata)

        return tmpdata

def read_data_alt(path="/n/home00/lblecha/dwarf_agn_data/",
                  fname='dwarf_sample_with_b13_extras.dat',
                  with_lx_only=True,has_b13_extras=True,
                  has_lx_a12=False):


    with open("%s/%s"%(path,fname),"r") as lines:
        colnames = np.genfromtxt(islice(lines,0,1),dtype='|S30')
        ncols = colnames.size
        print colnames.shape
        print colnames
    

    with open("%s/%s"%(path,fname),"r") as fp:

        fmtlist=np.concatenate((np.repeat('|S30',3),np.repeat(np.float,7),
                                np.array(['|S30']),np.repeat(np.float,14),
                                np.repeat('|S10',2),np.repeat(np.float,5))).flatten()
        fmtlist = fmtlist[:ncols]
        dtypes = np.dtype({'names':colnames,'formats':fmtlist})
        tmpdata = np.loadtxt(fp,dtype=dtypes,skiprows=1)
        #print tmpdata.size
        print tmpdata.shape
        print len(tmpdata[0])
        print tmpdata['optspec'].shape
        print tmpdata['optspec']
        #print tmpdata
        
        print "Loaded data for %d dwarf sources."%tmpdata.size
        if with_lx_only:
            name=tmpdata['name']
            Lx=tmpdata['Lx(14-195)']
            mask = np.zeros(len(Lx)).astype('bool')
            mask[Lx>0.0] = True
            if Lx[mask].size < Lx.size:
                print "cutting %d sources with negative Lx. names:"%(Lx[~mask].size,
                                                                     name[~mask])
                tmpdata = tmpdata[mask]
                print tmpdata.shape

        return tmpdata


def read_bass_data(path="/n/home00/lblecha/dwarf_agn_data/",
                   fname='dwarf_BASS_data.dat'):

    print "\nLoading BASS data..."
    with open("%s/%s"%(path,fname),"r") as lines:
        colnames = np.genfromtxt(islice(lines,0,1),dtype='|S30')
        ncols = colnames.size
        print "ncols=%d"%ncols
        print colnames
    
    with open("%s/%s"%(path,fname),"r") as fp:

        #BATindex BATname Name ra dec MBH_Hbeta MBH_vdisp MBH_vdisp_err FWHM_Halpha FWHM_Halpha_err Halpha_broad Halpha_broad_err Halpha Halpha_err Hbeta_broad Hbeta_broad_err corr_Hbeta Hbeta_err corr_OIII OIII_err corr_NII NII_err

        fmtlist=np.concatenate((np.array([np.int]),np.repeat('|S30',2),
                                np.repeat(np.float,19))).flatten()
        dtypes = np.dtype({'names':colnames,'formats':fmtlist})
        tmpdata = np.loadtxt(fp,dtype=dtypes,skiprows=1)
        print tmpdata.shape
        print len(tmpdata[0])
        
        print "Loaded BASS data for %d dwarf sources."%tmpdata.size

    return tmpdata
    

def write_ra_dec(path="/n/home00/lblecha/dwarf_agn_data/",
                 mrange=(0,15),exclude_nomass=False,fmt='deg'):

    if fmt not in ('deg','dms','hmsdms'):
    #if fmt not in ('deg'):
        print "Error: ra/dec format %s is not defined.\n"
        sys.exit()

    if mrange[0]<0 or mrange[1]>15:
        print "Error: invalid mass range specified: ", mrange
        sys.exit()

    catID, name, altname, ra, dec, Lx, Lbol, Bmag, AGNtype, z, optspec_str, bhmass_low, bhmass_hi, mstarK11, mstarother_low, mstarother_hi, vdisp_low, vdisp_low_err, vdisp_hi, vdisp_hi_err, pctAGNrband, dist, gal_gmag, gal_rmag, gminusr, hasHST, hasChandra = read_data(path=path)
    print "Read data for %d sources."%len(ra)

    if exclude_nomass:
        ix_m = np.where((mstarK11>=mrange[0])&(mstarK11<mrange[1]))[0]
        pstr="Retaining data for %d sources with 10^%.5g < M*(K11)/Msun <10^%.5g."
        nomass_str='_K11_only'
    else:
        ix_m = np.where((mstarK11==-1)|
                        ((mstarK11>=mrange[0])&(mstarK11<mrange[1])))[0]
        pstr="Retaining data for %d sources with no K11 mass or with 10^%.5g < M*(K11)/Msun <10^%.5g."
        nomass_str=''
    print pstr%(len(ix_m),mrange[0],mrange[1])
    ra = ra[ix_m]
    dec = dec[ix_m]
    
    if mrange[0]>0 or mrange[1]<15:
            mcut_str = '_M%.5g-%.5g'%(mrange[0],mrange[1])
    else: mcut_str = ''

    c = SkyCoord(ra=ra*u.degree, dec=dec*u.degree)
    if fmt=='dms':        
        cstr = c.to_string('dms')
    elif fmt=='hmsdms':
        cstr = c.to_string('hmsdms')
    else: cstr = c.to_string('decimal')
        
    fname = 'dwarf_sample_ra_dec_%s%s%s.dat'%(fmt,mcut_str,nomass_str)

    print "Writing data for %d sources."%len(ra)
    with open("%s/%s"%(path,fname),"w") as fpout:
        for coord in cstr:
            fpout.write('%s\n'%coord)

def print_subsample(path="/n/home00/lblecha/dwarf_agn_data/",
                    mkran=(-1,11),lxran=(38,48),zran=(0,0.1)):

    catID, name, altname, ra, dec, Lx, Lbol, Bmag, AGNtype, z, optspec_str, bhmass_low, bhmass_hi, mstarK11, mstarother_low, mstarother_hi, vdisp_low, vdisp_low_err, vdisp_hi, vdisp_hi_err, pctAGNrband, dist, gal_gmag, gal_rmag, gminusr, hasHST, hasChandra = read_data(path=path)
    print "Read data for %d sources."%len(ra)

    ix = np.where((mstarK11>=mkran[0])&(mstarK11<mkran[1])&
                  (Lx>=lxran[0])&(Lx<lxran[1])&
                  (z>=zran[0])&(z<zran[1]))[0]
    print len(ix)
    print mstarK11[catID=='058']
    print Lx[catID=='058']
    
    dL = np.array([cosmo.lum_dist(red) for red in z]) * 1e6*ac.PC
    Fbol = np.log10( 10**Lbol / ( 4*np.pi * dL*dL ) )

    print "min/max Lx:",min(Lx[ix]),max(Lx[ix])
    print "min/max dL:",min(dL[ix]/(1e6*ac.PC)),max(dL[ix]/(1e6*ac.PC))
    print "min/max Lbol:",min(Lbol[ix]),max(Lbol[ix])
    print "min/max Fbol:",min(Fbol[ix]),max(Fbol[ix])

    catID, name, altname, ra, dec, Lx, Lbol, Bmag, AGNtype, z, optspec_str, bhmass_low, bhmass_hi, mstarK11, mstarother_low, mstarother_hi, vdisp_low, vdisp_low_err, vdisp_hi, vdisp_hi_err, pctAGNrband, dist, gal_gmag, gal_rmag, gminusr, hasHST, hasChandra = lb.apply_cuts(ix, (catID, name, altname, ra, dec, Lx, Lbol, Bmag, AGNtype, z, optspec_str, bhmass_low, bhmass_hi, mstarK11, mstarother_low, mstarother_hi, vdisp_low, vdisp_low_err, vdisp_hi, vdisp_hi_err, pctAGNrband, dist, gal_gmag, gal_rmag, gminusr, hasHST, hasChandra))

    ixsort = ra.argsort(kind='mergesort')

    for i in range(ra.size):
        #print "%s %s %s %g %g %g %g"%(catID[ix[i]],name[ix[i]],altname[ix[i]],
        #                        Lx[ix[i]],Lbol[ix[i]],Fbol[ix[i]],mstarK11[ix[i]])

        #print "%g %g"%(ra[ixsort[i]],dec[ixsort[i]])
        print "%g %g"%(ra[i],dec[i])

def update_data(path="/n/home00/lblecha/dwarf_agn_data/",
                fname="dwarf_sample.dat",
                #new_fname="dwarf_sample_new.dat",
                update_lx_14_195=True,add_lx_15_55=False,
                add_b13_extra_data=True,bol_corr_x=15.0):

    if add_b13_extra_data and not update_lx_14_195: 
        print "keyword 'add_b13_extra_data requires update_lx_14_195=True."
        return -1
    if not update_lx_14_195 and not add_lx_15_55: return 0

    new_fname_ext=''
    if add_b13_extra_data: new_fname_ext=new_fname_ext+'_with_b13_extras'
    if add_lx_15_55: new_fname_ext=new_fname_ext+'_with_a12lx'
    new_fname = 'dwarf_sample%s.dat'%new_fname_ext

    catID, name, altname, ra, dec, Lx, Lbol, Bmag, AGNtype, z, optspec_str, bhmass_low, bhmass_hi, mstarK11, mstarother_low, mstarother_hi, vdisp_low, vdisp_low_err, vdisp_hi, vdisp_hi_err, pctAGNrband, dist, gal_gmag, gal_rmag, gminusr, hasHST, hasChandra = read_data(path=path,fname=fname)
    print "Read data for %d sources."%len(ra)

    ###first update Lx values so they *all* match baumgartner+13
    if add_b13_extra_data or update_lx_14_195:

        m_catID,catmatch_lx,catmatch_fx,catmatch_snr,catmatch_gamma = matchcat.batcat(cat='b13',max_sep=600,fname=fname,return_extra_data=True)
    
        catmatch_lbol = np.zeros((len(catmatch_lx))) - 1
        catmatch_lbol[catmatch_lx!=-1] = np.log10(bol_corr_x*10**catmatch_lx[catmatch_lx!=-1])
        diff=catmatch_lx-Lx
        Lx = catmatch_lx
        Lbol = catmatch_lbol

    if add_lx_15_55:
        m_catID,catmatch_lx_a12,catmatch_fx,catmatch_snr,catmatch_gamma = matchcat.batcat(cat='a12',max_sep=600,fname=fname,return_extra_data=True)

    new_fpath="%s/%s"%(path,new_fname)
    if os.path.exists(new_fpath):
        val=raw_input("Warning: file %s exists. Overwrite? [y/n] "%new_fpath)
        if val != 'y': return -1
            
    header = "catID name altname ra dec Lx(14-195) Lbol Bmag AGNtype z optspec bhmass_low bhmass_hi mstarK11 mstarother_low mstarother_hi vdisp_low vdisp_low_err vdisp_hi vdisp_hi_err pctAGNrband dist gmagcorr rmagcorr gminusr hasHST hasChandra "
    if add_b13_extra_data: header=header+"Fx SNR Gamma "
    if add_lx_15_55: header=header+"Lx(15-55) "
    with open(new_fpath,"w") as nfp:
        nfp.write("%s\n"%header)
        for i in range(len(catID)):
            pstr=(3*'%s '+2*'%.10g '+2*'%.4g '+3*'%g '+'%s '+14*'%g '+2*'%s ')%(catID[i], name[i], altname[i], ra[i], dec[i], Lx[i], Lbol[i], Bmag[i], AGNtype[i], z[i], optspec_str[i], bhmass_low[i], bhmass_hi[i], mstarK11[i], mstarother_low[i], mstarother_hi[i], vdisp_low[i], vdisp_low_err[i], vdisp_hi[i], vdisp_hi_err[i], pctAGNrband[i], dist[i], gal_gmag[i], gal_rmag[i], gminusr[i], hasHST[i], hasChandra[i])
            if add_b13_extra_data:
                pstr=pstr+(3*'%.4g ')%(catmatch_fx[i],catmatch_snr[i],
                                       catmatch_gamma[i])
            if add_lx_15_55: 
                pstr=pstr+('%.4g ')%(catmatch_lx_a12[i])
            nfp.write('%s\n'%pstr)

        


