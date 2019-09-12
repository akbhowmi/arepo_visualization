import pyfits, sys, string
import numpy as np
import lb_utils as lb
import process_dwarf_sample as pds
from copy import copy
from itertools import islice
import matplotlib.pyplot as plt

def join_recarrays(arrays):

    ### expects list of recarrays as input
    newdtype = np.sum([a.dtype.descr for a in arrays])
    newrecarray = np.empty(len(arrays[0]),dtype=newdtype)
    for a in arrays:
        for name in a.dtype.names:
            newrecarray[name] = a[name]
    return newrecarray

def matchpos(ra_1, dec_1, z_1, ra_2, dec_2, z_2,
             ignore_z=False, multi_match_fatal=True,
             max_sep_arcsec=10, max_dz=0.01):

    strz='' if ignore_z else ', dz=%g'%max_dz
    print "\nchecking for match w/n rad=%g\"%s of each source."%(max_sep_arcsec,strz)

    max_sep = max_sep_arcsec / 3600.0 ## everything in degrees

    ra_1[ra_1>180] = ra_1[ra_1>180]-360
    ra_2[ra_2>180] = ra_2[ra_2>180]-360

    mask_1 = np.zeros((ra_1.size)).astype('bool')
    mask_2 = np.zeros((ra_2.size)).astype('bool')

    idx2_matched_1 = -1*np.ones((ra_1.size)).astype('int')
    idx1_matched_2 = -1*np.ones((ra_2.size)).astype('int')
    sep_matched_1 = -1*np.ones((ra_1.size))
    dz_matched_1 = -1*np.ones((ra_1.size))
    

    # do an initial check for small projected separations:
    for i in np.arange(ra_2.size):
        ix = np.where( (np.abs(ra_2[i]-ra_1) < 10*max_sep) &
                       (np.abs(dec_2[i]-dec_1) < 10*max_sep) )[0]
        if len(ix)==0: continue
        sep = np.array([ np.sqrt((ra_2[i]-ra_1[ix[j]])**2 + 
                                 (dec_2[i]-dec_1[ix[j]])**2) 
                         for j in np.arange(len(ix)) ])
        if i in [0.,   149., 290., 301., 438., 520., 565., 657., 1036.]:
            print "sep[%d]=%g"%(i,sep*3600)

        if ignore_z:
            sub_ix = np.where(sep<max_sep)[0]
        else:
            dz = np.array([ np.abs(z_2[i]-z) for z in z_1[ix] ])
            if i in [0.,   149., 290., 301., 438., 520., 565., 657., 1036.]:
                print "dz[%d]=%g"%(i,dz)
            sub_ix = np.where((sep<max_sep)&(dz<max_dz))[0]


        if (len(sub_ix)>1) or mask_1[ix[sub_ix]]:
            pstr="\nWARNING!!\nMultiple catalog matches found w/n %g\"!\n"
            print pstr%(max_sep_arcsec)
            print "i=%d: ra_2, dec_2: %g %g"%(i,ra_2[i],dec_2[i])
            print "ra_1, dec_1, sep: "
            print ra_1[ix],dec_1[ix]
            print "sep: ",sep[sub_ix]*3600
            if not ignore_z:
                print "max dz=%g"%(max_dz)
                print "z_2: %g"%z_2[i]
                print "z_1: ",z_1[ix]
            if not multi_match_fatal:
                sub_ix = sub_ix[sep[sub_ix]==sep[sub_ix].min()]
                if len(sub_ix) != 1: 
                    print "Error: len(sub_ix)=",len(sub_ix)
                    return -1
                print "Retaining best match: sep, ra_1, dec_1, z_1="
                print sep[sub_ix]*3600,ra_1[ix[sub_ix]],dec_1[ix[sub_ix]],z_1[ix[sub_ix]]
                print "\n"
            else:
                return -1

        if len(sub_ix)==1:
            mask_1[ix[sub_ix]] = True
            mask_2[i] = True
            idx2_matched_1[ix[sub_ix]] = i
            idx1_matched_2[i] = ix[sub_ix]
            sep_matched_1[ix[sub_ix]] = sep[sub_ix]*3600
            if not ignore_z: dz_matched_1[ix[sub_ix]] = dz[sub_ix]

    print "%d of %d elements were matched in catalog 1."%(ra_1[mask_1].size,ra_1.size)
    print "min/max sep of matched elements [\"]: %g %g"%(sep_matched_1[sep_matched_1!=-1].min(),sep_matched_1.max())
    print "%d of %d elements were matched in catalog 2."%(ra_2[mask_2].size,ra_2.size)
    ix1 = np.where(sep_matched_1>300)[0]
    if len(ix1)>0:
        if ignore_z:
            print "sep1 > 5': (sep1, ra1, dec1):"
            print sep_matched_1[ix1],ra_1[ix1],dec_1[ix1]
        else:
            print "sep1 > 5': (sep1, dz1, ra1, dec1):"
            print sep_matched_1[ix1],dz_matched_1[ix1],ra_1[ix1],dec_1[ix1]

    #for i in range(len(ra_1[mask_1])):
    #    print '%d %g %g'%(i,ra_1[mask_1][i], dec_1[mask_1][i])

    return idx2_matched_1, idx1_matched_2


def galcat(cat='nsa',path="/n/home00/lblecha/dwarf_agn_data/",
           fname='dwarf_sample_all.dat',max_sep_arcsec=10):

    if cat=='nsa':
        catfile = '%s/nsa_v0_1_2.fits'%path
    elif cat=='vagc':
        catfile = '%s/kcorrect.nearest.petro.z0.00.fits'%path
    else:
        print "Error: choose either catalog 'nsa' or 'vagc'."
        return -1

    hdu=pyfits.open(catfile)

    cat_ra=hdu[1].data.field('RA')
    cat_dec=hdu[1].data.field('DEC')
    cat_mass=hdu[1].data.field('MASS')
    cat_z=hdu[1].data.field('Z')
    
    hdu.close()


    #catID, name, altname, ra, dec, Lx, Lbol, Bmag, AGNtype, z, optspec_str, bhmass_low, bhmass_hi, mstarK11, mstarother_low, mstarother_hi, vdisp_low, vdisp_low_err, vdisp_hi, vdisp_hi_err, pctAGNrband, dist, gal_gmag, gal_rmag, gminusr, hasHST, hasChandra

    tmpdata = pds.read_data(path=path,fname=fname)
    dwarf_catID=tmpdata[0]
    dwarf_name=tmpdata[1]
    dwarf_ra=tmpdata[3]
    dwarf_dec=tmpdata[4]
    dwarf_z=tmpdata[9]
    dwarf_mass=tmpdata[13]


    ixc_match,cidd_match = matchpos(cat_ra,cat_dec,cat_z,dwarf_catID,
                                    dwarf_ra,dwarf_dec,dwarf_z,max_sep_arcsec=max_sep_arcsec)

    nmatch = len(ixc_match)
    if nmatch==0:
        print "No matched of BAT dwarf AGN with %s catalog."%cat
        return 0

    print "Found %d matches of BAT dwarf AGN with %s catalog"%(nmatch,cat)

    m_cat_ra,m_cat_dec,m_cat_mass,m_cat_z = lb.apply_cuts(ixc_match, (cat_ra,cat_dec,cat_mass,cat_z))

    m_dwarf_catID,m_dwarf_name,m_dwarf_ra,m_dwarf_dec,m_dwarf_z,m_dwarf_mass = lb.apply_cuts(cidd_match, (dwarf_catID,dwarf_name,dwarf_ra,dwarf_dec,dwarf_z,dwarf_mass))

    print "Matched sources:"
    for k in range(nmatch):
        print (7*"%s ")%(m_dwarf_name[k],m_dwarf_ra[k],m_cat_ra[k],
                         m_dwarf_dec[k],m_cat_dec[k],
                         np.log10(m_cat_mass[k]),m_dwarf_mass[k])

    ix = np.where(np.log10(m_cat_mass)<10.0)[0]
    if len(ix)>0:
        print "\n%d Matched sources with M*(cat)<1e10:"%len(ix)
        for k in range(len(ix)):
            print (7*"%s ")%(m_dwarf_catID[ix[k]],
                             m_dwarf_ra[ix[k]],m_cat_ra[ix[k]],
                             m_dwarf_dec[ix[k]],m_cat_dec[ix[k]],
                             np.log10(m_cat_mass[ix[k]]),m_dwarf_mass[ix[k]])
            

    #fp=open('nyucat_matched_dwarfs.txt','w')
    #fp.write("Name           RA      Dec      z       M*(Koss)    M*(NYU)\n")
    #for i in range(dwarfra.size):
    #    print "Matched source %s:"%dwarfname[i]
    #	print "From Koss et al: ra,dec=%g,%g, z=%g, log mass=%g"%(dwarfra[i],dwarfdec[i],dwarfz[i],dwarfmass[i])
    #	print "From NYU cat: ra,dec=%g,%g, z=%g, log mass=%g"%(ra[match_index[i]],dec[match_index[i]],z[match_index[i]],np.log10(mass[match_index[i]]))
    #	print "%s %g %g %g %g %g"%(dwarfname[i],dwarfra[i],dwarfdec[i],dwarfz[i],dwarfmass[i],np.log10(mass[match_index[i]]))
    #	fp.write("%20s  %g  %g  %g  %g  %g\n"%(dwarfname[i],dwarfra[i],dwarfdec[i],dwarfz[i],dwarfmass[i],np.log10(mass[match_index[i]])))
    #fp.close()

    return 0

def batcat_xmatch(cat='a12',path="/n/home00/lblecha/dwarf_agn_data/",
                  fname='dwarf_sample_with_b13_extras.dat',max_sep_arcsec=600,
                  return_extra_data=False):

    if cat=='a12':
        catfile = '%s/Ajello_2012.fits'%path
        lstr='logLX'
        ffac=1.0e-11
    elif cat=='b13':
        catfile = '%s/Baumgartner_2013.fits'%path
        lstr='logL'
        ffac=1.0e-12
    else:        
        print "Error: choose either catalog 'a12' or 'b13'."
        return -1

    hdu=pyfits.open(catfile)

    cat_ra=hdu[1].data.field('RAJ2000')
    cat_dec=hdu[1].data.field('DEJ2000')
    cat_z=hdu[1].data.field('z')
    cat_flux=hdu[1].data.field('Flux') * ffac ## convert to erg/s/cm^2
    cat_lx=hdu[1].data.field(lstr) ## log erg/s
    cat_snr=hdu[1].data.field('S_N')
    cat_gamma=hdu[1].data.field('Gamma')
    #cat_cl=hdu[1].data.field('Cl') #(b13 only)

    ix_hasz = np.where((cat_z==cat_z)&(cat_z>0))[0]
    cat_ra,cat_dec,cat_z,cat_flux,cat_lx,cat_snr,cat_gamma = lb.apply_cuts(ix_hasz,(cat_ra,cat_dec,cat_z,cat_flux,cat_lx,cat_snr,cat_gamma))
    print "min/max cat_z:",cat_z.min(),cat_z.max()
    print "min/max cat_lx:",cat_lx.min(),cat_lx.max()

    hdu.close()
    
    tmpdata = pds.read_data(path=path,fname=fname)
    dwarf_catID=tmpdata[0]
    dwarf_name=tmpdata[1]
    dwarf_ra=tmpdata[3]
    dwarf_dec=tmpdata[4]
    dwarf_lx=tmpdata[5]
    dwarf_z=tmpdata[9]

    #ix=np.where(dwarf_z<0.002)[0]
    #print ["%g %g"%(dwarf_ra[ix[i]],dwarf_dec[ix[i]]) for i in range(len(ix))]
    #ixb=np.where(cat_z<0.002)[0]
    #print ["%g %g"%(cat_ra[ixb[i]],cat_dec[ixb[i]]) for i in range(len(ixb))]

    ixc_match,cidd_match = matchpos(cat_ra,cat_dec,cat_z,
                                    dwarf_ra,dwarf_dec,dwarf_z,
                                    max_sep_arcsec=max_sep_arcsec)
    nmatch = len(ixc_match)
    if nmatch==0:
        print "No matches of BAT dwarf AGN with %s catalog."%cat
        return 0

    print "Found %d matches of BAT dwarf AGN with %s catalog"%(nmatch,cat)

    tmparr = np.zeros((len(dwarf_z))) - 1
    dwarf_catmatch_lx = copy(tmparr)
    dwarf_catmatch_lx[cidd_match] = cat_lx[ixc_match]
    dwarf_catmatch_fx = copy(tmparr)
    dwarf_catmatch_fx[cidd_match] = cat_flux[ixc_match]
    dwarf_catmatch_snr = copy(tmparr)
    dwarf_catmatch_snr[cidd_match] = cat_snr[ixc_match]
    dwarf_catmatch_gamma = copy(tmparr)
    dwarf_catmatch_gamma[cidd_match] = cat_gamma[ixc_match]

    ### sanity check:
    dwarf_catmatch_ra = copy(tmparr)
    dwarf_catmatch_ra[cidd_match] = cat_ra[ixc_match]
    #print "dwarf_catmatch_ra[cidd_match], cat_ra[ixc_match]:"
    #print ["%g %g"%(dwarf_catmatch_ra[cidd_match[k]], cat_ra[ixc_match[k]]) for k in range(len(cidd_match))]
    print "len(dwarf_catmatch_ra[cidd_match])=%d, len(cat_ra[ixc_match])=%d"%(len(dwarf_catmatch_ra[cidd_match]),len(cat_ra[ixc_match]))
    print "this should be 0: %g"%(np.abs(dwarf_catmatch_ra[cidd_match]-cat_ra[ixc_match])).max()
    print "max ra diff: %g"%(np.abs(dwarf_catmatch_ra[cidd_match]-dwarf_ra[cidd_match])).max()

    m_cat_ra,m_cat_dec,m_cat_z,m_cat_flux,m_cat_lx,m_cat_snr,m_cat_gamma = lb.apply_cuts(ixc_match,(cat_ra,cat_dec,cat_z,cat_flux,cat_lx,cat_snr,cat_gamma))

    m_dwarf_catID,m_dwarf_name,m_dwarf_ra,m_dwarf_dec,m_dwarf_lx,m_dwarf_z = lb.apply_cuts(cidd_match, (dwarf_catID,dwarf_name,dwarf_ra,dwarf_dec,dwarf_lx,dwarf_z))

    if return_extra_data:
        print "min/max lx:",dwarf_catmatch_lx.min(),dwarf_catmatch_lx.max()
        return dwarf_catID,dwarf_catmatch_lx,dwarf_catmatch_fx,dwarf_catmatch_snr,dwarf_catmatch_gamma
    else:
        print "Matched sources:"
        for k in range(nmatch):
            #print (8*"%s ")%(m_dwarf_ra[k],m_cat_ra[k],
            #                 m_dwarf_dec[k],m_cat_dec[k],
            #                 m_dwarf_z[k],m_cat_z[k],
            #                 m_dwarf_lx[k],m_cat_lx[k])
            #print (4*"%0.3g  ")%(np.abs(m_dwarf_ra[k]-m_cat_ra[k]),
            #                     np.abs(m_dwarf_dec[k]-m_cat_dec[k]),
            #                     np.abs(m_dwarf_z[k]-m_cat_z[k]),
            #                     np.abs(m_dwarf_lx[k]-m_cat_lx[k]))
            #print "%s %.3g"%(m_dwarf_name[k],np.abs(m_dwarf_lx[k]-m_cat_lx[k]))
            print "%s  %g %g  %g"%(m_dwarf_name[k],m_dwarf_lx[k],m_cat_lx[k],
                                   (m_dwarf_lx[k]-m_cat_lx[k]))
        return (0,0)

def batcat(cat='b13',path="/n/home00/lblecha/dwarf_agn_data/",
                  max_sep_arcsec=600, return_extra_data=False):

    if cat=='a12':
        catfile = '%s/Ajello_2012.fits'%path
        fstr='Flux'
        ffac=1.0e-11
    elif cat=='b13':
        catfile = '%s/Baumgartner_2013.fits'%path
        fstr='Flux'
        ffac=1.0e-12
    elif cat=='b11':
        catfile = '%s/Burlon_2011.fits'%path
        fstr='Flux'
        ffac=1.0e-11
    elif cat=='t08':
        catfile = '%s/Tueller_2008.fits'%path
        fstr='fBAT'
        ffac=1.0e-11
    else:        
        print "Error: choose either catalog 'a12' or 'b13'."
        return -1

    print "\nOpening catalog file: %s\n"%catfile
    f=pyfits.open(catfile)
    d=f[1].data
    colnames = np.array(d.columns.names)
    print "column names:",colnames
    f.close()

    ### B13 column names: ###
    #['SWIFT', 'RAJ2000', 'DEJ2000', 'S_N', 'CName', 'OName', 'RACdeg',
    # 'DECdeg', 'Flux', 'e_Flux', 'E_Flux1', 'Cont', 'Gamma', 'e_Gamma',
    # 'E_Gamma1', 'chi2', 'z', 'u_z', 'logL', 'AS', 'Cl', 'Type', 'P',
    #  'SimbadName', 'NED'] 

    ### A12 column names: ###
    # ['SWIFT', 'RAJ2000', 'DEJ2000', 'ePos', 'Flux', 'S_N', 'A', 'ID',
    #  'u_ID', 'Type', 'z', 'Gamma', 'e_Gamma', 'logLX', 'XCat',
    #  'SimbadName', 'NED']

    ### B11 column names: ###
    # ['SWIFT', 'RAJ2000', 'DEJ2000', 'Flux', 'S_N', 'Name', 'Type', 'z',
    #  'Gamma', 'l_logNH', 'logNH', 'r_logNH', 'Xm', 'Simbad']

    ### T08 column names: ###
    # ['Seq', 'SWIFT', 'm_SWIFT', 'n_SWIFT', 'RAJ2000', 'DEJ2000', 'Flag',
    #  'SNR', 'fBAT', 'z', 'logL', 'logNH', 'Cplx', 'Type', 'Jmag',
    #  'fROSAT', 'SimbadName', 'RX', 'NED']

    print "Found %d catalog entries."%d.size
    #mask = (d['z']==d['z'])&(d['z']>0)
    #d = d[mask]
    #print "Retained %d catalog entries with valid redshifts."%d.size

    d[fstr] = ffac * d[fstr] ## convert to erg/s/cm^2

    print "min/max z:",d['z'].min(),d['z'].max()
    print "min/max flux:",d[fstr].min(),d[fstr].max()

    if cat=='b13':
        maskb = d['SimbadName']!='2MASX J15064412+0351444'
        print "Excluding source (%s) with incorrect ra,dec in catalog."%d['SimbadName'][~maskb]
        d = d[maskb]
    #ix= d['Name']=='2MASX J15064412+0351444'
    #print d['Name'][ix]
    #print d['RAJ2000'][ix]
    #print d['DEJ2000'][ix]

    return d

def kosscat(path="/n/home00/lblecha/dwarf_agn_data/",
            batcat_id='b13',remove_ucontam=True,PSr_cut=100, 
            include_b13_class2=True,use_BLflag=False, return_raw_cat_only=False, 
            write_file=False,load_nathan_masses=False):

    if batcat_id != 'b13' and not return_raw_cat_only:
        print "Error: invalid batcat_id=%s."%batcat_id
        print "  xmatch with Koss data currently requires BAT catalog 'b13'."
        
    kosscatfile='%s/Koss_2011_all.fits'%path
    
    f=pyfits.open(kosscatfile)
    d = f[1].data
    colnames = np.array(d.columns.names)
    print "processing data from %s..."%kosscatfile
    print "column names:",colnames

    #['Name', 'M', 'umag', 'gmag', 'rmag', 'imag', 'zmag', 'logM', 'PSr',
    # 'E_g-r', 'E_g-r_i', 'ucontam', 'contam', 'BAT', 'SimbadName', 'NED',
    # 'RA', 'DE', 'Name_kpno', 'Date_kpno', 'Type_kpno', 'z_kpno',
    # 'Dist_kpno', 'E_B-V_kpno', 'Air_kpno', 'PSF_kpno', 'Name_sdss',
    # 'Date_sdss', 'Type_sdss', 'z_sdss', 'Dist_sdss', 'E_B-V_sdss',
    # 'Air_sdss', 'PSF_sdss', 'BAT_sdss', 'SimbadName_sdss', 'NED_sdss',
    # 'RA_sdss', 'DE_sdss', 'Name_morph', 'Rp_morph', 'C_morph',
    # 'Class_morph', 'b_a_morph', 'BLflag_morph']
    
    if return_raw_cat_only:    return d

    print "%d table entries found in %s"%(d.size,kosscatfile)
    Nm10raw = d[d['logM']<10].size
    Nm97raw = d[d['logM']<9.7].size
    Nm95raw = d[d['logM']<9.5].size
    print "# with logM<10:",Nm10raw
    print "# with logM<9.7:",Nm97raw
    print "# with logM<9.5:",Nm95raw

    print "These don't have any matching data in the morphology or SDSS/KPNO tables: "
    print d['Name'][d['M']!='M']

    print "number with kpno z: " ,d[d['z_kpno']==d['z_kpno']].size
    print "number with sdss z: ", d[d['z_sdss']==d['z_sdss']].size
    print "number with both z's: ",d[(d['z_sdss']==d['z_sdss'])&
                                     (d['z_kpno']==d['z_kpno'])].size
    print "number with neither z: ",d[(d['z_sdss']!=d['z_sdss'])&
                                      (d['z_kpno']!=d['z_kpno'])].size
    print "names with neither z: ",d['Name'][(d['z_sdss']!=d['z_sdss'])&
                                             (d['z_kpno']!=d['z_kpno'])]
    print "masses with neither z: ",d['logM'][(d['z_sdss']!=d['z_sdss'])&
                                              (d['z_kpno']!=d['z_kpno'])]


    if load_nathan_masses:

        nd = nathan_masses()
        nathan_idxk,koss_idxn = matchpos(nd['CTPT_RA'],nd['CTPT_Dec'],nd['Best_z'],
                                         d['RA'],d['DE'],d['z_sdss'],
                                         max_sep_arcsec=600,ignore_z=True,
                                         multi_match_fatal=True)
        print d.size, nd.size
        mk=np.zeros(d.size).astype('bool')
        mk[(koss_idxn!=-1)]=True
        mn=np.zeros(nd.size).astype('bool')
        mn[(nathan_idxk!=-1)]=True
        print d[mk].size
        print nd[mn].size
        #print d['Name_sdss'][mk]
        #print d['z_sdss'][mk]
        #print d['Name_kpno'][mk]
        #print d['z_kpno'][mk]
        print d[d['logM']<10].size
        print d[mk][d['logM'][mk]<10].size
        print d[~mk][d['logM'][~mk]<10].size
        print d['Name'][~mk][d['logM'][~mk]<10]
        print d['Name_sdss'][~mk][d['logM'][~mk]<10]
        print d['Name_kpno'][~mk][d['logM'][~mk]<10]
        print nd['CTPT_Name'][~mn]
        print nd['CTPT_Name'][~mn].size
        print nd['CTPT_Name'][mn].size

        fig=plt.figure(figsize=(5,5))
        lims=(7.5,10.6)
        plt.xlim(lims)
        plt.ylim(lims)
        plt.xlabel('log M* (Koss 2011)')
        plt.ylabel('log M* (Nathan/NSA)')
        plt.plot([7.5,10.6],[7.5,10.6],'g-')
        plt.plot(d['logM'][mk],nd['MASS'][mn],'ko',ms=4)
        print len(d[~mk][d['logM'][~mk]<10])
        for i in range(len(d[~mk][d['logM'][~mk]<10])):
            print d[~mk][d['logM'][~mk]<10]['logM'][i]
            plt.plot(np.tile(d[~mk][d['logM'][~mk]<10]['logM'][i],2),[lims[0],lims[0]+0.25],'b-')
        print len(nd[~mn][nd['MASS'][~mn]<10])
        for i in range(len(nd[~mn][nd['MASS'][~mn]<10])):
            print nd[~mn][nd['MASS'][~mn]<10]['MASS'][i]
            plt.plot([lims[0],lims[0]+0.25],np.tile(nd[~mn][nd['MASS'][~mn]<10]['MASS'][i],2),'m-')
        fig.savefig('%s/mass_compare.pdf'%path)


        plt.clf()
        plt.cla()
        fig=plt.figure(figsize=(8,8))
        ax1=fig.add_subplot(221)
        ax1.plot(nd['L_BAT_bestz'][nd['MASS']<10],nd['MASS'][nd['MASS']<10],'ko',ms=4)
        ax1.set_xlabel('log L_BAT')
        ax2=fig.add_subplot(222)
        ax2.set_xlabel('log Lbol')
        ax2.plot(nd['log_L_bol_k8'][nd['MASS']<10],nd['MASS'][nd['MASS']<10],'ko',ms=4)
        ax3=fig.add_subplot(223)
        ax3.set_xlabel('W1-W2')
        ax3.plot(nd['w1'][nd['MASS']<10]-nd['w2'][nd['MASS']<10],nd['MASS'][nd['MASS']<10],'ko',ms=4)
        fig.subplots_adjust()
        fig.savefig('%s/nathan_lowmass_info.pdf'%path)

       ### column names: 
       # 'BAT_index' 'BAT_Name' 'CTPT_Name' 'Other_Name' 'CTPT_RA' 'CTPT_Dec'
       # 'NED_Type' 'NED_redshift' 'SNR' 'FLUX' 'FLUX_LO' 'FLUX_HI' 'CONTA'
       # 'Redshift_Ind_Dist' 'Best_z' 'Lum_dist' 'Best_dist' '70-mo_L14-195'
       # 'L_BAT_bestz' 'log_L_bol_k8' 'log_L_bol_high' 'log_L_bol_low' 'j_m_ext'
       # 'j_msig_ext' 'h_m_ext' 'h_msig_ext' 'k_m_ext' 'k_msig_ext' 'w1' 'w1e' 'w2'
       # 'w2e' 'w3' 'w3e' 'w4' 'w4e' 'uid' 'coeff' 'EBV' 'chi2' 'MASS' 



       # return



    tmpmask=((d['Type_sdss']!='')&(d['Type_kpno']!=''))
    if d[(d['Type_sdss']!=d['Type_kpno'])&(tmpmask)].size>0:
        print "Error: mismatch in Koss AGN types:"
        print d['Type_sdss'][tmpmask]
        print d['Type_kpno'][tmpmask]
        sys.exit()
    tmptype = d['Type_sdss']
    tmptype[tmptype==''] = d['Type_kpno'][tmptype=='']
    tmpmaskII = ( (tmptype!='Sy1')&(tmptype!='Sy1.2')&(tmptype!='Sy1.5') )
    tmpmaskI = ( (tmptype=='Sy1')|(tmptype=='Sy1.2')|(tmptype=='Sy1.5') )

    print "%d sources are TypeI, %d are TypeII"%(d[tmpmaskI].size,d[tmpmaskII].size)
    print "# Type I, Type II with logM<10: %g %g"%(d[(tmpmaskI)&(d['logM']<10)].size,d[(tmpmaskII)&(d['logM']<10)].size)
    print "# Type I, Type II with logM<9.7: %g %g"%(d[(tmpmaskI)&(d['logM']<9.7)].size,d[(tmpmaskII)&(d['logM']<9.7)].size)
    print "# Type I, Type II with logM<9.5: %g %g"%(d[(tmpmaskI)&(d['logM']<9.5)].size,d[(tmpmaskII)&(d['logM']<9.5)].size)

    print "\nsources with 'ucontam' flag (AGN>20% in u-band): "
    print "number=%d, frac of total=%g"%(d[d['ucontam']=='Y'].size,
                                         1.0*d[d['ucontam']=='Y'].size/d.size)
    print "N(logM<10)=%d, f(logM<10)=%g"%(d[(d['ucontam']=='Y')&(d['logM']<10)].size, 1.0*d[(d['ucontam']=='Y')&(d['logM']<10)].size/Nm10raw)
    print "N(logM<9.7)=%d, f(logM<9.7)=%g"%(d[(d['ucontam']=='Y')&(d['logM']<9.7)].size, 1.0*d[(d['ucontam']=='Y')&(d['logM']<9.7)].size/Nm97raw)
    print "N(logM<9.5)=%d, f(logM<9.5)=%g"%(d[(d['ucontam']=='Y')&(d['logM']<9.5)].size, 1.0*d[(d['ucontam']=='Y')&(d['logM']<9.5)].size/Nm95raw)

    print "number(typeII)=%d, frac of total typeII=%g"%(d[(d['ucontam']=='Y')&(tmpmaskII)].size, 1.0*d[(d['ucontam']=='Y')&(tmpmaskII)].size/d[tmpmaskII].size)
    print "N(typeII,logM<10)=%d, f(typeII,logM<10)=%g"%(d[(d['ucontam']=='Y')&(d['logM']<10)&(tmpmaskII)].size, 1.0*d[(d['ucontam']=='Y')&(d['logM']<10)&(tmpmaskII)].size/d[(d['logM']<10)&(tmpmaskII)].size)
    print "N(typeII,logM<9.7)=%d, f(typeII,logM<9.7)=%g"%(d[(d['ucontam']=='Y')&(d['logM']<9.7)&(tmpmaskII)].size, 1.0*d[(d['ucontam']=='Y')&(d['logM']<9.7)&(tmpmaskII)].size/d[(d['logM']<9.7)&(tmpmaskII)].size)
    print "N(typeII,logM<9.5)=%d, f(typeII,logM<9.5)=%g"%(d[(d['ucontam']=='Y')&(d['logM']<9.5)&(tmpmaskII)].size, 1.0*d[(d['ucontam']=='Y')&(d['logM']<9.5)&(tmpmaskII)].size/d[(d['logM']<9.5)&(tmpmaskII)].size)

    print "\nsources with 'contam' flag (AGN>35% in griz bands): "
    print "Ntot=%d, ftot=%g"%(d[d['contam']=='Y'].size,
                              1.0*d[d['contam']=='Y'].size/d.size)
    print "N(logM<10)=%d, f(logM<10)=%g"%(d[(d['contam']=='Y')&(d['logM']<10)].size, 1.0*d[(d['contam']=='Y')&(d['logM']<10)].size/Nm10raw)
    print "N(logM<9.7)=%d, f(logM<9.7)=%g"%(d[(d['contam']=='Y')&(d['logM']<9.7)].size, 1.0*d[(d['contam']=='Y')&(d['logM']<9.7)].size/Nm97raw)
    print "N(logM<9.5)=%d, f(logM<9.5)=%g"%(d[(d['contam']=='Y')&(d['logM']<9.5)].size, 1.0*d[(d['contam']=='Y')&(d['logM']<9.5)].size/Nm95raw)
    print "number(typeII)=%d, frac of total typeII=%g"%(d[(d['contam']=='Y')&(tmpmaskII)].size, 1.0*d[(d['contam']=='Y')&(tmpmaskII)].size/d[tmpmaskII].size)
    print "N(typeII,logM<10)=%d, f(typeII,logM<10)=%g"%(d[(d['contam']=='Y')&(d['logM']<10)&(tmpmaskII)].size, 1.0*d[(d['contam']=='Y')&(d['logM']<10)&(tmpmaskII)].size/d[(d['logM']<10)&(tmpmaskII)].size)
    print "N(typeII,logM<9.7)=%d, f(typeII,logM<9.7)=%g"%(d[(d['contam']=='Y')&(d['logM']<9.7)&(tmpmaskII)].size, 1.0*d[(d['contam']=='Y')&(d['logM']<9.7)&(tmpmaskII)].size/d[(d['logM']<9.7)&(tmpmaskII)].size)
    print "N(typeII,logM<9.5)=%d, f(typeII,logM<9.5)=%g"%(d[(d['contam']=='Y')&(d['logM']<9.5)&(tmpmaskII)].size, 1.0*d[(d['contam']=='Y')&(d['logM']<9.5)&(tmpmaskII)].size/d[(d['logM']<9.5)&(tmpmaskII)].size)

    print "\nsources with PSr>20 (AGN>20% in r-band): "
    print "Ntot=%d, ftot=%g"%(d[d['PSr']>20].size, 1.0*d[d['PSr']>20].size/d.size)
    print "N(logM<10)=%d, f(logM<10)=%g"%(d[(d['PSr']>20)&(d['logM']<10)].size, 1.0*d[(d['PSr']>20)&(d['logM']<10)].size/Nm10raw)
    print "N(logM<9.7)=%d, f(logM<9.7)=%g"%(d[(d['PSr']>20)&(d['logM']<9.7)].size, 1.0*d[(d['PSr']>20)&(d['logM']<9.7)].size/Nm97raw)
    print "N(logM<9.5)=%d, f(logM<9.5)=%g"%(d[(d['PSr']>20)&(d['logM']<9.5)].size, 1.0*d[(d['PSr']>20)&(d['logM']<9.5)].size/Nm95raw)

    print "\nsources with PSr>10 (AGN>10% in r-band): "
    print "Ntot=%d, ftot=%g"%(d[d['PSr']>10].size, 1.0*d[d['PSr']>10].size/d.size)
    for M in [10,9.7,9.5]:
        print "N(logM<%g)=%d, f(logM<%g)=%g"%(M,d[(d['PSr']>10)&(d['logM']<M)].size, M,1.0*d[(d['PSr']>10)&(d['logM']<M)].size/d[d['logM']<M].size)
    print "number(typeII)=%d, frac of total typeII=%g"%(d[(d['PSr']>10)&(tmpmaskII)].size, 1.0*d[(d['PSr']>10)&(tmpmaskII)].size/d[tmpmaskII].size)
    for M in [10,9.7,9.5]:
        print "N(typeII,logM<%g)=%d, f(typeII,logM<%g)=%g"%(M,d[(d['PSr']>10)&(d['logM']<M)&(tmpmaskII)].size, M,1.0*d[(d['PSr']>10)&(d['logM']<M)&(tmpmaskII)].size/d[(d['logM']<M)&(tmpmaskII)].size)



    #mask = (d['M']=='M')&(d['contam']!='Y')&(d['ucontam']!='Y')
    mask = ( (d['M']=='M')&(d['contam']!='Y')&
             ((d['z_kpno']==d['z_kpno'])|(d['z_sdss']==d['z_sdss'])) )
    pstr='griz contamination'
    if remove_ucontam:
        pstr='contamination'
        mask = ((mask)&(d['ucontam']!='Y'))
    if PSr_cut < 100:
        pstr=pstr+' (& PSr<%g)'%(PSr_cut)
        mask = ((mask)&(d['PSr']<PSr_cut))

    d = d[mask]
    print "\nRetained %d entries with observations, with z, and no %s."%(d.size,pstr)

    
    ### define a single redshift for each source
    ### (using the sdss redshift when both are present).
    z_koss = d['z_sdss']
    z_koss[z_koss!=z_koss] = d['z_kpno'][z_koss!=z_koss]
    if len(z_koss[z_koss!=z_koss]) > 0:
        print "error: still have nonzero entries with no redshift:"
        print "names with neither z: "
        print d['Name'][(d['z_sdss']!=d['z_sdss'])&(d['z_kpno']!=d['z_kpno'])]
        return -1


    ### Read data from BAT catalog (B13, specifically)
    bat_cat = batcat(cat=batcat_id)
    ### avoid duplicate key names:

    #print "koss:"
    #print d['Name'][(d['RA']<174)&(d['RA']>173)&(d['DE']>52)&(d['DE']<54)]
    #print d['Name'][(d['RA']<-7)&(d['RA']>-8)&(d['DE']>3)&(d['DE']<4)]
    #print d['Name'][(d['RA']>135)&(d['RA']<136)&(d['DE']>60)&(d['DE']<61)]
    #print "bat:"
    #print bat_cat['SimbadName'][(bat_cat['RAJ2000']<174)&(bat_cat['RAJ2000']>173)&
    #                      (bat_cat['DEJ2000']>52)&(bat_cat['DEJ2000']<54)]
    #print bat_cat['RACdeg'][(bat_cat['RAJ2000']>352)&(bat_cat['RAJ2000']<353)&
    #                      (bat_cat['DEJ2000']>3)&(bat_cat['DEJ2000']<4)]

    ### Cross-match with Koss catalog
    #koss_mask,bat_mask = matchpos(d['RA'],d['DE'],z_koss,bat_cat['RAJ2000'],
    koss_idxbat,bat_idxkoss = matchpos(d['RA'],d['DE'],z_koss,bat_cat['RAJ2000'],
                                       bat_cat['DEJ2000'],bat_cat['z'],
                                       max_sep_arcsec=600,ignore_z=False,
                                       multi_match_fatal=False)

    print "\nSources in Koss catalog unmatched in BAT catalog:"
    print 'Name:',d['Name'][koss_idxbat==-1]
    print 'logM:',d['logM'][koss_idxbat==-1]

    print "\nLowest-mass source names (logM<9.5):",d['Name'][d['logM']<9.5]

    d_new = join_recarrays([ koss_idxbat.astype([('koss_idxbat','>i4')]),
                             z_koss.astype([('z_koss','>f4')]), d ])
    bat_cat_new = join_recarrays([bat_idxkoss.astype([('bat_idxkoss','>i4')]), 
                                  bat_cat])
    bcnames = np.array(bat_cat_new.dtype.names)
    bcnames[bcnames=='SimbadName']='Simbad_BAT'
    bcnames[bcnames=='NED']='NED_BAT'
    bat_cat_new.dtype.names = tuple(bcnames)
    if d.size!=d_new.size or bat_cat.size!=bat_cat_new.size:
        print "Error: size mismatch b/t d,d_new or bat_cat,bat_cat_new."
        return -1

    dmatch = d_new[d_new['koss_idxbat']!=-1]
    bcmatch = bat_cat_new[bat_cat_new['bat_idxkoss']!=-1]
    bcmatch = bcmatch[np.argsort(bcmatch['bat_idxkoss'])]
    print "\ndmatch has %d entries (with %d cols)."%(dmatch.size,len(dmatch[0]))
    print "bcmatch has %d entries (with %d cols)."%(bcmatch.size,len(bcmatch[0]))
    if dmatch.size != bcmatch.size:
        print "Error: dmatch and bcmatch do not have same length."
        return -1

    #print dmatch['z_koss'].size
    #print bcmatch['z'].size
    #diff=(dmatch['z_koss']-bcmatch['z'])/bcmatch['z']
    #print "min/max/mean/med diff z: ",diff.min(),diff.max(),diff.mean(),np.median(diff)
    #print "min/max/mean/med abs(diff z): ",np.abs(diff).min(),np.abs(diff).max(),np.abs(diff).mean(),np.median(np.abs(diff))
    #print "abs diff z > 1%:"
    #print dmatch['Name'][np.abs(diff)>0.1].size
    #print dmatch['Name'][np.abs(diff)>0.1]
    
    #print np.argsort(bat_cat_new[bat_cat_new['bat_idxkoss']!=-1])
    #print np.argsort(bat_idxkoss[bat_idxkoss!=-1])
    #print np.argsort(bcmatch['bat_idxkoss'])
    #print np.sort(bcmatch['bat_idxkoss'])
    #print np.where(koss_idxbat!=-1)[0]

    alldata = join_recarrays([dmatch,bcmatch])
    print "alldata has %d entries (with %d cols).\n"%(alldata.size,len(alldata[0]))

    if batcat_id=='b13' and not include_b13_class2:
        c2mask = np.zeros((alldata.size)).astype('bool')
        c2mask[alldata['Cl']!=2] = True                                                         
        print "Excluding %d objects w/ B13 Class=2 ('Galaxy')"%alldata[~c2mask].size
        alldata = alldata[c2mask]
        print "Retained %d objects"%alldata.size

    #print alldata.dtype.names
    #print alldata['koss_idxbat']
    #print alldata['bat_idxkoss']
    #for j in range(alldata.size):
    #    print alldata['Name'][j], alldata['CName'][j]
    #    print alldata['RA'][j], alldata['RAJ2000'][j]


    if write_file:
        contamstr = 'ugriz' if remove_ucontam else 'griz'
        class2str = '' if include_b13_class2 else '_noClass2'
        with open(path+'/kosscat_no_%s_contam%s.dat'%(contamstr,class2str),"w") as fp:
            fp.write('Name  RA   DEC   logM   SwiftName   B13Class\n')
            for i in range(alldata.size):
                fp.write('%s %g %g %g %s %d\n'%(alldata['Name'][i],alldata['RA'][i],alldata['DE'][i],
                                                alldata['logM'][i],alldata['SWIFT'][i],alldata['Cl'][i]))
        return



    print "BL vs non-BL AGN (Koss definitions):"
    print "(Type=(Sy1,1.2,or1.5) & BLflag_morph=Y):"
    BLAGNstr = ['Sy1','Sy1.2','Sy1.5']
    #ix_noBL = [ j for j in range(alldata.size) 
    #            if (alldata['BLflag_morph'][j]!='Y'
    #                and alldata['Type_kpno'][j] not in BLAGNstr 
    #                and alldata['Type_sdss'][j] not in BLAGNstr) ]
    Type_koss = alldata['Type_sdss']
    Type_koss[Type_koss==''] = alldata['Type_kpno'][Type_koss=='']
    mask_noBL = ( (Type_koss!='Sy1')&(Type_koss!='Sy1.2')&(Type_koss!='Sy1.5') )
    mask_BL = ( (Type_koss=='Sy1')|(Type_koss=='Sy1.2')|(Type_koss=='Sy1.5') )


    #mask_noBL = ( (alldata['Type_kpno']!='Sy1')&(alldata['Type_sdss']!='Sy1')&
    #              (alldata['Type_kpno']!='Sy1.2')&(alldata['Type_sdss']!='Sy1.2')&
    #              (alldata['Type_kpno']!='Sy1.5')&(alldata['Type_sdss']!='Sy1.5') )
    #mask_BL = ( (alldata['Type_kpno']=='Sy1')|(alldata['Type_sdss']=='Sy1')|
    #            (alldata['Type_kpno']=='Sy1.2')|(alldata['Type_sdss']=='Sy1.2')|
    #            (alldata['Type_kpno']=='Sy1.5')|(alldata['Type_sdss']=='Sy1.5') )
    if use_BLflag:
        mask_noBL = ((mask_noBL)&(alldata['BLflag_morph']!='Y'))
        mask_BL = ((mask_BL)|(alldata['BLflag_morph']=='Y'))

    print "# noBL: %d"%alldata[mask_noBL].size
    print "# BL: %d"%alldata[mask_BL].size
    print "# noBL with BLflag_morph=='Y': %d"%alldata[(mask_noBL)&(alldata['BLflag_morph']=='Y')].size
    print "  Names:",alldata['Name'][(mask_noBL)&(alldata['BLflag_morph']=='Y')]
    print "  Type_kpno:",alldata['Type_kpno'][(mask_noBL)&(alldata['BLflag_morph']=='Y')]
    print "  Type_sdss:",alldata['Type_sdss'][(mask_noBL)&(alldata['BLflag_morph']=='Y')]
    print "# BL with BLflag_morph!='Y': %d"%alldata[(mask_BL)&(alldata['BLflag_morph']!='Y')].size
    name_noBL=alldata['Name'][mask_noBL]
    logM_noBL=alldata['logM'][mask_noBL]

    Nm10 = alldata[(alldata['logM']<10)].size
    Nm97 = alldata[(alldata['logM']<9.7)].size
    Nm95 = alldata[(alldata['logM']<9.5)].size

    print "# with logM<10: %d"%Nm10
    print "# with logM_noBL<10:",logM_noBL[logM_noBL<10].size
    print "# with logM<9.7: %d"%Nm97
    print "# with logM_noBL<9.7:",logM_noBL[logM_noBL<9.7].size
    print "# with logM<9.5: %d"%Nm95
    print "# with logM_noBL<9.5:",logM_noBL[logM_noBL<9.5].size

    print "\ntotal frac. w/ no BL: %g"%(1.0*alldata[mask_noBL].size/alldata.size)
    print "frac. logM<10 w/ no BL: %g"%(1.0*logM_noBL[logM_noBL<10].size/Nm10)
    print "frac. logM<9.7 w/ no BL: %g"%(1.0*logM_noBL[logM_noBL<9.7].size/Nm97)
    print "frac. logM<9.5 w/ no BL: %g16"%(1.0*logM_noBL[logM_noBL<9.5].size/Nm95)

    print "\nBL vs non-BL AGN (B13 definitions):"
    mask_cl4 = (alldata['Cl']==4)
    mask_cl5 = (alldata['Cl']==5)
    mask_cl45 = ((mask_cl4)|(mask_cl5))
    print "%d Class 2 (Galaxy)"%alldata[alldata['Cl']==2].size
    print "%d Class 4 (Sy1-1.5)"%alldata[mask_cl4].size
    print "%d Class 5 (Sy1.7-2.0)"%alldata[mask_cl5].size
    print "%d Class 6 (other)"%alldata[alldata['Cl']==6].size

    Nm10_cl45 = alldata[(alldata['logM']<10)&(mask_cl45)].size
    Nm97_cl45 = alldata[(alldata['logM']<9.7)&(mask_cl45)].size
    Nm95_cl45 = alldata[(alldata['logM']<9.5)&(mask_cl45)].size
    print "# with logM<10: %d"%Nm10_cl45
    print "# with Class5,logM<10:",alldata[(alldata['logM']<10)&(mask_cl5)].size
    print "# with logM<9.7: %d"%Nm97_cl45
    print "# with Class5,logM<9.7:",alldata[(alldata['logM']<9.7)&(mask_cl5)].size
    print "# with logM<9.5: %d"%Nm95_cl45
    print "# with Class5,logM<9.5:",alldata[(alldata['logM']<9.5)&(mask_cl5)].size

    print "\nfor Class 4 & 5 only (%d objects):"%alldata[mask_cl45].size
    print "total frac. class 5: %g"%(1.0*alldata[mask_cl5].size/alldata[mask_cl45].size)
    print "frac. logM<10, class 5: %g"%(1.0*alldata[(alldata['logM']<10)&(mask_cl5)].size/Nm10_cl45)
    print "frac. logM<9.7, class 5: %g"%(1.0*alldata[(alldata['logM']<9.7)&(mask_cl5)].size/Nm97_cl45)
    print "frac. logM<9.5, class 5: %g"%(1.0*alldata[(alldata['logM']<9.5)&(mask_cl5)].size/Nm95_cl45)

    print "\nsources with Koss vs B13 disagreement in AGN type:"
    print "# w/ Koss noBL, B13 Class 4: %d"%alldata[(mask_cl4)&(mask_noBL)].size
    print "# w/ Koss BL, B13 Class 5: %d"%alldata[(mask_cl5)&(mask_BL)].size
    print "# w/ Koss noBL, B13 Class 4, logM<10: %d"%alldata[(mask_cl4)&(mask_noBL)&(alldata['logM']<10)].size
    print "  Names: ", alldata['Name'][(mask_cl4)&(mask_noBL)&(alldata['logM']<10)]
    print "  Type_kpno: ", alldata['Type_kpno'][(mask_cl4)&(mask_noBL)&(alldata['logM']<10)]
    print "  Type_sdss: ", alldata['Type_sdss'][(mask_cl4)&(mask_noBL)&(alldata['logM']<10)]
    print "# w/ Koss BL, B13 Class 5, logM<10: %d"%alldata[(mask_cl5)&(mask_BL)&(alldata['logM']<10)].size
    print "  Names: ", alldata['Name'][(mask_cl5)&(mask_BL)&(alldata['logM']<10)]
    print "  Type_kpno: ", alldata['Type_kpno'][(mask_cl5)&(mask_BL)&(alldata['logM']<10)]
    print "  Type_sdss: ", alldata['Type_sdss'][(mask_cl5)&(mask_BL)&(alldata['logM']<10)]
    print "\n"

    f.close()
    return alldata
    #return

def nathan_masses(path="/n/home00/lblecha/dwarf_agn_data/", h=0.7,
                  return_raw_cat_only=False, write_file=False):

    nmassfile='%s/Nathan_masses/low_mass.fits'%path
    
    ### column names: 
    # 'BAT_index' 'BAT_Name' 'CTPT_Name' 'Other_Name' 'CTPT_RA' 'CTPT_Dec'
    # 'NED_Type' 'NED_redshift' 'SNR' 'FLUX' 'FLUX_LO' 'FLUX_HI' 'CONTA'
    # 'Redshift_Ind_Dist' 'Best_z' 'Lum_dist' 'Best_dist' '70-mo_L14-195'
    # 'L_BAT_bestz' 'log_L_bol_k8' 'log_L_bol_high' 'log_L_bol_low' 'j_m_ext'
    # 'j_msig_ext' 'h_m_ext' 'h_msig_ext' 'k_m_ext' 'k_msig_ext' 'w1' 'w1e' 'w2'
    # 'w2e' 'w3' 'w3e' 'w4' 'w4e' 'uid' 'coeff' 'EBV' 'chi2' 'MASS' 

    pyfits.info(nmassfile)
    f=pyfits.open(nmassfile)
    #print f[1].header
    colnames = np.array(f[1].columns.names)

    ### if we ever need these coeffs, need to figure out a different solution
    ### (like modifying the fits file directly). 'coeff' column contains len-4 lists,
    ### get errors when converting FITS_record to np.recarray.
    ix=np.where(colnames!='coeff')[0]
    colnames = colnames[ix]
    cols = [f[1].data.field(col) for col in colnames]
    f.close()
    d = np.rec.fromarrays(cols,names=list(colnames))
    print "processing data from %s..."%nmassfile
    print "column names:",colnames
    ## convert z values from string to float:
    dty = d.dtype.descr
    ix=np.where(colnames=='NED_redshift')[0][0]
    dty[ix] = (dty[ix][0],'>f8')
    ix=np.where(colnames=='Best_z')[0]
    dty[ix] = (dty[ix][0],'>f8')
    d = d.astype(dty)

    if return_raw_cat_only: return d

    print "%d table entries found in %s"%(d.size,nmassfile)

    ## get rid of spurious masses:
    mask_badmass = (d['MASS']>1.0e3)
    d = d[mask_badmass]
    print "%d entries after spurious masses removed."%(d.size)

    ## convert mass from Msun/h^2 to log Msun:
    d['MASS'] = np.log10(d['MASS'] / (h*h))
    print d['MASS']

    ix10 = (d['MASS']<10)
    ix97 = (d['MASS']<9.7)
    ix95 = (d['MASS']<9.5)

    Nm10raw = d[ix10].size
    Nm97raw = d[ix97].size
    Nm95raw = d[ix95].size
    print "\n# with logM<10:",Nm10raw
    print d['CTPT_Name'][ix10]
    print "\n# with logM<9.7:",Nm97raw
    print d['CTPT_Name'][ix97]
    print "\n# with logM<9.5:",Nm95raw
    print d['CTPT_Name'][ix95]

    return d



def koss_sdss(path="/n/home00/lblecha/dwarf_agn_data/",
              fix_errors_in_main_cat=False):

    tablefile='%s/Koss_2011_table2_SDSS.fits'%path
    
    f=pyfits.open(tablefile)
    d = f[1].data
    colnames = np.array(d.columns.names)
    print "column names:",colnames

    #['Name' 'Date' 'Type' 'z' 'Dist' 'E_B-V_' 'Air' 'PSF' 'BAT' 'SimbadName'
    # 'NED' '_RA' '_DE']

    dall = kosscat(return_raw_cat_only=True)
    all_mask,sdss_mask = matchpos(dall['RA'],dall['DE'],dall['z_sdss'],
                                   d['_RA'],d['_DE'],d['z'],
                                   max_sep_arcsec=60,ignore_z=True)

    print "Sources in sdss catalog unmatched in whole catalog:"
    print d['Name'][~sdss_mask]

    ##these are the main-catalog names:
    main_cat_names = ['MCG+10-17-061','Mrk 463E','NGC 0835','UGC 03995A','MRK 1044','Mrk 1210']
    ##corresponding sdss names:
    ##(phoenix galaxy also has wrong ra/dec in sdss cat)
    sdss_badnames = ['UGC 06732','Mrk 463','ARP 318','UGC 03995','Mrk 1044','Phoenix Galaxy']
    print [name for name in d['Name'] if name in sdss_badnames]

    if fix_errors_in_main_cat:
        
        for i,name in enumerate(sdss_badnames):
            ix_all = dall['Name']==main_cat_names[i]
            ix_sdss = d['Name']==name
            if (len(dall[ix_all])==1) and (len(d[ix_sdss])==1):                
                dall['Name_sdss'][ix_all] = main_cat_names[i]
                dall['Date_sdss'][ix_all] = d['Date'][ix_sdss]
                dall['Type_sdss'][ix_all] = d['Type'][ix_sdss]
                dall['z_sdss'][ix_all] = d['z'][ix_sdss]
                dall['Dist_sdss'][ix_all] = d['Dist'][ix_sdss]
                dall['E_B-V_sdss'][ix_all] = d['E_B-V_'][ix_sdss]
                dall['Air_sdss'][ix_all] = d['Air'][ix_sdss]
                dall['PSF_sdss'][ix_all] = d['PSF'][ix_sdss]
                dall['BAT_sdss'][ix_all] = d['BAT'][ix_sdss]
                dall['SimbadName_sdss'][ix_all] = d['SimbadName'][ix_sdss]
                dall['NED_sdss'][ix_all] = d['NED'][ix_sdss]
                if name != 'Phoenix Galaxy':
                    dall['RA_sdss'][ix_all] = d['_RA'][ix_sdss]
                    dall['DE_sdss'][ix_all] = d['_DE'][ix_sdss]

                print "updated values to write to main koss cat:"
                print dall[ix_all]

            else:
                print "error matching names. main_cat_namse[i],sdss_badnames[i]: ",main_cat_names[i],name
                print ix_all,ix_sdss
                return -1

        catfile='%s/Koss_2011_all.fits'%path
        f=pyfits.open(catfile,mode='update')
        pyfits.update(catfile,dall,ext=1)
        f.flush()
        f.close()
        
def koss_dwarf_xmatch(path="/n/home00/lblecha/dwarf_agn_data/",
                      fname='dwarf_sample_with_b13_extras.dat',
                      use_bass_data=False,bass_fname='dwarf_BASS_data.dat',
                      use_cat_z=True,use_BLflag=False,include_b13_class2=True,
                      logM_cut=10,PSr_cut=100,load_nathan_masses=True):

    kb = kosscat(batcat_id='b13',use_BLflag=use_BLflag,
                 include_b13_class2=include_b13_class2,PSr_cut=PSr_cut)
    #koss_z = kb['z'] if use_cat_z else kb['z_koss'] 
    koss_z = kb['z_koss'] 
    print "Loaded data for %d sources from koss catalog."%kb.size

    dw = pds.read_data_alt(path=path,fname=fname)
    print "Loaded data for %d sources from my dwarf catalog."%dw.size

    koss_idxd,dwarf_idxk = matchpos(kb['RA'],kb['DE'],koss_z,
                                    dw['ra'],dw['dec'],dw['z'],
                                    max_sep_arcsec=600,ignore_z=False,
                                    multi_match_fatal=False)

    mask = np.zeros(dw.size).astype('bool')
    mask[(dwarf_idxk!=-1)&(dw['mstarK11']<logM_cut)] = True
    #dwarf_data = [tup[dwarf_idxk!=-1] for tup in tmpdata]
    dw_orig = copy(dw)
    dw = dw[mask]

    print "\nRetained %d matched sources in dwarf sample (with logM<%g)."%(dw.size,logM_cut)

    pct_cut=10
    print "\n%g sources with pctAGNrband>%g:"%(dw[dw['pctAGNrband']>pct_cut].size,pct_cut)
    for i,pct in enumerate(dw['pctAGNrband']): 
        if pct>pct_cut:
            print "%s %g %g %g"%(dw['name'][i],dw['mstarK11'][i],
                                 dw['Lx(14-195)'][i],pct) 
    
    print "\nsources with Chandra data:"
    for i,flag in enumerate(dw['hasChandra']): 
        if flag.find('Y')!=-1:
            print "%s %g %g %g"%(dw['name'][i],dw['mstarK11'][i],
                                 dw['Lx(14-195)'][i],dw['pctAGNrband'][i]) 
    print "\nsources with HST data:"
    print [dw['name'][i] for i,flag in enumerate(dw['hasHST']) 
           if flag.find('Y')!=-1]
    print "\nvel disp:"
    print dw['vdisp_low']
    print dw['vdisp_hi']
    
    if load_nathan_masses:
        
        nd = nathan_masses()
        
        nathan_idxd,dwarf_idxn = matchpos(nd['CTPT_RA'],nd['CTPT_Dec'],nd['Best_z'],
                                          dw_orig['ra'],dw_orig['dec'],dw_orig['z'],
                                          max_sep_arcsec=600,ignore_z=False,
                                          multi_match_fatal=True)
        maskd = np.zeros(dw_orig.size).astype('bool')
        maskd[(dwarf_idxn!=-1)] = True
        masknd = np.zeros(nd.size).astype('bool')
        masknd[(nathan_idxd!=-1)] = True

        print "\n%d sources in orig dwarf sample."%(dw_orig.size)
        print "%d orig dwarf sample srcs have K11 mass<%g"%(dw_orig[dw_orig['mstarK11']<10].size,logM_cut)
        print "\n%d sources in koss-matched dwarf sample."%(dw.size)
        print "%d koss-matched dwarf sample srcs have K11 mass<%g"%(dw[dw['mstarK11']<10].size,logM_cut)
        print "%d orig dwarf sample sources matched with nathan_masses."%(dw_orig[maskd].size)
        print "%d matched objects have Nathan mass<%g."%(nd[masknd][nd['MASS'][masknd]<logM_cut].size,logM_cut)
        print "%d matched objects have K11 mass<%g."%(dw_orig[maskd][dw_orig['mstarK11'][maskd]<logM_cut].size,logM_cut)
        print "%d unmatched Nathan objects have Nathan mass<%g."%(nd[~masknd][nd['MASS'][~masknd]<logM_cut].size,logM_cut)
        print "%d unmatched Koss objects have K11 mass<%g."%(dw_orig[~maskd][dw_orig['mstarK11'][~maskd]<logM_cut].size,logM_cut)
        #print kb['Name_sdss'][mk]
        #print dw_orig['mstarK11']
        #print dw_orig['name']
        #print nd['CTPT_Name'][(nathan_idxd!=-1)]
        #print nd['MASS'][(nathan_idxd!=-1)]
        #print kb['Name_sdss'][(koss_idxd!=-1)]

        nathan_idxk,koss_idxn = matchpos(nd['CTPT_RA'],nd['CTPT_Dec'],nd['Best_z'],
                                         kb['RA'],kb['DE'],koss_z,
                                         max_sep_arcsec=600,ignore_z=False,
                                         multi_match_fatal=True)
        mk=np.zeros(kb.size).astype('bool')
        mk[(koss_idxn!=-1)]=True
        mn=np.zeros(nd.size).astype('bool')
        mn[(nathan_idxk!=-1)]=True
        print "%d matched objects have Nathan mass < 10."%nd[mn][nd['MASS'][mn]<10].size
        print "%d matched objects have K11 mass < 10."%kb[mk][kb['logM'][mk]<10].size
        print "%d unmatched Nathan objects have Nathan mass < 10."%nd[~mn][nd['MASS'][~mn]<10].size
        print "%d unmatched Koss objects have K11 mass < 10."%kb[~mk][kb['logM'][~mk]<10].size
        #print kb['Name_sdss'][mk]
        #print kb['z_sdss'][mk]
        #print kb['Name_kpno'][mk]
        #print kb['z_kpno'][mk]


    if use_bass_data:
        bd = pds.read_bass_data(path=path,fname=bass_fname)
        print "Loaded BASS data for %d sources from my dwarf catalog."%bd.size

        #BATindex BATname Name ra dec MBH_Hbeta MBH_vdisp MBH_vdisp_err FWHM_Halpha FWHM_Halpha_err Halpha_broad Halpha_broad_err Halpha Halpha_err Hbeta_broad Hbeta_broad_err corr_Hbeta Hbeta_err corr_OIII OIII_err corr_NII NII_err
        NII_Ha = np.repeat(np.nan,bd.size)
        OIII_Hb = np.repeat(np.nan,bd.size)
        ix = np.where((bd['corr_NII']>0)&(bd['Halpha']>0))[0]
        #print len(ix)
        NII_Ha[ix] = np.log10(bd['corr_NII'][ix]/bd['Halpha'][ix])
        ix = np.where((bd['corr_OIII']>0)&(bd['corr_Hbeta']>0))[0]
        #print len(ix)
        OIII_Hb[ix] = np.log10(bd['corr_OIII'][ix]/bd['corr_Hbeta'][ix])
        print "\nOf %d sources,"%bd.size
        print "%d are missing [NII] and/or Halpha."%(NII_Ha[NII_Ha!=NII_Ha].size)
        print "%d are missing [OIII] and/or Hbeta."%(OIII_Hb[OIII_Hb!=OIII_Hb].size)
        print "%d are missing both line ratios."%(OIII_Hb[(OIII_Hb!=OIII_Hb)&
                                                          (NII_Ha!=NII_Ha)].size)

        ixnb = np.in1d(bd['BATindex'],nd['BAT_index'])
        #print bd['BATindex'][ixnb]
        print "%d objects have BASS data and Nathan masses."%bd['BATindex'][ixnb].size
        ixbn = np.in1d(nd['BAT_index'],bd['BATindex'])
        #print nd['BAT_index'][ixbn]
        #print nd['MASS'][ixbn]
        ixnb_mcut = np.in1d(bd['BATindex'],nd['BAT_index'][nd['MASS']<10])
        print "%d objects with nathan M<10 have BASS data & Nathan masses."%bd['BATindex'][ixnb_mcut].size
        
        dw_idxb,bd_idxd = matchpos(dw['ra'],dw['dec'],[0],
                                   bd['ra'],bd['dec'],[0],
                                   max_sep_arcsec=60,ignore_z=True,
                                   multi_match_fatal=True)

        #bd_mask = np.zeros(dw.size).astype('bool')
        #bd_mask[dw_idxb!=-1] = True
        #print "\n%d dwarf sources (with logM<%g) matched with BASS data."%(dw[bd_mask].size,logM_cut)
        bd_mask = np.zeros(bd.size).astype('bool')
        bd_mask[bd_idxd!=-1] = True
        print "\n%d dwarf sources (with logM<%g) matched with BASS data."%(bd[bd_mask].size,logM_cut)

        ## Kewley et al. 2006 BPT
        x1 = np.arange(-1.28,0.0,0.05)
        y1 = 0.61/(x1-0.05)+1.3
        x2 = np.arange(-2.2,0.4,0.05)
        y2 = 0.61/(x2-0.47)+1.19

        fig = plt.figure(figsize=(5,5))
        plt.clf()
        ax = fig.add_subplot(111)
        ax.set_xlim(-2.2,0.8)
        ax.set_ylim(-1.2,1.5)
        ax.set_xlabel("log([NII]/Halpha)")
        ax.set_ylabel("log([OIII]/Hbeta)")
        allp,=ax.plot(NII_Ha,OIII_Hb,'o',mec='b',fillstyle='none')
        lowM,=ax.plot(NII_Ha[bd_mask],OIII_Hb[bd_mask],'o',mec='b',color='b')
        allnb,=ax.plot(NII_Ha[ixnb],OIII_Hb[ixnb],'^',mec='m',fillstyle='none',ms=4)
        nb_lowM,=ax.plot(NII_Ha[ixnb_mcut],OIII_Hb[ixnb_mcut],'^',mec='m',color='m',ms=4)
        ax.plot(x1,y1,'--')
        ax.plot(x2,y2)
        ax.legend((allp,lowM,allnb,nb_lowM),('parent sample','Koss M<10','Nathan M/h<10','Nathan M<10'),
                  loc='lower left',fontsize=9)
        fig.savefig(path+'/bass_BPT_quickplot.pdf')
        
        bd = bd[bd_mask]
        NII_Ha = NII_Ha[bd_mask]
        OIII_Hb = OIII_Hb[bd_mask]

        print "\nOf %d matched sources,"%bd.size
        print "%d are missing [NII] and/or Halpha."%(NII_Ha[NII_Ha!=NII_Ha].size)
        print "%d are missing [OIII] and/or Hbeta."%(OIII_Hb[OIII_Hb!=OIII_Hb].size)
        print "%d are missing both line ratios."%(OIII_Hb[(OIII_Hb!=OIII_Hb)&
                                                          (NII_Ha!=NII_Ha)].size)

        
        ## dw colnames: catID name altname ra dec Lx(14-195) Lbol Bmag AGNtype z optspec bhmass_low bhmass_hi mstarK11 mstarother_low mstarother_hi vdisp_low vdisp_low_err vdisp_hi vdisp_hi_err pctAGNrband dist gmagcorr rmagcorr gminusr hasHST hasChandra Fx SNR Gamma
        ## not done, not working yet:
        #ix_hbeta = np.where(bd[bd_mask]['MBH_Hbeta']>0)[0]
        #ix_vdisp = np.where(bd[bd_mask]['MBH_vdisp']>0)[0]
        #ix_dw_hbeta = np.where(dw
        print bd['MBH_Hbeta'].size
        print dw[dw_idxb!=-1]['Lbol'].size
        fEdd_Hbeta = 1.26e38*10**(bd['MBH_Hbeta'])/10.0**(dw[dw_idxb!=-1]['Lbol'])
        fEdd_vdisp = 1.26e38*10**(bd['MBH_Hbeta'])/10.0**(dw[dw_idxb!=-1]['Lbol'])

        fig = plt.figure(figsize=(5,8))
        plt.clf()
        ax1 = fig.add_subplot(121)
        ax1.plot(bd['MBH_Hbeta'],dw[dw_idxb!=-1]['mstarK11'],'ro')
        ax1.plot(bd['MBH_vdisp'],dw[dw_idxb!=-1]['mstarK11'],'ko')
        ax2 = fig.add_subplot(122)
        ax2.plot(bd['MBH_Hbeta'],fEdd_Hbeta,'ro')
        ax2.plot(bd['MBH_vdisp'],fEdd_vdisp,'ko')
        fig.subplots_adjust()
        fig.savefig(path+'/bass_mbhmstar_quickplot.pdf')
        

def nathan_b13_xmatch(path="/n/home00/lblecha/dwarf_agn_data/",
                      max_sep_arcsec=600):

    catfile = '%s/Baumgartner_2013.fits'%path
    lstr='logL'
    ffac=1.0e-12

    print "\nOpening catalog file: %s\n"%catfile
    f=pyfits.open(catfile)
    bd=f[1].data
    colnames = np.array(bd.columns.names)
    print "column names:",colnames
    f.close()

    ### Note: column 'BAT_index' in Nathan's file corresponds to the array index in the B13 catalog,
    ### but indexed from 1 instead of 0. i.e., array index 0 in b13 has 'BAT_index'=1.
    nd = nathan_masses()

    print nd['MASS'][nd['MASS']<10]
    



