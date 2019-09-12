import pyfits
import numpy as np


def update(path='/n/home00/lblecha/dwarf_agn_data/'):

    catfile='%s/Koss_2011_all.fits'%path
    
    f=pyfits.open(catfile,mode='update')
    #f=pyfits.open(catfile)
    d = f[1].data
    hdr = f[1].header
    colnames = np.array(d.columns.names)
    print "column names:",colnames
    

    newcolnames = np.array(['Name', 'M', 'umag', 'gmag', 'rmag', 'imag',
                            'zmag', 'logM', 'PSr', 'E_g-r', 'E_g-r_i', 
                            'ucontam', 'contam', 'BAT', 'SimbadName', 
                            'NED', 'RA', 'DE', 'Name_kpno', 'Date_kpno', 
                            'Type_kpno', 'z_kpno', 'Dist_kpno', 'E_B-V_kpno', 
                            'Air_kpno', 'PSF_kpno', 'Name_sdss', 'Date_sdss', 
                            'Type_sdss', 'z_sdss', 'Dist_sdss', 'E_B-V_sdss', 
                            'Air_sdss', 'PSF_sdss', 'BAT_sdss', 
                            'SimbadName_sdss', 'NED_sdss', 'RA_sdss', 
                            'DE_sdss', 'Name_morph', 'Rp_morph', 'C_morph', 
                            'Class_morph', 'b_a_morph', 'BLflag_morph'],
                           dtype='|S15')

    for i in np.arange(len(colnames)):
        print "%d %s TTYPE%d %s"%(i,colnames[i],(i+1),newcolnames[i])
        print hdr['TTYPE%d'%(i+1)]
        hdr.set('TTYPE%d'%(i+1),newcolnames[i])
        print hdr['TTYPE%d'%(i+1)]
        

    print "Finished updating header keys."
    f.close()
    ##pyfits.update(catfile,d,hdr,ext=1)
    ##f.flush()
    ##print f[1].header
    ##f.writeto('%s/Koss_2011_all_new.fits'%path)
    ##f.close()
