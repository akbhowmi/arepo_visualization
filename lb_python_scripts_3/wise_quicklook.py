import pyfits, glob, re
import numpy as np
import broadband as bb


msun = 1.989e33
yr = 3.15e7
c = 3.0e10
radeff=0.1

## Jarrett et al. 2011
W1_Vega_to_AB = 2.699
W2_Vega_to_AB = 3.339
W3_Vega_to_AB = 5.174
W4_Vega_to_AB = 6.620

def colors(path='/oasis/projects/nsf/hvd115/lblecha',
           subdir='q0.5_fg0.3_sunruns/test_fiducial',mycalc=False):
    #path='/oasis/projects/nsf/hvd115/lblecha/test_sunruns/test_whole_SED_lowres'
    #dir_arr = ['snap_123','snap_441','snap_646']
    #bb_arr = ['broadband_123.fits','broadband_0441.fits','broadband_646.fits']
    #mc_arr = ['mcrx_123.fits','mcrx_0441.fits','mcrx_646.fits']
    #tmrg = 1.56
    #time_arr = [1.221,1.539,1.724]
    #sfr_arr = [0.4,8.5,0.4]
    #mdot_arr = [0.0017,0.045,0.018]

    mcfiles = glob.glob(path+'/'+subdir+'/mcrx*fits')   
    snaps = np.sort(np.array([np.int(re.search('(?<=mcrx_)\d+',f).group(0))
                              for f in mcfiles]))

    if len(snaps)==0:
        print "no snapshots found in %s/%s."%(path,subdir)
        return -1

    print "\n%s/%s:"%(path,subdir)
    print "using the following snapshots:",snaps

    for snap in snaps:

        fbb = pyfits.open('%s/%s/broadband_%s.fits'%(path,subdir,str(snap).zfill(3)))
        
        w1data = fbb[17].data[110]
        lambda_eff_w1 = w1data[2]
        w2data = fbb[17].data[111]
        lambda_eff_w2 = w2data[2]

        #using 'scatter' values only
        w1lum = np.array(w1data[4:][::2][:7])*lambda_eff_w1
        w1magAB = np.array(w1data[5:][::2][:7])
        w2lum = np.array(w2data[4:][::2][:7])*lambda_eff_w2
        w2magAB = np.array(w2data[5:][::2][:7])
        
        w1magVega = w1magAB - W1_Vega_to_AB
        w2magVega = w2magAB - W2_Vega_to_AB
        
        w1minusw2_Vega = w1magVega - w2magVega

        fbb.close()

        fmc = pyfits.open('%s/%s/mcrx_%s.fits'%(path,subdir,str(snap).zfill(3)))
        
        lbol_tot = fmc[48].header['L_bol_grid']
        lbol_cam = np.array([fmc[48].header['L_SCAT%d'%i] for i in range(7)])
            
        if mycalc:
            lam = fmc[5].data.field('lambda')
            llam_cam = np.array([fmc[48].data.field('L_lambda_out%d'%i) 
                                 for i in range(7)])
            my_w1magAB = np.array([ bb.ABMag('WISE-W1',lam,llam_cam[i],sed_lam_units='m')
                                    for i in range(7) ])
            my_w2magAB = np.array([ bb.ABMag('WISE-W2',lam,llam_cam[i],sed_lam_units='m')
                                    for i in range(7) ])

            my_w1magVega = my_w1magAB - W1_Vega_to_AB
            my_w2magVega = my_w2magAB - W2_Vega_to_AB
            
            my_w1minusw2_Vega = my_w1magVega - my_w2magVega
            
        fmc.close()

        print "\nsnap %d:"%snap
        print "[value avg. over viewing angle (min, max)]"
        print "L_W1 [W] = %.3g (%.3g, %.3g)"%(w1lum.mean(),w1lum.min(),w1lum.max())
        print "L_W2 [W] = %.3g (%.3g, %.3g)"%(w2lum.mean(),w2lum.min(),w2lum.max())
        print "W1 (AB) = %.3g (%.3g, %.3g)"%(w1magAB.mean(),w1magAB.min(),w1magAB.max())
        print "W2 (AB) = %.3g (%.3g, %.3g)"%(w2magAB.mean(),w2magAB.min(),w2magAB.max())
        print "W1-W2 (Vega) = %.3g (%.3g, %.3g)"%(w1minusw2_Vega.mean(),w1minusw2_Vega.min(),
                                                  w1minusw2_Vega.max())
        if mycalc:
            print "my W1 (AB) = %.3g (%.3g, %.3g)"%(my_w1magAB.mean(),my_w1magAB.min(),my_w1magAB.max())
            print "my W2 (AB) = %.3g (%.3g, %.3g)"%(my_w2magAB.mean(),my_w2magAB.min(),my_w2magAB.max())
            print "my W1-W2 (Vega) = %.3g (%.3g, %.3g)"%(my_w1minusw2_Vega.mean(),my_w1minusw2_Vega.min(),
                                                         my_w1minusw2_Vega.max())
        print "Lbol,cam [W] = %.3g (%.3g, %.3g)"%(lbol_cam.mean(),lbol_cam.min(),lbol_cam.max())
        print "Lbol,tot [W] = %.4g"%lbol_tot

    return 0
