import numpy as np
import pyfits,sys
from scipy import interpolate

abmag_zp = 3.631e-23 ## W/m^2/Hz 
clight = 3e8 ## everything in m

def convert_lam_units(lam,unit):
    if unit == 'm':
        pass
    elif unit=='angstrom':
        lam = lam * 1.0e-10
    elif unit=='micron':
        lam = lam * 1.0e-6
    elif unit=='cm':
        lam = lam * 0.01
    else: raise Exception("Invalid lambda units specified: %s.\nPlease enter 'angstrom', 'micron', 'cm', or 'm'."%sed_lam_units)

    return lam

def convert_lum_units(lum,unit):
    if unit == 'W':
        pass
    elif unit == 'ergs':
        lum = lum * 1.0e-7
    elif unit == 'Lsun':
        lum = lum * 3.83e26
    else: raise Exception("Invalid units specified: %s.\nPlease enter 'W', 'ergs', or 'Lsun'."%unit)

    return lum

def slice_lam_Llam(lam,Llam,filter_lam):
    ix_cut = np.where((lam>=filter_lam.min())&
                      (lam<=filter_lam.max()))[0]

    return lam[ix_cut],Llam[ix_cut]


def Leff(filter_name,sed_lambda,sed_Llambda,
         sed_lam_units='m',sed_lum_units='W',
         filter_path='/home/lblecha/sundata/filters',return_ABMag=False):

    # lambda in m
    sed_lambda = convert_lam_units(sed_lambda, sed_lam_units)

    if sed_Llambda.max() != sed_Llambda.max():
        print "Warning: before conversion, sed_Llambda.max()=%g"%sed_Llambda.max()
        #print ['%g %g'%(sed_lambda[i],sed_Llambda[i]) for i in range(len(sed_lambda))]
        #sys.exit()

    ## luminosity in W
    sed_Llambda = convert_lum_units(sed_Llambda, sed_lum_units)

    with open('%s/%s.res'%(filter_path,filter_name),'r') as fpfilter:
        filter_lambda,filter_resp_from_file=np.loadtxt(fpfilter,unpack=True)
        filter_lambda = filter_lambda *1.0e-10 ## convert Angstrom to m        

    sed_lambda,sed_Llambda = slice_lam_Llam(sed_lambda,sed_Llambda,filter_lambda)

    filter_resp_interp = interpolate.interp1d(filter_lambda, filter_resp_from_file)
    filter_resp = filter_resp_interp(sed_lambda)
    sed_in_band = filter_resp * sed_Llambda

    band_Leff = np.trapz(sed_in_band,sed_lambda)
    #if band_Leff != band_Leff:
    #    print "Error: band_Leff=%g"%band_Leff
    #    #print "filter_lambda:",filter_lambda
    #    #print "filter_resp_from_file:",filter_resp_from_file
    #    print "sed_lambda:",sed_lambda
    #    print "sed_in_band:",sed_in_band
    #    print "sed_Llambda:",sed_Llambda
    #    print "filter_resp:",filter_resp
    #    #sys.exit()

    if return_ABMag:
        ewidth_lambda = np.trapz(filter_resp, sed_lambda)
        #band_Llambda_eff = band_Leff/ewidth_lambda

        lambda_eff = np.trapz(sed_lambda*filter_resp,sed_lambda)/ewidth_lambda
        ewidth_nu = clight*np.trapz(filter_resp/(sed_lambda**2),sed_lambda)

        ## calculate flux at 10pc, converting pc to m. 
        ## 5 + 85.19 = -2.5*log10((1/10**2) * (pc/m)**2/(4pi))
        band_Mag = -2.5 * np.log10(band_Leff/(ewidth_nu*abmag_zp)) + 5 + 85.19
        return band_Leff,band_Mag
    else:
        return band_Leff


def ABMag(filter_name,sed_lambda,sed_Llambda,
          sed_lam_units='micron',sed_lum_units='W',
          filter_path='/home/lblecha/sundata/filters'):

    # lambda in m
    sed_lambda = convert_lam_units(sed_lambda, sed_lam_units)

    ## luminosity in W
    sed_Llambda = convert_lum_units(sed_Llambda, sed_lum_units)

    with open('%s/%s.res'%(filter_path,filter_name),'r') as fpfilter:
        filter_lambda,filter_resp_from_file=np.loadtxt(fpfilter,unpack=True)
        filter_lambda = filter_lambda *1.0e-10 ## convert Angstrom to m        

    sed_lambda,sed_Llambda = slice_lam_Llam(sed_lambda,sed_Llambda,filter_lambda)

    filter_resp_interp = interpolate.interp1d(filter_lambda, filter_resp_from_file)
    filter_resp = filter_resp_interp(sed_lambda)

    band_Leff = Leff(filter_name,sed_lambda,sed_Llambda,sed_lam_units=sed_lam_units,
                     sed_lum_units=sed_lum_units,filter_path=filter_path)

    ewidth_lambda = np.trapz(filter_resp, sed_lambda)
    #band_Llambda_eff = band_Leff/ewidth_lambda

    lambda_eff = np.trapz(sed_lambda*filter_resp,sed_lambda)/ewidth_lambda
    ewidth_nu = clight*np.trapz(filter_resp/(sed_lambda**2),sed_lambda)

    ## calculate flux at 10pc, converting pc to m. 
    ## 5 + 85.19 = -2.5*log10((1/10**2) * (pc/m)**2/(4pi))
    band_Mag = -2.5 * np.log10(band_Leff/(ewidth_nu*abmag_zp)) + 5 + 85.19

    #print "%s: lambda_eff=%g, ewidth_lambda=%g, ewidth_nu=%g"%(filter_name,lambda_eff,ewidth_lambda,ewidth_nu)
    #print "    band_Leff=%g, band_Mag=%g"%(band_Leff,band_Mag)

    return band_Mag
