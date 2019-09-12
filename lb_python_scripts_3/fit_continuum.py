from numpy import *


def fit_continuum(vlambda, vflux, nloops, sigma_parameter, window_width):

    #print "start of fit_continuum: vflux min/max:",min(vflux),max(vflux)

    vcontinuum_level = array(vflux)
    vlam = array(vlambda)
    nlam = vlam.shape[0]
    vfinal_sigma = zeros((nlam,nloops),dtype=float64)

    for i in range(nloops):

        for j in range(nlam):
            
            jmin = (value_locate(vlam,vlam[j]-window_width))[0]
            jmax = (value_locate(vlam,vlam[j]+window_width))[0]

            #print "in fit: i,j=",i,j
            #print "window = ",vlam[j]-window_width,vlam[j]+window_width
            #print "jmin,jmax=",jmin,jmax

            if jmin == -1: jmin = 0
            if jmax == nlam: jmax -= 1
            
            vflux_cut = vcontinuum_level[jmin : jmax]
            cut_sigma = sqrt(var(vflux_cut))
            cut_median = median(vflux_cut)        

            #print "jmin,jmax=",jmin,jmax
            #print "cut_sigma, cut_median=",cut_sigma,cut_median
            
            if (vcontinuum_level[j] - cut_median) > 2.0*sigma_parameter*cut_sigma:
                vcontinuum_level[j] = cut_median
                #print "reducing bin %d to cut_median=%f"%(j,cut_median)
            if (vcontinuum_level[j] - cut_median) <= -1.0*sigma_parameter*cut_sigma:
                vcontinuum_level[j] = cut_median
                #print "increasing bin %d to cut_median=%f"%(j,cut_median)
                    
            vfinal_sigma[j, i] = cut_sigma

    #print "end of fit_continuum: vflux min/max:",min(vflux),max(vflux)


    return vcontinuum_level, vfinal_sigma


def value_locate(refx,x):
    """
    Python version of IDL's useful VALUE_LOCATE procedure.
    """
    
    refx=array(refx)
    x=atleast_1d(x)
    loc=zeros(len(x),dtype='int')
    
    for i in xrange(len(x)):
        ix=x[i]
        ind=((refx-ix) <= 0).nonzero()[0]
        if len(ind) == 0: loc[i]=-1
        else: loc[i]=ind[-1]

    return loc
