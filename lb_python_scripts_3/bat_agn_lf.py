import numpy as np
import matplotlib as mpl
#mpl.use('agg')
import matplotlib.pyplot as plt
from copy import copy
import sys, pyfits, timeit
import PyCos as pc
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
from astropy import constants as const
import scipy.optimize as optimize
from scipy import stats, interpolate
import astro_constants as ac
import lb_utils as lb
import process_dwarf_sample as pds
import match_catalog as matchcat

def unit_comoving_vol(z, cosmo):

    D_H = const.c.to('km/s').value/cosmo.H0.value

    #return  D_H * cosmo.D_M(0.0,z)**2 / np.sqrt(omega_m*(1+z)**3 + omega_l) 
    return D_H * cosmo.comoving_transverse_distance(z).value**2 * cosmo.inv_efunc(z)

def get_zmax(zmax, *data):

    cosmo,fx_lim,lglx,gamma = data
    fac = np.sqrt( 10**lglx / (4*np.pi*fx_lim) ) / (1e6*ac.PC)    
    #print fx_lim, lglx, gamma
    return ( fac * (1+zmax)**(1-0.5*gamma) - 
             cosmo.luminosity_distance(zmax).value )

def print_array_attr(arr,name):
    arrstr = "min/max/mean/med %s: %g %g %g %g"
    print arrstr%(name,arr.min(),arr.max(),
                  arr.mean(),np.median(arr))


def mcrab_to_cgs(emin,emax):
    ## 1 mCrab in the 15-55 keV band = 1.27e-11 erg/cm^2/s
    ## 1 mCrab in the 14-195 keV band = 2.386e-11 erg/cm^2/s (B12 Eqn6)
    ## Crab spctrum: F(E) = 10.17 * E^-2.15 (photons/cm^2/s/keV)

    ## NOTE:
    ## flim(5sig) = 0.99 mCrab (t / 1 Ms)^-0.5 (Tueller+10, based on 22mo data)
    ## flim(5sig) = 1.18 mCrab (t / 1 Ms)^-0.5 (Baumgartner+12, " " 70mo data)


    if emax <= emin or emin < 14:
        print "Invalid energy range: (%g,%g) keV."%(emin,emax)
        return -1

    return 1e-3 * 1.6021773e-9 * (10.17/0.15) * ( emin**-0.15 - emax**-0.15 )

def shift_band(emin,emax,emin_new,emax_new,gamma):
    
    if emin<14 or emin_new<14 or emax<emin or emax_new<emin_new:
        print "invalid energy range(s) entered:"
        print "emin1,emax1=%g,%g; emin2,emax2=%g,%g"%(emin,emax,emin_new,
                                                      emax_new)

    fac2 = ( (np.log(emax_new)-np.log(emin_new))/
             (np.log(emax)-np.log(emin)) )
    fac = np.zeros((gamma.size)) + fac2
    ix = np.where(gamma!=2)[0]    
    if len(ix)>0:
        fac[ix] = ( (emax_new**(2-gamma[ix])-emin_new**(2-gamma[ix]))/
                    (emax**(2-gamma[ix])-emin**(2-gamma[ix])) )

    print_array_attr(fac,'fac')
    return fac


def calc_Vmax(cosmo, lglx, flux_tab, area_tab, fx_lim=5.7e-12, 
              zarr_min=1.0e-5, zarr_max=0.4, dz=1.0e-5,
              gamma=2, snr_limited=False):

    start_time=timeit.default_timer()
    if not snr_limited:
        ## "new" method, uses sky coverage vs flux with a constant 
        ## limiting flux for the survey
        Vmax = np.zeros(len(lglx))
        zmax = copy(Vmax)
        zarr = np.arange(zarr_min,zarr_max,dz)
        interpA = interpolate.interp1d(flux_tab, area_tab, kind='cubic')
        dV_dOmegadz = unit_comoving_vol(zarr,cosmo)
        dLfac_zarr = 4*np.pi*(cosmo.luminosity_distance(zarr).to(u.cm).value)**2 
        print "calculating Vmax for %d sources."%len(lglx)
        for src,l in enumerate(lglx):
            if src % 20 == 0: print "source %d of %d"%(src,len(lglx))
            A = np.zeros(len(zarr))
            mask = copy(A).astype('bool')
            #fx_zarr = 10**l / dLfac_zarr
            fx_zarr = 10**l / dLfac_zarr * (1+zarr)**(2.0-gamma[src])
            #print "l:",l
            #print "min/max fx_zarr=%g,%g."%(fx_zarr.min(),fx_zarr.max())
            #if len(fx_zarr[mask])==0 or len(fx_zarr[mask])==len(fx_zarr):
            if fx_zarr.min()>fx_lim or fx_zarr.max()<fx_lim:
                print "Error: fx_lim outside range of fx_zarr."
                print "lglx = %g, fx_lim=%g"%(l,fx_lim)
                print "min/max fx_zarr=%g,%g."%(fx_zarr.min(),fx_zarr.max())
                return -1
            mask[(fx_zarr>=fx_lim)] = True
            zmax[src] = zarr[mask].max()
            #zmax[src] = zarr[mask][-1]
            #print len(fx_zarr),len(fx_zarr[mask])
            #print zmax[src]
            A[mask] = np.array([interpA(fx) for fx in fx_zarr[mask]])
            Vmax[src] = np.sum(dV_dOmegadz[mask] * A[mask] * dz)

    else:
        ## "old" method using same coverage area for whole survey, but limiting 
        ## flux depends on SNR of each source. Equivalent to above method if a 
        ## constant area is used in that calculation. 

        print "min/max lglx:",lglx.min(),lglx.max() 
        print "min/max fx_lim=",fx_lim.min(),fx_lim.max()
        print "min/max gamma:",gamma.min(),gamma.max()
        zmax = np.array([ optimize.fsolve(get_zmax,np.float64(1.5),
                                          args=(cosmo,fx_lim[i],lglx[i],gamma[i])) 
                          for i in range(lglx.size) ]).reshape(lglx.size)
        print "min/max zmax:",zmax.min(),zmax.max() 
 
        ### note: unlike the old cosmo routines (PyCos), comoving volume
        ### is defined in astropy.cosmology as the total comoving volume,
        ### not the volume per unit area
        Vmax = cosmo.comoving_volume(zmax).value
        print "min/max Vmax:",Vmax.min(),Vmax.max()
    end_time = timeit.default_timer()        
    print "calc_Vmax run time: %g s"%(end_time-start_time)
    #return zmax, Vmax
    return Vmax


def load_cat(cosmo, path="/n/home00/lblecha/dwarf_agn_data/",
             cat='a12',include_b13_class2=False):

    ## BAT AGN catalog notes:
    ## 'b13' denotes Baumgartner et al. 2008; BAT flux is 14-195keV [1e-12 erg/s/cm^2]. Catalog includes *all* BAT sources, so we must filter out the AGN.
    ## 'a12' denotes Ajello et al. 2012; BAT flux is 15-55keV [1e-11 erg/s/cm^2]
    ## 'b11' denotes Burlon et al. 2011; BAT flux is 15-55keV [1e-11 erg/s/cm^2]
    ## 't08' denotes Tueller et al. 2008; BAT flux is 14-195keV [1e-11 erg/s/cm^2]

    fname = {'b13':'Baumgartner_2013.fits', 'a12':'Ajello_2012.fits', 
             'b11':'Burlon_2011.fits', 't08':'Tueller_2008.fits'}
    if cat not in fname:
        print "Error: cat %s is not defined."%cat
        return 0

    snrstr = {'b13':'S_N', 'a12':'S_N', 'b11':'S_N', 't08':'SNR'}
    Fstr = {'b13':'Flux', 'a12':'Flux', 'b11':'Flux', 't08':'fBAT'}
    Ffac = {'b13':1e-12, 'a12':1e-11, 'b11':1e-11, 't08':1e-11}
    f=pyfits.open('%s/%s'%(path,fname[cat]))

    cat_flux=f[1].data.field(Fstr[cat]) * Ffac[cat] ## convert to erg/s/cm^2
    cat_ra=f[1].data.field('RAJ2000')
    cat_dec=f[1].data.field('DEJ2000')
    cat_z=f[1].data.field('z')
    cat_snr=f[1].data.field(snrstr[cat])
    cat_swiftname=f[1].data.field('SWIFT')
    if cat=='b13':
        cat_gamma=f[1].data.field('Gamma')
        cat_nh=np.zeros((cat_z.size)) ## not defined for this sample
        cat_isCT=np.zeros((cat_z.size),dtype=bool)
        ## we will calculate from flux
        #cat_lglx=f[1].data.field('logL').astype('float64') ## log erg/s
        cat_type=f[1].data.field('Cl')
        print "\nRead %d total BAT sources from b13 catalog."%cat_z.size
        ### class 2: 'Galaxy' (extended objects w/ no firm evidence that they host an AGN)
        ### class 4: Sy1 (Sy1.0-1.5)
        ### class 5: Sy2 (Sy1.7-2.0)
        ### class 6: other AGN
        classarr = [2,4,5,6] if include_b13_class2 else [4,5,6]
        ix_agn = [j for j in range(len(cat_type)) if cat_type[j] in classarr]
        print "retained cat_type values: ",np.unique(cat_type[ix_agn])
        print "%d Cl=2, %d Cl=4, %d Cl=5, %d Cl=6."%(cat_type[cat_type==2].size,cat_type[cat_type==4].size,cat_type[cat_type==5].size,cat_type[cat_type==6].size)
        #print "unique 'Type' strings (%d):"%len(np.unique(f[1].data.field('Type')[ix_agn]))
        #print np.unique(f[1].data.field('Type')[ix_agn])        
        data = lb.apply_cuts(ix_agn, (cat_ra,cat_dec,cat_z,cat_flux,cat_snr, 
                                       cat_gamma,cat_nh,cat_isCT,cat_swiftname) )
        cat_ra,cat_dec,cat_z,cat_flux,cat_snr,cat_gamma,cat_nh,cat_isCT,cat_swiftname=data
        ## sources named Jxxx+xxxB,C,D have ditto marks for snr value, which print as nan's. Ditto marks, I tell you! Need to pick up the 'A' snr value instead.
        ix_bad_snr = np.where(cat_snr!=cat_snr)[0]
        while len(ix_bad_snr)>0:
            cat_snr[ix_bad_snr] = cat_snr[ix_bad_snr-1]
            ix_bad_snr = np.where(cat_snr!=cat_snr)[0]
        print "Retained %d BAT AGN from b13 catalog.\n"%cat_z.size
    elif cat=='a12':
        cat_gamma=f[1].data.field('Gamma')
        cat_nh=np.zeros((cat_z.size)) ## not defined for this sample        
        cat_isCT=np.zeros((cat_z.size),dtype=bool)
        cat_isCT[f[1].data.field('A')=='A'] = True
        ## we will calculate from flux
        cat_lglx=f[1].data.field('logLX').astype('float64') ## log erg/s
        cat_type=f[1].data.field('Type')
        print "Read %d total BAT sources from a12 catalog."%cat_z.size
        ix_agn = [j for j in range(len(cat_type)) 
                  if cat_type[j] not in ['BLAZAR','RG']]
        data = lb.apply_cuts(ix_agn, (cat_ra,cat_dec,cat_z,cat_flux,cat_snr, 
                                       cat_gamma,cat_nh,cat_isCT,cat_swiftname,cat_lglx) )
        cat_ra,cat_dec,cat_z,cat_flux,cat_snr,cat_gamma,cat_nh,cat_isCT,cat_swiftname,cat_lglx=data
        print "retained cat_type values: ",np.unique(cat_type[ix_agn])
        print "Retained %d BAT AGN from a12 catalog.\n"%cat_z.size
    elif cat=='b11':
        cat_gamma=f[1].data.field('Gamma')
        cat_nh=f[1].data.field('logNH')
        cat_isCT=np.zeros((cat_z.size),dtype=bool)
        cat_isCT[cat_nh>=24] = True
    elif cat=='t08':
        cat_gamma=2.0+np.zeros((cat_z.size)) ## not defined for this sample
        cat_nh=f[1].data.field('logNH')
        cat_isCT=np.zeros((cat_z.size),dtype=bool)
        cat_isCT[cat_nh>=24] = True
        ## we will calculate from flux
        #cat_lglx=f[1].data.field('logL').astype('float64') ## log erg/s
        cat_type=f[1].data.field('Type')
        print "Read %d total BAT sources from t08 catalog."%cat_z.size
        print "Cutting sources with S/N<4.8, |b|<15deg, or type BL Lac/Blazar."
        ## sample contains some objects with S/N<4.8! 
        ## we're also supposed to eliminate objects at |b|<15deg...
        c = SkyCoord(ra=cat_ra*u.degree, dec=cat_dec*u.degree)
        cat_b = c.galactic.b.degree
        ix_cut = [j for j in range(len(cat_snr)) 
                  if cat_snr[j]>=4.8 and np.abs(cat_b[j])>=15 
                  and cat_type[j] not in ['BL Lac','Blazar']]
        data = lb.apply_cuts(ix_cut, (cat_ra,cat_dec,cat_z,cat_flux,
                                      cat_snr,cat_gamma,cat_nh,
                                      cat_isCT,cat_swiftname) )
        cat_ra,cat_dec,cat_z,cat_flux,cat_snr,cat_gamma,cat_nh,cat_isCT,cat_swiftname=data
        print "Retained %d BAT AGN from t08 catalog.\n"%cat_z.size

    f.close()

    ## check for 'nan' nh and set to 0:
    if len(cat_nh[cat_nh!=cat_nh])>0:
        cat_nh[cat_nh!=cat_nh]=0
    ## eliminate sources w/ no redshift:
    #ix_hasz = np.where((cat_z==cat_z)&(cat_z>0))[0]
    ix_hasz = np.where((cat_z==cat_z)&(cat_z>0))[0]
    tmp_lglx = np.zeros((cat_z.size))
    if len(ix_hasz) > 0: 
        dL = cosmo.luminosity_distance(cat_z[ix_hasz]).to('cm').value
        ## includes k-correction, cf. B11 Eqn 2 [in log erg/s]:
        tmp_lglx[ix_hasz]= np.log10( 4*np.pi * dL**2 * cat_flux[ix_hasz] / 
                                     (1+cat_z[ix_hasz])**(2-cat_gamma[ix_hasz]) )
        if cat=='a12':
            print_array_attr(np.abs((tmp_lglx[ix_hasz]-cat_lglx[ix_hasz])/cat_lglx[ix_hasz]),'diff lglx')
        cat_lglx = copy(tmp_lglx)
        data = lb.apply_cuts(ix_hasz, (cat_ra,cat_dec,cat_z,cat_flux,cat_snr, 
                                       cat_gamma,cat_nh,cat_isCT,cat_swiftname,cat_lglx) )
        return data
    else: return [-1]
                

def lum_func(cosmo, binedge, lglx, Vmax, flux_tab, area_tab, 
             fx_lim=None, gamma=2, use_snr_limited_vmax=False):

    binwidth = binedge[1:]-binedge[:-1] 
    bins = binedge[:-1] + 0.5*binwidth
    Lhist,tmp = np.histogram(lglx, bins=binedge)
    L_Phi,tmp = np.histogram(lglx, bins=binedge, weights=1.0/Vmax)
    L_Phi = L_Phi / binwidth
    Phi_err,tmp = np.histogram(lglx, bins=binedge, weights=1.0/Vmax**2)
    Phi_err = np.sqrt(Phi_err) / binwidth
    print "lx binedges: "
    print binedge
    print "L_Phi: "
    print L_Phi
    print "Phi_err: "
    print Phi_err
    print "Lhist: "
    print Lhist
    print "# w/ lglx<41:",lglx[lglx<41].size
    print lglx[lglx<41]
    print Vmax[lglx<41]
    ix0 = np.where(L_Phi==0)[0]
    if len(ix0)>0:
        tmp = np.zeros(len(ix0))
        print "calculating Vmax for LF upper limits."
        Vmax_ul = calc_Vmax(cosmo, bins[ix0], flux_tab, area_tab,
                            zarr_max=4.0, dz=1.0e-4, fx_lim=fx_lim, 
                            gamma=gamma, snr_limited=use_snr_limited_vmax)
        Phi_err[ix0] = 1.0 / Vmax_ul / binwidth[ix0]

    return Lhist, L_Phi, Phi_err

def cum_lum_func(cosmo, binedge, lglx, Vmax):


    #l_hist,tmp = np.histogram(lglx, bins=binedge)
    l_phi,tmp = np.histogram(lglx, bins=binedge, weights=1.0/Vmax)
    cum_l_phi = np.cumsum(l_phi[::-1])[::-1]
    cum_l_phi[cum_l_phi==0] = 1.0e-20

    #print "in function cum_lum_func():"
    #print "binedge: ",binedge
    #print "l_phi: ",l_phi
    #print "cum_l_phi: ",cum_l_phi

    return cum_l_phi

def cat_lf_fitfunc(cat='a12',lglxmin=39,lglxmax=47,dlglx=0.25,
                   logvals=True,include_CT=True):

    n=(lglxmax-lglxmin)/.25
    lglx = np.arange(n)*dlglx + lglxmin
    lx = np.array([10**l for l in lglx])
    lxarr = lglx if logvals else lx/1.0e44

    if cat=='b13': return lxarr, np.zeros((n))

    if not include_CT: cat=cat+'noCT'

    A =     {'a12':113.1e-7, 'a12noCT':122.4e-7, 'b11':1.53e-5, 't08':1.8e-5}
    Lstar = {'a12':0.51e44,  'a12noCT':0.48e44,  'b11':0.53e44, 't08':10**43.85}
    gam1 =  {'a12':0.79,     'a12noCT':0.72,     'b11':0.74,    't08':0.84}
    gam2 =  {'a12':2.39,     'a12noCT':2.37,     'b11':2.60,    't08':2.55}

    fitfunc = ( A[cat] / np.log(10) / 
                ( (lx/Lstar[cat])**gam1[cat] + 
                  (lx/Lstar[cat])**gam2[cat] ) )
    if logvals: fitfunc = np.log10(fitfunc)

    return lxarr,fitfunc

#def kstest(cat, binedge, dwarf_lf, cat_lf, min_lx=40.0):
def kstest(lglxA, lglxB, VmaxA, VmaxB, min_lglx=40.0, 
           binwidth=0.005,manual=False,make_plot=False,
           path="/n/home00/lblecha/dwarf_agn_data/",lblA='A',lblB='B'):

    ixA = np.where(lglxA>=min_lglx)[0]
    ixB = np.where(lglxB>=min_lglx)[0]    
    NobjA=len(ixA)
    NobjB=len(ixB)
    print "min_lglx = %g"%min_lglx
    print "NobjA, NobjB = %d, %d"%(NobjA,NobjB)
    if NobjA==0 or NobjB==0:
        print "array A or B has size 0 for min_lglx=%g."%(NobjA,NobjB,min_lglx)
        return 0

    if manual:
        ## large-N approximation for critical D-value;
        ## multiply significance level by this factor:
        Nfac = np.sqrt((NobjA+NobjB)/(1.0*NobjA*NobjB))
        a = {0.1: 1.22, 0.05: 1.36, 0.025: 1.48, 0.01: 1.63}

        lhistA=np.array([lglxA[ixA].size])
        lhistB=np.array([lglxB[ixB].size])
        while lhistA.max()>1 or lhistB.max()>1:
            binwidth = 0.25*binwidth
            binedge = np.arange(min_lglx,47.5,binwidth)
            lhistA,tmp = np.histogram(lglxA[ixA],bins=binedge)
            lhistB,tmp = np.histogram(lglxB[ixB],bins=binedge)
        wtd_lhistA,tmp = np.histogram(lglxA[ixA],bins=binedge,weights=1/VmaxA[ixA])
        wtd_lhistB,tmp = np.histogram(lglxB[ixB],bins=binedge,weights=1/VmaxB[ixB])
        cdfA = 1.0*np.cumsum(lhistA)/lhistA.sum()
        cdfB = 1.0*np.cumsum(lhistB)/lhistB.sum()
        wtd_cdfA = 1.0*np.cumsum(wtd_lhistA)/wtd_lhistA.sum()
        wtd_cdfB = 1.0*np.cumsum(wtd_lhistB)/wtd_lhistB.sum()
        #wtd_cdfA = 1.0*np.cumsum(wtd_lhistA[::-1])[::-1]/wtd_lhistA.sum()
        #wtd_cdfB = 1.0*np.cumsum(wtd_lhistB[::-1])[::-1]/wtd_lhistB.sum()
        #print "wtd_lhistA.sum()=",wtd_lhistA.sum()
        #print "1/vmaxA sum: ",np.sum(1/VmaxA[ixA])
        #print "wtd_lhistB.sum()=",wtd_lhistB.sum()
        #print "1/vmaxB sum: ",np.sum(1/VmaxB[ixB])
        D_man = np.abs(cdfB-cdfA).max()
        bins = binedge[1:]-0.5*binwidth
        lx_D_man = bins[np.where(np.abs(cdfB-cdfA)==D_man)[0]]
        D_man_wtd = np.abs(wtd_cdfB-wtd_cdfA).max()
        lx_D_man_wtd = bins[np.where(np.abs(wtd_cdfB-wtd_cdfA)==D_man_wtd)[0]]

        if make_plot:
            fig = plt.figure(figsize=(6,4))
            plt.clf()
            plt.cla()
            #plt.yscale('log')
            plt.xlim(min_lglx-0.25,46.5)
            plt.ylim(1e-3,1.1)
            plt.plot(binedge[:-1],cdfA,'k')
            plt.plot(binedge[:-1],cdfB,'g')
            plt.plot(binedge[:-1],wtd_cdfA,'r:')
            plt.plot(binedge[:-1],wtd_cdfB,'b:')
            fig.subplots_adjust(left=0.15,right=0.95,bottom=0.15,top=0.9)
            plotname="%s/test_cumlf_lmin%g_%s_%s.eps"%(path,min_lglx,lblA,lblB)
            fig.savefig(plotname)
            #plt.show()

        print "P<0.01 requires D >~ %g"%(a[0.01]*Nfac)
        print "binwidth = %g, manual D = %g at lglx=%g"%(binwidth,D_man,np.median(lx_D_man))
        print "manual weighted D = %g at lglx=%g"%(D_man_wtd,np.median(lx_D_man_wtd))
        

    D,P = stats.ks_2samp(lglxA[ixA], lglxB[ixB])
    print "2-sample KS test: D = %g, P = %g"%(D,P)

    return D,P

def plot_lf(ax,lglx_binedge,lphi,phierr,color='k',logvals=True,
            xleg=0,yleg=0,leg_str='',xlabel=True,**kwargs):

    mask = np.ones(len(lphi),dtype=bool)
    ix0 = np.where(lphi==0)[0]
    if len(ix0)>0: mask[ix0]=False

    lglx_binwidth = lglx_binedge[1:]-lglx_binedge[:-1] 
    lglx_bins = lglx_binedge[:-1]+0.5*lglx_binwidth 
    xll = (10**lglx_bins-10**lglx_binedge[:-1])/1e44
    xul = (10**lglx_binedge[1:]-10**lglx_bins)/1e44

    x = lglx_bins if logvals else (10**lglx_bins)/1e44
    xerr = (0.5*lglx_binwidth[mask],)*2 if logvals else (xll[mask],xul[mask])

    y = np.log10(lphi[mask]) if logvals else lphi[mask]
    yerr = phierr[mask]/np.log(10)/lphi[mask] if logvals else phierr[mask]
    #print "phierr=",phierr[mask]
    #print "lphi=",lphi[mask]

    ax.errorbar(x[mask], y, yerr=yerr, xerr=xerr, color=color,
                markeredgecolor=color,fmt='o')

    if len(ix0)>0:
        y0 = np.log10(phierr[~mask]) if logvals else phierr[~mask]
        y0err = 0.75 if logvals else 10**(np.log10(phierr[~mask])-.1)
        #x0err = (0.5*lglx_binwidth[~mask],)*2 if logvals else (xll[~mask],xul[~mask])
        #print "y0,y0err=",y0,y0err
        ax.errorbar(x[~mask],y0,yerr=y0err,color=color,
                    uplims=True,ls='None',elinewidth=1)

    if logvals:
        if xlabel: plt.xlabel(r'log L$_{\rm X}$ [erg s$^{-1}]$')
        plt.ylabel(r'log L$_{\rm X}$ $\Phi$(L$_{\rm X}$) [Mpc$^{-3}]$')
    else: 
        if xlabel: plt.xlabel(r'L$_{\rm X}$ [10$^{44}$ erg s$^{-1}$]')
        plt.ylabel(r'L$_{\rm X}$ $\Phi$(L$_{\rm X}$) [Mpc$^{-3}]$')

    ax.plot([xleg],[yleg],color=color,markeredgecolor=color,marker='o')
    plt.text(xleg+0.2,yleg-0.1,leg_str,fontsize=11)

def plot_lf_ratio(ax,lglx_binedge,lphi1,lphi2,color='k',fillstyle='full',
                  lbl1='1',lbl2='2',logvals=True,xlabel=True,**kwargs):
                                           
    if len(lphi1)!=len(lphi2) or len(lphi1)!=len(lglx_binedge)-1:                               
        print "Error: size mismatch for lglx_binedge,lphi,lphi2."                               
        print len(lglx_binedge),len(lphi1),len(lphi2)                                           
        sys.exit()                                                                              
                                                                                                
    mask = np.ones(len(lphi1),dtype=bool)                                                       
    ix0 = np.where((lphi1==0)|(lphi2==0))[0]                                                    
    if len(ix0)>0: mask[ix0]=False                                                              
                                                                                                
    lglx_binwidth = lglx_binedge[1:]-lglx_binedge[:-1]                                          
    lglx_bins = lglx_binedge[:-1]+0.5*lglx_binwidth                                             
    xll = (10**lglx_bins-10**lglx_binedge[:-1])/1e44                                            
    xul = (10**lglx_binedge[1:]-10**lglx_bins)/1e44                                             
                                                                                                
    x = lglx_bins if logvals else (10**lglx_bins)/1e44                                          
    xerr = (0.5*lglx_binwidth[mask],)*2 if logvals else (xll[mask],xul[mask])                   
    y = np.log10(lphi1[mask]/lphi2[mask]) if logvals else (lphi1[mask]/lphi2[mask])             
                                                                                                
    ax.errorbar(x[mask], y, xerr=xerr, color=color,
                markeredgecolor=color, fillstyle=fillstyle, fmt='o')   
                                                                                                
    if logvals:                                                                                 
        if xlabel: plt.xlabel(r'log L$_{\rm X}$ [erg s$^{-1}]$')
        plt.ylabel(r'log $\Phi_{%s}$(L$_{\rm X}$)/$\Phi_{%s}$(L$_{\rm X}$)'%(lbl1,lbl2))        
    else:                                                                                       
        if xlabel: plt.xlabel(r'L$_{\rm X}$ [10$^{44}$ erg s$^{-1}$]')
        plt.ylabel(r'$\Phi_{\rm %s}$(L$_{\rm X}$)/$\Phi_{\rm %s}$(L$_{\rm X}$)'%(lbl1,lbl2))


def plot_cum_lf_ratio(ax,lglx_binedge,cumlf1,cumlf2,color='k',lbl1='1',lbl2='2',
                      linestyle='-',logvals=True,xlabel=True,**kwargs):

    if len(cumlf1)!=len(cumlf2) or len(cumlf1)!=len(lglx_binedge)-1:
        print "Error: size mismatch for lglx_binedge,cumlf1,cumlf2."
        print len(lglx_binedge),len(cumlf1),len(cumlf2)
        sys.exit()

    mask = np.ones(len(cumlf1),dtype=bool)
    ix0 = np.where((cumlf1<=1e-20)|(cumlf2<=1e-20))[0]
    if len(ix0)>0: mask[ix0]=False

    x = lglx_binedge[:-1] if logvals else (10**lglx_binedge[:-1])/1e44
    y = np.log10(cumlf1[mask]/cumlf2[mask]) if logvals else (cumlf1[mask]/cumlf2[mask])

    ax.plot(x[mask], y, color=color, linestyle=linestyle, drawstyle='steps-post')

    if logvals:
        if xlabel: plt.xlabel(r'log L$_{\rm X}$ [erg s$^{-1}]$')
        plt.ylabel(r'log $\Phi_{\rm %s}$(>L$_{\rm X}$)/$\Phi_{\rm %s}$(>L$_{\rm X}$)'%(lbl1,lbl2))
    else: 
        if xlabel: plt.xlabel(r'L$_{\rm X}$ [10$^{44}$ erg s$^{-1}$]')
        plt.ylabel(r'$\Phi_{\rm %s}$(>L$_{\rm X}$)/$\Phi_{\rm %s}$(>L$_{\rm X}$)'%(lbl1,lbl2))


def load_sky_coverage(cat, path="/n/home00/lblecha/dwarf_agn_data/"):

    area_filename = {'b13':'B13_fig10_data.txt','a12':'A12_fig1_data.txt'}
    ## b13 sensitivity curve is in mCrab (see eqn 6)
    flux_tab_units = {'b13':2.386e-11,'a12':1.0e-12}
    area_tab_units = {'b13':4*np.pi, 'a12':(np.pi/180.)**2}

    flux, area = np.loadtxt('%s/%s'%(path,area_filename[cat]),unpack=True)

    ## convert units and add a very bright point to flux/area table
    ## to account for high-lum sources in z-iteration:
    flux = np.append(flux, 1.0e12) * flux_tab_units[cat]
    area = np.append(area, area[-1]) * area_tab_units[cat]

    return flux, area

def calc(path="/n/home00/lblecha/dwarf_agn_data/", cat='b13', 
         sep_agn_type=False, logM_cut=10, use_cat_z=True,
         include_b13_class2=False, use_BLflag=False,
         write_src_list=False, dwarfs_in_cat_only=True,
         use_snr_limited_vmax=False, compare_manual_KS=False,
         plot_logvals=True,cosmopar=(0,)):

    ### completeness is a real issue for the comparison with tueller 08, and we're still ignoring small corrections for more recent samples. will need to ask richard how to deal with this. 
    ### giving up on the de-absoprtion correction for the CT agn -- no analytic formula
    ## and try making binned data from model as in a12, for comparison/extra sanity check
    
    cat_name = {'b13':'Baumgartner+2013', 'a12':'Ajello+2012', 
                'b11':'Burlon+2011', 't08':'Tueller+2008'}
    cat_snr_cut = {'b13':4.8, 'a12':5.0, 'b11':5.0, 't08':4.8}
    cat_emin = {'b13':14, 'a12':15, 'b11':15, 't08':14}
    cat_emax = {'b13':195, 'a12':55, 'b11':55, 't08':195}
    ## limiting all-sky flux for a12: 1mCrab (15-55keV)
    ## median 5sigma sensitivity for b13: 0.43mCrab (14-195keV)
    ## 90%-coverage 5sigma sensitivity for b13: 0.56mCrab (14-195keV)
    ## min/max 4.8sigma sensitivity for b13: ~0.29/0.67mCrab (14-195keV)
    cat_allsky_fx_lim = {'b13':1.6e-11,'a12':1.27e-11,'b11':7.3e-12,'t08':6e-11}
    cat_abs_fx_lim = {'b13':6.9e-12,'a12':6e-12,'b11':7.3e-12,'t08':3e-11}
    cat_cosmo = {'b13':(0.3,0.7),'a12':(0.27,0.71),
                 'b11':(0.27,0.7),'t08':(0.3,0.7)}

    if cat not in ('b13','a12'): 
        print "Error: sky coverage only defined for B13 and A12 catalogs."
        return
    #if cat not in ('b13','a12','b11','t08'):
    #    print "Error: cat %s is not defined."
    #    return -1

    ## default cosmology is wmap9.
    ## if catalog is defined, set cosmopar='wmap9' to override catalog cosmology
    ## can also manually set cosmopar to len-2 tuple of (Om0, h0)
    if len(cosmopar)==2:
        cosmo = FlatLambdaCDM(Om0=cosmopar[0], H0=100*cosmopar[1])
    elif cat and cosmopar!='wmap9':
        cosmo = FlatLambdaCDM(Om0=cat_cosmo[cat][0], H0=100*cat_cosmo[cat][1])
    else:
        from astropy.cosmology import WMAP9
        cosmo = WMAP9
    #cosmo = pc.Cosmology(omega_m,omega_l,0.0,-1.0,h)


    flux_tab,area_tab = load_sky_coverage('b13')
    #print flux_tab
    #print area_tab
    

    ## testing new Vmax calculation:
    ##ztmp,Vtmp = calc_Vmax(cosmo,[41.,43.,45.],flux_tab,area_tab)
    #for elem in np.arange(41.0,46.0):
    #    ztmp,Vtmp = calc_Vmax(cosmo,[elem],flux_tab,area_tab)
    #    ztmp2,Vtmp2 = calc_Vmax_old(cosmo,np.array([elem]),np.array([5.7e-12]),np.array([2.0]),np.array([area_tab[-1]]))
    #    print Vtmp/Vtmp2
    #    print ztmp/ztmp2
    #return

    ### Load (xmatched) Koss data and BAT catalog
    kb = matchcat.kosscat(batcat_id='b13',use_BLflag=use_BLflag,
                          include_b13_class2=include_b13_class2)

    #print np.unique(kb['Type'])
    #print len(kb[kb['Type']=='AGN'])
    #print len(kb[kb['Type']=='double AGN'])
    #print len(kb[kb['Type']=='XBONG'])
    #print len(kb[kb['Type']=='Galaxy'])
    #print len(kb[kb['Type']=='Galaxy Pair'])
    #print len(kb[kb['Type']=='Galaxy, XBONG'])
    #print len(kb[kb['Type']=='Sy2'])
    #print len(kb[kb['Type']=='Sy1'])
    #print np.unique(kb['Type_kpno'])
    #print np.unique(kb['Type_sdss'])
    #print len(kb[(kb['Type_kpno']=='')&(kb['Type_sdss']=='')])
    #print len(kb[(kb['Type_kpno']!='')&(kb['Type_sdss']!='')])
    print "%d Cl=2, %d Cl=4, %d Cl=5, %d Cl=6."%(kb[kb['Cl']==2].size,kb[kb['Cl']==4].size,
                                                 kb[kb['Cl']==5].size,kb[kb['Cl']==6].size)

    ## this is the B13 catalog z:
    redshift = kb['z'] if use_cat_z else kb['z_koss']
    dL = cosmo.luminosity_distance(redshift).to('cm').value
    lglx = np.log10( 4*np.pi * dL**2 * kb['Flux'] / 
                     (1+redshift)**(2-kb['Gamma']) )
    #BLAGNstr = ['Sy1','Sy1.2','Sy1.5']
    tmpmask=((kb['Type_sdss']!='')&(kb['Type_kpno']!=''))
    if kb[(kb['Type_sdss']!=kb['Type_kpno'])&(tmpmask)].size>0:
        print "Error: mismatch in Koss AGN types:"
        print kb['Type_sdss'][tmpmask]
        print kb['Type_kpno'][tmpmask]
        sys.exit()
    Type_koss = kb['Type_sdss']
    Type_koss[Type_koss==''] = kb['Type_kpno'][Type_koss=='']
    mask_noBL = ( (Type_koss!='Sy1')&(Type_koss!='Sy1.2')&(Type_koss!='Sy1.5') )
    mask_BL = ( (Type_koss=='Sy1')|(Type_koss=='Sy1.2')|(Type_koss=='Sy1.5') )
    if use_BLflag:
        mask_noBL = ((mask_noBL)&(kb['BLflag_morph']!='Y'))
        mask_BL = ((mask_BL)|(kb['BLflag_morph']=='Y'))
    print "%d BL and %d non-BL AGN."%(Type_koss[mask_BL].size,Type_koss[mask_noBL].size)
    #print "%d BL and %d non-BL AGN w/ logM<10."%(kb[~ix2][kb['logM']<10].size,
    #                                             kb[ix2][kb['logM']<10].size)
    print_array_attr(kb['logM'],'logM')
    print_array_attr(kb['Flux'],'flux')
    print_array_attr(lglx,'lglx')
    print_array_attr(kb['S_N'],'snr')
    print "\n"


    if use_snr_limited_vmax:
        fx_lim = cat_snr_cut['b13'] / kb['S_N'] * kb['Flux']
        print_array_attr(fx_lim,'fx_lim')
    else:
        fx_lim = cat_abs_fx_lim['b13']


    print "sources with lglx>43, 9.5<mstar<9.7 (name, lglx, mstar, %AGN):"
    ix=np.where((lglx>43)&(kb['logM']>=9.5)&(kb['logM']<9.7)&(kb['logM']>0))[0]
    for i in ix:
        print "%s  %g %g  %g %s"%(kb['Name'][i],lglx[i],kb['logM'][i],
                                  kb['PSr'][i],Type_koss[i])
    print "sources with lglx>43, mstar<9.5 (name, lglx, mstar, %AGN, type):"
    ix=np.where((lglx>43)&(kb['logM']<9.5)&(kb['logM']>0))[0]
    for i in ix:
        print "%s %g %g %g %s"%(kb['Name'][i],lglx[i],kb['logM'][i],
                                kb['PSr'][i],Type_koss[i])
    print "\nsources with %AGN>15%:"
    ix=np.where(kb['PSr']>15.0)[0]
    for i in ix:
        print "%s %g %g %g %s"%(kb['Name'][i],lglx[i],kb['logM'][i],
                                kb['PSr'][i],Type_koss[i])
    print "\n number of type 2 sources: %d"%Type_koss[mask_noBL].size

    #return

    ## Now calculate LF using Vmax method:
        
    lglx_binedge = np.arange(40.5,47,1.0)
    #lglx_binedge = np.arange(39.5,47,0.75)
    #lglx_binedge = np.arange(39.5,47,0.3)

    lglx_binedge_fine = np.arange(40.5,47,0.05)

    print "calculating Vmax for matched dwarf sample."
    Vmax = calc_Vmax(cosmo, lglx, flux_tab, area_tab, fx_lim=fx_lim,
                     gamma=kb['Gamma'], snr_limited=use_snr_limited_vmax)
    Lhist,L_Phi,Phi_err = lum_func(cosmo, lglx_binedge, lglx, Vmax, flux_tab, 
                                   area_tab, fx_lim=fx_lim, gamma=kb['Gamma'],
                                   use_snr_limited_vmax=use_snr_limited_vmax)
    Phi_tot = np.sum( 1.0 / Vmax ) ## in Mpc^-3    
    print_array_attr(lglx,'lglx')
    print_array_attr(Vmax,'Vmax [Mpc^3]')
    print "\nTotal space density of dwarf BAT AGN [Mpc-3]: %g"%Phi_tot
    cum_L_Phi_fine = cum_lum_func(cosmo,lglx_binedge_fine,lglx,Vmax)

    #print "Sources in each lglx bin:"
    #for i in range(lglx_bins.size):
        #print "lglx_bins[%d]=%g"%(i,lglx_bins[i])
        ##print "Vmax_bins[%d]=%g"%(i,Vmax_bins[i])
        #print "lglx=",lglx[(lglx>=lglx_binedge[i])&(lglx<lglx_binedge[i+1])]
        #print "1/Vmax=",1.0/Vmax[(lglx>=lglx_binedge[i])&(lglx<lglx_binedge[i+1])]
        ##print "zmax=",zmax[(lglx>=lglx_binedge[i])&(lglx<lglx_binedge[i+1])]
        #print "L_Phi[%d]=%g"%(i,L_Phi[i])
        ##print "Phi_binavg[%d]=%g"%(i,1.0*Lhist[i]/Vmax_bins[i])
        #print "sources:"
        #print name[(lglx>=lglx_binedge[i])&(lglx<lglx_binedge[i+1])]

    if logM_cut:
        Mmask = np.zeros((kb.size)).astype('bool')
        Mmask[kb['logM']<logM_cut] = True
        print "calculating LF for %d sources w/ logM<%g."%(kb[Mmask].size,logM_cut)
        fx_lim_Mmask=fx_lim if np.isscalar(fx_lim) else fx_lim[Mmask]
        Vmax_Mcut = calc_Vmax(cosmo, lglx[Mmask], flux_tab, area_tab, fx_lim=fx_lim_Mmask,
                              gamma=kb['Gamma'][Mmask], snr_limited=use_snr_limited_vmax)
        Phi_tot_Mcut = np.sum( 1.0 / Vmax_Mcut ) ## in Mpc^-3    
        #Lhist_Mcut,L_Phi_Mcut,Phi_Mcut_err = lum_func(cosmo, lglx_binedge, lglx[Mmask], 
        LF_Mcut = lum_func(cosmo, lglx_binedge, lglx[Mmask], Vmax_Mcut, flux_tab, 
                           area_tab, fx_lim=fx_lim_Mmask, gamma=kb['Gamma'],
                           use_snr_limited_vmax=use_snr_limited_vmax)
        cum_L_Phi_fine_Mcut = cum_lum_func(cosmo,lglx_binedge_fine,lglx[Mmask],Vmax_Mcut)

        plt.plot(kb['logM'][Mmask],lglx[Mmask],'ko')
        print kb['logM'][Mmask]
        print kb['Name'][Mmask]
        for i in range(kb[Mmask].size):
            print "%s %g %g"%(kb['Name'][Mmask][i],kb['RA'][Mmask][i],kb['DE'][Mmask][i])
            #print "%g %g"%(kb['RA'][Mmask][i],kb['DE'][Mmask][i])

    if sep_agn_type:
        Phi_tot_type1 = np.sum( 1.0 / Vmax[mask_BL] )
        Phi_tot_type2 = np.sum( 1.0 / Vmax[mask_noBL] )
        print "%d BL AGN in sample."%Type_koss[mask_BL].size
        print "Total BL AGN space density [Mpc-3]: %g"%Phi_tot_type1
        print "%d non-BL AGN in sample."%Type_koss[mask_noBL].size
        print "Total non-BL AGN space density [Mpc-3]: %g"%Phi_tot_type2

        fx_lim_type1=fx_lim if np.isscalar(fx_lim) else fx_lim[mask_BL]
        fx_lim_type2=fx_lim if np.isscalar(fx_lim) else fx_lim[mask_noBL]
        LF_type1=lum_func(cosmo, lglx_binedge,lglx[mask_BL],Vmax[mask_BL],
                          flux_tab,area_tab,fx_lim=fx_lim_type1,gamma=kb['Gamma'],
                          use_snr_limited_vmax=use_snr_limited_vmax)
        LF_type2=lum_func(cosmo, lglx_binedge,lglx[mask_noBL],Vmax[mask_noBL],
                          flux_tab,area_tab,fx_lim=fx_lim_type2,gamma=kb['Gamma'],
                          use_snr_limited_vmax=use_snr_limited_vmax)

    ### Get external catalog data
    catdata = load_cat(cosmo, path=path,cat=cat,
                       include_b13_class2=include_b13_class2)
    if catdata:
        cat_ra,cat_dec,cat_z,cat_flux,cat_snr,cat_gamma,cat_nh,cat_isCT,cat_swiftname,cat_lglx = catdata
        print_array_attr(cat_z,'cat_z')
        print_array_attr(cat_nh,'cat_nh')
        #print "CT nh: ",cat_nh[cat_nh>=24]
        
        cat_flux_tab,cat_area_tab = load_sky_coverage(cat)

        if use_snr_limited_vmax:
            cat_fx_lim = cat_snr_cut[cat] / cat_snr * cat_flux
        else:
            cat_fx_lim = cat_abs_fx_lim[cat]

        print_array_attr(cat_gamma,'cat_gamma')
        print_array_attr(cat_snr,'cat_snr')
        print_array_attr(cat_lglx,'cat_lglx')
        print_array_attr(cat_flux,'cat_flux')
    
        print "from koss data:"
        print 'NGC4395:'
        print "lglx:",lglx[lglx<41]
        print "flux:",kb['Flux'][lglx<41]
        print "z: ",redshift[lglx<41]
        print kb['Name'][lglx<41]
        print "from cat data:"
        print "lglx:",cat_lglx[cat_lglx<41]
        print "flux:",cat_flux[cat_lglx<41]
        print "z:",cat_z[cat_lglx<41]

        print "calculating Vmax and LF for cat data (%d objects)..."%len(cat_lglx)
        cat_Vmax = calc_Vmax(cosmo, cat_lglx, cat_flux_tab, cat_area_tab, 
                             zarr_max=0.6,fx_lim=cat_fx_lim, gamma=cat_gamma, 
                             snr_limited=use_snr_limited_vmax)
        LF_cat = lum_func(cosmo, lglx_binedge, cat_lglx, cat_Vmax, cat_flux_tab, cat_area_tab,
                          fx_lim=cat_fx_lim, gamma=cat_gamma, use_snr_limited_vmax=use_snr_limited_vmax)
        cat_Phi_tot = np.sum( 1.0/cat_Vmax ) ## in Mpc^-3    
        print "\nTotal space density of cat objects [Mpc-3]: %g"%cat_Phi_tot
        cat_Phi_CT_tot = np.sum( 1.0/cat_Vmax[cat_isCT] ) ## in Mpc^-3    
        print "Total space density of CT cat objects [Mpc-3]: %g"%cat_Phi_CT_tot
        n_CT = len(cat_isCT[cat_isCT])
        print "number of CT objects: %d\n"%n_CT
        cum_L_Phi_fine_cat = cum_lum_func(cosmo,lglx_binedge_fine,cat_lglx,cat_Vmax)

        if n_CT>0:
            cat_fx_lim_CT=cat_fx_lim if np.isscalar(cat_fx_lim) else cat_fx_lim[cat_isCT]
            LF_cat_CT = lum_func(cosmo,lglx_binedge,cat_lglx[cat_isCT],cat_Vmax[cat_isCT], 
                                 flux_tab, area_tab,fx_lim=cat_fx_lim_CT,gamma=cat_gamma[cat_isCT],
                                 use_snr_limited_vmax=use_snr_limited_vmax)

        print "calculating LF for cat data at z<0.05 only (%d objects)..."%len(cat_lglx[cat_z<0.05])
        cat_fx_lim_zcut = cat_fx_lim if np.isscalar(cat_fx_lim) else cat_fx_lim[cat_z<0.05]
        LF_cat_zcut = lum_func(cosmo, lglx_binedge, cat_lglx[cat_z<0.05], cat_Vmax[cat_z<0.05], 
                               cat_flux_tab, cat_area_tab, fx_lim=cat_fx_lim_zcut, 
                               gamma=cat_gamma[cat_z<0.05], use_snr_limited_vmax=use_snr_limited_vmax)
        cum_L_Phi_fine_cat_zcut = cum_lum_func(cosmo,lglx_binedge_fine,
                                               cat_lglx[cat_z<0.05],cat_Vmax[cat_z<0.05])

        ### LF fit formulae from catalog paper(s):
        lxarr,cat_lffit = cat_lf_fitfunc(cat,logvals=plot_logvals)
        if cat=='a12': 
            lxarr_noct,cat_lffit_noct = cat_lf_fitfunc(cat,logvals=plot_logvals,
                                                       include_CT=False)

        ## K-S test:
        print "\nK-S test (Koss vs B13):"
        lglx_ks = [40.5,41.5,42.5]
        D_kb_b13,P_kb_b13=np.array(zip(*[kstest(lglx,cat_lglx,Vmax,cat_Vmax,min_lglx=l,
                                                manual=compare_manual_KS,lblA='Koss',lblB='B13') 
                                         for l in lglx_ks ]) )
        #Dtmp,Ptmp=kstest(lglx,cat_lglx,Vmax,cat_Vmax,min_lglx=lglx_ks[0],
        #                 manual=compare_manual_KS,lblA='Koss',lblB='B13') 
        print "\nK-S test (Koss vs B13, z<0.05 only):"
        tmp=np.array(zip(*[kstest(lglx,cat_lglx[cat_z<0.05],Vmax,cat_Vmax[cat_z<0.05],
                                  min_lglx=l,manual=compare_manual_KS,lblA='Koss',lblB='B13') 
                           for l in lglx_ks ]) )
        D_kb_b13zcut,P_kb_b13zcut = tmp
        print "Dstat array (z<0.05 only): ",D_kb_b13zcut
        print "Pval array (z<0.05 only): ",P_kb_b13zcut
        #return 0

        if logM_cut:
            print "\nK-S test (Koss Mcut vs B13):"
            tmp=np.array(zip(*[kstest(lglx[Mmask],cat_lglx,Vmax_Mcut,cat_Vmax,min_lglx=l,
                                      manual=compare_manual_KS,lblA='Mcut',lblB='B13') 
                               for l in lglx_ks]))
            D_Mcut_b13,P_Mcut_b13 = tmp
            #Dtmp,Ptmp=kstest(lglx[Mmask],cat_lglx,Vmax_Mcut,cat_Vmax,min_lglx=lglx_ks[0],
            #                 manual=compare_manual_KS,lblA='Mcut',lblB='B13') 
            print "\nK-S test (Koss Mcut vs B13, z<0.05 only):"
            tmp=np.array(zip(*[kstest(lglx[Mmask],cat_lglx[cat_z<0.05],Vmax_Mcut,cat_Vmax[cat_z<0.05],
                                      min_lglx=l,manual=compare_manual_KS,lblA='Mcut',lblB='B13')
                                                 for l in lglx_ks]))
            D_Mcut_b13zcut,P_Mcut_b13zcut = tmp
            print "Dstat array (z<0.05 only): ",D_Mcut_b13
            print "Pval array (z<0.05 only): ",P_Mcut_b13

            print "\nK-S test (Koss Mcut vs Koss total):"
            tmp=np.array(zip(*[kstest(lglx[Mmask],lglx,Vmax_Mcut,Vmax,min_lglx=l,
                                      manual=compare_manual_KS,lblA='Mcut',lblB='Koss') 
                               for l in lglx_ks ]) )
            D_kb_Mcut,P_kb_Mcut = tmp
            print "Dstat array: ",D_kb_Mcut
            print "Pval array: ",P_kb_Mcut

    extra=''
    if logM_cut: extra='_mcut%g'%logM_cut
    if use_snr_limited_vmax: extra=extra+'_snrVmax'
    fmtstr=(6*"%g "+"\n")
    with open("%s/kstest_%scat%s.dat"%(path,cat,extra),"w") as fp:
        for j in range(len(lglx_ks)):
            fp.write(fmtstr%(lglx_ks[j],P_kb_b13[j],P_kb_b13zcut[j],
                             P_Mcut_b13[j],P_Mcut_b13zcut[j],P_kb_Mcut[j]))

    ### Make LF plot: ###
    plt.clf()
    plt.cla()
    plt.close('all')
    fig = plt.figure(figsize=(5.5,8))
    kwargs = dict(logvals=plot_logvals)
    x0,x1 = 40.2,46.3
    y0,y1 = -11.2,-1.7
    xlim = (x0,x1) if plot_logvals else (10**x0/1e44,10**x1/1e44)
    ylim = (y0,y1) if plot_logvals else (10**y0,10**y1) 
    xleg=40.5
    yleg0=-7.5

    ax1 = fig.add_subplot(311)
    plt.xlim(xlim)
    plt.ylim(ylim)
    if not plot_logvals: plt.xscale('log')
    if not plot_logvals: plt.yscale('log')
    # LF_[..] vars are len-3 tuples: Lhist, L_Phi, Phi_err
    if catdata:
        plot_lf(ax1,lglx_binedge,LF_cat[1],LF_cat[2],color='k',xleg=xleg,yleg=yleg0-0.8,
                leg_str='All BAT AGN (N=%d)'%cat_lglx.size,xlabel=False,**kwargs)
        ax1.plot(lxarr, cat_lffit, 'k')
        if cat=='a12': ax1.plot(lxarr_noct, cat_lffit_noct, 'k', linestyle='dashed')
    #if sep_agn_type:
    #    plot_lf(ax1,lglx_binedge,LF_type1[1],LF_type1[2],color='b',xlabel=False,**kwargs)
    #    plot_lf(ax1,lglx_binedge,LF_type2[1],LF_type2[2],color='r',xlabel=False,**kwargs)
    #else:
    plot_lf(ax1,lglx_binedge,L_Phi,Phi_err,color='m',xleg=xleg,yleg=yleg0-1.6,
            leg_str='Koss sample (N=%d)'%lglx.size,xlabel=False,**kwargs)
    plt.text(xleg,yleg0-2.5,'Lx>%g:  K-S P=%.2g'%(lglx_ks[0],P_kb_b13[0]),fontsize=11)
    plt.text(xleg,yleg0-3.2,'Lx>%g:  K-S P=%.2g'%(lglx_ks[1],P_kb_b13[1]),fontsize=11)

    ax2 = fig.add_subplot(312)
    plt.xlim(xlim)
    plt.ylim(ylim)
    if not plot_logvals: plt.xscale('log')
    if not plot_logvals: plt.yscale('log')
    if catdata:
        plot_lf(ax2,lglx_binedge,LF_cat[1],LF_cat[2],color='k',xleg=xleg,yleg=yleg0-0.8,
                leg_str='All BAT AGN (N=%d)'%cat_lglx.size,xlabel=False,**kwargs)
        ax2.plot(lxarr, cat_lffit, 'k')
        if cat=='a12': ax2.plot(lxarr_noct, cat_lffit_noct, 'k', linestyle='dashed')

    if logM_cut:
        plot_lf(ax2,lglx_binedge,LF_Mcut[1],LF_Mcut[2],color='c',xleg=xleg,yleg=yleg0-1.6,
                leg_str='Koss, logM<%g (N=%d)'%(logM_cut,lglx[Mmask].size),xlabel=False,**kwargs)
        plt.text(xleg,yleg0-2.5,'Lx>%g:  K-S P=%.2g'%(lglx_ks[0],P_Mcut_b13[0]),fontsize=11)
        plt.text(xleg,yleg0-3.2,'Lx>%g:  K-S P=%.2g'%(lglx_ks[1],P_Mcut_b13[1]),fontsize=11)

    ax3 = fig.add_subplot(313)
    plt.xlim(xlim)
    plt.ylim(ylim)
    if not plot_logvals: plt.xscale('log')
    if not plot_logvals: plt.yscale('log')
    if catdata:
        plot_lf(ax3,lglx_binedge,LF_cat[1],LF_cat[2],color='lightgray',xleg=xleg,yleg=yleg0+0.1,
                leg_str='All BAT AGN (N=%d)'%cat_lglx.size,**kwargs)
        ax3.plot(lxarr, cat_lffit, 'k')
        if cat=='a12': ax3.plot(lxarr_noct, cat_lffit_noct, 'k', linestyle='dashed')
    plot_lf(ax3,lglx_binedge,L_Phi,Phi_err,color='m',xleg=xleg,yleg=yleg0-0.8,
            leg_str='Koss sample (N=%d)'%lglx.size,**kwargs)
    if logM_cut: 
        plot_lf(ax3,lglx_binedge,LF_Mcut[1],LF_Mcut[2],color='c',xleg=xleg,yleg=yleg0-1.6,
                leg_str='Koss, logM<%g (N=%d)'%(logM_cut,lglx[Mmask].size),**kwargs)
        plt.text(xleg,yleg0-2.5,'Lx>%g:  K-S P=%.2g'%(lglx_ks[0],P_kb_Mcut[0]),fontsize=11)
        plt.text(xleg,yleg0-3.2,'Lx>%g:  K-S P=%.2g'%(lglx_ks[1],P_kb_Mcut[1]),fontsize=11)


    extra = {'b13':'70-month','a12':'60-month','b11':'36-month','t08':'9-month'}
    titlestr='L(%d-%dkeV) BAT AGN LF (%s %s data)'%(cat_emin[cat],cat_emax[cat],
                                                    cat_name[cat],extra[cat])
    fig.suptitle(titlestr)
    #fig.subplots_adjust(left=0.15,right=0.95,bottom=0.15,top=0.9)
    fig.subplots_adjust(left=0.16,right=0.92,bottom=0.07,top=0.95,hspace=0.18)

    print '\nsaving LF plot..'
    extra=''
    if logM_cut: extra='_mcut%g'%logM_cut
    if use_snr_limited_vmax: extra=extra+'_snrVmax'
    if plot_logvals: extra=extra+'_log'
    #plotname="%s/bat_agn_lf_%scat%s.eps"%(path,cat,extra)
    plotname="%s/lf_%scat%s.eps"%(path,cat,extra)
    fig.savefig(plotname)


    ### Make LF plot for z<0.05 only: ###
    plt.clf()
    plt.cla()
    plt.close('all')
    fig = plt.figure(figsize=(5.5,8))
    kwargs = dict(logvals=plot_logvals)
    xlim = (x0,x1) if plot_logvals else (10**x0/1e44,10**x1/1e44)
    ylim = (y0,y1) if plot_logvals else (10**y0,10**y1) 

    ax1 = fig.add_subplot(311)
    plt.xlim(xlim)
    plt.ylim(ylim)
    if not plot_logvals: plt.xscale('log')
    if not plot_logvals: plt.yscale('log')
    # LF_[..] vars are len-3 tuples: Lhist, L_Phi, Phi_err
    if catdata:
        plot_lf(ax1,lglx_binedge,LF_cat_zcut[1],LF_cat_zcut[2],color='k',xleg=xleg,yleg=yleg0-0.8,
                leg_str='z<0.05 BAT AGN (N=%d)'%cat_lglx.size,xlabel=False,**kwargs)
    plot_lf(ax1,lglx_binedge,L_Phi,Phi_err,color='m',xleg=xleg,yleg=yleg0-1.6,
            leg_str='Koss sample (N=%d)'%lglx.size,xlabel=False,**kwargs)
    plt.text(xleg,yleg0-2.5,'Lx>%g:  K-S P=%.2g'%(lglx_ks[0],P_kb_b13zcut[0]),fontsize=11)
    plt.text(xleg,yleg0-3.2,'Lx>%g:  K-S P=%.2g'%(lglx_ks[1],P_kb_b13zcut[1]),fontsize=11)

    ax2 = fig.add_subplot(312)
    plt.xlim(xlim)
    plt.ylim(ylim)
    if not plot_logvals: plt.xscale('log')
    if not plot_logvals: plt.yscale('log')
    if catdata:
        plot_lf(ax2,lglx_binedge,LF_cat_zcut[1],LF_cat_zcut[2],color='k',xleg=xleg,yleg=yleg0-0.8,
                leg_str='z<0.05 BAT AGN (N=%d)'%cat_lglx[cat_z<0.05].size,xlabel=False,**kwargs)
    if logM_cut:
        plot_lf(ax2,lglx_binedge,LF_Mcut[1],LF_Mcut[2],color='c',xleg=xleg,yleg=yleg0-1.6,
                leg_str='Koss, logM<%g (N=%d)'%(logM_cut,lglx[Mmask].size),xlabel=False,**kwargs)
        plt.text(xleg,yleg0-2.5,'Lx>%g:  K-S P=%.2g'%(lglx_ks[0],P_Mcut_b13zcut[0]),fontsize=11)
        plt.text(xleg,yleg0-3.2,'Lx>%g:  K-S P=%.2g'%(lglx_ks[1],P_Mcut_b13zcut[1]),fontsize=11)

    ax3 = fig.add_subplot(313)
    plt.xlim(xlim)
    plt.ylim(ylim)
    if not plot_logvals: plt.xscale('log')
    if not plot_logvals: plt.yscale('log')
    if catdata:
        plot_lf(ax3,lglx_binedge,LF_cat_zcut[1],LF_cat_zcut[2],color='lightgray',xleg=xleg,
                yleg=yleg0+0.1,leg_str='z<0.05 BAT AGN (N=%d)'%cat_lglx[cat_z<0.05].size,**kwargs)
    plot_lf(ax3,lglx_binedge,L_Phi,Phi_err,color='m',xleg=xleg,yleg=yleg0-0.8,
            leg_str='Koss sample (N=%d)'%lglx.size,**kwargs)
    if logM_cut: 
        plot_lf(ax3,lglx_binedge,LF_Mcut[1],LF_Mcut[2],color='c',xleg=xleg,yleg=yleg0-1.6,
                leg_str='Koss, logM<%g (N=%d)'%(logM_cut,lglx[Mmask].size),**kwargs)
        plt.text(xleg,yleg0-2.5,'Lx>%g:  K-S P=%.2g'%(lglx_ks[0],P_kb_Mcut[0]),fontsize=11)
        plt.text(xleg,yleg0-3.2,'Lx>%g:  K-S P=%.2g'%(lglx_ks[1],P_kb_Mcut[1]),fontsize=11)


    extra = {'b13':'70-month','a12':'60-month','b11':'36-month','t08':'9-month'}
    titlestr='L(%d-%dkeV) BAT AGN LF (%s %s data)'%(cat_emin[cat],cat_emax[cat],
                                                    cat_name[cat],extra[cat])
    fig.suptitle(titlestr)
    fig.subplots_adjust(left=0.16,right=0.92,bottom=0.07,top=0.95,hspace=0.18)

    print '\nsaving LF zcut plot..'
    extra=''
    if logM_cut: extra='_mcut%g'%logM_cut
    if use_snr_limited_vmax: extra=extra+'_snrVmax'
    if plot_logvals: extra=extra+'_log'
    #plotname="%s/bat_agn_lf_%scat%s.eps"%(path,cat,extra)
    plotname="%s/lf_zcut_%scat%s.eps"%(path,cat,extra)
    fig.savefig(plotname)


    ### Make LF ratio plot: ###
    plt.clf()
    plt.cla()
    plt.close()
    fig = plt.figure(figsize=(5.5,8))
    kwargs = dict(logvals=plot_logvals)
    y0,y1 = -2.5,0.2
    xlim = (x0,x1) if plot_logvals else (10**x0/1e44,10**x1/1e44)
    ylim = (y0,y1) if plot_logvals else (10**y0,10**y1) 

    ax1 = fig.add_subplot(311)
    plt.xlim(xlim)
    plt.ylim(ylim)
    if not plot_logvals: plt.xscale('log')
    if not plot_logvals: plt.yscale('log')
    lbl2='B+13' if cat=='b13' else cat
    # LF_[..] vars are len-3 tuples: Lhist, L_Phi, Phi_err
    if catdata: 
        plot_lf_ratio(ax1,lglx_binedge,L_Phi,LF_cat[1],color='m',
                      lbl1='Koss',lbl2=lbl2,xlabel=False,**kwargs)        
        plot_lf_ratio(ax1,lglx_binedge,L_Phi,LF_cat_zcut[1],color='m',
                      fillstyle='none',lbl1='Koss',lbl2=lbl2,xlabel=False,**kwargs)
        ax1.plot([40.6],[-2.05],'mo',mec='m')
        ax1.plot([40.6],[-2.25],'mo',mec='m',fillstyle='none')
        plt.text(40.8,-2.1,'All BAT AGN')
        plt.text(40.8,-2.3,'z<0.05')

    ax2 = fig.add_subplot(312)
    plt.xlim(xlim)
    plt.ylim(ylim)
    if not plot_logvals: plt.xscale('log')
    if not plot_logvals: plt.yscale('log')
    if logM_cut and catdata: 
        plot_lf_ratio(ax2,lglx_binedge,LF_Mcut[1],LF_cat[1],color='c',
                      mec='c',lbl1='Mcut',lbl2=lbl2,xlabel=False,**kwargs)
        plot_lf_ratio(ax2,lglx_binedge,LF_Mcut[1],LF_cat_zcut[1],color='c',
                      fillstyle='none',lbl1='Mcut',lbl2=lbl2,xlabel=False,**kwargs)
        ax2.plot([40.6],[-2.05],'co',mec='c')
        ax2.plot([40.6],[-2.25],'co',mec='c',fillstyle='none')
        plt.text(40.8,-2.1,'All BAT AGN')
        plt.text(40.8,-2.3,'z<0.05')

    ax3 = fig.add_subplot(313)
    plt.xlim(xlim)
    plt.ylim(ylim)
    if not plot_logvals: plt.xscale('log')
    if not plot_logvals: plt.yscale('log')
    if logM_cut: plot_lf_ratio(ax3,lglx_binedge,LF_Mcut[1],L_Phi,
                               color='g',lbl1='Mcut',lbl2='Koss',**kwargs)
    
    extra = {'b13':'70-month','a12':'60-month','b11':'36-month','t08':'9-month'}
    titlestr='L(%d-%dkeV) BAT AGN LF (%s %s data)'%(cat_emin[cat],cat_emax[cat],
                                                    cat_name[cat],extra[cat])
    fig.suptitle(titlestr)
    #fig.subplots_adjust(left=0.15,right=0.95,bottom=0.15,top=0.9)
    #fig.subplots_adjust(left=0.15,right=0.95,bottom=0.1,top=0.92,hspace=0.3)
    fig.subplots_adjust(left=0.16,right=0.92,bottom=0.07,top=0.95,hspace=0.18)

    print '\nsaving LF ratio plot..'
    extra=''
    if logM_cut: extra='_mcut%g'%logM_cut
    if use_snr_limited_vmax: extra=extra+'_snrVmax'
    if plot_logvals: extra=extra+'_log'
    plotname="%s/lf_ratio_%scat%s.eps"%(path,cat,extra)
    fig.savefig(plotname)

    ### Make cum LF plot: ###
    plt.clf()
    plt.cla()
    plt.close('all')
    fig = plt.figure(figsize=(6,4))
    kwargs = dict(logvals=plot_logvals)
    y0,y1 = -11.2,-1.7
    xlim = (x0,x1) if plot_logvals else (10**x0/1e44,10**x1/1e44)
    ylim = (y0,y1) if plot_logvals else (10**y0,10**y1) 
    xleg=40.5
    yleg0=-7.5
    ax1 = fig.add_subplot(111)
    plt.xlim(xlim)
    plt.ylim(ylim)
    if not plot_logvals: plt.xscale('log')
    if not plot_logvals: plt.yscale('log')

    if catdata:
        ax1.plot(lglx_binedge_fine[:-1],np.log10(cum_L_Phi_fine_cat),
                 color='gray',linewidth=2,drawstyle='steps-post')
        ax1.plot(lglx_binedge_fine[:-1],np.log10(cum_L_Phi_fine_cat_zcut),
                 color='k',drawstyle='steps-post')
    ax1.plot(lglx_binedge_fine[:-1],np.log10(cum_L_Phi_fine),color='m',drawstyle='steps-post')
    if logM_cut:
        ax1.plot(lglx_binedge_fine[:-1],np.log10(cum_L_Phi_fine_Mcut),
                 color='c',drawstyle='steps-post')
    if plot_logvals:
        plt.xlabel(r'log L$_{\rm X}$ [erg s$^{-1}]$')
        plt.ylabel(r'log $\Phi$(>L$_{\rm X}$) [Mpc$^{-3}]$')
    else: 
        plt.xlabel(r'L$_{\rm X}$ [10$^{44}$ erg s$^{-1}$]')
        plt.ylabel(r'$\Phi$(>L$_{\rm X}$) [Mpc$^{-3}]$')

    xlbl0=44.1
    ylbl0=-2.5
    ax1.plot([xlbl0,xlbl0+0.3],[ylbl0,ylbl0],'gray')
    ax1.plot([xlbl0,xlbl0+0.3],[ylbl0-0.5,ylbl0-0.5],'k')
    ax1.plot([xlbl0,xlbl0+0.3],[ylbl0-1,ylbl0-1],'m')
    ax1.plot([xlbl0,xlbl0+0.3],[ylbl0-1.5,ylbl0-1.5],'c')
    plt.text(xlbl0+0.4,ylbl0-0.1,'All BAT AGN',fontsize=11)
    plt.text(xlbl0+0.4,ylbl0-0.6,'z<0.05',fontsize=11)
    plt.text(xlbl0+0.4,ylbl0-1.1,'Koss sample',fontsize=11)
    plt.text(xlbl0+0.4,ylbl0-1.6,'Koss (M<%g)'%logM_cut,fontsize=11)


    extra = {'b13':'70-mo.','a12':'60-mo.','b11':'36-mo.','t08':'9-mo.'}
    titlestr='L(%d-%dkeV) BAT AGN cumulative LF (%s data)'%(cat_emin[cat],cat_emax[cat],
                                                               extra[cat])
    fig.suptitle(titlestr)
    fig.subplots_adjust(left=0.15,right=0.95,bottom=0.15,top=0.9)
    #fig.subplots_adjust(left=0.16,right=0.92,bottom=0.07,top=0.95,hspace=0.18)

    print '\nsaving cum LF plot..'
    extra=''
    if logM_cut: extra='_mcut%g'%logM_cut
    if use_snr_limited_vmax: extra=extra+'_snrVmax'
    if plot_logvals: extra=extra+'_log'
    #plotname="%s/bat_agn_cumlf_%scat%s.eps"%(path,cat,extra)
    plotname="%s/cumlf_%scat%s.eps"%(path,cat,extra)
    fig.savefig(plotname)

    ### Make cum LF ratio plot: ###
    plt.clf()
    plt.cla()
    plt.close()
    fig = plt.figure(figsize=(5.5,8))
    kwargs = dict(logvals=plot_logvals)
    y0,y1 = -2.5,0.2
    xlim = (x0,x1) if plot_logvals else (10**x0/1e44,10**x1/1e44)
    ylim = (y0,y1) if plot_logvals else (10**y0,10**y1) 

    ax1 = fig.add_subplot(311)
    plt.xlim(xlim)
    plt.ylim(ylim)
    if not plot_logvals: plt.xscale('log')
    if not plot_logvals: plt.yscale('log')
    if catdata:
        # LF_[..] vars are len-3 tuples: Lhist, L_Phi, Phi_err
        plot_cum_lf_ratio(ax1,lglx_binedge_fine,cum_L_Phi_fine,cum_L_Phi_fine_cat,
                          color='m',lbl1='Koss',lbl2=lbl2,xlabel=False,**kwargs)
        plot_cum_lf_ratio(ax1,lglx_binedge_fine,cum_L_Phi_fine,cum_L_Phi_fine_cat_zcut,
                          color='m',linestyle=':',lbl1='Koss',lbl2=lbl2,xlabel=False,**kwargs)
        ax1.plot([40.5,40.8],[-2.05,-2.05],'m')
        ax1.plot([40.5,40.8],[-2.25,-2.25],'m:')
        plt.text(40.9,-2.1,'All BAT AGN')
        plt.text(40.9,-2.3,'z<0.05')
        
    ax2 = fig.add_subplot(312)
    plt.xlim(xlim)
    plt.ylim(ylim)
    if not plot_logvals: plt.xscale('log')
    if not plot_logvals: plt.yscale('log')
    lbl2='B+13' if cat=='b13' else cat
    if catdata and logM_cut:
        plot_cum_lf_ratio(ax2,lglx_binedge_fine,cum_L_Phi_fine_Mcut,cum_L_Phi_fine_cat,
                          color='c',lbl1='Mcut',lbl2=lbl2,xlabel=False,**kwargs)
        plot_cum_lf_ratio(ax2,lglx_binedge_fine,cum_L_Phi_fine_Mcut,cum_L_Phi_fine_cat_zcut,
                          color='c',linestyle=':',lbl1='Mcut',lbl2=lbl2,xlabel=False,**kwargs)
        ax2.plot([40.5,40.8],[-2.05,-2.05],'c')
        ax2.plot([40.5,40.8],[-2.25,-2.25],'c:')
        plt.text(40.9,-2.1,'All BAT AGN')
        plt.text(40.9,-2.3,'z<0.05')

    ax3 = fig.add_subplot(313)
    plt.xlim(xlim)
    plt.ylim(ylim)
    if not plot_logvals: plt.xscale('log')
    if not plot_logvals: plt.yscale('log')
    if logM_cut:
        plot_cum_lf_ratio(ax3,lglx_binedge_fine,cum_L_Phi_fine_Mcut,cum_L_Phi_fine,
                          color='g',lbl1='Mcut',lbl2='Koss',**kwargs)
    
    extra = {'b13':'70-mo.','a12':'60-mo.','b11':'36-mo.','t08':'9-mo.'}
    titlestr='L(%d-%dkeV) BAT AGN LF (%s data)'%(cat_emin[cat],cat_emax[cat],
                                                 extra[cat])
    fig.suptitle(titlestr)
    #fig.subplots_adjust(left=0.15,right=0.95,bottom=0.15,top=0.9)
    #fig.subplots_adjust(left=0.15,right=0.95,bottom=0.1,top=0.92,hspace=0.3)
    fig.subplots_adjust(left=0.16,right=0.92,bottom=0.07,top=0.95,hspace=0.18)

    print '\nsaving cum LF ratio plot..'
    extra=''
    if logM_cut: extra='_mcut%g'%logM_cut
    if use_snr_limited_vmax: extra=extra+'_snrVmax'
    if plot_logvals: extra=extra+'_log'
    plotname="%s/cumlf_ratio_%scat%s.eps"%(path,cat,extra)
    fig.savefig(plotname)


    ### Make redshift distribution plot: ###
    print '\nmax cat_z=%g, koss z=%g, mcut koss z=%g'%(cat_z.max(),redshift.max(),redshift[Mmask].max())
    print "2-sample KS test for redshift distributions:"
    D,P = stats.ks_2samp(redshift,cat_z[cat_z<0.05])
    print "Koss vs B13(z<0.05): D = %g, P = %g"%(D,P)
    D,P = stats.ks_2samp(redshift[Mmask],cat_z[cat_z<0.05])
    print "Mcut vs B13(z<0.05): D = %g, P = %g"%(D,P)
    D,P = stats.ks_2samp(redshift,cat_z)
    print "Koss vs B13: D = %g, P = %g"%(D,P)
    D,P = stats.ks_2samp(redshift[Mmask],cat_z)
    print "Mcut vs B13: D = %g, P = %g"%(D,P)
    D,P = stats.ks_2samp(redshift[Mmask],redshift)
    print "Mcut vs Koss: D = %g, P = %g\n"%(D,P)

    plt.clf()
    plt.cla()
    plt.close()
    fig = plt.figure(figsize=(5,5))
    ax1 = fig.add_subplot(211)
    #plt.xlabel('redshift')
    plt.ylabel('N')
    kwargs = dict(bins=30,range=(0.0,0.2),histtype='step',normed=False,cumulative=False)
    if catdata: ax1.hist(cat_z,color='gray',linewidth=2.0,**kwargs)
    ax1.hist(redshift,color='m',**kwargs)
    if logM_cut: ax1.hist(redshift[Mmask],color='c',linestyle='dotted',**kwargs)

    ax2 = fig.add_subplot(212)
    plt.xlabel('redshift')
    plt.ylabel('normalized CDF')
    kwargs = dict(bins=100,range=(0.0,0.2),histtype='step',normed=True,cumulative=True)
    if catdata: ax2.hist(cat_z,color='gray',linewidth=2.0,**kwargs)
    if catdata: ax2.hist(cat_z[cat_z<0.05],color='k',linewidth=0.7,**kwargs)
    ax2.hist(redshift,color='m',**kwargs)
    if logM_cut: ax2.hist(redshift[Mmask],color='c',linestyle='dotted',**kwargs)

    xlbl0=0.12
    ylbl0=0.42
    ax2.plot([xlbl0,xlbl0+0.015],[ylbl0,ylbl0],'gray',linewidth=2)
    ax2.plot([xlbl0,xlbl0+0.015],[ylbl0-0.1,ylbl0-0.1],'k')
    ax2.plot([xlbl0,xlbl0+0.015],[ylbl0-0.2,ylbl0-0.2],'m')
    ax2.plot([xlbl0,xlbl0+0.015],[ylbl0-0.3,ylbl0-0.3],'c:')
    plt.text(xlbl0+0.02,ylbl0-0.02,'All BAT AGN',fontsize=11)
    plt.text(xlbl0+0.02,ylbl0-0.12,'z<0.05',fontsize=11)
    plt.text(xlbl0+0.02,ylbl0-0.22,'Koss sample',fontsize=11)
    plt.text(xlbl0+0.02,ylbl0-0.32,'Koss (M<%g)'%logM_cut,fontsize=11)


    print '\nsaving z distribution plot..'
    extra=''
    if logM_cut: extra='_mcut%g'%logM_cut
    fig.savefig('%s/zhist_%scat%s.eps'%(path,cat,extra))




