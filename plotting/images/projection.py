import numpy as np
import ctypes
import physicalmodels.stellarproperties.stellar_luminosities as star_props
import physicalmodels.attenuation.attenuate as atten
import util.cast



##
## routine to use 'raytrace_projection_compute' to make mock stellar images,
##   treating the starlight as sources and accounting for gas extinction via ray-tracing
##
## important to set KAPPA_UNITS appropriately: code loads opacity (kappa) for
##   the bands of interest in cgs (cm^2/g), must be converted to match units of input
##   mass and size. the default it to assume gadget units (M=10^10 M_sun, l=kpc)
##
def stellar_raytrace(stellar_x, stellar_y, stellar_z, \
                     stellar_mass, stellar_age, stellar_metallicity, stellar_hsml, \
                     gas_x, gas_y, gas_z, gas_mass, gas_metallicity, gas_hsml, \
                     xrange=0, yrange=0, zrange=0, pixels=720,
                     KAPPA_UNITS=2.08854068444, \
                     IMF_CHABRIER=1, IMF_SALPETER=0 , \
                     ADD_BASE_METALLICITY=0.0, ADD_BASE_AGE=0.0,\
                     BAND_IDS=[9,10,11], \
                     unit_length_mpc=False, \
                     **kwargs ):
    
    Nbands=len(np.array(BAND_IDS)); Nstars=len(np.array(stellar_mass)); Ngas=len(np.array(gas_mass));
    if (Nbands != 3): print "stellar_raytrace needs 3 bands, you gave",Nbands; return -1,-1,-1,-1;
    ## check if stellar metallicity is a matrix
    if (len(stellar_metallicity.shape)>1): stellar_metallicity=stellar_metallicity[:,0];
    if (len(gas_metallicity.shape)>1): gas_metallicity=gas_metallicity[:,0];
    
    ## get opacities and luminosities at frequencies we need:
    stellar_metallicity[stellar_metallicity>0] += ADD_BASE_METALLICITY;
    gas_metallicity[gas_metallicity>0] += ADD_BASE_METALLICITY;
    stellar_age += ADD_BASE_AGE;
    kappa=np.zeros([Nbands]); lums=np.zeros([Nbands,Nstars]);
    for i_band in range(Nbands):
        nu_eff = star_props.colors_table(np.array([1.0]),np.array([1.0]), \
                                   BAND_ID=BAND_IDS[i_band],RETURN_NU_EFF=1);
        kappa[i_band] = atten.opacity_per_solar_metallicity(nu_eff);
        l_m_ssp = star_props.colors_table( stellar_age, stellar_metallicity/0.02, \
                                    BAND_ID=BAND_IDS[i_band], CHABRIER_IMF=IMF_CHABRIER, SALPETER_IMF=IMF_SALPETER, CRUDE=1, \
                                    UNITS_SOLAR_IN_BAND=1); ## this is such that solar-type colors appear white
        l_m_ssp[l_m_ssp >= 300.] = 300. ## just to prevent crazy values here
        l_m_ssp[l_m_ssp <= 0.] = 0. ## just to prevent crazy values here
        lums[i_band,:] = stellar_mass * l_m_ssp
    
    gas_lum=np.zeros(Ngas); ## gas has no 'source term' for this calculation
    stellar_mass_attenuation = np.zeros(Nstars); ## stars have no 'attenuation term'
    gas_mass_metal = gas_mass * (gas_metallicity/0.02);
    if unit_length_mpc:
        kappa *= KAPPA_UNITS / 1000.0 / 1000.0;
    else:
        kappa *= KAPPA_UNITS;


    ## combine the relevant arrays so it can all be fed into the ray-tracing
    x=np.concatenate([stellar_x,gas_x]); y=np.concatenate([stellar_y,gas_y]); z=np.concatenate([stellar_z,gas_z]);
    mass=np.concatenate([stellar_mass_attenuation,gas_mass_metal]);
    hsml=np.concatenate([stellar_hsml,gas_hsml]);
    wt1=np.concatenate([lums[0,:],gas_lum]); wt2=np.concatenate([lums[1,:],gas_lum]); wt3=np.concatenate([lums[2,:],gas_lum]);
    k1=kappa[0]; k2=kappa[1]; k3=kappa[2];



    return raytrace_projection_compute(x,y,z,hsml,mass,wt1,wt2,wt3,k1,k2,k3,\
                                   xrange=xrange,yrange=yrange,zrange=zrange,pixels=pixels,TRIM_PARTICLES=1);


##
## routine to use 'raytrace_projection_compute' to make mock gas images,
##   with three color channels for different temperature ranges
##
def gas_raytrace_temperature( TEMPERATURE_CUTS, \
                             gas_x, gas_y, gas_z, gas_temperature, gas_mass, gas_hsml, \
                             xrange=0, yrange=0, zrange=0, pixels=720,
                             KAPPA_UNITS=2.08854068444, kernel_width=0.05, use_log_t=1 , \
                             isosurfaces=0 , add_temperature_weights=0,
                            **kwargs):

    import scipy.special

    wtfn = gas_mass #* np.sqrt(gas_temperature/1.0e4)
    if (add_temperature_weights==1): wtfn *= np.sqrt(1. + gas_temperature/1.0e4)
    # weighting by sqrt(temp) makes the contributions more similar by temperature bins

    tcuts=TEMPERATURE_CUTS;
    tval=gas_temperature
    if(use_log_t==1):
        tcuts=np.log10(tcuts);
        tval=np.log10(tval);

    if np.array(kernel_width).size > 1:
        w = kernel_width
    else:
        w = kernel_width + np.zeros(3)
    ## continuous smoothing with gaussians for the temperature:
    if (isosurfaces==1):
        wt1 = np.exp(-(tval-tcuts[0])*(tval-tcuts[0])/(2.*w[0]*w[0]));
        wt2 = np.exp(-(tval-tcuts[1])*(tval-tcuts[1])/(2.*w[1]*w[1]));
        wt3 = np.exp(-(tval-tcuts[2])*(tval-tcuts[2])/(2.*w[2]*w[2]));
    else: ## isosurfaces==0, so do total in integral ranges set by temperature_cuts
        wt1 = 0.5*(1.0-scipy.special.erf((tval-tcuts[0])/(np.sqrt(2.)*w[0])));
        wt3 = 0.5*(1.0-scipy.special.erf((tcuts[1]-tval)/(np.sqrt(2.)*w[1])));
        wt2 = 1.-wt1-wt3; wt2[wt2<0.]=0.;

    wt1*= wtfn; wt2*=wtfn; wt3*=wtfn;
    kappa = 200. * (1.+np.zeros((3)));
    kappa *= KAPPA_UNITS;
    print 'KAPPA == ',kappa

    print "gas_temperature:"
    print gas_temperature

    print "tval:"
    print np.min(tval), np.max(tval)
    print tval

    print "tcuts:"
    print tcuts

    print "wt1/2/3:"
    print np.min(wt1), np.max(wt1), np.median(wt1)
    print np.min(wt2), np.max(wt2), np.median(wt2)
    print np.min(wt3), np.max(wt3), np.median(wt3)
    print wt1
    print wt2
    print wt3

    print "xrange/yrange/zrange:"
    print xrange
    print yrange
    print zrange

#    sys.exit()

    return raytrace_projection_compute(gas_x,gas_y,gas_z,gas_hsml,gas_mass,\
                                       wt1,wt2,wt3,kappa[0],kappa[1],kappa[2],\
                                       xrange=xrange,yrange=yrange,zrange=zrange,pixels=pixels,TRIM_PARTICLES=1);






##
##  Wrapper for raytrace_rgb, program which does a simply line-of-sight projection
##    with multi-color source and self-extinction along the sightline: here called
##    from python in its most general form: from the c-code itself:
##
##  int raytrace_rgb(
##    int N_xy, // number of input particles/positions
##    float *x, float *y, // positions (assumed already sorted in z)
##    float *hsml, // smoothing lengths for each
##    float *Mass, // total weight for 'extinction' part of calculation
##    float *wt1, float *wt2, float *wt3, // weights for 'luminosities'
##    float KAPPA1, float KAPPA2, float KAPPA3, // opacities for each channel
##    float Xmin, float Xmax, float Ymin, float Ymax, // boundaries of output grid
##    int Xpixels, int Ypixels, // dimensions of grid
##    float *OUT0, float *OUT1, float *OUT2, float*OUT3 ) // output vectors with final weights
##
def raytrace_projection_compute( x, y, z, hsml, mass, wt1, wt2, wt3, \
                                kappa_1, kappa_2, kappa_3, xrange=0, yrange=0, zrange=0, pixels=720, \
                                TRIM_PARTICLES=1, **kwargs ):
    



    ## define bounaries
    if(util.cast.checklen(xrange)<=1): xrange=[np.min(x),np.max(x)];
    if(util.cast.checklen(yrange)<=1): yrange=[np.min(y),np.max(y)];
    if(util.cast.checklen(zrange)<=1): zrange=[np.min(z),np.max(z)];
    xr=xrange; yr=yrange; zr=zrange;
    print xr, yr, zr
    print xr[1], xr[0]
    x00=0.5*(xr[1]+xr[0]); y00=0.5*(yr[1]+yr[0]); z00=0.5*(zr[1]+zr[0]);
    tolfac = 1.0e10;
    if (TRIM_PARTICLES==1): tolfac = 0.05;
    
    ## clip to particles inside those
    xlen=0.5*(xr[1]-xr[0]); ylen=0.5*(yr[1]-yr[0]); zlen=0.5*(zr[1]-zr[0]);
    x-=x00; y-=y00; z-=z00; dx=xlen*(1.+tolfac*2.); dy=ylen*(1.+tolfac*2.); dz=zlen*(1.+tolfac*2.);
    ok=util.cast.ok_scan(x,xmax=dx) & util.cast.ok_scan(y,xmax=dy) & util.cast.ok_scan(z,xmax=dz) & \
        util.cast.ok_scan(hsml,pos=1) & util.cast.ok_scan(mass+wt1+wt2+wt3,pos=1)
    #& ok_scan(mass) & ok_scan(wt1) & ok_scan(wt2) & ok_scan(wt3);
    x=x[ok]; y=y[ok]; z=z[ok]; hsml=hsml[ok]; mass=mass[ok]; wt1=wt1[ok]; wt2=wt2[ok]; wt3=wt3[ok];
    N_p=util.cast.checklen(x); xmin=-xlen; xmax=xlen; ymin=-ylen; ymax=ylen;
    if(N_p<=1):
        print ' UH-OH: EXPECT ERROR NOW, there are no valid source/gas particles to send!'; return -1,-1,-1,-1;

    ## now sort these in z (this is critical!)
    s=np.argsort(z);
    x=x[s]; y=y[s]; z=z[s]; hsml=hsml[s]; mass=mass[s]; wt1=wt1[s]; wt2=wt2[s]; wt3=wt3[s];
    ## cast new copies to ensure the correct formatting when fed to the c-routine:
    x=util.cast.fcor(x); y=util.cast.fcor(y); z=util.cast.fcor(z); hsml=util.cast.fcor(hsml); mass=util.cast.fcor(mass); wt1=util.cast.fcor(wt1); wt2=util.cast.fcor(wt2); wt3=util.cast.fcor(wt3);

    ## load the routine we need
    exec_call=util.dir.c_routines_dir()+'/RayTrace_RGB/raytrace_rgb.so'
    routine=ctypes.cdll[exec_call];
    
    ## cast the variables to store the results
    aspect_ratio=ylen/xlen; Xpixels=util.cast.int_round(pixels); Ypixels=util.cast.int_round(aspect_ratio*np.float(Xpixels));
    
    print "some diags in raytrace projection "
    print "Xpixels = {:d}".format(Xpixels)
    print "Ypixels = {:d}".format(Ypixels)
    print "aspect_ratio = {:f}".format(aspect_ratio)


    N_pixels=Xpixels*Ypixels; out_cast=ctypes.c_float*N_pixels;
    out_0=out_cast(); out_1=out_cast(); out_2=out_cast(); out_3=out_cast();
    
    ## main call to the calculation routine
    routine.raytrace_rgb( ctypes.c_int(N_p), \
                         util.cast.vfloat(x),   util.cast.vfloat(y),    util.cast.vfloat(hsml), util.cast.vfloat(mass), \
                         util.cast.vfloat(wt1), util.cast.vfloat(wt2),  util.cast.vfloat(wt3), \
                         ctypes.c_float(kappa_1), ctypes.c_float(kappa_2), ctypes.c_float(kappa_3), \
                         ctypes.c_float(xmin), ctypes.c_float(xmax), ctypes.c_float(ymin), ctypes.c_float(ymax), \
                         ctypes.c_int(Xpixels), ctypes.c_int(Ypixels), \
                         ctypes.byref(out_0), ctypes.byref(out_1), ctypes.byref(out_2), ctypes.byref(out_3) );
                         
    ## now put the output arrays into a useful format
    out_0 = np.copy(np.ctypeslib.as_array(out_0));
    out_1 = np.copy(np.ctypeslib.as_array(out_1));
    out_2 = np.copy(np.ctypeslib.as_array(out_2));
    out_3 = np.copy(np.ctypeslib.as_array(out_3));
    out_0 = out_0.reshape([Xpixels,Ypixels]);
    out_1 = out_1.reshape([Xpixels,Ypixels]);
    out_2 = out_2.reshape([Xpixels,Ypixels]);
    out_3 = out_3.reshape([Xpixels,Ypixels]);

    return out_0, out_1, out_2, out_3;
