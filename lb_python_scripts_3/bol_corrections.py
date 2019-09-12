import numpy as np
import sys
import scipy.optimize as opt
import astro_constants as ac


def calc_bc_from_lbol(bc_model,lbol):
    
    log_lbol_12 = np.log10(lbol/ac.LSUN)-12.0

    if bc_model>6:
        # constant BC value
        bc = bc_model
    elif bc_model==1:
        # Marconi et al. 2004 B-band formula for BC:
        bc = 10**( 0.80 - 0.067*log_lbol_12 + 0.017*log_lbol_12*log_lbol_12 - 
                   0.0023*log_lbol_12*log_lbol_12*log_lbol_12)
    elif bc_model==2:
        # Marconi et al. 2004 0.5-2keV formula for BC: 
        bc = 10**( 1.65 + 0.22*log_lbol_12 + 0.012*log_lbol_12*log_lbol_12 - 
                   0.0015*log_lbol_12*log_lbol_12*log_lbol_12)
    elif bc_model==3:
        # Marconi et al. 2004 2-10keV formula for BC:
        bc = 10**( 1.54 + 0.24*log_lbol_12 + 0.012*log_lbol_12*log_lbol_12 - 
                   0.0015*log_lbol_12*log_lbol_12*log_lbol_12)
    elif bc_model==4:
        # Hopkins et al. 2007 B-band formula for BC:
        bc = 6.25*(lbol/(1.0e10*ac.LSUN))**-0.37 + 9.00*(lbol/(1.0e10*ac.LSUN))**-0.012
    elif bc_model==5:
        # Hopkins et al. 2007 0.5-2keV formula for BC:
        bc = 17.87*(lbol/(1.0e10*ac.LSUN))**0.28 + 10.03*(lbol/(1.0e10*ac.LSUN))**-0.02
    elif bc_model==6:
        # Hopkins et al. 2007 2-10keV formula for BC:
        bc = 10.83*(lbol/(1.0e10*ac.LSUN))**0.28 + 6.08*(lbol/(1.0e10*ac.LSUN))**-0.02
    else:
        print "Error: bc model %d undefined."%bc_model
        sys.exit()

    return bc


def invbc_m04(loglbol, *data):
    # loglband must be in units log10(lband/LSUN)
    a0,a1,a2,a3,loglband = data
    return (a0 + a1*(loglbol-12.0) + a2*(loglbol-12.0)**2 + 
            a3*(loglbol-12.0)**3 - loglbol + loglband)

def invbc_h07(lbol, *data):
    # lband must be in units of 1e10 LSUN
    c1,k1,c2,k2,lband = data
    return c1*(lbol**k1) + c2*(lbol**k2) - lbol/lband


def calc_bc_from_lband(bc_model,lband):
    
    if bc_model>6:
        # const value for bc
        bc = bc_model

    elif bc_model==1:
        func_args = (0.80,-0.067,0.017,-0.0023,np.log10(lband/ac.LSUN))
        loglbol = opt.fsolve(invbc_m04,np.log10(10*lband),args=func_args)
        bc = (10**loglbol)*ac.LSUN/lband

    elif bc_model==2:
        func_args = (1.65,0.22,0.012,-0.0015,np.log10(lband/ac.LSUN))
        loglbol = opt.fsolve(invbc_m04,np.log10(10*lband),args=func_args)
        bc = (10**loglbol)*ac.LSUN/lband

    elif bc_model==3:
        func_args = (1.54,0.24,0.012,-0.0015,np.log10(lband/ac.LSUN))
        loglbol = opt.fsolve(invbc_m04,np.log10(10*lband),args=func_args)
        bc = (10**loglbol)*ac.LSUN/lband

    elif bc_model==4:
        func_args = (6.25,-0.37,9.00,-0.012,lband/(1.0e10*ac.LSUN))
        lbol = opt.fsolve(invbc_h07,10*lband/(1.0e10*ac.LSUN),args=func_args)
        bc = lbol*1.0e10*ac.LSUN/lband

    elif bc_model==5:
        func_args = (17.87,0.28,10.03,-0.02,lband/(1.0e10*ac.LSUN))
        lbol = opt.fsolve(invbc_h07,10*lband/(1.0e10*ac.LSUN),args=func_args)
        bc = lbol*1.0e10*ac.LSUN/lband

    elif bc_model==6:
        func_args = (10.83,0.28,6.08,-0.02,lband/(1.0e10*ac.LSUN))
        lbol = opt.fsolve(invbc_h07,10*lband/(1.0e10*ac.LSUN),args=func_args)
        bc = lbol*1.0e10*ac.LSUN/lband

    else: 
        print "Error: bc_model %d undefined."%bc_model
        sys.exit()

    return bc
