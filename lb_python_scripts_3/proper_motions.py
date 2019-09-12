import numpy as np
import matplotlib.pyplot as plt


#import PyCos as pc
from astropy import units as u
from astropy.coordinates import SkyCoord
#from astropy.cosmology import FlatLambdaCDM as cosmo
from astropy.cosmology import WMAP9 as cosmo
from astropy import constants as const


def vel_kms_to_masyr(vel = 1000.0, z=1.0):

    vel = vel * u.kilometer / u.second
    
    vel_mas_per_year = (vel.to(u.kpc/u.yr) *
                        cosmo.arcsec_per_kpc_proper(z).to(u.mas/u.kpc) )
    if z<0.1:
        d=cosmo.angular_diameter_distance(z).to(u.kpc)
        loc_vel_mas_per_year = ((vel.to(u.kpc/u.yr)/d)*(u.rad)).to(u.mas/u.yr)

    #table 6 says 0.7 mas is max resolution, for 93GHz band

    print("D_A (z=%g) = "%z, cosmo.angular_diameter_distance(z))
    print("D_L (z=%g) = "%z, cosmo.luminosity_distance(z))
    print("arcsec/kpc (z=%g) ="%z, cosmo.arcsec_per_kpc_proper(z))
    print("mas/kpc (z=%g) ="%z, cosmo.arcsec_per_kpc_proper(z).to(u.mas/u.kpc))
    print("physical res for 10 mas (z=%g) = "%z, 10*u.mas/(cosmo.arcsec_per_kpc_proper(z).to(u.mas/u.pc)))
    print("vel [km/s] = %g"%vel.to_value())
    print("vel [km/s] =",vel)
    print("vel [pc/Myr] =",vel.to(u.pc/u.Myr))
    print(vel)
    print("vel [kpc/yr] =",vel.to(u.kpc/u.yr))
    print("vel [mas/yr] = ",vel_mas_per_year)
    if z<0.1: print("loc vel [mas/yr] = ",loc_vel_mas_per_year)
    print("motion in 10 yr lifetime = ", vel_mas_per_year*10*u.yr)

    return vel_mas_per_year

