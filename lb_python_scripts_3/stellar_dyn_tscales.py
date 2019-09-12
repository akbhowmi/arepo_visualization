import numpy as np
import matplotlib.pyplot as plt
import matplotlib

G = 6.67259D-8
MSUN = 1.989D33
PC = 3.0857D18
YR = 3.1556936D7
C = 2.99792458D10
MP = 1.6726231D-24
KB = 1.380658D-16
SB = 5.67051D-5
KES = 0.4
THOMSON = 6.65245D-25
LSUN = 3.83D33



def t_gasdrag_subsonic(r):

    rho_g = foo
    csound = foo

    tscale = 3.0/(4*np.pi) * csound^3 / (G * rho_g * mstar) / YR

    return tscale


def t_gasdrag_supersonic(mach, r, t_subsonic):

    fgeom = foo
    coulomb = foo
    
    tscale = mach^3 * fgeom / coulomb * t_subsonic

    return tscale

