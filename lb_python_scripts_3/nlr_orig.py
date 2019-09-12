import pyfits, numpy,pdb, sys
from numpy import *

def read_hbeta_file(file, snap):
    print "Reading file %s"%file
    sys.stdout.flush()
    f=open(file,'r')

    while(True):
        l=f.readline().split()
        n=int(l[-1])
        d=[]
        for i in range(n):
            d.append(f.readline().split())

        if int(l[0])==snap:
            break

    dd=[]
    for x in d:
        dd.append([float(xx) for xx in x])
    hbf=numpy.array(dd)

    return hbf

def make_hbeta_sed(lum, lam):
    """lum is hbeta lum in Lsun."""
    
    hbetai=numpy.where(lam>4861e-10)[0][0]+1

    sed=zeros(lam.shape[0])
    deltalam=0.5*(lam[hbetai+1]-lam[hbetai-1])
    
    sed[hbetai] = 3.83e26*lum/deltalam

    return sed

def add_at_bottom(original, add):
    s=(original.shape[0]+add.shape[0],original.shape[1])
    assert original.shape[1]==add.shape[1]
    new=numpy.empty(shape=s)
    new[:original.shape[0],:]=original
    new[original.shape[0]:,:]=add
    return new
    
def add_hbeta_sources(velfile, otherfile, sfrhistfile, snapnum):
    """Reads Laura's velfile and otherfile and gets the Hbeta
    luminosity for the NLR particles. Then adds these particles as
    pure Hbeta sources to the specified sfrhist file."""
    
    velf=read_hbeta_file(velfile, snapnum)
    otherf=read_hbeta_file(otherfile, snapnum)

    hbetalum=velf[:,-1]
    print "Total Hbeta luminosity: %g Lsun = %g W"%(sum(hbetalum),sum(hbetalum)*3.83e26)

    print "Reading data for existing sources"
    sys.stdout.flush()
    sf=pyfits.open(sfrhistfile)

    phdu=sf['PARTICLEDATA']

    nold=phdu.data.shape[0]
    nadd = velf.shape[0]

    lam=sf['LAMBDA'].data.field('lambda')
    #pdb.set_trace()
    sed=phdu.data.field('L_lambda')
    pos=phdu.data.field('position')
    vel=phdu.data.field('velocity')
    age=phdu.data.field('age')

    # figure out which BH is which. Laura's file contains the pos of
    # the particle in ref to the BH, so by subtracting those vectors
    # we get - delta pos of the black holes. We want the position in
    # relation to the black hole that's the first column in her files.
    bhpart=where(age!=age)[0]

    if bhpart.shape[0]==1:
        # if we only have one then it's easy
        bhpos=pos[bhpart[0],:]
        # check that the other bh has junk
        assert otherf[0,3]==-1
    else:
        bhdeltapos=pos[bhpart[0],:]-pos[bhpart[1],:]
        lbhdeltapos=otherf[0,0:3]-otherf[0,3:6]
        #if these vectors are parallel we have the BHs switched
        if dot(bhdeltapos, lbhdeltapos)>0:
            bhpos=pos[bhpart[1],:]
        else:
            bhpos=pos[bhpart[0],:]

    print "Generating new PARTICLEDATA table"
    sys.stdout.flush()

    # create new table hdu. this also copies the data
    hdu = pyfits.new_table(phdu.columns, header=phdu.header, nrows=nold+nadd)

    hbetased = array([make_hbeta_sed(hbetalum[i],lam) for i in range(nadd)])
    print "Integrated Hbeta L_lambda: %g W/m"%(sum(hbetased))
                      
    # append new data
    hdu.data.field('position')[nold:] = otherf[:,0:3]+bhpos
    hdu.data.field('velocity')[nold:] = velf[:,1:4]*1.0226903e-09
    hdu.data.field('L_lambda')[nold:] = log10(hbetased+1e-300)
    hdu.data.field('L_bol')[nold:] = otherf[:,-1]*3.83e26

    # calculate particle size from density and mass
    h=sf['gadget'].header['HubbleParam']
    mtomsun=sf['gadget'].header['UnitMass_in_g']*5.0273993e-34/h
    ltokpc=sf['gadget'].header['UnitLength_in_cm']*3.2407793e-22/h
    
    mass=otherf[:,8]*mtomsun
    density = otherf[:,6]*mtomsun/ltokpc**3
    
    hdu.data.field('radius')[nold:] = (3*mass/(4*pi*density))**(1./3)

    # fudge up new data
    hdu.data.field('mass')[nold:] = 0
    hdu.data.field('metallicity')[nold:] = 0
    hdu.data.field('formation_time')[nold:] = nan
    hdu.data.field('parent_ID')[nold:] = 0
    hdu.data.field('age')[nold:] = nan
    hdu.data.field('creation_mass')[nold:] = 0.
    sf.close

    print "Updating file %s"%sfrhistfile
    sys.stdout.flush()

    pyfits.update(sfrhistfile, hdu.data, hdu.header, 'PARTICLEDATA')
    

    
    
    

    

    
