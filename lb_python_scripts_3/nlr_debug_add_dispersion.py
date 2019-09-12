import pyfits, numpy,pdb, sys, gc
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

    f.close()

    return hbf

def make_hbeta_sed(lum, lam):
    """lum is hbeta lum in Lsun."""
    
    hbetai=numpy.where(lam>4861e-10)[0][0]+1

    sed=zeros(lam.shape[0])
    deltalam=0.5*(lam[hbetai+1]-lam[hbetai-1])
    
    sed[hbetai] = 3.83e26*lum/deltalam

    return sed

def make_hbeta_sed_DEBUG_add_dispersion(lum, lam, csnd):
    """lum is hbeta lum in Lsun."""
    
    hbetai=numpy.where(lam>4861e-10)[0][0]+1

    sed=zeros(lam.shape[0])
    #deltalam=0.5*(lam[hbetai+1]-lam[hbetai-1])

    deltalam = zeros(lam.shape[0])
    deltalam[0] = 0.5*(lam[2]-lam[0])
    for i in range(1,lam.shape[0]-1):
        deltalam[i] = 0.5*(lam[i+1]-lam[i-1])
    deltalam[0] = deltalam[1]
    deltalam[-1] = deltalam[-2]

    #dispersion is a gaussian with FWHM = dlambda = lambda_0*v/c
    sigma = lam[hbetai] * (0.5*csnd/2.99792458e+5) / (2.0*sqrt(2*log(2)))
    #print 'sigma = ',sigma
    sed = 3.83e26*lum * (1.0/(sqrt(2*pi)*sigma))*exp(-(lam-lam[hbetai])**2/(2.0*sigma**2))
    #print 'min/max sed=',min(sed),max(sed)
    
    return sed

def add_at_bottom(original, add):
    s=(original.shape[0]+add.shape[0],original.shape[1])
    assert original.shape[1]==add.shape[1]
    new=numpy.empty(shape=s)
    new[:original.shape[0],:]=original
    new[original.shape[0]:,:]=add
    return new

def memory_usage():
    """Memory usage of the current process in kilobytes."""
    status = None
    result = {'peak': 0, 'rss': 0}
    try:
        # This will only work on systems with a /proc file system
        # (like Linux).
        status = open('/proc/self/status')
        for line in status:
            parts = line.split()
            key = parts[0][2:-1].lower()
            if key in result:
                result[key] = int(parts[1])
    finally:
        if status is not None:
            status.close()
    return result
    
def add_hbeta_sources_DEBUG_ADD_DISP(dir, model_fpath, sfrhist_fpath, sfrhist_fname, snapnum):
    """Reads Laura's velfile and otherfile and gets the Hbeta
    luminosity for the NLR particles. Then adds these particles as
    pure Hbeta sources to the specified sfrhist file."""

    print "Current/peak memory:", memory_usage()
    gc.get_referrers()

    velfile='LHb_data_vel.dat'
    otherfile='LHb_data_other.dat'

    vpath=dir+'/'+model_fpath+'/'+velfile
    opath=dir+'/'+model_fpath+'/'+otherfile
    
    print dir, model_fpath, velfile, otherfile
    print vpath
    print opath
    
    velf=read_hbeta_file(vpath, snapnum)
    otherf=read_hbeta_file(opath, snapnum)

    hbetalum=velf[:,-1]
    print "Total Hbeta luminosity: %g Lsun = %g W"%(sum(hbetalum),sum(hbetalum)*3.83e26)

    csnd=velf[:,-2]  # local sound speed around NL particle, in km/s

    sfrhistfile = dir+'/'+sfrhist_fpath+'/'+sfrhist_fname

    print "Current/peak memory:", memory_usage()

    print "Reading data for existing sources"
    sys.stdout.flush()
    sf=pyfits.open(sfrhistfile,memmap=True)

    print "Current/peak memory:", memory_usage()

    phdu=sf['PARTICLEDATA']

    print "getting npart..."
    sys.stdout.flush()

    nold=phdu.data.shape[0]
    nadd = velf.shape[0]

    print "getting particle data..."
    sys.stdout.flush()

    lam=sf['LAMBDA'].data.field('lambda')
    #pdb.set_trace()
    sed=phdu.data.field('L_lambda')
    pos=phdu.data.field('position')
    vel=phdu.data.field('velocity')
    age=phdu.data.field('age')


    print "getting BH info..."
    sys.stdout.flush()

    # figure out which BH is which. Laura's file contains the pos of
    # the particle in ref to the BH, so by subtracting those vectors
    # we get - delta pos of the black holes. We want the position in
    # relation to the black hole that's the first column in her files.
    # (Laura's note: last particles in unmodified fits file are BHs, and
    # have age=nan.
    bhpart=where(age!=age)[0]

    if bhpart.shape[0]==1:
        # if we only have one then it's easy
        bhpos=pos[bhpart[0],:]
        # check that the other bh has junk
        assert otherf[0,3]==-1
        nbh = 1
    else:
        bhdeltapos=pos[bhpart[0],:]-pos[bhpart[1],:]
        lbhdeltapos=otherf[0,0:3]-otherf[0,3:6]
        #if these vectors are parallel we have the BHs switched
        nbh = 2
        if dot(bhdeltapos, lbhdeltapos)>0:
            bhpos=pos[bhpart[1],:]
        else:
            bhpos=pos[bhpart[0],:]
        

    print "Generating new PARTICLEDATA table"
    sys.stdout.flush()

    # create new table hdu. this also copies the data
    hdu = pyfits.new_table(phdu.columns, header=phdu.header, nrows=nold+nadd)

    print "ID: ",shape(hdu.data.field('ID'))
    print "position: ",shape(hdu.data.field('position'))
    print "mass: ",shape(hdu.data.field('mass'))

    print "nold-nbh=",nold-nbh
    print "nold+nadd-nbh=",nold+nadd-nbh
    
    # move BH data to end of table
    hdu.data.field('ID')[nold+nadd-nbh:] = hdu.data.field('ID')[nold-nbh:nold]
    hdu.data.field('position')[nold+nadd-nbh:] = hdu.data.field('position')[nold-nbh:nold]
    hdu.data.field('velocity')[nold+nadd-nbh:] = hdu.data.field('velocity')[nold-nbh:nold]
    hdu.data.field('L_lambda')[nold+nadd-nbh:] = hdu.data.field('L_lambda')[nold-nbh:nold]
    hdu.data.field('L_bol')[nold+nadd-nbh:] = hdu.data.field('L_bol')[nold-nbh:nold]
    hdu.data.field('radius')[nold+nadd-nbh:] = hdu.data.field('radius')[nold-nbh:nold]
    hdu.data.field('mass')[nold+nadd-nbh:] = hdu.data.field('mass')[nold-nbh:nold]
    hdu.data.field('metallicity')[nold+nadd-nbh:] = hdu.data.field('metallicity')[nold-nbh:nold]
    hdu.data.field('formation_time')[nold+nadd-nbh:] = hdu.data.field('formation_time')[nold-nbh:nold]
    hdu.data.field('parent_ID')[nold+nadd-nbh:] = hdu.data.field('parent_ID')[nold-nbh:nold]
    hdu.data.field('age')[nold+nadd-nbh:] = hdu.data.field('age')[nold-nbh:nold]
    hdu.data.field('creation_mass')[nold+nadd-nbh:] = hdu.data.field('creation_mass')[nold-nbh:nold]
 
    print "ID: ",shape(hdu.data.field('ID'))
    print "position: ",shape(hdu.data.field('position'))
    print "mass: ",shape(hdu.data.field('mass'))


    hbetased = array([make_hbeta_sed_DEBUG_add_dispersion(hbetalum[i],lam,csnd[i]) for i in range(nadd)])
    print "Integrated Hbeta L_lambda: %g W/m"%(sum(hbetased))

    print "hbetased:",hbetased.shape
    ised_hbeta = array([sum(hbetased[:,i]) for i in range(lam.size)])
    print "ised_hbeta:",ised_hbeta.shape
    
    # add new data before BH particles
    hdu.data.field('position')[nold-nbh:nold+nadd-nbh] = otherf[:,0:3]+bhpos
    hdu.data.field('velocity')[nold-nbh:nold+nadd-nbh] = velf[:,1:4]*1.0226903e-09
    #hdu.data.field('L_lambda')[nold-nbh:nold+nadd-nbh] = log10(hbetased+1e-30)
    #TEMPORARY, DEBUG ONLY:
    #hdu.data.field('L_lambda')[nold-nbh:nold+nadd-nbh] = log10(hbetased+1.0e+30)
    hdu.data.field('L_lambda')[nold-nbh:nold+nadd-nbh] = log10(hbetased+1.0e-5)
    #hdu.data.field('L_lambda')[nold-nbh:nold+nadd-nbh] = log10(hbetased+1.0e+10)
    #hdu.data.field('L_lambda')[nold-nbh:nold+nadd-nbh] = log10(hbetased+1.0e+20)
    hdu.data.field('L_bol')[nold-nbh:nold+nadd-nbh] = velf[:,-1]*3.83e26


    # calculate particle size from density and mass
    h=sf['gadget'].header['HubbleParam']
    mtomsun=sf['gadget'].header['UnitMass_in_g']*5.0273993e-34/h
    ltokpc=sf['gadget'].header['UnitLength_in_cm']*3.2407793e-22/h
    
    mass=otherf[:,8]*mtomsun
    density = otherf[:,6]*mtomsun/ltokpc**3

    ### try to get rid of some stuff we're done with:
    print "before del: Current/peak memory:", memory_usage()
    del velf
    del otherf    
    del hbetased
    print "after del: Current/peak memory:", memory_usage()

    hdu.data.field('radius')[nold-nbh:nold+nadd-nbh] = (3*mass/(4*pi*density))**(1./3)

    # fudge up new data

    #hdu.data.field('ID')[nold-nbh:nold+nadd-nbh] = 0 #0 only for NLs
    # new crazy scheme: assign fake, continguous IDs to the NL particles. 
    hdu.data.field('ID')[nold-nbh:nold+nadd-nbh] = hdu.data.field('ID')[nold-nbh-1]+range(1,nadd+1)     #0 only for NLs
    hdu.data.field('mass')[nold-nbh:nold+nadd-nbh] = 0     #0 only for NLs
    #hdu.data.field('mass')[nold-nbh:nold+nadd-nbh] = mass     #0 only for NLs
    hdu.data.field('metallicity')[nold-nbh:nold+nadd-nbh] = 0     #0 for BHs
    hdu.data.field('formation_time')[nold-nbh:nold+nadd-nbh] = 0     #0, not nan, for BHs
    hdu.data.field('parent_ID')[nold-nbh:nold+nadd-nbh] = 0     #0 for BHs
    #hdu.data.field('parent_ID')[nold-nbh:nold+nadd-nbh] = hdu.data.field('ID')[nold-nbh-1]+range(1,nadd+1)     #0 for BHs
    #hdu.data.field('age')[nold-nbh:nold+nadd-nbh] = nan     #nan for BHs
    hdu.data.field('age')[nold-nbh:nold+nadd-nbh] = 0     #nan for BHs
    hdu.data.field('creation_mass')[nold-nbh:nold+nadd-nbh] = 0.    #0 for BHs

    # update the integrated sed
    ihdu = sf['INTEGRATED_QUANTITIES']
    ised = ihdu.data.field('L_lambda')
    print "lg_ised",ised.shape,ised.dtype
    print "lg_ised[101]=",ised[101]
    print "ised_hbeta[101]=",ised_hbeta[101]    
    #ised_new = ised+ised_hbeta+1e-30
    # TEMPORARY; DEBUG ONLY:
    #ised_new = ised+ised_hbeta+1e+30
    ised_new = ised+ised_hbeta+1e-5
    #ised_new = ised+ised_hbeta+1e+10
    #ised_new = ised+ised_hbeta+1e+20
    print "lg_ised_new",ised_new.shape
    print "lg_ised_new[101]", ised_new[101]
    ihdu.data.field('L_lambda')[:] = ised_new

    #sf.close()

    print "Current/peak memory:", memory_usage()

    # add keyword with NLR particle info
    hdu.header.update('NNLR',nadd,'number of NLR particles added',after='LOGFLUX')
    hdu.header.update('NLRDIR',model_fpath,'NLR data subdirectory',after='LOGFLUX')

    print "Updating PARTICLEDATA hdu in file %s"%sfrhistfile
    sys.stdout.flush()

    print "Current/peak memory:", memory_usage()

    pyfits.update(sfrhistfile, hdu.data, hdu.header, 'PARTICLEDATA')

    del hdu.data
    
    print "Updating INTEGRATED_QUANTITIES hdu in file %s"%sfrhistfile
    sys.stdout.flush()

    print "Current/peak memory:", memory_usage()

    pyfits.update(sfrhistfile, ihdu.data, ihdu.header, 'INTEGRATED_QUANTITIES')

    del ihdu.data
    
    print "Current/peak memory:", memory_usage()

    sf.close()

    print "Current/peak memory:", memory_usage()

    gc.collect()
    gc.collect()
    
    print "Current/peak memory:", memory_usage()
    gc.get_referrers()

    return
    

    


if __name__ == "__main__":
    add_hbeta_sources_DEBUG_ADD_DISP(sys.argv[1],
                                     sys.argv[2], 
                                     sys.argv[3],
                                     sys.argv[4],
                                     long(sys.argv[5]))
