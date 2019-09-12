import numpy as np
import h5py
import sys, os, glob, re
import astro_constants as ac
from copy import copy

UnitMass_in_g = 1.989e+43 
UnitTime_in_s = 3.08568e+16 
UnitTime_in_Gyr = UnitTime_in_s / (1.0e9*ac.YR)
UnitVelocity_in_cm_per_s = 100000 
UnitDensity_in_cgs = 6.76991e-22 
UnitEnergy_in_cgs = 1.989e+53 


class snapshot_header:

    def __init__(self, *args, **kwargs):
        if len(args)!=1:
            print "Error: class snapshot_header takes exactly 1 argument (filename)."
            sys.stdout.flush()
            sys.exit()

        fname = args[0]

        if os.path.exists(fname):
            filename = fname
        elif os.path.exists(fname+'.hdf5'):
            filename = fname+'.hdf5'
        else: 
            print "Error: file %s not found."%fname
            sys.stdout.flush()
            sys.exit()
            
        f=h5py.File(filename,'r')

        hdr=f['Header'].attrs 
        ##print hdr.items()
        self.npart = hdr['NumPart_Total']
        self.time = hdr['Time']
        self.h = hdr['HubbleParam']

        f.close()


class snapshot_bhs:

    def calc_mdotedd(self, mbh, reff_edd=0.1):
        return mbh * 4*np.pi*ac.G*ac.MP*ac.YR / (ac.THOMSON*reff_edd*ac.C) 

    def calc_reff(self, mdot, mdotedd, eddval=0.1, riaf=True):

        reff = eddval + np.zeros((mdot.size))
        if riaf:
            fac = mdot/(0.01*mdotedd)
            reff[fac<1] =  reff[fac<1] * fac[fac<1]
        return reff

    def __init__(self, *args, **kwargs):

        if len(args)!=1:
            print "Error: class snapshot_bhs takes exactly 1 argument (filename)."
            sys.stdout.flush()
            sys.exit()

        fname = args[0]

        riaf = kwargs.get('riaf',True)
        reff_edd = kwargs.get('reff_edd',0.1)        

        if os.path.exists(fname):
            filename = fname
        elif os.path.exists(fname+'.hdf5'):
            filename = fname+'.hdf5'
        else: 
            print "Error: file %s not found."%fname
            sys.stdout.flush()
            sys.exit()            

        f = h5py.File(filename,'r')

        self.bh = f['PartType5']
    
        hubble = f['Header'].attrs['HubbleParam'] 

        self.pos = np.array(self.bh['Coordinates'])/hubble
        self.mdot = np.array(self.bh['BH_Mdot'])*(UnitMass_in_g/ac.MSUN/UnitTime_in_s*ac.YR)
        self.mass = np.array(self.bh['BH_Mass'])/hubble*UnitMass_in_g/ac.MSUN
        self.dynmass = np.array(self.bh['Masses'])/hubble*UnitMass_in_g/ac.MSUN
        self.ids = np.array(self.bh['ParticleIDs']).astype('uint64')

        self.mdot_edd = self.calc_mdotedd(self.mass, reff_edd=reff_edd)
        self.rad_eff = self.calc_reff(self.mdot,self.mdot_edd, riaf=riaf, eddval=reff_edd)

        self.lbol = self.rad_eff * self.mdot*ac.MSUN/ac.YR * ac.C * ac.C 
        self.lbol_edd = reff_edd * self.mdot_edd*ac.MSUN/ac.YR * ac.C * ac.C

        #print " "
        #for grp in f.keys()[1:]:        
        #    print grp
        #    print f[grp].keys()
        #    print " "

        f.close()


class all_snapshot_bhs:

    def __init__(self, *args, **kwargs):

        if len(args)!=1:
            #print "Error: class snapshot_bhs takes exactly 1 argument (array of snapshot filenames)."
            print "Error: class snapshot_bhs takes exactly 1 argument (snapshot directory path)."
            sys.stdout.flush()
            sys.exit()
        path = args[0]

        if not os.path.exists(path):
            print "Error: file %s not found."%path
            sys.stdout.flush()
            sys.exit()            

        tmin = kwargs.get('tmin',0.0)
        tmax = kwargs.get('tmax',15.0)
        tmax_postmrg = kwargs.get('tmax_postmrg',-1.0)

        self.snapfiles = np.array(glob.glob(path+'/snapshot*hdf5'))
        self.snaps = np.array([ np.int(re.findall('(\d+)',s.strip('hdf5'))[-1])
                                for s in self.snapfiles ])
        self.nsnaps = len(self.snaps)
        if self.nsnaps<=0:
            print "Error: no snapshot files passed to class all_snapshot_bhs."
            sys.stdout.flush()
            sys.exit()            
        ix_sort = self.snaps.argsort(kind='mergesort')
        self.snaps = self.snaps[ix_sort]
        self.snapfiles = self.snapfiles[ix_sort]

        self.time = np.repeat(np.nan, self.nsnaps)
        self.nbh = copy(self.time)

        self.pos = np.repeat(np.nan,2*3*self.nsnaps).reshape(2,3,self.nsnaps)
        self.mdot = np.repeat(np.nan,2*self.nsnaps).reshape(2,self.nsnaps)
        self.mass = copy(self.mdot)
        self.dynmass = copy(self.mdot)
        self.ids = np.zeros((2,self.nsnaps)).astype('uint64')

        self.mdot_edd = copy(self.mdot)
        self.rad_eff = copy(self.mdot)
        self.lbol = copy(self.mdot)
        self.lbol_edd = copy(self.mdot)        

        for (i,snfile) in enumerate(self.snapfiles):

            head = snapshot_header(snfile)        
            self.time[i] = head.time*UnitTime_in_Gyr/head.h
            #nbh = head.npart[5]
            #print "snap: %d, time: %g Gyr, nbh: %d"%(snaps[i],t,nbh)
            #assert 0<nbh<=2, "Error: found %d BHs in snap."%nbh                

            #print "snap %i"%i
            bhs = snapshot_bhs(snfile)
            self.nbh[i] = len(bhs.mdot)

            self.pos[:self.nbh[i],:,i] = bhs.pos
            self.mdot[:self.nbh[i],i] = bhs.mdot
            self.mass[:self.nbh[i],i] = bhs.mass
            self.dynmass[:self.nbh[i],i] = bhs.dynmass
            self.ids[:self.nbh[i],i] = bhs.ids

            self.mdot_edd[:self.nbh[i],i] = bhs.mdot_edd
            self.rad_eff[:self.nbh[i],i] = bhs.rad_eff
            self.lbol[:self.nbh[i],i] = bhs.lbol
            self.lbol_edd[:self.nbh[i],i] = bhs.lbol_edd
            
        self.tmrg = self.time[self.nbh==1].min() if self.snaps[self.nbh==1].size<self.snaps.size else -1
        if tmax_postmrg >= 0 and self.tmrg >= 0:
            mask = ((self.time>=tmin)&(self.time<=tmax)&(self.time<=self.tmrg+tmax_postmrg))
        else:
            mask = ((self.time>=tmin)&(self.time<=tmax))
            if tmax_postmrg >=0: print "Warning: No BH merger found. Ignoring keyword tmax_postmrg."
        if self.snaps[mask].size < self.snaps.size:
            self.snapfiles = self.snapfiles[mask]
            self.snaps = self.time[mask]
            self.time = self.time[mask]
            self.nbh = self.nbh[mask]
            self.pos = self.pos[:,:,mask]
            self.mdot = self.mdot[:,mask]
            self.mass = self.mass[:,mask]
            self.dynmass = self.dynmass[:,mask]
            self.ids = self.ids[:,mask]
            self.mdot_edd = self.mdot_edd[:,mask]
            self.rad_eff = self.rad_eff[:,mask]
            self.lbol = self.lbol[:,mask]
            self.lbol_edd = self.lbol_edd[:,mask]
        
        
def read_header(path,snapnum):

    head = snapshot_header(path+'/snapshot_%.3d.hdf5'%snapnum)

    print "\nTime: %g [Gyr/h] (%g Gyr)"%(head.time*UnitTime_in_Gyr,
                                         head.time*UnitTime_in_Gyr/head.h)
    print "NumPart_Total:",head.npart
    

def snap(path,snapnum):

    head = snapshot_header(path+'/snapshot_%.3d.hdf5'%snapnum)

    print "\nTime: %g [Gyr/h] (= %g Gyr)"%(head.time,head.time/head.h)
    print "nBH:",head.npart[5]
    assert 0<head.npart[5]<=2, "Error: found %d BHs in snap."%head.npart[5]

    bhs = snapshot_bhs(path+'/snapshot_%.3d.hdf5'%snapnum)
    
    pos = np.array(bh['Coordinates'])
    mdot = np.array(bh['BH_Mdot'])
    mass = np.array(bh['BH_Mass'])
    dynmass = np.array(bh['Masses'])
    ids = np.array(bh['ParticleIDs'])

    print pos.shape, mdot.shape, mass.shape, dynmass.shape, ids.shape
    return time, ids, mdot, mass, dynmass, pos
    #return np.array(f['PartType5']['BH_Mass'])/h*1.0e10


def write_all(path,return_data=False):

    #get file list
    #snapfiles = np.array(glob.glob(path+'/snapshot*hdf5'))
    #print snapfiles.shape
    #snaps = np.array([ np.int(re.findall('(\d+)',s.strip('hdf5'))[-1])
    #                   for s in snapfiles ])
    #print snapfiles.size, snaps.size
    #ix_sort = snaps.argsort(kind='mergesort')
    #snaps = snaps[ix_sort]
    #snapfiles = snapfiles[ix_sort]


    allbhs = all_snapshot_bhs(path)

    with open(path+'/bh_snap_info.txt','w') as fp:

        for i,snap in enumerate(allbhs.snaps):

            #head = snapshot_header(snfile)        
            #t = head.time*UnitTime_in_Gyr/head.h
            #nbh = head.npart[5]
            ##print "snap: %d, time: %g Gyr, nbh: %d"%(snaps[i],t,nbh)
            #assert 0<nbh<=2, "Error: found %d BHs in snap."%nbh

            pstr1 = " ".join([str(p) for p in allbhs.pos[0,:,i]])
            pstr2 = " ".join([str(p) for p in allbhs.pos[1,:,i]])
            fp.write(("%d %g %d "+2*"%s %g %g "+"\n")%(snap, allbhs.time[i], allbhs.nbh[i],
                                                       pstr1, allbhs.lbol[0,i], allbhs.mass[0,i], 
                                                       pstr2, allbhs.lbol[1,i], allbhs.mass[1,i]))
                
    #print allbhs.mass.shape
    #print allbhs.mass

    if return_data: return allbhs

