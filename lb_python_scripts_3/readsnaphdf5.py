import numpy as np
import h5py


def readheader(path,snapnum):

    fname=path+"/snapshot_%.3d.hdf5"%snapnum
    f=h5py.File(fname,'r')
    
    hdr=f['Header'].attrs 
    #print hdr.items()
    time=hdr['Time']
    h=hdr['HubbleParam']
    npart=hdr['NumPart_Total']

    print "\nTime: %g [Gyr/h] (= %g Gyr)"%(time,time/h)
    print "NumPart_Total:",npart
    
    #print " "
    #for grp in f.keys()[1:]:        
    #    print grp
    #    print f[grp].keys()
    #    print " "

    #return time,npart

def bh_mass(path,snapnum):

    fname=path+"/snapshot_%.3d.hdf5"%snapnum
    f=h5py.File(fname,'r')
    
    hdr=f['Header'].attrs 
    #print hdr.items()
    time=hdr['Time']
    h=hdr['HubbleParam']
    npart=hdr['NumPart_Total']


    print "\nTime: %g [Gyr/h] (= %g Gyr)"%(time,time/h)
    print "nBH:",npart[5]
    
    return np.array(f['PartType5']['BH_Mass'])/h*1.0e10


