import numpy as np
import h5py


def read(path,snapnum):

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
    for grp in f.keys():        
    #for grp in f.keys()[1:]:        
        print "grp: ",grp
        print f[grp].keys()
        print " "
        if grp=='PartType5':
            print np.array(f[grp]['BH_Mass'])/.7*1.0e10

    #return time,npart
