import matplotlib.pyplot as plt
import matplotlib
from numpy import *
import pyfits
import sys

def multicam(snap, ncam):
    """ Plot Hbeta profiles for all cameras,
        with & without NLR particles. """

    #file='../sunruns/mcrx_%d.hdf5.fits'%snap
    #nlrfile='mcrx_%d.hdf5.fits'%snap
    file='../nonlr_data/mcrx_0%d.hdf5.fits'%snap
    nlrfile='mcrx_0%d.hdf5.fits'%snap
    #file='nonlr_data/mcrx_0%d.hdf5.fits'%snap
    #nlrfile='nlr_add_disp_snap%d/mcrx_0%d.hdf5.fits'%{snap,snap}
    
    nof=pyfits.open(file)
    nlf=pyfits.open(nlrfile)    
    lam=nof['lambda'].data.field('lambda')
    nllam=nlf['lambda'].data.field('lambda')

    subp0=331
    fig = plt.figure(figsize=(8,8))
    plt.rc("font",size=10)
    ax1 = fig.add_subplot(subp0)

    #cam=0
    #ins=(nof['CAMERA%d-NONSCATTER'%cam].data.sum(axis=2)).sum(axis=1)
    #ax1.plot(lam, ins)
    
    for cam in range(0,ncam):        
        print "cam=%d"%cam
        i=nof['CAMERA%d'%cam].data.sum(axis=2).sum(axis=1)
        nli=nlf['CAMERA%d'%cam].data.sum(axis=2).sum(axis=1)
        
        ins=nof['CAMERA%d-NONSCATTER'%cam].data.sum(axis=2).sum(axis=1)
        nlins=nlf['CAMERA%d-NONSCATTER'%cam].data.sum(axis=2).sum(axis=1)


        #print "i=",i
        #print "nli=",nli
        
        ax=fig.add_subplot(subp0+cam,sharex=ax1,sharey=ax1) 
        ax.plot(lam, ins,'k')
        ax.plot(nllam, nlins,'b')
        #ax=fig.add_subplot(subp0+1,sharex=ax1,sharey=ax1) 
        ax.plot(lam, i,'r')
        ax.plot(nllam, nli,'m')
            

    ax.set_xlim(4.84e-7,4.88e-7)
    #ax.set_xlim(4.83e-7,4.89e-7)
    #ax.set_ylim(0,8.5e6)
    ax.set_ylim(0.25,3e6)
    #ax.set_ylim(0.25,1e7)
    #ax.set_ylim(0,1e6)
    ax.set_xticks([4.84e-7,4.86e-7,4.88e-7])
    ax.set_yticks([1e6,3e6,5e6,7e6])
    ax.set_xlabel('lambda')
    ax.set_ylabel('flux')
    ax.ticklabel_format(style='sci',scilimits=(-3,3),axis='both')
    plt.subplots_adjust(wspace=0.2,hspace=0.3)
    plt.suptitle(file)
    #plt.show()
    fig.savefig('hbeta_prof_lowlum_%d.pdf'%snap)

        
    nof.close()
    nlf.close()


if __name__ == "__main__":
    multicam(int(sys.argv[1]),int(sys.argv[2]))
