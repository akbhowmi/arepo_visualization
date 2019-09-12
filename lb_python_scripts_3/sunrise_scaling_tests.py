import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


def plot(maindir='/oasis/projects/nsf/hvd115/lblecha/q0.5_fg0.3_allrx10_sunruns/',
         subdir='hires_kin/scaling_tests',sfrhist_file='sfrhist_scaling_test_snap170.txt',
         mcrx_file='mcrx_scaling_test_snap170.txt',title=''):

    sfrhist_path="%s/%s/%s"%(maindir,subdir,sfrhist_file)
    with open(sfrhist_path,'r') as sf:
        sn_sf,nc_sf,maxrss_sf,traw_sf,tcpu_sf,tio_sf,telap_sf = np.loadtxt(sf,unpack=True,dtype='float64')
        maxrss_sf = maxrss_sf/1.0e6 # kB to GB
        traw_sf = traw_sf/3600.0 # s to h
        tcpu_sf = tcpu_sf/3600.0 # s to h
        tio_sf = traw_sf/3600.0 # s to h
        telap_sf = telap_sf/3600.0 # s to h

    mcrx_path="%s/%s/%s"%(maindir,subdir,mcrx_file)
    with open(mcrx_path,'r') as mf:
        sn_mc,nc_mc,maxrss_mc,traw_mc,tcpu_mc,tio_mc,telap_mc = np.loadtxt(mf,unpack=True,dtype='float64')
        maxrss_mc = maxrss_mc/1.0e6 # kB to GB
        traw_mc = traw_mc/3600.0 # s to h
        tcpu_mc = tcpu_mc/3600.0 # s to h
        tio_mc = traw_mc/3600.0 # s to h
        telap_mc = telap_mc/3600.0 # s to h

    if 'hires_kin' in subdir:
        if not title: title='hi-res partial SED, w/ kinematics'
        xmax=33
    elif 'late' in subdir:
        if not title: title='low-res full SED, no kinematics'
        xmax=25

    #plt.close()
    #plt.clf()
    #plt.cla()
    fig = plt.figure(figsize=(6,4))

    ax1 = fig.add_subplot(221)
    ax1.set_yscale('log')
    #ax1.set_xlabel('# cores')
    ax1.set_ylabel('Total time [CPUh]')
    ax1.set_title('sfrhist')
    ax1.set_xlim(0,xmax)
    plt.plot(nc_sf,traw_sf,'ko-')
    #plt.plot(nc_sf,tcpu_sf,'bo-')
    #plt.plot(nc_sf,tio_sf,'go-')
    #plt.plot(nc_sf,telap_sf,'ro-')

    #ax2 = fig.add_subplot(222)
    
    #ax2.set_xlabel('# cores')
    #ax2.set_ylabel('Max memory usage [GB]')
    #plt.plot(nc_sf,maxrss_sf)


    ax2 = fig.add_subplot(222)
    ax2.set_yscale('log')
    #ax2.set_xlabel('# cores')
    #ax2.set_ylabel('Total time [CPUh]')
    ax2.set_title('mcrx')
    ax2.set_xlim(7,xmax)
    plt.plot(nc_mc,traw_mc,'ko-')
    #plt.plot(nc_mc,tcpu_mc,'bo-')
    #plt.plot(nc_mc,tio_mc,'go-')
    #plt.plot(nc_mc,telap_mc,'ro-')

    ax3 = fig.add_subplot(223)
    ax3.set_yscale('log')
    ax3.set_xlabel('# cores')
    ax3.set_ylabel('Elapsed time [h]')
    ax3.set_xlim(0,xmax)
    #ax3.set_title('sfrhist')
    #plt.plot(nc_sf,traw_sf,'ko-')
    #plt.plot(nc_sf,tcpu_sf,'bo-')
    #plt.plot(nc_sf,tio_sf,'go-')
    plt.plot(nc_sf,telap_sf,'ro-')

    ax4 = fig.add_subplot(224)
    ax4.set_yscale('log')
    ax4.set_xlabel('# cores')
    ax4.set_xlim(7,xmax)
    #ax4.set_ylabel('Elapsed time [h]')
    #ax4.set_title('mcrx')
    #plt.plot(nc_mc,traw_mc,'ko-')
    #plt.plot(nc_mc,tcpu_mc,'bo-')
    #plt.plot(nc_mc,tio_mc,'go-')
    plt.plot(nc_mc,telap_mc,'ro-')

    #ax4 = fig.add_subplot(224)
    
    #ax4.set_xlabel('# cores')
    #ax4.set_ylabel('Max memory usage [GB]')
    #plt.plot(nc_mc,maxrss_mc)

    fig.subplots_adjust(hspace=0.3,wspace=0.25,bottom=0.15,right=0.95,top=0.87)
    fig.suptitle(title)
    fig.savefig("%s/%s/sfrhist_mcrx_scaling.pdf"%(maindir,subdir))
    plt.close()
    plt.clf()
    plt.cla()
