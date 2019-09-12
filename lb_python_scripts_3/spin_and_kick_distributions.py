import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from copy import copy
import sys
import PyCos as pc
import scipy.optimize as opt
import ctypes

def readfile(fname):

    try:
        fp=open(fname,"r")
    except IOError:
        print "Error opening file %s."%fname
        sys.exit()

    #dtypes = np.dtype({'names':(['f{}'.format(i) for i in range(8)]),
    #                   'formats':[np.float,np.float,np.float,np.float,
    #                              np.float,np.float,np.float,np.float]})
    tmpdata = np.loadtxt(fp, unpack=True, dtype=np.float)

    return tmpdata

def plot(dir='/n/home00/lblecha/recoil_distributions/n1e3/',pubversion=True):

    title=''
    #file_arr = ["amag9_arand_lgqrand_params.dat",
    #            "aallrand_lgqrand_params.dat",
    #            "amag9_arand5_lgqrand_params.dat",
    #            "hot_lgqrand_params.dat",
    #            "cold_lgqrand_params.dat",
    #            "hot_qcosmo_params.dat",
    #            "cold_qcosmo_params.dat",
    #            "hot_qrand_params.dat",
    #            "cold_qrand_params.dat"]

    #plotname = 'spin_kick_distributions_hot_cold.eps'
    #file_arr = ["hot_lgqrand_params.dat",
    #            "cold_lgqrand_params.dat",
    #            "hot_qcosmo_params.dat",
    #            "cold_qcosmo_params.dat",
    #            "hot_qrand_params.dat",
    #            "cold_qrand_params.dat"]
    #c_arr = ['m','c','r','b','orange','g']

    #file_arr = ["hotpri_coldsec_qcosmo_params.dat",
    #            "coldpri_hotsec_qcosmo_params.dat",
    #            "hotpri_coldsec_qrand_params.dat",
    #            "coldpri_hotsec_qrand_params.dat"]    
    #c_arr = ['r','b','orange','cyan']
    #title = 'red/orange = hot primary, cold secondary\nblue/cyan=cold primary, hot secondary'

    plotname = 'spin_kick_distributions.eps'
    if not pubversion: plotname = 'spin_kick_distributions_with_los.eps'
    #file_arr = ["amag9_arand_lgqrand_params.dat",
    #            "aallrand_lgqrand_params.dat",
    #            "hot_lgqrand_params.dat",
    #            "cold_lgqrand_params.dat",               
    #            "amag9_arand5_lgqrand_params.dat"]
    file_arr = ["amag9_arand_lgqrand_params.dat",
                "amrg_thetarand_lgqrand_params.dat",
                "hot_lgqrand_params.dat",
                "cold_lgqrand_params.dat",               
                "amag9_arand5_lgqrand_params.dat"]
    c_arr = ['k','g','r','b','orange','c']

    #plotname = 'spin_kick_distributions_with_b08.eps'
    #file_arr = ["amag9_arand_lgqrand_params.dat",
    #            "aallrand_lgqrand_params.dat",
    #            "hot_lgqrand_params.dat",
    #            "cold_lgqrand_params.dat",               
    #            "amag9_arand5_lgqrand_params.dat",
    #            "b08_amag9_arand_lgqrand_params.dat"]

    #c_arr = ['k','g','r','b','orange','c']
    #ls_arr = ['solid','dashed','dashed','solid','solid']

    plt.clf()
    plt.cla()
    plt.close('all')
    
    #fig = plt.figure(figsize=(8,8))
    fig = plt.figure(figsize=(4.5,6.5))
    #kwargs = dict(histtype='step',normed=True)
    kwargs = dict(histtype='step',bins=40,normed=True)
    mpl.rcParams['font.size'] = 10
    mpl.rcParams['lines.linewidth'] = 1.5
    mpl.rcParams['axes.linewidth'] = 1.5

    #ax1 = fig.add_subplot(2,2,1)
    ax1 = fig.add_subplot(3,1,1)
    plt.xlim(0,1)
    plt.ylim(0,5)
    plt.yticks([1,2,3,4])
    plt.xlabel(r'$a$',fontsize=12)
    plt.ylabel(r'$P (a)$',fontsize=12)
    #plt.text(0.05,4.4,'random',color='g',fontsize=10)
    plt.text(0.05,4.4,'dry',color='g',fontsize=10)
    plt.text(0.05,4.0,r'$a=0.9$',color='k',fontsize=10)
    #plt.text(0.05,4.0,'high',color='k',fontsize=10)
    plt.text(0.05,3.6,'hot',color='r',fontsize=10)
    plt.text(0.05,3.2,'cold',color='b',fontsize=10)

    #ax2 = fig.add_subplot(2,2,2)
    ax2 = fig.add_subplot(3,1,2)
    plt.xlim(0,90)
    #plt.xlim(0,180)
    #plt.ylim(1.0e-5,1)
    plt.yticks([0.05,0.1,0.15,0.2,0.25])
    plt.ylim(0,0.17)
    plt.xlabel(r'$\theta_{1,2}$ [deg]',fontsize=12)
    plt.ylabel(r'$P (\theta_{1,2}$)',fontsize=12)
    plt.text(68,0.15,'random',color='k',fontsize=10)
    plt.text(68,0.135,r'$\Delta\theta<5^{\circ}$',color='orange',fontsize=10)
    plt.text(68,0.12,'hot',color='r',fontsize=10)
    plt.text(68,0.105,'cold',color='b',fontsize=10)

    #ax3 = fig.add_subplot(2,2,3)
    #plt.xlim(-2,0)
    #plt.ylim(0,0.6)
    #plt.xlabel('log(q)')

    #ax4 = fig.add_subplot(2,2,4)
    ax4 = fig.add_subplot(3,1,3)
    plt.xlim(-0.1,3.7)
    #plt.ylim(5.0e-3,1)
    #plt.xlim(0,4500)
    #plt.ylim(1.0e-6,0.1)
    plt.yticks([0.2,0.4,0.6,0.8,1.0,1.2])
    plt.ylim(0,1.1)
    plt.xlabel(r'$v_k$ [log km s$^{-1}$]',fontsize=12)
    plt.ylabel(r'$P(v_k$)',fontsize=12)
    if pubversion:
        plt.text(0.1,0.95,'random-high',color='k',fontsize=10)
        #plt.text(0.1,0.86,'random-random',color='g',fontsize=10)
        plt.text(0.1,0.86,'random-dry',color='g',fontsize=10)
        plt.text(0.1,0.77,'hot',color='r',fontsize=10)
        plt.text(0.1,0.68,'cold',color='b',fontsize=10)
        plt.text(0.1,0.59,'5deg-high',color='orange',fontsize=10)
        if len(file_arr)==6:
            plt.text(0.1,0.5,'B08 random-high',color='c',fontsize=10)

    for i,file in enumerate(file_arr):

        kwargs.update(color=c_arr[i])

        tmpdata = readfile(dir+"/"+file)
        q, a2, a1, costh2, costh1, cosTheta12, phidiff, vk, vklos = tmpdata
        nmrg = len(q)
        print "Read %d lines from file %s."%(nmrg,file)
        print "vkmin=%g, vkmax=%g."%(vk.min(),vk.max())
        if q.min()<=0.0 or vk.min()<=0.0:
            print "%d mergers had q=0."%(q[q==0].size)
            print "%d mergers had vk=0."%(vk[vk==0].size)
            q[q==0] = 1.0e-10
            vk[vk==0] = 1.0e-10

        print "model %s:"%file
        print "Fraction of kicks > 100 km/s: %g"%(vk[vk>100.0].size/(1.0*nmrg))
        print "Fraction of kicks 100-200 km/s: %g"%(vk[(vk>100.0)&(vk<=200.0)].size/(1.0*nmrg))
        print "Fraction of kicks > 500 km/s: %g"%(vk[vk>500.0].size/(1.0*nmrg))
        print "For q<0.3, fraction of kicks > 100 km/s: %g"%(vk[(vk>100.0)&(q<0.3)].size/(1.0*nmrg))
        print "For q<0.3, raction of kicks 100-200 km/s: %g"%(vk[(vk>100.0)&(vk<=200.0)&
                                                                 (q<0.3)].size/(1.0*nmrg))
        print "For q<0.3, fraction of kicks > 500 km/s: %g\n"%(vk[(vk>500.0)&
                                                                  (q<0.3)].size/(1.0*nmrg))

        print "%s %g %g"%("f(vk),f(vkLOS) < 100 km/s:",
                           vk[vk<100.0].size/(1.0*nmrg),
                           vklos[vklos<100.0].size/(1.0*nmrg))
        print "%s %g %g"%("f(vk),f(vkLOS) 100-200 km/s:",
                          vk[(vk>100.0)&(vk<=200.0)].size/(1.0*nmrg),
                          vklos[(vklos>100.0)&(vklos<=200.0)].size/(1.0*nmrg))
        print "%s %g %g"%("f(vk),f(vkLOS) 200-300 km/s:",
                          vk[(vk>200.0)&(vk<=300.0)].size/(1.0*nmrg),
                          vklos[(vklos>200.0)&(vklos<=300.0)].size/(1.0*nmrg))
        print "%s %g %g"%("f(vk),f(vkLOS) 300-400 km/s:",
                          vk[(vk>300.0)&(vk<=400.0)].size/(1.0*nmrg),
                          vklos[(vklos>300.0)&(vklos<=400.0)].size/(1.0*nmrg))
        print "%s %g %g"%("f(vk),f(vkLOS) 400-500 km/s:",
                          vk[(vk>400.0)&(vk<=500.0)].size/(1.0*nmrg),
                          vklos[(vklos>400.0)&(vklos<=500.0)].size/(1.0*nmrg))
        print "%s %g %g"%("f(vk),f(vkLOS) 500-1000 km/s:",
                          vk[(vk>500.0)&(vk<=1000.0)].size/(1.0*nmrg),
                          vklos[(vklos>500.0)&(vklos<=1000.0)].size/(1.0*nmrg))
        print "%s %g %g"%("f(vk),f(vkLOS) 1000-1500 km/s:",
                          vk[(vk>1000.0)&(vk<=1500.0)].size/(1.0*nmrg),
                          vklos[(vklos>1000.0)&(vklos<=1500.0)].size/(1.0*nmrg))
        print "%s %g %g"%("f(vk),f(vkLOS) 1500-2000 km/s:",
                          vk[(vk>1500.0)&(vk<=2000.0)].size/(1.0*nmrg),
                          vklos[(vklos>1500.0)&(vklos<=2000.0)].size/(1.0*nmrg))
        print "%s %g %g"%("f(vk),f(vkLOS) 2000-2500 km/s:",
                          vk[(vk>2000.0)&(vk<=2500.0)].size/(1.0*nmrg),
                          vklos[(vklos>2000.0)&(vklos<=2500.0)].size/(1.0*nmrg))
        print "%s %g %g"%("f(vk),f(vkLOS) 2500-3000 km/s:",
                          vk[(vk>2500.0)&(vk<=3000.0)].size/(1.0*nmrg),
                          vklos[(vklos>2500.0)&(vklos<=3000.0)].size/(1.0*nmrg))
        print "%s %g %g"%("f(vk),f(vkLOS) 3000-3500 km/s:",
                          vk[(vk>3000.0)&(vk<=3500.0)].size/(1.0*nmrg),
                          vklos[(vklos>3000.0)&(vklos<=3500.0)].size/(1.0*nmrg))
        print "%s %g %g"%("f(vk),f(vkLOS) 3500-4000 km/s:\n",
                          vk[(vk>3500.0)&(vk<=4000.0)].size/(1.0*nmrg),
                          vklos[(vklos>3500.0)&(vklos<=400.0)].size/(1.0*nmrg))

        #nbins=100 if i==0 else 30
        nbins=50
        a2hist,a2_binedg = np.histogram(a2,bins=nbins,normed=True)
        hfbinw = 0.5*(a2_binedg[1]-a2_binedg[0])
        #ls='dashed' if i==4 else 'solid'
        if (not pubversion) or (i > 0 and i < 4):
            #ax1.plot(a2_binedg[:-1]+hfbinw,a2hist,color=c_arr[i])
            ax1.plot(a2_binedg+hfbinw,np.append(a2hist,0),color=c_arr[i])
        ax1.plot([0.9,0.9],[0,5],color=c_arr[0])
        #ax1.hist(a2,**kwargs)
        #ax1.hist(a1,color=c_arr[i],linestyle='dashed',**kwargs)
        
        anglehist,angle_binedg = np.histogram((180/np.pi)*np.arccos(np.abs(cosTheta12)),bins=50,normed=True)
        hfbinw = 0.5*(angle_binedg[1]-angle_binedg[0])
        #ls='dashed' if i==1 else 'solid'
        if (not pubversion) or (i != 1 and i != 5):
            ax2.plot(angle_binedg+hfbinw,np.append(anglehist,anglehist[-1]),
                     color=c_arr[i])
        #ax2.hist((180/np.pi)*np.arccos(np.abs(cosTheta12)),**kwargs)
        #ax2.hist((180/np.pi)*np.arccos(cosTheta12),**kwargs)
        #ax2.hist((180/np.pi)*np.arccos(cosTheta12),log=True,**kwargs)

        #lgqhist,lgq_binedg = np.histogram(np.log10(q),bins=10,normed=True)
        #hfbinw = 0.5*(lgq_binedg[1]-lgq_binedg[0])
        #ax3.plot(lgq_binedg[:-1]+hfbinw,lgqhist,color=c_arr[i])
        ##ax3.hist(np.log10(q),**kwargs)

        lgvkhist,lgvk_binedg = np.histogram(np.log10(vk),bins=50,normed=True)
        hfbinw = 0.5*(lgvk_binedg[1]-lgvk_binedg[0])
        ax4.plot(lgvk_binedg[:-1]+hfbinw,lgvkhist,color=c_arr[i])
        if not pubversion:
            lgvkloshist, lgvklos_binedg = np.histogram(np.log10(vklos),bins=70,normed=True)
            hfbinw = 0.5*(lgvklos_binedg[1]-lgvklos_binedg[0])
            ax4.plot(lgvklos_binedg[:-1]+hfbinw,lgvkloshist,color=c_arr[i],linestyle='dotted')
        #ax4.hist(np.log10(vk),**kwargs)
        #ax4.hist(vk,log=True,**kwargs)
        #ax4.hist(np.log10(vk),log=True,**kwargs)

    #fig.suptitle(title)
    fig.subplots_adjust(left=0.16,bottom=0.08,top=0.98,right=0.94,hspace=0.33)
    #fig.tight_layout()
    fig.savefig(dir+'/'+plotname)
    plt.clf()
    plt.cla()
    plt.close('all')

    print "Finished making plot."
