import numpy as np
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import PyCos as pc
import sys
import merging_bh_info as mbi
from scipy import sparse
from copy import copy
import subprocess
import astro_constants as ac
import bol_corrections as bc



def plot(dir='/n/hernquistfs2/lblecha/illustris_data/',
         mrgdir='recoils_v0515_gas',obsdir='hst',
         run='Illustris-1', zmin=0, zmax=6, 
         min_proghalo_npart=300,min_progstar_npart=20,
         min_bh_mass=1.0e6,min_tothalo_mass=0.0,
         min_totstar_hfmass=0.0,):

    #recdir_arr = ['random_spins','dryrandom_spins']
    #fname='kick_distributions_random.eps'
    #recdir_arr = ['fgassf0.1_spins','fgassf0.1_dryrandom_spins']
    #fname='kick_distributions_hybrid.eps'
    recdir_arr = ['hot_spins','cold_spins','5deg_spins']
    fname='kick_distributions_aligned.eps'

    fig = plt.figure(figsize=(5,8))

    ax1=fig.add_subplot(2,1,1)
    plt.xlim(1.2,3.8)
    plt.ylim(0.7,5e3)
    plt.yscale('log')

    ax2=fig.add_subplot(2,1,2)
    plt.xlim(1.2,3.8)
    #plt.ylim(0.7,2e4)
    plt.ylim(0,1.05)
    #plt.yscale('log')
    
    #kwargs = dict(drawstyle='steps-mid')
    kwargs = dict(drawstyle='default')

    carr = ['k','m','g']
    icarr = ['b','orange','r']

    for i,recdir in enumerate(recdir_arr):

        mrgpath = dir+'/'+run+'/'+mrgdir+'/'+obsdir
        recpath = mrgpath+'/'+recdir

        mrgdata,progdata=mbi.load_mrg_and_prog_data(mrgpath=mrgpath,run=run,sort_z=False,min_proghalo_npart=min_proghalo_npart,min_progstar_npart=min_progstar_npart,min_bh_mass=min_bh_mass,min_tothalo_mass=min_tothalo_mass,min_totstar_hfmass=min_totstar_hfmass,make_cuts=True)

        imrg,task,time,redshift,id1,bhmass1,id2,bhmass2,bhmdot1,bhlrho1,bhlc1 = mrgdata

        k_imrg,k_z,k_mbh,k_q,k_mhalo,k_vkick = mbi.load_kickdata(recpath,imrg=imrg,sort_z=False,make_cuts=True)

        i_imrg,i_z,i_jmax,i_status,i_tstop,i_rstop,i_phistop,i_zstop,i_vrstop,i_vphistop,i_vzstop,i_lbol_final,i_lbol_const,i_lbol_lim,i_dR_lim,i_dv_lim,i_thetaLOS,i_toffAGN_dr,i_toffAGN_dv,i_toffAGN_const_dr,i_toffAGN_const_dv,i_toffAGN_drLOS,i_toffAGN_dvLOS,i_toffAGN_const_drLOS,i_toffAGN_const_dvLOS = mbi.load_intdata(recpath,imrg=imrg,k_imrg=k_imrg,sort_z=False,make_cuts=True)
        
        ix_i = [n for n,elem in enumerate(k_imrg) if elem in i_imrg]
        i_vkick = k_vkick[ix_i]

        vkhist,binedg = np.histogram(np.log10(k_vkick),range=(-4,3.8),bins=80)
        #vkhist[vkhist==0] = 1.0e-10
        ivkhist,tmp = np.histogram(np.log10(i_vkick),bins=binedg)
        #ivkhist[ivkhist==0] = 1.0e-10
        ivkhist_offdrLOS,tmp = np.histogram(np.log10(i_vkick[i_toffAGN_drLOS>0]),
                                            bins=binedg,
                                            weights=i_toffAGN_drLOS[i_toffAGN_drLOS>0])
        #ivkhist_offdrLOS[ivkhist_offdrLOS==0] = 1.0e-10
        ivkhist_offdvLOS,tmp = np.histogram(np.log10(i_vkick[i_toffAGN_dvLOS>0]),
                                            bins=binedg,
                                            weights=i_toffAGN_dvLOS[i_toffAGN_dvLOS>0])
        #ivkhist_offdvLOS[ivkhist_offdvLOS==0] = 1.0e-10
        binsize=binedg[1]-binedg[0]
        bins = binedg[:-1]+0.5*binsize


        if i_toffAGN_dr.max()>0:
            print i_vkick[i_toffAGN_drLOS>0].min()
            diff=np.abs(1.0*np.cumsum(ivkhist_offdrLOS)/np.sum(ivkhist_offdrLOS)-0.5)
            ix_median = np.where(diff==diff.min())[0]
            print 10**bins[ix_median]
        if i_toffAGN_dv.max()>0:
            print i_vkick[i_toffAGN_dvLOS>0].min()
            diff=np.abs(1.0*np.cumsum(ivkhist_offdvLOS)/np.sum(ivkhist_offdvLOS)-0.5)
            ix_median = np.where(diff==diff.min())[0]
            print 10**bins[ix_median]
 
        ax1.plot(bins,vkhist,color=carr[i],linewidth=1.5,**kwargs)
        ax1.plot(bins,ivkhist,
                 color=icarr[i],linestyle='dashed',linewidth=1.5,**kwargs)
        ax1.plot(bins,ivkhist_offdrLOS/i_toffAGN_drLOS.max(),
                 color=icarr[i],linestyle='dotted',linewidth=1.5,**kwargs)
        ax1.plot(bins,ivkhist_offdvLOS/i_toffAGN_dvLOS.max(),
                 color=icarr[i],linestyle='dashdot',linewidth=1.5,**kwargs)

        #print 1.0*np.cumsum(ivkhist)/np.sum(vkhist)
        ax2.plot(bins,1.0*np.cumsum(vkhist)/np.sum(vkhist),
                 color=carr[i],linewidth=1.5,**kwargs)
        ax2.plot(bins,1.0*np.cumsum(ivkhist)/np.sum(ivkhist),
                 color=icarr[i],linestyle='dashed',linewidth=1.5,**kwargs)
        ax2.plot(bins,1.0*np.cumsum(ivkhist_offdrLOS)/np.sum(i_toffAGN_drLOS),
                 color=icarr[i],linestyle='dotted',linewidth=1.5,**kwargs)
        ax2.plot(bins,1.0*np.cumsum(ivkhist_offdvLOS)/np.sum(i_toffAGN_dvLOS),
                 color=icarr[i],linestyle='dashdot',linewidth=1.5,**kwargs)

        print "\ntotal mergers in vkhist for %s: %d"%(recdir,np.sum(vkhist))
        print "total mergers in ivkhist for %s: %d\n"%(recdir,np.sum(ivkhist))
        

    fig.savefig(mrgpath+'/'+fname)
    plt.clf()
    plt.cla()
    plt.close()
