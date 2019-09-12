#! /n/sw/python-2.7.1/bin/python

import nlr_debug_add_dispersion
import numpy as np

### WARNING! ###
### The program called here hemmorhages memory in a yet-unidentified manner.
### Use of this wrapper is strongly discouraged until this is fixed.
### Instead, use the add_nlr.bsub bash script to submit snaps individually.
################

def add_nlr():

    #nlrdir_arr=['test_Urng_nmax6_mc4_ffa1']
    #for i in range(len(nlrdir_arr)):
    #    nlrdir_arr[i] = '/'.join(['nlr_model_data',nlrdir_arr[i]])        

    #maindir='/n/scratch2/hernquist_lab/lblecha/nlr_sims_v112110/q1_fg0.1_allrx10_nomrg/hisnapres/'
    #rundir_arr  = ['sunsnaps/add_nlr_nmax6_mc4_ffa1']
    #snap_arr = np.array([413,414,420,421,422,424,425,426,431,433,435,436,439,441,443,446,447,448,452])

    #maindir='/n/scratch2/hernquist_lab/lblecha/nlr_sims_v112110/q1_fg0.1_allrx10_nomrg/hisnapres/'
    #nlrdir_arr=['test_Urng_nmax6_mc200_ffa1']
    #rundir_arr  = ['sunsnaps/add_nlr_nmax6_mc200_ffa1']
    ##snap_arr = np.array([413,414,420,421,422,424,425,426,431,433,435,436,439,441,443,446,447,448,452])



    #maindir='/n/scratch2/hernquist_lab/lblecha/nlr_sims_v112110/q0.5_fg0.04_allrx10_nomrg/hisnapres/'
    ##nlrdir_arr=['test_Urng_nmax6_mc200_ffa1']
    #nlrdir_arr=['test_Urng_nmax6_mc4_ffa1']
    ##snap_arr = np.array([272, 274, 275, 276, 277, 278, 279, 281, 282, 283, 
    ##                     284, 285, 286, 288, 291, 292, 293, 294, 303, 342,	
    ##                     346, 351, 355, 360, 362, 363, 369, 371, 376, 377,
    ##                     382, 384, 385, 386, 389])
    #snap_arr = np.array([275, 276, 277, 278, 279, 281, 282, 283, 
    #                     284, 285, 286, 288, 291, 292, 293, 294, 303, 342,	
    #                     346, 351, 355, 360, 362, 363, 369, 371, 376, 377,
    #                     382, 384, 385, 386, 389])
    ##rundir_arr  = ['sunruns_add_nlr_nmax6_mc200_ffa1']
    #rundir_arr  = ['sunruns_add_nlr_nmax6_mc4_ffa1']

    #maindir='/n/scratch2/hernquist_lab/lblecha/nlr_sims_v112110/q0.5_fg0.1_b0.2add_HRbh_mdyn5e-4_allrx10_nomrg/hisnapres/'
    ##nlrdir_arr=['test_Urng_nmax6_mc200_ffa1']
    #nlrdir_arr=['test_Urng_nmax6_mc4_ffa1']
    ##snap_arr = np.array([273,274,514,580,599,602,604,605,606,611,614,615,616,617,618,623,624])
    #snap_arr = np.array([580,599,602,604,605,606,611,614,615,616,617,618,623,624])
    ##rundir_arr  = ['sunruns_add_nlr_nmax6_mc200_ffa1']
    #rundir_arr  = ['sunruns_add_nlr_nmax6_mc4_ffa1']

    #maindir='/n/scratch2/hernquist_lab/lblecha/nlr_sims_v112110/q0.5_fg0.3_b0.3add_HRbh_mdyn5e-4_allrx10_nomrg/hisnapres/'
    #nlrdir_arr=['test_Urng_nmax6_mc200_ffa1']
    #snap_arr = np.array([394,395,396,397,399,400,402,404,407,410,411,419,420,421])
    #rundir_arr  = ['sunruns_add_nlr_nmax6_mc200_ffa1']

    maindir='/n/nss2a/hernquist_scratch1/lblecha/nlr_sims_v112110/q0.5_fg0.1_allrx10_nomrg/hisnapres'
    #nlrdir_arr=['test_Urng_nmax6_mc200_ffa1']
    nlrdir_arr=['test_Urng_nmax6_mc4_ffa1']
    snap_arr = np.array([299, 300, 301, 303, 304, 305, 306, 307, 308, 309,
                         310, 311, 312, 314, 317, 318, 319, 320, 321, 322,
                         323, 324, 325, 326, 336, 390])
    #rundir_arr = ['sunruns_add_nlr_nmax6_mc200_ffa1']
    rundir_arr = ['sunruns_add_nlr_nmax6_mc4_ffa1']

    sfrhist_files = np.array(['sfrhist_'+str(snap_arr[j])+'.hdf5.fits' for j in range(snap_arr.size)])
    #sfrhist_files = np.array(['sfrhist_0'+str(snap_arr[j])+'.hdf5.fits' for j in range(snap_arr.size)])
    print nlrdir_arr
    print rundir_arr
    print sfrhist_files
    
    #subdir_arr = np.array(['test_Urng_nmax5_ffa0.1', 'test_Urng_nmax6_mc200_ffa1', 'test_Urng_nmax7_ffa0.2', 'test_Urng_nmax5_ffa0.2', 'test_Urng_nmax6_mc2.5_ffa1', 'test_Urng_nmax8_ffa0.2', 'test_Urng_nmax6_ffa0.1', 'test_Urng_nmax6_mc3_ffa1', 'test_Urng_nmax6_mc4_ffa1'])
    #subdir_arr = np.array(['test_Urng_nmax5_ffa0.1', 'test_Urng_nmax6_mc200_ffa1', 'test_Urng_nmax7_ffa0.2', 'test_Urng_nmax5_ffa0.2', 'test_Urng_nmax6_mc2.5_ffa1', 'test_Urng_nmax8_ffa0.2', 'test_Urng_nmax6_ffa0.1', 'test_Urng_nmax6_mc3_ffa1', 'test_Urng_nmax6_ffa0.2', 'test_Urng_nmax6_mc4_ffa1'])


    for i in range(len(rundir_arr)):
        print rundir_arr[i]
        for j in range(snap_arr.size):
            print "snap",snap_arr[j]
            nlr_debug_add_dispersion.add_hbeta_sources_DEBUG_ADD_DISP(maindir, nlrdir_arr[i],
                                                                      rundir_arr[i], sfrhist_files[j],
                                                                      snap_arr[j])
            


if __name__ == "__main__":
    add_nlr()
