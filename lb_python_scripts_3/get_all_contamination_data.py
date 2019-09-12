import numpy as np
import column_density as cd
import sys

#print "calculating for z=0"
#sys.stdout.flush()
#cd.merger_phases_agnx0compare(subdir_arr='fid',sfrhist_w_grid=False,plot_wedge_allsim=True,valtype='median',
#                              tmax_postmrg=0.1,alt_wedge=True,errtype='iqr',wise_fluxlim_type='',z=0.0,
#                              contam_file_only=True)


for limtype in ('','SN3','SN10'):

    if limtype=='': 
        print "calculations with no wise lum cut:"
    else: 
        print "calculations for wise lum cut at %s limit"%limtype
    sys.stdout.flush()

    #for z in np.arange(0.1,1.1,0.1):
    for z in np.arange(1.0,1.1,0.1):

        print "calculating for z=%g"%z
        sys.stdout.flush()
        cd.merger_phases_agnx0compare(subdir_arr='fid',sfrhist_w_grid=False,plot_wedge_allsim=True,valtype='median',
                                      tmax_postmrg=0.1,alt_wedge=True,errtype='iqr',wise_fluxlim_type=limtype,
                                      z=z,contam_file_only=True)


    #### all fid sims complete for 1.5, 2.0, & 2.5
    ##for z in np.arange(1.5,3.0,0.5):

    ##    print "calculating for z=%g"%z
    ##    sys.stdout.flush()
    ##    cd.merger_phases_agnx0compare(subdir_arr='fid',sfrhist_w_grid=False,plot_wedge_allsim=True,valtype='median',
    ##                                  tmax_postmrg=0.1,alt_wedge=True,errtype='iqr',wise_fluxlim_type=limtype,
    ##                                  z=z, contam_file_only=True)
