import numpy as np
import column_density as cd

def all_z(maindir='/oasis/projects/nsf/hvd115/lblecha',subdir_arr='fid',
          min_z=0.1,max_z=0.5,dz=0.1,agnscale=False):

    basepath_arr = cd.define_sunrise_path_arr(basepath=maindir,name=subdir_arr,
                                              sfrhist_w_grid=False,return_full_path=False)
    basepath_arr = [bp.split('/')[0] for bp in basepath_arr]
    print basepath_arr

    for i_path,bp in enumerate(basepath_arr):

        tmax=4.2 if 'q0.2_fg0.3_sunruns' in bp else 2.95
        for z in np.arange(min_z,max_z+dz,dz):
        
            print bp+'/z%.1f'%z
            print bp+'/z%.1f_agnx0'%z
            cd.color_color_compare(subdirA=bp+'/z%.1f'%z,
                                   subdirB=bp+'/z%.1f_agnx0'%z,
                                   tmax=tmax,pubstyle=False,agnscale=agnscale)


def all_z0_tscale(maindir='/oasis/projects/nsf/hvd115/lblecha',subdir_arr='fid',pubstyle=True):
    
          
    path_arr = cd.define_sunrise_path_arr(basepath=maindir,name=subdir_arr,
                                          sfrhist_w_grid=False,return_full_path=False)

    for i_path,p in enumerate(path_arr):
        
        tmax=4.2 if 'q0.2_fg0.3_sunruns' in p else 2.95
        cd.color_color_compare(subdirA=p, subdirB=p+'_agnx0',
                               tmax=tmax,pubstyle=pubstyle,agnscale=False)

def all_z0_agnscale(maindir='/oasis/projects/nsf/hvd115/lblecha',subdir_arr='fid',pubstyle=True):
    
          
    path_arr = cd.define_sunrise_path_arr(basepath=maindir,name=subdir_arr,
                                          sfrhist_w_grid=False,return_full_path=False)

    for i_path,p in enumerate(path_arr):
        
        cd.color_color_compare(subdirA=p, subdirB=p+'_agnx0',pubstyle=pubstyle)


def min_w2w3(maindir='/oasis/projects/nsf/hvd115/lblecha',subdir_arr='fid',
             z0=True,min_z=0.1,max_z=0.9,dz=0.1,total_only=True):

    if z0:
        path_arr = cd.define_sunrise_path_arr(basepath=maindir,name=subdir_arr,
                                              sfrhist_w_grid=False,return_full_path=False)
    else:
        basepath_arr = cd.define_sunrise_path_arr(basepath=maindir,name=subdir_arr,
                                                  sfrhist_w_grid=False,return_full_path=False)
        basepath_arr = [bp.split('/')[0] for bp in basepath_arr]
        path_arr=[]
        for i_path,bp in enumerate(basepath_arr):
            path_arr = path_arr + [bp+'/z%.1f'%z for z in np.arange(min_z,max_z+dz,dz)]

    min_w2w3=1e6
    min_w2w3_above=1e6
    for i_path,p in enumerate(path_arr):
        time,snap,w1w2,w2w3=cd.load_wise(maindir=maindir,subdir=p)
        if not total_only: 
            print "\n"+p
            print "min w2w3: %g"%w2w3.min()
        if w1w2.max()>0.5:
            if not total_only:print "min w2w3[w1w2>0.5]: %g"%w2w3[w1w2>0.5].min()
            if w2w3[w1w2>0.5].min()<min_w2w3_above: min_w2w3_above=w2w3[w1w2>0.5].min()
        else: 
            if not total_only:print "no w1w2>0.5."
        if w2w3.min()<min_w2w3: min_w2w3=w2w3.min()

    print "\nmin w2w3 for all sims: %g"%min_w2w3
    print "min w2w3[w1w2>0.5] for all sims: %g"%min_w2w3_above
