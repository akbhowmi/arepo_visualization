import numpy as np
import wise_agn as wa
import sys
from column_density import define_sunrise_path_arr

##def all_z_fid(maindir='/oasis/projects/nsf/hvd115/lblecha',subdir_arr='fid',
##          min_z=0.1,max_z=0.5,dz=0.1):
#def all_z_fid(maindir='/oasis/projects/nsf/hvd115/lblecha',subdir_arr='fid',
#              min_z=0.6,max_z=1.0,dz=0.1):
def all_z_fid(maindir='/oasis/projects/nsf/hvd115/lblecha',subdir_arr='fid',
              min_z=0.1,max_z=1.0,dz=0.1):

    basepath_arr = define_sunrise_path_arr(basepath=maindir,name=subdir_arr,
                                              sfrhist_w_grid=False,return_full_path=False)
    basepath_arr = [bp.split('/')[0] for bp in basepath_arr]
    print basepath_arr

    for i_path,bp in enumerate(basepath_arr):

        for z in np.arange(min_z,max_z+dz,dz):
        
            print bp+'/z%.1f'%z
            wa.wise_colors(path=maindir,subdir=bp+'/z%.1f'%z,write_llambda=True,input_type='bb')

def all_z_fid_agnx0(maindir='/oasis/projects/nsf/hvd115/lblecha',subdir_arr='fid',
                    min_z=0.1,max_z=1.0,dz=0.1):
    
    basepath_arr = define_sunrise_path_arr(basepath=maindir,name=subdir_arr,
                                              sfrhist_w_grid=False,return_full_path=False)
    basepath_arr = [bp.split('/')[0] for bp in basepath_arr]
    print basepath_arr

    for i_path,bp in enumerate(basepath_arr):

        for z in np.arange(min_z,max_z+dz,dz):
        
            print bp+'/z%.1f_agnx0'%z
            wa.wise_colors(path=maindir,subdir=bp+'/z%.1f_agnx0'%z,write_llambda=True,input_type='bb')


def onesim_z(maindir='/oasis/projects/nsf/hvd115/lblecha',subdir_arr='q1_fg0.3_sunruns',
             min_z=1.5,max_z=4.5,dz=0.5):
    
    basepath = define_sunrise_path_arr(basepath=maindir,name=subdir_arr,
                                           sfrhist_w_grid=False,return_full_path=False)
    print basepath

    for z in np.arange(min_z,max_z+dz,dz):
            wa.wise_colors(path=maindir,subdir=basepath[0]+'/z%.1f'%z,write_llambda=True,input_type='bb')

def onesim_z_agnx0(maindir='/oasis/projects/nsf/hvd115/lblecha',subdir_arr='q1_fg0.3_sunruns',
                    min_z=1.5,max_z=4.5,dz=0.5):
    
    basepath = define_sunrise_path_arr(basepath=maindir,name=subdir_arr,
                                           sfrhist_w_grid=False,return_full_path=False)
    print basepath

    for z in np.arange(min_z,max_z+dz,dz):
            wa.wise_colors(path=maindir,subdir=basepath[0]+'/z%.1f_agnx0'%z,write_llambda=True,input_type='bb')


def all_z0_fid(maindir='/oasis/projects/nsf/hvd115/lblecha',subdir_arr='fid',pubstyle=True):
    
          
    path_arr = define_sunrise_path_arr(basepath=maindir,name=subdir_arr,
                                       sfrhist_w_grid=False,return_full_path=False)

    for i_path,p in enumerate(path_arr):
        
        print p
        wa.wise_colors(path=maindir,subdir=p,write_llambda=True,input_type='bb')

def all_z0_fid_agnx0(maindir='/oasis/projects/nsf/hvd115/lblecha',subdir_arr='fid_agnx0',pubstyle=True):
    
          
    path_arr = define_sunrise_path_arr(basepath=maindir,name=subdir_arr,
                                       sfrhist_w_grid=False,return_full_path=False)

    for i_path,p in enumerate(path_arr):
        
        print p
        wa.wise_colors(path=maindir,subdir=p,write_llambda=True,input_type='bb')


def all_z0_hires(maindir='/oasis/projects/nsf/hvd115/lblecha',subdir_arr='hires',pubstyle=True):
    
          
    path_arr = define_sunrise_path_arr(basepath=maindir,name=subdir_arr,
                                       sfrhist_w_grid=False,return_full_path=False)

    for i_path,p in enumerate(path_arr):
        
        print p
        #wa.wise_colors(path=maindir,subdir=p,write_lums=True,input_type='mcrx')
        wa.wise_colors(path=maindir,subdir=p,write_llambda=True,input_type='bb')


if __name__=='__main__':

    func=sys.argv[1]
    if len(sys.argv)==3:
        sim=sys.argv[2]
        if func=='all_z_fid': all_z_fid(subdir_arr=[sim],min_z=0.1,max_z=1.0,dz=0.1)
        if func=='all_z_fid_agnx0': all_z_fid_agnx0(subdir_arr=[sim],min_z=0.1,max_z=1.0,dz=0.1)
    else:
        if func=='all_z_fid': all_z_fid(min_z=1.5,max_z=1.5,dz=0.1)
    #if func=='all_z_fid': all_z_fid(min_z=1.5,max_z=1.5,dz=0.5)
    #if func=='all_z_fid' and sim=='q1_fg0.3_sunruns': all_z_fid(subdir_arr=[sim],min_z=0.7,max_z=1.0,dz=0.1)
    #if func=='all_z_fid' and sim!='q1_fg0.3_sunruns': all_z_fid(subdir_arr=[sim],min_z=0.6,max_z=1.0,dz=0.1)
    #if func=='all_z_fid': all_z_fid()
    if func=='all_z0_fid': all_z0_fid()
    if func=='all_z0_fid_agnx0': all_z0_fid_agnx0()
    if func=='all_z0_hires': all_z0_hires()
    if func=='onesim_z': onesim_z()
    if func=='onesim_z_agnx0': onesim_z_agnx0()
