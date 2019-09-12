import column_density as cd


fidres_grid_res=0.136 ## 0.136 kpc = 0.034 h^-1 kpc / 0.7 * 2.8 (when grav smoothing kernal becomes newtonian)
hires_grid_res=0.064 ## 0.064 kpc = 0.016 h^-1 kpc / 0.7 * 2.8 (when grav smoothing kernal becomes newtonian)

fid_arr = cd.define_sunrise_path_arr(name='fid',return_full_path=False)
hires_arr = cd.define_sunrise_path_arr(name='hires',return_full_path=False)
arr = fid_arr+hires_arr
##arr=hires_arr
#print "\nmaking single-sim 'true Ngas' plots for these runs:"
#print arr

##for subdir in arr:
##    print "\n",subdir
##    cd.merger_phases(subdir_arr=subdir,sfrhist_w_grid=True,single_sim_plots=False,
##                     equal_bins=False,ngas_only=True,tmax_postmrg=0.25,grid_nh_res=0)

#for subdir in fid_arr:
#    print "\n",subdir
#    cd.merger_phases(subdir_arr=subdir,sfrhist_w_grid=True,single_sim_plots=False,
#                     equal_bins=False,ngas_only=True,tmax_postmrg=0.25,grid_nh_res=fidres_grid_res)

#for subdir in hires_arr:
#    print "\n",subdir
#    cd.merger_phases(subdir_arr=subdir,sfrhist_w_grid=True,single_sim_plots=False,
#                     equal_bins=False,ngas_only=True,tmax_postmrg=0.25,grid_nh_res=fidres_grid_res)
#    print "\n",subdir
#    cd.merger_phases(subdir_arr=subdir,sfrhist_w_grid=True,single_sim_plots=False,
#                     equal_bins=False,ngas_only=True,tmax_postmrg=0.25,grid_nh_res=hires_grid_res)


name_arr = ['fid','gasrich','gaspoor','fid_major','gasrich_major','fid_compare_to_hires','hires']
fidres_name_arr = ['fid','gasrich','gaspoor','fid_major','gasrich_major','fid_compare_to_hires']
hires_name_arr = ['hires']
print "\nmaking multi-sim 'true Ngas' plots for these run types:"
print name_arr

##for name in name_arr:
##    print "\n%s"%name
#    #cd.merger_phases(subdir_arr=name,sfrhist_w_grid=True,single_sim_plots=False,
#    #                 equal_bins=False,ngas_only=True,tmax_postmrg=0.25,grid_nh_res=0)

#for name in fidres_name_arr:
#    cd.merger_phases(subdir_arr=name,sfrhist_w_grid=True,single_sim_plots=False,
#                     equal_bins=False,ngas_only=True,tmax_postmrg=0.25,grid_nh_res=fidres_grid_res)

name='fid'
cd.merger_phases(subdir_arr=name,sfrhist_w_grid=True,single_sim_plots=False,
                 equal_bins=False,ngas_only=True,tmax_postmrg=0.25,grid_nh_res=fidres_grid_res)

#for name in hires_name_arr:
#    cd.merger_phases(subdir_arr=name,sfrhist_w_grid=True,single_sim_plots=False,
#                     equal_bins=False,ngas_only=True,tmax_postmrg=0.25,grid_nh_res=fidres_grid_res)
#    cd.merger_phases(subdir_arr=name,sfrhist_w_grid=True,single_sim_plots=False,
#                     equal_bins=False,ngas_only=True,tmax_postmrg=0.25,grid_nh_res=hires_grid_res)




