
## used to make the plot q05fg03rx10_sed_compare_snaps_agnx0x1_tot_cam3.pdf
quick_seds.compare_manual(rundir='q0.5_fg0.3_allrx10_sunruns',
                          compare_subdirs=['all','all_agnx0'],compare_fbase='',
                          snap_arr=[176,191,205,219],sed_type='tot',ylog=True,
                          lambda_units='microns',xlim=(3e-7,1e-4),flux_units='nuFnu',
                          ylim=(3e35,8e38),oplot_wise=True,single_cam=3,
                          single_snap_plots=False,snap_titles=False,
                          snap_labels=['a','b','c','d'])

## lowres analogue:
quick_seds.compare_manual(rundir='q0.5_fg0.3_sunruns',
                          compare_subdirs=['new_fiducial_all','new_fiducial_all_agnx0'],
                          compare_fbase='',snap_arr=[180,191,205,219],sed_type='tot',
                          ylog=True,lambda_units='microns',xlim=(3e-7,1e-4),
                          flux_units='nuFnu',ylim=(2e35,8e38),oplot_wise=True,
                          single_cam=3,single_snap_plots=False,snap_titles=False,
                          snap_labels=['a_lowres','b_lowres','c_lowres','d_lowres'])
