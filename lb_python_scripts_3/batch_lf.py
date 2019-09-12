import bat_agn_lf

#for mcut in [10.0,9.7,9.5]:
for mcut in [9.7,9.5]:
    bat_agn_lf.calc(use_snr_limited_vmax=True,logM_cut=mcut)
    bat_agn_lf.calc(use_snr_limited_vmax=True,logM_cut=mcut,plot_logvals=False)
    bat_agn_lf.calc(logM_cut=mcut)
    bat_agn_lf.calc(logM_cut=mcut,plot_logvals=False)

bat_agn_lf.calc(use_snr_limited_vmax=True,cat='a12')
bat_agn_lf.calc(use_snr_limited_vmax=True,cat='a12',plot_logvals=False)
bat_agn_lf.calc(cat='a12')
bat_agn_lf.calc(cat='a12',plot_logvals=False)

