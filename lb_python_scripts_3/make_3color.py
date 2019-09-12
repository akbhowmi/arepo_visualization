import uvoir_img, glob, re
import numpy as np
import sys, os, datetime

def multisnaps(path='/n/hernquistfs2/lblecha/sunrise_images/',subdir='A1A0rx10',
               filter_dir='sdss/',fname='gri',scale=(),cam_list=[3],snaps=[],
               plbls=[],plbls2=[],label_only=False,plotscale=False):


    if fname not in ('gri','ugz','w123'):
        print('Error: filter set %s is undefined.'%fname)
        sys.exit()

    if ( (fname in ('gri','ugz') and 'sdss' not in filter_dir) or 
         (fname=='w123' and 'wise' not in filter_dir) ):
        print('Error: filter_dir %s does not correspond to fname %s.'%(filter_dir,fname))
        sys.exit()

    if not snaps:
        bbfiles = glob.glob(path+subdir+'/broadband*fits')
        snaps = np.sort(np.array([np.int(re.search('(?<=broadband_)\d+',f).group(0)) 
                                  for f in bbfiles]))
    print('Making 3-color %s images for the following cameras:')
    print(cam_list)
    print('And the following snaps:')
    print(snaps)

    if not plbls:
        plbls=np.tile([''],len(snaps))
    #    plbls=np.array([str(s) for s in snaps])
    if not plbls2:
        plbls2=np.tile([''],len(snaps))


    ## good for ugz (esp. scale 'a' & 'b')
    myscale_aa=(1.25,0.92,1.25)
    myscale_a=(1.25,0.9,1.25)
    myscale_b=(1.15,0.8,1.15)
    myscale_c=(1.85,1.5,1.85)
    myscale_d=(1.55,1.2,1.55)

    ## good for gri
    myscale_e=(1.7,1.45,1.32)
    ###myscale_e=(1.55,1.3,1.2)
    ###myscale_e=(1.65,1.35,1.25)
    ###myscale_e=(1.5,1.5,1.5)
    ###myscale_e=(1.04,1.0,0.94)
    ###myscale_e=(1.8,1.55,1.45)
    ###myscale_e=(2.0,1.75,1.65)
    ##myscale_f=(1.75,1.45,1.2)
    myscale_f=(1.25,0.96,0.8)
    ##myscale_f=(1.5,1.5,1.5)
    
    ## good for wise w123
    myscale_wa=(100.0,100.0,30.0)
    

    filters = { 'gri':('g_SDSS','r_SDSS','i_SDSS'), 
                'ugz':('u_SDSS','g_SDSS','z_SDSS'),
                'w123':('WISE-W1','WISE-W2','WISE-W3') }
    if not scale:
        scale = { 'gri':myscale_f,
                  'ugz':myscale_a,
                  'w123':myscale_wa }

    print scale[fname][0]
    print scale[fname][1]
    print scale[fname][2]
    with open('%s/%s/make_3color.out'%(path,subdir),'a') as fout:
        fout.write('\n\n%s\n'%str(datetime.datetime.now()))
        fout.write('In make_3color.multisnaps():')
        fout.write('Making %s images with scale (%g,%g,%g) '%(fname,scale[fname][0],
                                                               scale[fname][1],scale[fname][2]))
        fout.write('\nfor the following cameras: ')
        for c in cam_list: fout.write("%d "%c)
        fout.write('and the following snaps:\n')
        for s in snaps: fout.write("%d "%s)


    for i,snap in enumerate(snaps):    
        print(i, snap, plbls[i], plbls2[i])
        print(i, snap)
        for jcam in cam_list:        
            uvoir_img.make_broadband(path=path+subdir,snaps=[snap],cam=jcam,
                                     filter_dir=filter_dir,filters=filters[fname],fname=fname,
                                     scale=scale[fname],plotlabel=plbls[i],plotlabel2=plbls2[i],
                                     label_only=label_only,plotscale=plotscale)



if __name__ == "__main__":

    multisnaps()




#snaps=np.append(5,np.arange(10,100,10))
#filter_dir='wise/'
##snaps=np.array([176,191,205,219])
#plbls=np.array(['a','b','c','d'])

#subdir='A1E0rx10'
#snaps=np.arange(199,218)
#plbls=np.array([str(s) for s in snaps])

#subdir='A1A0'
#subdir='A0A0'

#filter_dir=''
#subdir='A1A0rx10_old'

#filters=('g_SDSS','r_SDSS','i_SDSS')
#fname='gri'
#filters=('u_SDSS','g_SDSS','z_SDSS')
#fname='ugz'
#filters=('WISE-W1','WISE-W2','WISE-W3')
#fname='w123'


#snaps=[191] 
    
#uvoir_img.make_broadband(path=path+subdir,snaps=snaps,filter_dir=filter_dir,filters=filters,fname=fname)

#uvoir_img.make_broadband(path=path+subdir,snaps=snaps,filter_dir=filter_dir,filters=filters,fname=fname,
#                         scale=myscale_aa,extra='aa',plotlabel='b',plotscale=True)
#uvoir_img.make_broadband(path=path+subdir,snaps=snaps,filter_dir=filter_dir,filters=filters,fname=fname,
#                         scale=myscale_a,extra='a')
#uvoir_img.make_broadband(path=path+subdir,snaps=snaps,filter_dir=filter_dir,filters=filters,fname=fname,
#                         scale=myscale_b,extra='b')
#uvoir_img.make_broadband(path=path+subdir,snaps=snaps,filter_dir=filter_dir,filters=filters,fname=fname,
#                         scale=myscale_c,extra='c')
#uvoir_img.make_broadband(path=path+subdir,snaps=snaps,filter_dir=filter_dir,filters=filters,fname=fname,
#                         scale=myscale_d,extra='d')

#uvoir_img.make_broadband(path=path+subdir,snaps=snaps,filter_dir=filter_dir,filters=filters,fname=fname,
#                         scale=myscale_e,extra='e')
#uvoir_img.make_broadband(path=path+subdir,snaps=snaps,filter_dir=filter_dir,filters=filters,fname=fname,
#                         scale=myscale_f,extra='f')

#uvoir_img.make_broadband(path=path+subdir,snaps=snaps,filter_dir=filter_dir,filters=filters,fname=fname,
#                         scale=myscale_wa,extra='wa')
#uvoir_img.make_broadband(path=path+subdir,snaps=snaps,filter_dir=filter_dir,filters=filters,fname=fname)


## for A1A0rx10 pubstyle images:
##snaps=np.arange(200,230)
###snaps=np.arange(180,200)
###plbls=np.repeat('',snaps.size)
##snaps=np.array([176,191,206,220])

#filters=('u_SDSS','g_SDSS','z_SDSS')
#fname='ugz'
#for i,snap in enumerate(snaps):
#    print(i, snap, plbls[i])
#    uvoir_img.make_broadband(path=path+subdir,snaps=[snap],cam=2,
#                             filter_dir=filter_dir,filters=filters,fname=fname,
#                             scale=myscale_aa,extra='aa',plotlabel=plbls[i],plotscale=True)
#    uvoir_img.make_broadband(path=path+subdir,snaps=[snap],
#                             filter_dir=filter_dir,filters=filters,fname=fname,
#                             scale=myscale_aa,extra='aa',plotlabel=plbls[i],plotscale=True)
#    uvoir_img.make_broadband(path=path+subdir,snaps=[snap],cam=4,
#                             filter_dir=filter_dir,filters=filters,fname=fname,
#                             scale=myscale_aa,extra='aa',plotlabel=plbls[i],plotscale=True)


#original images:
#/oasis/projects/nsf/hvd115/lblecha/q0.5_fg0.3_allrx10_sunruns/old/whole_sed_lowres_newcode    
#170: t=1.574gyr
#320: t=1.721
#470: t=1.868
#620: t=2.014
# final choices for myscale:
#sdss_img.make_broadband(snaps=[170],scale=(1.25,0.9,1.25))               
#sdss_img.make_broadband(snaps=[320],scale=(1.15,0.8,1.15))
#sdss_img.make_broadband(snaps=[470],scale=(1.85,1.5,1.85))
#sdss_img.make_broadband(snaps=[620],scale=(1.55,1.2,1.55))


#closest corresponding snaps in new q0.5_fg0.3_allrx10_sunruns:
#160 - 1.564
#176 - 1.721
#191 - 1.868
#206 - 2.014

#other snaps used:
#205 -
#219 -
#220 - 2.15
#230 - 2.25
#240 - 2.35
#250 - 2.44
#260 - 2.54
#270 - 2.64

#tmrg = 2.11




