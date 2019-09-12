##  projection plotting for particle data

import numpy as np
import ctypes
import colorsys
import math
import plotting.colors as colors
#import plotting.colors

import util
import os 


import matplotlib
import visualization.colors as viscolors

def contour_makepic( x, y, z, hsml, weight, 
        weight2=0, weight3=0, 
        xlen = 1, 
        pixels = 128, set_aspect_ratio = 1.0, 
        set_maxden = 1.0e-1, 			## (gadget units, 10^10 msun/kpc^2 = 10^4 msun/pc^2)
        set_dynrng = 1.0e4, 
        set_percent_maxden = 0, set_percent_minden = 0 ):
        
    ## set some basic initial values
    xypixels=pixels; xpixels=xypixels; ypixels = np.around(float(xpixels)*set_aspect_ratio).astype(int);
    ylen=xlen*set_aspect_ratio;
    ma = set_maxden; mi = ma / set_dynrng;
    xmin=-xlen; xmax=xlen; ymin=-ylen; ymax=ylen;

    ## load the routine we need
    exec_call=util.dir.c_routines_dir()+'/SmoothedProjPFH/allnsmooth.so'
    smooth_routine=ctypes.cdll[exec_call];

    ## make sure values to be passed are in the right format
    N=checklen(x); x=fcor(x); y=fcor(y); M1=fcor(weight); M2=fcor(weight2); M3=fcor(weight3); H=fcor(hsml)
    xpixels=np.int(xpixels); ypixels=np.int(ypixels)
    ## check for whether the optional extra weights are set
    NM=1; 
    if(checklen(M2)==checklen(M1)): 
        NM=2; 
        if(checklen(M3)==checklen(M1)):
            NM=3;
        else:
            M3=np.copy(M1);
    else:
        M2=np.copy(M1);
        M3=np.copy(M1);
    ## initialize the output vector to recieve the results
    XYpix=xpixels*ypixels; MAP=ctypes.c_float*XYpix; MAP1=MAP(); MAP2=MAP(); MAP3=MAP();
    ## main call to the imaging routine
    smooth_routine.project_and_smooth( \
        ctypes.c_int(N), \
        vfloat(x), vfloat(y), vfloat(H), \
        ctypes.c_int(NM), \
        vfloat(M1), vfloat(M2), vfloat(M3), \
        cfloat(xmin), cfloat(xmax), cfloat(ymin), cfloat(ymax), \
        ctypes.c_int(xpixels), ctypes.c_int(ypixels), \
        ctypes.byref(MAP1), ctypes.byref(MAP2), ctypes.byref(MAP3) );
    ## now put the output arrays into a useful format 
    MassMap1=np.ctypeslib.as_array(MAP1).reshape([xpixels,ypixels]);
    MassMap2=np.ctypeslib.as_array(MAP2).reshape([xpixels,ypixels]);
    MassMap3=np.ctypeslib.as_array(MAP3).reshape([xpixels,ypixels]);
 
    # set boundaries and do some clipping
    MassMap = np.copy(MassMap1);
    print "MassMap : max: ", np.max(MassMap), "   min: ", np.min(MassMap)
    if (set_percent_maxden !=0) or (set_percent_minden !=0):
        print 'percent max/min = ',set_percent_maxden,set_percent_minden
        Msort=np.sort(MassMap);
        if (set_percent_maxden != 0): ma=Msort[set_percent_maxden*float(checklen(MassMap)-1)];
        mi=ma/set_dynrng;
        if (set_percent_minden != 0): mi=Msort[set_percent_minden*float(checklen(MassMap)-1)];
        if (set_percent_maxden == 0): ma=mi*set_dynrng;
        ok=(Msort > 0.) & (np.isnan(Msort)==False)
        if (mi <= 0) or (np.isnan(mi)): mi=np.min(Msort[ok]);
        if (ma <= 0) or (np.isnan(ma)): ma=np.max(Msort[ok]);
    print "Clipping at   ma= ", ma, " mi= ", mi
    MassMap[MassMap < mi]=mi; MassMap[MassMap > ma]=ma;

    # now set colors
    cols=255. # number of colors
    Pic = (np.log(MassMap/mi)/np.log(ma/mi)) * (cols-3.) + 2.
    Pic[Pic > 255.]=255.; Pic[Pic < 1]=1.;
    backgrd = np.where((Pic<=2) | (np.isnan(Pic)))
    
    if (NM>1):
        return MassMap1,MassMap2,MassMap3, Pic
    else:
        return MassMap, Pic


def make_threeband_image_process_bandmaps(r,g,b, \
                                          dont_make_image=0, maxden=0, dynrange=0, pixels=720, \
                                          color_scheme_nasa=1, color_scheme_sdss=0 , \
                                          filterset = ['r','g','b'], **kwargs ):
    
    ## now clip the maps and determine saturation levels
    cmap_m=np.zeros((util.cast.checklen(r[:,0]),util.cast.checklen(r[0,:]),3),dtype='f');
    cmap_m[:,:,0]=r; cmap_m[:,:,1]=g; cmap_m[:,:,2]=b;
    if (dont_make_image==1): return cmap_m;
    
    if (maxden<=0):
        f_saturated=0.005 ## fraction of pixels that should be saturated
        x0=util.cast.int_round( f_saturated * (np.float(util.cast.checklen(r)) - 1.) );
        for rgb_v in [r,g,b]:
            rgbm=sorted(np.reshape(rgb_v,rgb_v.size),reverse=True)       #single_vec_sorted(rgb_v,reverse=True);
            if(rgbm[x0]>maxden): maxden=rgbm[x0]
        print "NO maxden VALUE WAS FOUND.  SETTING TO {:16.8f}".format( maxden )




    if (dynrange<=0):
        f_zeroed=0.1 ## fraction of pixels that should be black
        x0=util.cast.int_round( f_zeroed * (np.float(util.cast.checklen(r)) - 1.) ); minden=np.max(r);
        
        rgbm=sorted(np.reshape(r+g+b,r.size),reverse=False)      #single_vec_sorted(r+g+b,reverse=False);
        if(rgbm[x0]<minden): minden=rgbm[x0];
        #for rgb_v in [r,g,b]:
        #    rgbm=single_vec_sorted(rgb_v,reverse=False);
        #    if(rgbm[x0]<minden): minden=rgbm[x0];
        if (minden<=0):
            minden = np.min(np.concatenate((r[r>0.],g[g>0.],b[b>0.])));
        dynrange = maxden/minden;

        print "NO dynrange VALUE WAS FOUND.  SETTING TO {:16.8f}".format( dynrange )



    ## now do the color processing on the maps
    maxnorm=maxden; minnorm=(maxden/dynrange);
    print 'maxnorm == ',maxnorm,' dynrange == ',dynrange,' minnorm == ',minnorm;
    
    i = (r+g+b)/3.
    f_i = np.log10(i/minnorm) / np.log10(maxnorm/minnorm);
    f_i[i>=maxnorm]=1.; f_i[i<=minnorm]=0.
    
    if (color_scheme_sdss==1):
        q=9.; alpha=0.3;
        f_i = np.arcsinh( alpha * q * (i/minnorm) ) / q;
        wt=f_i/i; r*=wt; g*=wt; b*=wt;
    if (color_scheme_nasa==1):
        r = np.log(r/minnorm) / np.log(maxnorm/minnorm);
        g = np.log(g/minnorm) / np.log(maxnorm/minnorm);
        b = np.log(b/minnorm) / np.log(maxnorm/minnorm);
    
    ## rescale to saturation limit
    bad=(i<=0.); maxrgb=0.;
    if (util.cast.checklen(i[bad])>0): r[bad]=0.; g[bad]=0.; b[bad]=0.;
    f_saturated=0.0004  ## fraction of pixels that should be saturated
    f_saturated=0.0001  ## fraction of pixels that should be saturated
    x0=util.cast.int_round( f_saturated * (np.float(util.cast.checklen(r)) - 1.) );
    for rgb_v in [r,g,b]:
        rgbm=sorted(np.reshape(rgb_v,rgb_v.size),reverse=True)      #single_vec_sorted(rgb_v,reverse=True);
        if(rgbm[x0]>maxrgb): maxrgb=rgbm[x0]
    if (maxrgb > 1.): r/=maxrgb; g/=maxrgb; b/=maxrgb;
    ## rescale to 256-colors to clip the extremes (rescales back to 0-1):
    max_c=255; min_c=2;
    r=clip_256(r,max=max_c,min=min_c);
    g=clip_256(g,max=max_c,min=min_c);
    b=clip_256(b,max=max_c,min=min_c);
    image24=np.zeros((util.cast.checklen(r[:,0]),util.cast.checklen(r[0,:]),3),dtype='f');
    image24[:,:,0]=r; image24[:,:,1]=g; image24[:,:,2]=b;
    
    ## ok have r, g, b -- really just three re-scaled maps. no reason they
    ##   have to map to r, g, b colors: use the filter set given to map them ::
    image24_new = 0.*image24
    colors.load_my_custom_color_tables();

    for i in [0,1,2]:
        im=image24[:,:,i]
        if filterset[i]=='r': image24_new[:,:,0] = im
        if filterset[i]=='g': image24_new[:,:,1] = im
        if filterset[i]=='b': image24_new[:,:,2] = im
        if (filterset[i] != 'r') & (filterset[i] != 'g') & (filterset[i] != 'b'):
            my_cmap = matplotlib.cm.get_cmap(filterset[i])
            rgb_im = my_cmap(im)
            image24_new += rgb_im[:,:,0:3] ## dropping the alpha channel here!
    image24 = image24_new
    
    return image24, cmap_m; ## return both processed image and massmap



def layer_band_images(ims, maps):
    print "trying to layer band images..."
    matplotlib.pyplot.imshow(0.*ims[:,:,0]); ## zero out background

    viscolors.load_my_custom_color_tables();
    print "have we made it here?"
    nx=ims[:,0,0].size; ny=ims[0,:,0].size;
    im_new = np.zeros((nx,ny,3))
    map_cum = np.zeros((nx,ny))

    map_sum = 0.*map_cum
    for i in range(ims.shape[2]):
        map_sum += maps[:,:,i]
    map_min = 0.5 * np.min(map_sum[map_sum > 0])

    for i in range(ims.shape[2]):
        im = ims[:,:,i]
        map = maps[:,:,i]

        #im_0=im/np.max(im); # if want more saturated images
        ####alpha_im=maps[:,:,i]/map_sum; ## deprected
        cm=pick_custom_cmap(i);
        my_cmap=matplotlib.cm.get_cmap(cm); ## load cmap
        #rgb_im=my_cmap(im_0); ## get rgba values of image as mapped by this cmap

        rgb_im = my_cmap(im)[:,:,0:3]
        for j in [0,1,2]:
            im_new[:,:,j] = (map_cum*im_new[:,:,j] + map*rgb_im[:,:,j]) / (map_min + map_cum + map)
        map_cum += map

        #rgb_im[:,:,3]=alpha_im; ## replace alpha channel for this image
        #rgb_im[:,:,3]=0.*alpha_im+1.;
        #matplotlib.pyplot.imshow(rgb_im); ## plot it, with appropriate (new) alpha channel

    return im_new



def include_lighting( image24, massmap ):
    light = viscolors.CustomLightSource(azdeg=0,altdeg=65)
    if (len(massmap.shape)>2):
        ## do some clipping to regulate the lighting:
        elevation = massmap.sum(axis=2)
        minden = maxden / dynrange
        print " "
        print "In lighting routine the max/min/mean elevations are {:.2f}/{:.2f}/{:.2f}".format( elevation.min(), elevation.max(), elevation.mean() )
        print "In lighting routine the min/max dens are {:.2f}/{:.2f}".format( minden, maxden )
        print " "
        elevation = (elevation - minden) / (maxden - minden)
        elevation[elevation < 0.] = 0.
        elevation[elevation > 1.] = 1.
        elevation *= maxden
        grad_max = maxden / 5.
        grad_max = maxden / 6.
        #image24_lit = light.shade_rgb(image24, massmap.sum(axis=2))
        image24_lit = light.shade_rgb(image24, elevation, vmin=-grad_max, vmax=grad_max)
    else:
        image24_lit = light.shade(massmap, matplotlib.cm.get_cmap('hot'))   #reversed args          # ptorrey -- important
    image24 = image24_lit


def fcor(x):
    return np.array(x,dtype='f',ndmin=1)
def vfloat(x):
    return x.ctypes.data_as(ctypes.POINTER(ctypes.c_float));
def cfloat(x):
    return ctypes.c_float(x);
def checklen(x):
    return len(np.array(x,ndmin=1));

def clip_256(x,max=255,min=2):
    x = x*(max-min) + min;
    x[x >= max]=max;
    x[x <= min]=min;
    x[np.isnan(x)]=min;
    x /= 256.;
    return x;



def pick_custom_cmap(i):
    cm='hot'
    if(i==0): cm='heat_purple'
    if(i==1): cm='heat_green'
    if(i==2): cm='heat_blue'
    if(i==3): cm='heat_yellow'
    if(i==4): cm='heat_red'
    if(i==5): cm='heat_orange'
    if(i==6): cm='heat_redyellow'
    if(i==7): cm='pink' # whites out fairly quickly
    if(i==8): cm='bone' # pretty nice for white-ish (X-ray like colors)
    if(i==9): cm='copper' # orange-y
    if(i==10): cm='gray' # basic greyscale

    if(i==11): cm='spring'
    if(i==12): cm='summer'
    if(i==13): cm='winter'
    if(i==14): cm='autumn'
    if(i==15): cm='gist_earth'
    if(i==16): cm='Blues_r'
    if(i==17): cm='Greens_r'
    if(i==18): cm='Oranges_r'
    if(i==19): cm='Purples_r'
    if(i==20): cm='RdPu_r'
    if(i==21): cm='Reds_r'


    if(i==0): cm='heat_orange'
    if(i==0): cm='heat_redyellow'
    if(i==1): cm='heat_green'
    if(i==2): cm='heat_purple'
    return cm;


