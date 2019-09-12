import make_3color
import glob,sys,re,pyfits
import numpy as np

path='/n/hernquistfs2/lblecha/sunrise_images'
print(path)
#subdir='A1A0rx10/gal2_center/zoom_movie/cam2/cropped'
#subdir='A1A0rx10/gal2_center/zoom2_movie/cam2/'
subdir='A1A0rx10/gal2_center/zoom2/'

#files = glob.glob("%s/%s/sdss*jpg"%(path,subdir))
#print(files)
#files.sort()
#print(files)
## this won't work, don't have all the snaps in this file. 
##tmpdata = np.loadtxt(path+'/bh_sep.txt',unpack=True)
##snaps = [int(tmp) for tmp in tmpdata[0]]
#fsnaps = [np.int(re.search('(?<=gri_)\d+',f).group(0))
#          for f in files]

#need to load bb files and do this instead:
snaptime = []
#fbnames=glob.glob(path+"/A1A0rx10/gal2_center/broadband*fits")
fbnames=glob.glob(path+"/A1A0rx10/gal2_center/zoom2/broadband*fits")
fbnames.sort()
#print(fbnames)
fsnaps = [np.int(re.search('(?<=broadband_)\d+',f).group(0))
          for f in fbnames]
#fst=glob.glob(path+"/A1A0rx10/gal2_center/snaptimes.txt")
fst=glob.glob(path+"/A1A0rx10/gal2_center/zoom2/snaptimes.txt")
if fst:
    snaptime=np.loadtxt(fst[0],unpack=True)
else:
    print("reading snaptimes from broadband files and creating file snaptimes.txt...")
    #fout=open(path+"/A1A0rx10/gal2_center/snaptimes.txt","w")
    fout=open(path+"/A1A0rx10/gal2_center/zoom2/snaptimes.txt","w")
    for fb in fbnames:
        bb = pyfits.open(fb)
    #for s in fsnaps:
    #    bb = pyfits.open(path+"/A1A0rx10/gal2_center/broadband_%.3d"%s)
        t = bb['SFRHIST'].header['SNAPTIME']/1.0e9
        fout.write("%g\n"%t)
        snaptime = snaptime + [t]
        #print fb,snaptime
        bb.close()
    fout.close()
    
#print(snaptime)
print("processing these snaps:")
print(fsnaps)
if len(fsnaps) != len(snaptime):
    print("Error: mismatch between broadband and jpg files.")
    sys.exit()

#print(snaps[mask])
#mask=np.in1d(snaps,fsnaps)
#print(len(mask),len(snaps),len(fsnaps))
#snaps=snaps[mask]
#time=time[mask]
#time = list(tmpdata[1])

timelbls = ['t=%.2f Gyr'%t for t in snaptime]
imgtitle = ['Merger simulation']*len(snaptime)
##print(snaps)
##print(files)
##print(time)
#print(timelbls)
#print(lbl_line2[0],len(lbl_line2))
#sys.exit()

### test
#make_3color.multisnaps(subdir=subdir,cam_list=[2],snaps=[5],
#                       plbls=imgtitle,plbls2=timelbls,plotscale=False)
##                       plbls=timelbls,plbls2=lbl_line2,plotscale=False)
make_3color.multisnaps(subdir=subdir,cam_list=[2],snaps=fsnaps,
                       plbls=imgtitle,plbls2=timelbls,plotscale=False)
