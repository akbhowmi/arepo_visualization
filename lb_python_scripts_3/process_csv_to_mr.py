import numpy as np
import re
from collections import OrderedDict
import matplotlib.pyplot as plt

def loadfile():

    fname = 'dwarfinfointableform_rev2_All_Dwarf_Data.csv'

    chandra_list = ['2MASXJ03534246+3714077','IRAS04124-0803','2MASXJ05054575-2351139','ESO506-G027','Mrk477','NGC5995','NGC4102','NGC4138','IC2461','MCG+10-17-061']

    with open(fname,"r") as fp:
        
        data = []
        source_count = 0
        for line in fp:            
            cols = line.strip().split(',')
            if len(cols) != 24:
                print "Skipping line with invalid format:"
                print line
                continue
            if cols[0] == 'Name' or cols[0] == '' or re.match('Sources',cols[0]):
                print "Skipping line with invalid format:"
                print line
                continue
            cols = [cols[i].replace(' ','') if (i!=9 and i!=23 and i>0)
                    else cols[i] for i in range(24)]
            source = OrderedDict()
            #source = {}
            source['catID'] = '%0.3d'%source_count
            source_count += 1
            if re.match('.*\(.*\)',cols[0]):
                tmp = re.split('\(',cols[0])
                if len(tmp) != 2: 
                    print "Error: could not separate alt name."
                    return -1
                source['name'] = safeassign(re.sub('\s','',tmp[0]),mytype=np.str)
                source['altname'] = safeassign(re.sub('\s','',tmp[1].strip().strip(')')),mytype=np.str)
            else:
                source['name'] = safeassign(re.sub('\s','',cols[0]),mytype=np.str)
                source['altname'] = 'None'
            source['ra'] =  safeassign(cols[1])
            source['dec'] =  safeassign(cols[2])
            source['Lx'] = safeassign(cols[3])
            source['Lbol'] = safeassign(cols[4])
            source['Bmag'] = safeassign(cols[5])
            source['AGNtype'] = safeassign(cols[6])
            source['z'] = safeassign(cols[7])
            source['optspec'] = safeassign(cols[8],mytype=np.str)
            source['morphnotes'] = safeassign(cols[9],mytype=np.str)
            if re.match('\d',cols[10]):
                tmp = re.split(';',cols[10])
                mbh_tmparr=[safeassign(elem) for elem in tmp]
                source['bhmass_low'] = min(mbh_tmparr)
                source['bhmass_hi'] = max(mbh_tmparr)
            else:
                source['bhmass_low'] = -1.0
                source['bhmass_hi'] = -1.0
            #source['LEdd'] = cols[11]
            #source['fEdd'] = cols[12]
            source['mstarK11'] = safeassign(cols[13])
            if re.match('\d',cols[14]):
                tmp = re.split(';',cols[14])
                ms_tmparr=[safeassign(elem) for elem in tmp]
                source['mstarother_low'] = min(ms_tmparr)
                source['mstarother_hi'] = max(ms_tmparr)
            else:
                source['mstarother_low'] = -1.0
                source['mstarother_hi'] = -1.0
            if re.match('\d',cols[15]):
                tmp = re.split(';',cols[15])
                vdisp_tmparr=[]
                for elem in tmp:
                    if re.match('.*\+\/\-',elem):
                        tmp2 = re.split('\+\/\-',elem)
                        if len(tmp2) != 2:
                            print "Error separating vel disp errs: ",elem                            
                            return -1
                        vdisp_tmparr.append((safeassign(tmp2[0]),safeassign(tmp2[1])))
                    else:
                        vdisp_tmparr.append((safeassign(elem),-1.0))
                if len(tmp) == 1:
                    source['vdisp_low'] = vdisp_tmparr[0][0]
                    source['vdisp_low_err'] = vdisp_tmparr[0][1]
                    source['vdisp_hi'] = source['vdisp_low']
                    source['vdisp_hi_err'] = source['vdisp_low_err']
                else:
                    source['vdisp_low'] = min([elem[0] for i,elem in enumerate(vdisp_tmparr)])
                    source['vdisp_low_err'] = ([v[1] for i,v in enumerate(vdisp_tmparr)
                                                if v[0]==source['vdisp_low']])[0]
                    source['vdisp_hi'] = max([elem[0] for i,elem in enumerate(vdisp_tmparr)])
                    source['vdisp_hi_err'] = ([v[1] for i,v in enumerate(vdisp_tmparr)
                                               if v[0]==source['vdisp_hi']])[0]
            else: 
                source['vdisp_low'] = -1.0
                source['vdisp_low_err'] = -1.0
                source['vdisp_hi'] = -1.0
                source['vdisp_hi_err'] = -1.0
            source['pctAGNrband'] = safeassign(cols[16].strip().strip('%'))
            source['dist'] = safeassign(cols[17])
            source['gmagcorr'] = safeassign(cols[18])
            source['rmagcorr'] = safeassign(cols[19])
            source['gminusr'] = safeassign(cols[20])
            source['hasHST'] = safeassign(cols[21],mytype=np.str)
            source['hasChandra'] = safeassign(cols[22],mytype=np.str)
            source['comments'] = safeassign(cols[23],mytype=np.str)
            data.append(source)

    print "number of sources read = %d"%(len(data))
    print "should match this: %d"%(source_count+1)
    print "\n1st source:",data[0]
    print "\n source names:",([src['name'] for src in data ])
    #print "\n source Bmag:",([src['Bmag'] for src in data ])
    #print "\n source mstarK11:",([src['mstarK11'] for src in data ])

    ## perform check:
    if len([src for src in data if ((src['Bmag']>-19.5) and (src['Bmag']!=-1)) or
            ((src['mstarK11']<10.0) and (src['mstarK11']!=-1))]) != len(data):
        print "Error in matching initial selection criteria."
        print "Problem sources:"
        print [src['name'] for src in data if ( ((src['Bmag']==-1) or (src['Bmag']<=-19.5)) and 
                                                ((src['mstarK11']>=10.0) or (src['mstarK11']==-1)) )]
        return

    with open('dwarf_sample_comments.txt',"w") as fpcomm:
        fpcomm.write('catID   Comments\n')
        for src in data:
            fpcomm.write('%s     %s\n'%(src['catID'],src['comments']))
    with open('dwarf_sample_morph_notes.txt',"w") as fpmorph:
        fpmorph.write('catID   Morphology notes\n')
        for src in data:
            fpmorph.write('%s     %s\n'%(src['catID'],src['morphnotes']))

    with open('dwarf_sample_all.dat',"w") as fpall:
        fpall.write(' '.join([key for key in data[0]
                    if (key!='comments') and (key!='morphnotes')])+"\n")
        for src in data:
            fpall.write(' '.join([str(src[key]) for key in src
                                  if (key!='comments') and (key!='morphnotes')]))
            fpall.write("\n")


    print "\n all NGC objects: ",[src['name'] for src in data if re.match('NGC.*',src['name']) ]

    print "\n all obj with vdisp: ",[src['name'] for src in data if (((src['vdisp_low']!=-1) or (src['vdisp_low']!=-1)) and src['mstarK11']<10.0)]
    print [src['vdisp_low'] for src in data if ((src['vdisp_low']!=-1) or (src['vdisp_low']!=-1))]
    print [src['vdisp_hi'] for src in data if ((src['vdisp_low']!=-1) or (src['vdisp_low']!=-1))]

    print "\n# sources with mstarK11<9.7:",len([src for src in data if ((src['mstarK11']<9.7) and (src['mstarK11']!=-1))])
    print "\n# sources with mstarK11<9.7 and vdisp_hi<100:",len([src for src in data if ((src['mstarK11']<9.7) and (src['mstarK11']!=-1) and (src['vdisp_hi']<100))])
    print [src['name'] for src in data if ((src['mstarK11']<9.7) and (src['mstarK11']!=-1) and (src['vdisp_hi']<100))]
    print "\n# sources with no mstarK11 and Bmag>-19.5:",len([src for src in data if ((src['mstarK11']==-1) and (src['Bmag']>-19.5) and (src['Bmag']!=-1))])
    print "\n# sources with no mstarK11 and Bmag>-19.5 and vdisp_hi<100:",len([src for src in data if ((src['mstarK11']==-1) and (src['Bmag']>-19.5) and (src['Bmag']!=-1) and (src['vdisp_hi']<100))])
    print [src['name'] for src in data if ((src['mstarK11']==-1) and (src['Bmag']>-19.5) and (src['Bmag']!=-1) and (src['vdisp_hi']<100))]


    print "\n",[(src['name'],src['mstarother_low'],src['mstarother_hi']) for src in data if ((src['mstarK11']==-1) and ((src['mstarother_low']!=-1) or (src['mstarother_hi']!=-1)))]

    print "\nSources excluded based on vel disp alone (i.e., need to be sure of the vel disp values!)\n"
    for src in data:
        if ( (((src['mstarK11']<9.7) and (src['mstarK11']!=-1)) or
              ((src['Bmag']>-19.5) and (src['Bmag']!=-1) and
               ((src['mstarK11']<9.7) or (src['mstarK11']==-1)))) and
             (src['vdisp_hi']>=100.0)):
            print src['name'],src['Bmag'],src['mstarK11'],src['vdisp_low'],src['vdisp_hi']


    print "\nMost concerning sources (no K11 masses OR vel disp:):\n"
    for src in data:
        if ((src['mstarK11']==-1) and 
            (src['Bmag']>-19.5) and (src['Bmag']!=-1) and
            (src['vdisp_hi']==-1) and (src['vdisp_low']==-1)):
            print src['name'],src['Bmag']

    print "\nSources with no K11 masses that do have vel disp (double check these too!):\n"
    for src in data:
        if ((src['mstarK11']==-1) and 
            (src['Bmag']>-19.5) and (src['Bmag']!=-1) and
            (src['vdisp_hi']!=-1) and (src['vdisp_hi']<100)):
            print src['name'],src['Bmag'], src['vdisp_low'],src['vdisp_hi']

    print "\nSources with K11 masses that don't have vel disp (double check literature):\n"
    for src in data:
        if ((src['mstarK11']!=-1) and (src['mstarK11']<9.7) and  
            (src['vdisp_hi']==-1)):
            print src['name'],src['mstarK11'],src['Bmag']


def load_BASS_file(path='/n/home00/lblecha/dwarf_agn_data'):

    fname = '%s/dwarf_BASS_data.csv'%path

    with open(fname,"r") as fp:
        
        data = []
        source_count = 0
        for line in fp:            
            cols = line.strip().split(',')
            if source_count==0: "column names: ",cols
            if len(cols) != 22 or cols[0]=='' or cols[0]=='BAT index':
                print "Skipping line with invalid format:"
                print line
                continue

            source = OrderedDict()
            #source = {}
            source['BATindex'] = safeassign(cols[0],mytype=np.int)
            source['BATname'] = safeassign(re.sub('\s','',cols[1]),mytype=np.str)
            source['Name'] = safeassign(re.sub('\s','',cols[2]),mytype=np.str)
            source['ra'] =  safeassign(cols[3])
            source['dec'] =  safeassign(cols[4])
            source['MBH_Hbeta'] = safeassign(cols[5])
            source['MBH_vdisp'] = safeassign(cols[6])
            source['MBH_vdisp_err'] = safeassign(cols[7])
            source['FWHM_Halpha'] = safeassign(cols[8])
            source['FWHM_Halpha_err'] = safeassign(cols[9])
            source['Halpha_broad'] = safeassign(cols[10])
            source['Halpha_broad_err'] = safeassign(cols[11])
            source['Halpha'] = safeassign(cols[12])
            source['Halpha_err'] = safeassign(cols[13])
            source['Hbeta_broad'] = safeassign(cols[14])
            source['Hbeta_broad_err'] = safeassign(cols[15])
            source['corr_Hbeta'] = safeassign(cols[16])*1.0e-15
            source['Hbeta_err'] = safeassign(cols[17])
            source['corr_OIII'] = safeassign(cols[18])*1.0e-15
            source['OIII_err'] = safeassign(cols[19])
            source['corr_NII'] = safeassign(cols[20])*1.0e-15
            source['NII_err'] = safeassign(cols[21])
            data.append(source)
            source_count += 1

    print "number of sources read = %d"%(len(data))
    print "should match this: %d"%(source_count)
    print "\n1st source:",data[0]
    print "\n source names:",([src['Name'] for src in data ])
    #print "\n narrow Halpha:",(["%g "%(src['Halpha']) for src in data ])
    #print "\n NII:",(["%g "%src['corr_NII'] for src in data ])
    #print "\n narrow Hbeta:",(["%g "%(src['corr_Hbeta']) for src in data ])
    #print "\n OIII:",(["%g "%src['corr_OIII'] for src in data ])
    
    NII_Ha = np.repeat(np.nan,source_count)
    OIII_Hb = np.repeat(np.nan,source_count)
    for i,src in enumerate(data):
        if (src['corr_NII']>0 and src['Halpha']>0):
            NII_Ha[i] = np.log10(src['corr_NII']/src['Halpha']) 
        if (src['corr_OIII']>0 and src['corr_Hbeta']>0):
            OIII_Hb[i] = np.log10(src['corr_OIII']/src['corr_Hbeta']) 

    print "\nOf %d sources,"%source_count
    print "%d are missing [NII] and/or Halpha."%(NII_Ha[NII_Ha!=NII_Ha].size)
    print "%d are missing [OIII] and/or Hbeta."%(OIII_Hb[OIII_Hb!=OIII_Hb].size)
    print "%d are missing both line ratios."%(OIII_Hb[(OIII_Hb!=OIII_Hb)&
                                                      (NII_Ha!=NII_Ha)].size)
    ## Kewley et al. 2006 BPT
    x1 = np.arange(-1.28,0.0,0.05)
    y1 = 0.61/(x1-0.05)+1.3 
    x2 = np.arange(-2.2,0.4,0.05)
    y2 = 0.61/(x2-0.47)+1.19
    #y3 = 0.61/(x-0.05)+1.19

    plt.clf()
    plt.xlim(-2.2,0.8)
    plt.ylim(-1.2,1.5)
    plt.xlabel("log([NII]/Halpha)")
    plt.ylabel("log([OIII]/Hbeta)")
    plt.plot(NII_Ha,OIII_Hb,'o')
    plt.plot(x1,y1,'--')
    plt.plot(x2,y2)

    with open('%s/dwarf_BASS_data.dat'%path,"w") as fpall:
        fpall.write(' '.join([key for key in data[0]])+'\n')
        for src in data:
            fpall.write(' '.join([str(src[key]) for key in src]))
            fpall.write("\n")


def safeassign(val,mytype=np.float,default=-1,exception=(AssertionError,ValueError)):

    if mytype==np.str and default==-1: default='None'
    try:
        assert val != ''
        return mytype(val)
    except exception:
        return mytype(default)
