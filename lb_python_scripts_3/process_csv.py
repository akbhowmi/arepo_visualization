import numpy as np
import re
from collections import OrderedDict

def loadfile():

    fname = 'dwarfinfointableform_rev2_All_Dwarf_Data.csv'

    chandra_list = ['2MASXJ03534246+3714077','IRAS04124-0803','2MASXJ05054575-2351139','ESO506-G027','Mrk477','NGC5995','NGC4102','NGC4138','IC2461','MCG+10-17-061']

    with open(fname,"r") as fp:
        
        data = []
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
            source = OrderedDict()
            #source = {}
            if re.match('.*\(.*\)',cols[0]):
                tmp = re.split('\(',cols[0])
                if len(tmp) != 2: 
                    print "Error: could not separate alt name."
                    return
                source['name'] = re.sub('\s','',tmp[0])
                source['altname'] = re.sub('\s','',tmp[1].strip().strip(')'))
            else:
                source['name'] = re.sub('\s','',cols[0])
                source['altname'] = ''
            source['ra'] =  safeassign(cols[1])
            source['dec'] =  safeassign(cols[2])
            source['Lx'] = safeassign(cols[3])
            source['Lbol'] = safeassign(cols[4])
            source['Bmag'] = safeassign(cols[5])
            source['AGNtype'] = safeassign(cols[6])
            source['z'] = safeassign(cols[7])
            source['optspec'] = cols[8]
            source['morphnotes'] = cols[9]
            if re.match('\d',cols[10]):
                tmp = re.split(';',cols[10])
                mbh_tmparr=[safeassign(elem) for elem in tmp]
                source['bhmass_low'] = min(mbh_tmparr)
                source['bhmass_hi'] = max(mbh_tmparr)
            else:
                source['bhmass_low'] = -1
                source['bhmass_hi'] = -1
            #source['LEdd'] = cols[11]
            #source['fEdd'] = cols[12]
            source['mstarK11'] = safeassign(cols[13])
            if re.match('\d',cols[14]):
                tmp = re.split(';',cols[14])
                ms_tmparr=[safeassign(elem) for elem in tmp]
                source['mstarother_low'] = min(ms_tmparr)
                source['mstarother_hi'] = max(ms_tmparr)
            else:
                source['mstarother_low'] = -1
                source['mstarother_hi'] = -1
            if re.match('\d',cols[15]):
                tmp = re.split(';',cols[15])
                vdisp_tmparr=[]
                for elem in tmp:
                    if re.match('.*\+\/\-',elem):
                        tmp2 = re.split('\+\/\-',elem)
                        if len(tmp2) != 2:
                            print "Error separating vel disp errs: ",elem                            
                            return
                        vdisp_tmparr.append((safeassign(tmp2[0]),safeassign(tmp2[1])))
                    else:
                        vdisp_tmparr.append((safeassign(elem),-1))
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
                source['vdisp_low'] = -1
                source['vdisp_low_err'] = -1
                source['vdisp_hi'] = -1
                source['vdisp_hi_err'] = -1
            source['pctAGNrband'] = safeassign(cols[16].strip().strip('%'))
            source['dist'] = safeassign(cols[17])
            source['gmagcorr'] = safeassign(cols[18])
            source['rmagcorr'] = safeassign(cols[19])
            source['gminusr'] = safeassign(cols[20])
            source['hasHST'] = cols[21]
            source['hasChandra'] = cols[22]
            source['comments'] = cols[23]
            data.append(source)

    print "number of sources read = "+str(len(data))
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

    with open('dwarf_sample_all_withcomments.txt',"w") as fpallc:
        fpallc.write(', '.join([key for key in data[0]])+"\n")
        for src in data:
            fpallc.write(', '.join([str(src[key]) for key in src]))
            fpallc.write("\n")
    with open('dwarf_sample_all.txt',"w") as fpall:
        fpall.write(', '.join([key for key in data[0] 
                               if (key!='comments') and (key!='morphnotes')])+"\n")
        for src in data:
            fpall.write(', '.join([str(src[key]) for key in src
                                   if (key!='comments') and (key!='morphnotes')]))
            fpall.write("\n")

    strictsample1_data =  [src for src in data if (src['Bmag']>-19.5) and
                           (src['mstarK11']<10.0) and ((src['Bmag']!=-1) or (src['mstarK11']!=-1))]    
    print "\nstrict sample 1 sources (sample size="+str(len(strictsample1_data))+"):\n"
    print [src['name'] for src in strictsample1_data]
    print [src['Bmag'] for src in strictsample1_data]
    print [src['mstarK11'] for src in strictsample1_data]
    print [src['mstarother_low'] for src in strictsample1_data]
    print [src['mstarother_hi'] for src in strictsample1_data]

    sample2_data =  [src for src in data if ((src['Bmag']>-19.5) and (src['Bmag']!=-1)) or
                     ((src['mstarK11']<9.7) and (src['mstarK11']!=-1))]    
    print "\nsample 2 sources (sample size="+str(len(sample2_data))+"):\n"
    print [src['name'] for src in sample2_data]
    print [src['Bmag'] for src in sample2_data]
    print [src['mstarK11'] for src in sample2_data]
    print [src['mstarother_low'] for src in sample2_data]
    print [src['mstarother_hi'] for src in sample2_data]

 
    sample3_data =  [src for src in data if ( ( ((src['Bmag']>-19.5) and (src['Bmag']!=-1)) or
                                                ((src['mstarK11']<10.0) and (src['mstarK11']!=-1)) ) and 
                                              (src['vdisp_hi']<100.0) )]    
    print "\nsample 3 sources (sample size="+str(len(sample3_data))+"):\n"
    print [src['name'] for src in sample3_data]
    print [src['Bmag'] for src in sample3_data]
    print [src['mstarK11'] for src in sample3_data]
    print [src['mstarother_low'] for src in sample3_data]
    print [src['mstarother_hi'] for src in sample3_data]

    sample4_data =  [src for src in data if ( ( ((src['Bmag']>-19.5) and (src['Bmag']!=-1)) or
                                                ((src['mstarK11']<9.7) and (src['mstarK11']!=-1)) ) and 
                                              (src['vdisp_hi']<100.0) )]    
    print "\nsample 4 sources (sample size="+str(len(sample4_data))+"):\n"
    print [src['name'] for src in sample4_data]
    print [src['Bmag'] for src in sample4_data]
    print [src['mstarK11'] for src in sample4_data]
    print [src['mstarother_low'] for src in sample4_data]
    print [src['mstarother_hi'] for src in sample4_data]

    strictsample4_data =  [src for src in data if ( (src['Bmag']>-19.5) and (src['mstarK11']<9.7) and 
                                                    ((src['Bmag']!=-1) and (src['mstarK11']!=-1)) and 
                                                    (src['vdisp_hi']<100.0) )]    
    print "\n'Strict' sample 4: mstarK11<9.7 *and* Bmag>-19.5, plus exclude *all* with vdisp_hi<100."
    print "(sample size="+str(len(strictsample4_data))+"):\n"
    print [src['name'] for src in strictsample4_data]
    print [src['Bmag'] for src in strictsample4_data]
    print [src['mstarK11'] for src in strictsample4_data]
    print [src['mstarother_low'] for src in strictsample4_data]
    print [src['mstarother_hi'] for src in strictsample4_data]

    print "\n with newly-approved Chandra observations:"
    print [(src['name'],src['Bmag'],src['mstarK11'],src['mstarother_low'],src['mstarother_hi'],src['vdisp_hi']) for src in data if src['name'] in chandra_list]
    

    sample5_data =  [src for src in data if ( ((src['mstarK11']<10.0) and (src['mstarK11']!=-1)) or
                                              ((src['Bmag']>-19.5) and (src['Bmag']!=-1) and
                                               ((src['mstarK11']<10.0) or (src['mstarK11']==-1))) )]
    print "\nSample 5: mstarK11<10 or Bmag>-19.5, plus exclude *all* with mstarK11>10."
    print "(sample size="+str(len(sample5_data))+"):\n"
    print [src['name'] for src in sample5_data]
    print [src['Bmag'] for src in sample5_data]
    print [src['mstarK11'] for src in sample5_data]
    print [src['vdisp_hi'] for src in sample5_data]
    print [src['mstarother_low'] for src in sample5_data]
    print [src['mstarother_hi'] for src in sample5_data]

    print "\nexcluded due to K11 mass > 10.0:"
    print [src['name'] for src in data if src['mstarK11']>10.0]
    print [src['Bmag'] for src in data if src['mstarK11']>10.0]
    print [src['mstarK11'] for src in data if src['mstarK11']>10.0]
    print [src['vdisp_hi'] for src in data if src['mstarK11']>10.0]
    print [src['mstarother_low'] for src in data if src['mstarK11']>10.0]
    print [src['mstarother_hi'] for src in data if src['mstarK11']>10.0]


    with open('dwarf_sample_mstar10_bmag-19.5_withcomments.txt',"w") as fp5c:
        fp5c.write(', '.join([key for key in sample5_data[0]])+"\n")
        for src in sample5_data:
            fp5c.write(', '.join([str(src[key]) for key in src]))
            fp5c.write("\n")
    with open('dwarf_sample_mstar10_bmag-19.5.txt',"w") as fp5:
        fp5.write(', '.join([key for key in sample5_data[0] 
                             if (key!='comments') and (key!='morphnotes')])+"\n")
        for src in sample5_data:
            fp5.write(', '.join([str(src[key]) for key in src 
                                 if (key!='comments') and (key!='morphnotes')]))
            fp5.write("\n")


    sample6_data =  [src for src in data if ( (((src['mstarK11']<10.0) and (src['mstarK11']!=-1)) or
                                               ((src['Bmag']>-19.5) and (src['Bmag']!=-1) and
                                                ((src['mstarK11']<10.0) or (src['mstarK11']==-1)))) and
                                             (src['vdisp_hi']<100.0) )]
    print "\nSample 6: mstarK11<10 or Bmag>-19.5, plus exclude *all* with mstarK11>10 or vdisp_hi>=100."
    print "(sample size="+str(len(sample6_data))+"):\n"
    print [src['name'] for src in sample6_data]
    print [src['Bmag'] for src in sample6_data]
    print [src['mstarK11'] for src in sample6_data]
    print [src['vdisp_low'] for src in sample6_data]
    print [src['vdisp_hi'] for src in sample6_data]
    print [src['mstarother_low'] for src in sample6_data]
    print [src['mstarother_hi'] for src in sample6_data]

    sample7_data =  [src for src in data if ( (((src['mstarK11']<9.7) and (src['mstarK11']!=-1)) or
                                               ((src['Bmag']>-19.5) and (src['Bmag']!=-1) and
                                                ((src['mstarK11']<9.7) or (src['mstarK11']==-1)))) )]
    print "\nSample 7: mstarK11<9.7 or Bmag>-19.5, plus exclude *all* with mstarK11>9.7."
    print "(sample size="+str(len(sample7_data))+"):\n"
    print [src['name'] for src in sample7_data]
    print [src['Bmag'] for src in sample7_data]
    print [src['mstarK11'] for src in sample7_data]
    print [src['vdisp_low'] for src in sample7_data]
    print [src['vdisp_hi'] for src in sample7_data]
    print [src['mstarother_low'] for src in sample7_data]
    print [src['mstarother_hi'] for src in sample7_data]

    with open('dwarf_sample_mstar9.7_bmag-19.5_withcomments.txt',"w") as fp7c:
        fp7c.write(', '.join([key for key in sample7_data[0]])+"\n")
        for src in sample7_data:
            fp7c.write(', '.join([str(src[key]) for key in src]))
            fp7c.write("\n")
    with open('dwarf_sample_mstar9.7_bmag-19.5.txt',"w") as fp7:
        fp7.write(', '.join([key for key in sample7_data[0] 
                             if (key!='comments') and (key!='morphnotes')])+"\n")
        for src in sample7_data:
            fp7.write(', '.join([str(src[key]) for key in src 
                                 if (key!='comments') and (key!='morphnotes')]))
            fp7.write("\n")

    #sample8_data =  [src for src in data if ( (((src['mstarK11']<9.7) and (src['mstarK11']!=-1)) or
    #                                           ((src['Bmag']>-19.5) and (src['Bmag']!=-1) and
    #                                            ((src['mstarK11']<9.7) or (src['mstarK11']==-1)))) and
    #                                          (src['vdisp_hi']<100.0) and (src['mstarother_low']<9.7))]
    sample8_data =  [src for src in data if ( (((src['mstarK11']<9.7) and (src['mstarK11']!=-1)) or
                                               ((src['Bmag']>-19.5) and (src['Bmag']!=-1) and
                                                ((src['mstarK11']<9.7) or (src['mstarK11']==-1)))) and
                                              (src['vdisp_hi']<100.0))]
    print "\nSample 8: mstarK11<9.7 or Bmag>-19.5, plus exclude *all* with mstarK11>9.7 or vdisp_hi>=100."
    print "(sample size="+str(len(sample8_data))+"):\n"
    print [src['name'] for src in sample8_data]
    print [src['Bmag'] for src in sample8_data]
    print [src['mstarK11'] for src in sample8_data]
    print [src['vdisp_low'] for src in sample8_data]
    print [src['vdisp_hi'] for src in sample8_data]
    print [src['mstarother_low'] for src in sample8_data]
    print [src['mstarother_hi'] for src in sample8_data]
    print [src['bhmass_low'] for src in sample8_data]
    print [src['bhmass_hi'] for src in sample8_data]
    print [(src['AGNtype'],src['morphnotes']) for src in sample8_data]

    with open('dwarf_sample_mstar9.7_bmag-19.5_vdisp100_withcomments.txt',"w") as fp8c:
        fp8c.write(', '.join([key for key in sample8_data[0]])+"\n")
        for src in sample8_data:
            fp8c.write(', '.join([str(src[key]) for key in src]))
            fp8c.write("\n")
    with open('dwarf_sample_mstar9.7_bmag-19.5_vdisp100.txt',"w") as fp8:
        fp8.write(', '.join([key for key in sample8_data[0] 
                             if (key!='comments') and (key!='morphnotes')])+"\n")
        for src in sample8_data:
            fp8.write(', '.join([str(src[key]) for key in src 
                                 if (key!='comments') and (key!='morphnotes')]))
            fp8.write("\n")
    with open('lowmass_srclist_mstar9.7_bmag-19.5_vdisp100.txt',"w") as fp8:
        for src in sample8_data:
            fp8.write(' '.join([str(src[key]) for key in src 
                                 if (key=='name') or (key=='ra') or (key=='dec')]))
            fp8.write("\n")


    sample9_data =  [src for src in data if ( (((src['mstarK11']<9.5) and (src['mstarK11']!=-1)) or
                                               ((src['Bmag']>-19.0) and (src['Bmag']!=-1) and
                                                ((src['mstarK11']<9.5) or (src['mstarK11']==-1)))) )]
    print "\nSample 9: mstarK11<9.5 or Bmag>-19.0, plus exclude *all* with mstarK11>9.5."
    print "(sample size="+str(len(sample9_data))+"):\n"
    print [src['name'] for src in sample9_data]
    print [src['Bmag'] for src in sample9_data]
    print [src['mstarK11'] for src in sample9_data]
    print [src['vdisp_low'] for src in sample9_data]
    print [src['vdisp_hi'] for src in sample9_data]
    print [src['mstarother_low'] for src in sample9_data]
    print [src['mstarother_hi'] for src in sample9_data]
    print [src['bhmass_low'] for src in sample9_data]
    print [src['bhmass_hi'] for src in sample9_data]
    print [(src['AGNtype'],src['morphnotes']) for src in sample9_data]

    with open('dwarf_sample_mstar9.5_bmag-19.0_withcomments.txt',"w") as fp9c:
        fp9c.write(', '.join([key for key in sample9_data[0]])+"\n")
        for src in sample9_data:
            fp9c.write(', '.join([str(src[key]) for key in src]))
            fp9c.write("\n")
    with open('dwarf_sample_mstar9.5_bmag-19.0.txt',"w") as fp9:
        fp9.write(', '.join([key for key in sample9_data[0] 
                             if (key!='comments') and (key!='morphnotes')])+"\n")
        for src in sample9_data:
            fp9.write(', '.join([str(src[key]) for key in src 
                                 if (key!='comments') and (key!='morphnotes')]))
            fp9.write("\n")
    with open('dwarf_srclist_mstar9.5_bmag-19.0.txt',"w") as fp9:
        for src in sample9_data:
            fp9.write(' '.join([str(src[key]) for key in src 
                                if (key=='name') or (key=='ra') or (key=='dec')]))
            fp9.write("\n")


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



def safeassign(val,mytype=np.float,default=-1,exception=ValueError):

    try:
        return mytype(val)
    except exception:
        return default
