import numpy as np
import matplotlib.pyplot as plt
import matplotlib

def testplot(kicksize,tmax):

    if kicksize == 'mid' and tmax == 1.0e+9:
        
        n2=np.array([])
        t2=np.array([])
        r2=np.array([])
        v2=np.array([])
        with open('testint_midvk_eps1e-2.dat','r') as f2:
            for line in f2:
                line = line.strip() ## get rid of newline
                columns = line.split()
                n2 = np.append(n2,int(columns[0]))
                t2 = np.append(t2,float(columns[2]))
                r2 = np.append(r2,np.abs(float(columns[3])))
                v2 = np.append(v2,float(columns[4]))

        n3=np.array([])
        t3=np.array([])
        r3=np.array([])
        v3=np.array([])
        with open('testint_midvk_eps1e-3.dat','r') as f3:
            for line in f3:
                line = line.strip() ## get rid of newline
                columns = line.split()
                n3 = np.append(n3,int(columns[0]))
                t3 = np.append(t3,float(columns[2]))
                r3 = np.append(r3,np.abs(float(columns[3])))
                v3 = np.append(v3,float(columns[4]))

        n4=np.array([])
        t4=np.array([])
        r4=np.array([])
        v4=np.array([])
        with open('testint_midvk_eps1e-4.dat','r') as f4:
            for line in f4:
                line = line.strip() ## get rid of newline
                columns = line.split()
                n4 = np.append(n4,int(columns[0]))
                t4 = np.append(t4,float(columns[1]))
                r4 = np.append(r4,np.abs(float(columns[2])))
                v4 = np.append(v4,float(columns[3]))

        n5=np.array([])
        t5=np.array([])
        r5=np.array([])
        v5=np.array([])
        with open('testint_midvk_eps1e-5.dat','r') as f5:
            for line in f5:
                line = line.strip() ## get rid of newline
                columns = line.split()
                n5 = np.append(n5,int(columns[0]))
                t5 = np.append(t5,float(columns[2]))
                r5 = np.append(r5,np.abs(float(columns[3])))        
                v5 = np.append(v5,float(columns[4]))

        n6=np.array([])
        t6=np.array([])
        r6=np.array([])
        v6=np.array([])
        with open('testint_midvk_eps1e-6.dat','r') as f6:
            for line in f6:
                line = line.strip() ## get rid of newline
                columns = line.split()
                n6 = np.append(n6,int(columns[0]))
                t6 = np.append(t6,float(columns[1]))
                r6 = np.append(r6,np.abs(float(columns[2])))
                v6 = np.append(v6,float(columns[3]))

        nrel6=np.array([])
        trel6=np.array([])
        rrel6=np.array([])
        vrel6=np.array([])
        with open('testint_midvk_epsrel1e-6.dat','r') as frel6:
            for line in frel6:
                line = line.strip() ## get rid of newline
                columns = line.split()
                nrel6 = np.append(nrel6,int(columns[0]))
                trel6 = np.append(trel6,float(columns[4]))
                rrel6 = np.append(rrel6,np.abs(float(columns[5])))
                vrel6 = np.append(vrel6,float(columns[6]))

        nrel7=np.array([])
        trel7=np.array([])
        rrel7=np.array([])
        vrel7=np.array([])
        with open('testint_midvk_epsrel1e-7.dat','r') as frel7:
            for line in frel7:
                line = line.strip() ## get rid of newline
                columns = line.split()
                nrel7 = np.append(nrel7,int(columns[0]))
                trel7 = np.append(trel7,float(columns[4]))
                rrel7 = np.append(rrel7,np.abs(float(columns[5])))
                vrel7 = np.append(vrel7,float(columns[6]))

        nrel8=np.array([])
        trel8=np.array([])
        rrel8=np.array([])
        vrel8=np.array([])
        with open('testint_midvk_epsrel1e-8.dat','r') as frel8:
            for line in frel8:
                line = line.strip() ## get rid of newline
                columns = line.split()
                nrel8 = np.append(nrel8,int(columns[0]))
                trel8 = np.append(trel8,float(columns[4]))
                rrel8 = np.append(rrel8,np.abs(float(columns[5])))
                vrel8 = np.append(vrel8,float(columns[6]))

        nrel9=np.array([])
        trel9=np.array([])
        rrel9=np.array([])
        vrel9=np.array([])
        with open('testint_midvk_epsrel1e-9.dat','r') as frel9:
            for line in frel9:
                line = line.strip() ## get rid of newline
                columns = line.split()
                nrel9 = np.append(nrel9,int(columns[0]))
                trel9 = np.append(trel9,float(columns[4]))
                rrel9 = np.append(rrel9,np.abs(float(columns[5])))
                vrel9 = np.append(vrel9,float(columns[6]))

        nrel10=np.array([])
        trel10=np.array([])
        rrel10=np.array([])
        vrel10=np.array([])
        with open('testint_midvk_epsrel1e-10.dat','r') as frel10:
            for line in frel10:
                line = line.strip() ## get rid of newline
                columns = line.split()
                nrel10 = np.append(nrel10,int(columns[0]))
                trel10 = np.append(trel10,float(columns[4]))
                rrel10 = np.append(rrel10,np.abs(float(columns[5])))
                vrel10 = np.append(vrel10,float(columns[6]))


    elif kicksize == 'mid' and tmax == 1.0e+10:
        
        ntH3=np.array([])
        ttH3=np.array([])
        rtH3=np.array([])
        vtH3=np.array([])
        with open('testint_midvk_eps1e-3_tH.dat','r') as ftH3:
            for line in ftH3:
                line = line.strip() ## get rid of newline
                columns = line.split()
                ntH3 = np.append(ntH3,int(columns[0]))
                ttH3 = np.append(ttH3,float(columns[4]))
                rtH3 = np.append(rtH3,np.abs(float(columns[5])))
                vtH3 = np.append(vtH3,float(columns[6]))

        ntH4=np.array([])
        ttH4=np.array([])
        rtH4=np.array([])
        vtH4=np.array([])
        with open('testint_midvk_eps1e-4_tH.dat','r') as ftH4:
            for line in ftH4:
                line = line.strip() ## get rid of newline
                columns = line.split()
                ntH4 = np.append(ntH4,int(columns[0]))
                ttH4 = np.append(ttH4,float(columns[4]))
                rtH4 = np.append(rtH4,np.abs(float(columns[5])))
                vtH4 = np.append(vtH4,float(columns[6]))

        ntHrel8=np.array([])
        ttHrel8=np.array([])
        rtHrel8=np.array([])
        vtHrel8=np.array([])
        with open('testint_midvk_epsrel1e-8_tH.dat','r') as ftHrel8:
            for line in ftHrel8:
                line = line.strip() ## get rid of newline
                columns = line.split()
                ntHrel8 = np.append(ntHrel8,int(columns[0]))
                ttHrel8 = np.append(ttHrel8,float(columns[4]))
                rtHrel8 = np.append(rtHrel8,np.abs(float(columns[5])))
                vtHrel8 = np.append(vtHrel8,float(columns[6]))

        ntHrel10=np.array([])
        ttHrel10=np.array([])
        rtHrel10=np.array([])
        vtHrel10=np.array([])
        with open('testint_midvk_epsrel1e-10_tH.dat','r') as ftHrel10:
            for line in ftHrel10:
                line = line.strip() ## get rid of newline
                columns = line.split()
                ntHrel10 = np.append(ntHrel10,int(columns[0]))
                ttHrel10 = np.append(ttHrel10,float(columns[4]))
                rtHrel10 = np.append(rtHrel10,np.abs(float(columns[5])))
                vtHrel10 = np.append(vtHrel10,float(columns[6]))

        ntHrel12=np.array([])
        ttHrel12=np.array([])
        rtHrel12=np.array([])
        vtHrel12=np.array([])
        with open('testint_midvk_epsrel1e-12_tH.dat','r') as ftHrel12:
            for line in ftHrel12:
                line = line.strip() ## get rid of newline
                columns = line.split()
                ntHrel12 = np.append(ntHrel12,int(columns[0]))
                ttHrel12 = np.append(ttHrel12,float(columns[4]))
                rtHrel12 = np.append(rtHrel12,np.abs(float(columns[5])))
                vtHrel12 = np.append(vtHrel12,float(columns[6]))

        ntHrel12a0=np.array([])
        ttHrel12a0=np.array([])
        rtHrel12a0=np.array([])
        vtHrel12a0=np.array([])
        with open('testint_midvk_epsrel1e-12_epsabs0_tH.dat','r') as ftHrel12a0:
            for line in ftHrel12a0:
                line = line.strip() ## get rid of newline
                columns = line.split()
                ntHrel12a0 = np.append(ntHrel12a0,int(columns[0]))
                ttHrel12a0 = np.append(ttHrel12a0,float(columns[4]))
                rtHrel12a0 = np.append(rtHrel12a0,np.abs(float(columns[5])))
                vtHrel12a0 = np.append(vtHrel12a0,float(columns[6]))


    elif kicksize == 'tiny':
        
        ntiny4=np.array([])
        ttiny4=np.array([])
        rtiny4=np.array([])
        vtiny4=np.array([])
        with open('testint_tinyvk_eps1e-4.dat','r') as ftiny4:
            for line in ftiny4:
                line = line.strip() ## get rid of newline
                columns = line.split()
                ntiny4 = np.append(ntiny4,int(columns[0]))
                ttiny4 = np.append(ttiny4,float(columns[4]))
                rtiny4 = np.append(rtiny4,np.abs(float(columns[5])))
                vtiny4 = np.append(vtiny4,float(columns[6]))

        ntinyrel12=np.array([])
        ttinyrel12=np.array([])
        rtinyrel12=np.array([])
        vtinyrel12=np.array([])
        with open('testint_tinyvk_epsrel1e-12.dat','r') as ftinyrel12:
            for line in ftinyrel12:
                line = line.strip() ## get rid of newline
                columns = line.split()
                ntinyrel12 = np.append(ntinyrel12,int(columns[0]))
                ttinyrel12 = np.append(ttinyrel12,float(columns[4]))
                rtinyrel12 = np.append(rtinyrel12,np.abs(float(columns[5])))
                vtinyrel12 = np.append(vtinyrel12,float(columns[6]))

    elif kicksize == 'df':

        ndf=np.array([])
        tdf=np.array([])
        rdf=np.array([])
        vdf=np.array([])
        with open('testint_df.dat','r') as fdf:
        #with open('testint_df_Hq.dat','r') as fdf:
            for line in fdf:
                line = line.strip() ## get rid of newline
                columns = line.split()
                ndf = np.append(ndf,int(columns[0]))
                tdf = np.append(tdf,float(columns[1]))
                rdf = np.append(rdf,np.abs(float(columns[2])))
                vdf = np.append(vdf,float(columns[3]))

        tf=np.array([])
        rf=np.array([])
        vf=np.array([])
        agrav=np.array([])
        adf=np.array([])

        # with open('forces.dat','r') as fforce:
        ##with open('forces_Hq.dat','r') as fforce:
        #count=0
        #for line in fforce:
        #    count += 1
        #    if (count % 10) == 0:
        #        line = line.strip() ## get rid of newline
        #        columns = line.split()
        #        tf = np.append(tf,float(columns[0]))
        #rf = np.append(rf,float(columns[1]))
        #vf = np.append(vf,np.abs(float(columns[2])))
        #        agrav = np.append(agrav,float(columns[3]))
        #        adf = np.append(adf,float(columns[4]))
        #        if (count % 1000) == 0:
        #            print "read ",count," lines"
                
    else:

        print "Invalid values of kicksize or tmax."
        return 1
    

    plt.clf()

######### Radius plot #########
    
    ax1 = plt.subplot(211)
    #ax1 = plt.subplot(411)

    if kicksize == 'mid' and tmax == 1.0e+9:

##     plt.plot(t2,r2,'y')
##     plt.plot(t4,r4,'r')
##     plt.plot(t5,r5,'k')
##     plt.plot(t6,r6,'g')
##     plt.plot(trel6,rrel6,'b')
##     plt.yscale('log')

        plt.plot(t5,np.abs((r5-r4)/r5),'r')
        plt.plot(t5,np.abs((r5-r3)/r5),'c')
        plt.plot(t5,np.abs((r5-r2)/r5),'y')
        plt.plot(t5,np.abs((r5-rrel6)/r5),'m')
        plt.plot(t5,np.abs((r5-rrel8)/r5),'g')
        plt.plot(t5,np.abs((r5-rrel10)/r5),'b')
        plt.yscale('log')

    if kicksize == 'mid' and tmax == 1.0e+10:

        plt.plot(ttH4,np.abs((rtH4-rtH3)/rtH4),'r')
        plt.plot(ttHrel8,np.abs((rtH4-rtHrel8)/rtH4),'m:')
        plt.plot(ttHrel10,np.abs((rtH4-rtHrel10)/rtH4),'g')
        plt.plot(ttHrel12,np.abs((rtH4-rtHrel12)/rtH4),'b')
        plt.plot(ttHrel12a0,np.abs((rtH4-rtHrel12a0)/rtH4),'y')
        plt.yscale('log')

    if kicksize == 'tiny':

        #plt.plot(ttiny4,rtiny4,'k')
        #plt.plot(ttinyrel12,rtinyrel12,'r')

        plt.plot(ttinyrel12,np.abs((rtinyrel12-rtiny4)/rtinyrel12),'r')
        plt.yscale('log')

        plt.xlim(0.0,1.0e+7)

    if kicksize == 'df':

	rinfl = 2.51776
        plt.plot(tdf,rdf,'k')
        #plt.xlim(1.0e6,1.0e+9)
	plt.plot([tdf[0],tdf[len(tdf)-1]],[rinfl,rinfl],'b')
        plt.yscale('log')
        #plt.xscale('log')
        #plt.xlim(0.0,4.0e7)

######### Velocity plot #########

    ax2 = plt.subplot(212)
    #ax2 = plt.subplot(412)

    if kicksize == 'mid' and tmax == 1.0e+9:

##     plt.plot(t2,v2,'y')
##     plt.plot(t4,v4,'r')
##     plt.plot(t5,v5,'k')
##     plt.plot(t6,v6,'g')
##     plt.plot(trel6,vrel6,'b')
##     plt.yscale('linear')

        plt.plot(t5,np.abs((v5-v4)/v5),'r')
        plt.plot(t5,np.abs((v5-v3)/v5),'c')
        plt.plot(t5,np.abs((v5-v2)/v5),'y')
        plt.plot(t5,np.abs((v5-vrel6)/v5),'m')
        plt.plot(t5,np.abs((v5-vrel8)/v5),'g')
        plt.plot(t5,np.abs((v5-vrel10)/v5),'b')
        plt.yscale('log')

    if kicksize == 'mid' and tmax == 1.0e+10:

        plt.plot(ttH4,np.abs((vtH4-vtH3)/vtH4),'r')
        plt.plot(ttHrel8,np.abs((vtH4-vtHrel8)/vtH4),'m:')
        plt.plot(ttHrel10,np.abs((vtH4-vtHrel10)/vtH4),'g')
        plt.plot(ttHrel12,np.abs((vtH4-vtHrel12)/vtH4),'b')
        plt.plot(ttHrel12a0,np.abs((vtH4-vtHrel12a0)/vtH4),'y')
        plt.yscale('log')

    if kicksize == 'tiny':

        #plt.plot(ttiny4,vtiny4,'k')
        #plt.plot(ttinyrel12,vtinyrel12,'r')

        plt.plot(ttinyrel12,np.abs((vtinyrel12-vtiny4)/vtinyrel12),'r')
        plt.yscale('log')

        plt.xlim(0.0,1.0e+7)

    if kicksize == 'df':

        plt.plot(tdf,vdf,'k')
        #plt.xlim(1.0e6,1.0e+9)
        plt.yscale('linear')
        #plt.xscale('log')
        #plt.xlim(0.0,4.0e7)


    #ax3 = plt.subplot(413)

    #if kicksize == 'df':

        #plt.plot(tf, agrav,'b')
        ##plt.plot(tf, adf,'g')
        ##plt.plot(tf, np.abs(agrav))
        ##plt.plot(tf, np.abs(adf))
        ##plt.ylim(1.0e-13,1.0e-6)
        ##plt.yscale('log')
        ##plt.xlim(1.9e8,2e8)

    #ax4 = plt.subplot(414)

    #if kicksize == 'df':

        #plt.plot(tf, adf,'g')
        ##plt.plot(tf, np.abs(agrav))
        ##plt.plot(tf, np.abs(adf))
        ##plt.xlim(3.8e8,4.4e8)
        ##plt.ylim(1.0e-13,1.0e-6)
        ##plt.yscale('log')
        ##plt.xlim(1.9e8,2e8)

    plt.show()

    return 0
