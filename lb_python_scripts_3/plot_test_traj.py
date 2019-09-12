import numpy as np
import matplotlib.pyplot as plt
import sys

jstep = np.array([]).astype('int')
time = np.array([])
r = np.array([])
rphi = np.array([])
z = np.array([])
vr = np.array([])
vphi = np.array([])
vz = np.array([])

fnstr = '_vm08compare_h+b_c5.5_2mbh_mhcorr_trunc_1.0vesc'
fname = 'alltraj'+fnstr+'.dat'
#pname = 'traj'+fnstr+'.eps'
pname = 'traj'+fnstr+'.png'

try:
    fp = open(fname,"r")
except IOError:
    print "Error opening file "+fname    
else:
    linct = 0
    for line in fp:
        linarr = line.strip().split()
        if (linct % 1000) == 0:
            print str(linct)+" lines read."
        #if (linct > 50000):
        #    break
        if len(linarr) != 8:
            print "Invalid line read from file."
            sys.exit()
        linct += 1
        jstep = np.append(jstep,int(linarr[0]))        
        time = np.append(time,np.float64(linarr[1]))
        r = np.append(r,np.float64(linarr[2]))
        rphi = np.append(rphi,np.float64(linarr[3]))
        z = np.append(z,np.float64(linarr[4]))
        vr = np.append(vr,np.float64(linarr[5]))
        vphi = np.append(vphi,np.float64(linarr[6]))
        vz = np.append(vz,np.float64(linarr[7]))

print str(linct)+" lines read from file "+fname
fp.close()

phi = rphi/r
R = np.sqrt(r*r + z*z)
print "min/max r: ",r.min(),r.max()
print "min/max phi: ",phi.min(),phi.max()
print "min/max z: ",z.min(),z.max()
print "min/max R: ",R.min(),R.max()
print "min/max vr: ",vr.min(),vr.max()
print "min/max vphi: ",vphi.min(),vphi.max()
print "min/max vz: ",vz.min(),vz.max()
plt.clf()
fig = plt.figure(figsize=(8,3.5))

ax1 = fig.add_subplot(1,2,1)
#ax1.xlabel('')
plt.xlim(-r.max(),r.max())
plt.ylim(-r.max(),r.max())
#if linct > 75000:
#    ax1.plot(r[-50000:]*np.cos(phi[-50000:]), r[-50000:]*np.sin(phi[-50000:]),'k-')
#else:
ax1.plot(r*np.cos(phi), r*np.sin(phi),'k-')
#ax1.plot(r[0:50000]*np.cos(phi[0:50000]), r[0:50000]*np.sin(phi[0:50000]),'k-')
#ax1.plot(r[50001:100000]*np.cos(phi[50001:100000]), r[50001:100000]*np.sin(phi[50001:100000]),'b-')
#ax1.plot(r[100001:]*np.cos(phi[100001:]), r[100001:]*np.sin(phi[100001:]),'r-')
#ax1.title('')

ax2 = fig.add_subplot(1,2,2)
#if linct > 75000:
#    ax2.plot(time[-50000:], r[-50000:],'k-')
#else:
#     ax2.plot(time, r,'k-')
ax2.set_xscale('log')
ax2.set_yscale('log')
#plt.xlim(1.0e6,5.0e9)
plt.xlim(1.0e6,1.4e11)
plt.ylim(100,2.0e7)
plt.xlabel('time [yr]')
plt.ylabel('R [pc]')
ax2.plot(time, R,'k-')

#fig.savefig('./'+pname)
fig.savefig('./'+pname,dpi=80)
plt.close()
