import numpy as np
import matplotlib.pyplot as plt
import sys

lgr1 = np.array([])
mass1 = np.array([])
rho1 = np.array([])
sigma1 = np.array([])

lgr2 = np.array([])
mass2 = np.array([])
rho2 = np.array([])
sigma2 = np.array([])

lgr3 = np.array([])
mass3 = np.array([])
rho3 = np.array([])
sigma3 = np.array([])


fnstr1 = '_halo_only_hqmod_rsoft1_nosigflat'
fnstr2 = '_halo_only_hqmod_rsoft1'
fnstr3 = '_halo_only_hqmod_rsoft2'
#fnstr3 = '_halo_only_hqmod_rsoft5'
fname1 = 'prof_info'+fnstr1+'.dat'
fname2 = 'prof_info'+fnstr2+'.dat'
fname3 = 'prof_info'+fnstr3+'.dat'
pname = 'prof_halo_only_hqmod_rsoft_with_nosigflat.png'

try:
    fp = open(fname1,"r")
except IOError:
    print "Error opening file "+fname1
else:
    linct = 0
    for line in fp:
        linarr = line.strip().split()
        if len(linarr) != 4:
            print "Invalid line read from file."
            sys.exit()
        linct += 1
        lgr1 = np.append(lgr1,np.float64(linarr[0]))
        mass1 = np.append(mass1,np.float64(linarr[1]))
        rho1 = np.append(rho1,np.float64(linarr[2]))
        sigma1 = np.append(sigma1,np.float64(linarr[3]))
 
print str(linct)+" lines read from file "+fname1
fp.close()

try:
    fp = open(fname2,"r")
except IOError:
    print "Error opening file "+fname2
else:
    linct = 0
    for line in fp:
        linarr = line.strip().split()
        if len(linarr) != 4:
            print "Invalid line read from file."
            sys.exit()
        linct += 1
        lgr2 = np.append(lgr2,np.float64(linarr[0]))
        mass2 = np.append(mass2,np.float64(linarr[1]))
        rho2 = np.append(rho2,np.float64(linarr[2]))
        sigma2 = np.append(sigma2,np.float64(linarr[3]))
 
print str(linct)+" lines read from file "+fname2
fp.close()

try:
    fp = open(fname3,"r")
except IOError:
    print "Error opening file "+fname3
else:
    linct = 0
    for line in fp:
        linarr = line.strip().split()
        if len(linarr) != 4:
            print "Invalid line read from file."
            sys.exit()
        linct += 1
        lgr3 = np.append(lgr3,np.float64(linarr[0]))
        mass3 = np.append(mass3,np.float64(linarr[1]))
        rho3 = np.append(rho3,np.float64(linarr[2]))
        sigma3 = np.append(sigma3,np.float64(linarr[3]))
 
print str(linct)+" lines read from file "+fname3
fp.close()

print "min/max lgr1: ",lgr1.min(),lgr1.max()
print "min/max enclosed mass1: ",mass1.min(),mass1.max()
print "min/max rho1: ",rho1.min(),rho1.max()
print "min/max sigma1: ",sigma1.min(),sigma1.max()

plt.clf()
fig = plt.figure(figsize=(4,8))

ax1 = fig.add_subplot(3,1,1)
#ax1.xlabel('')
ax1.set_yscale('log')
ax1.plot(lgr1, mass1,'k-')
ax1.plot(lgr2, mass2,'b--')
ax1.plot(lgr3, mass3,'r:')

ax2 = fig.add_subplot(3,1,2)
ax2.set_yscale('log')
ax2.plot(lgr1, rho1,'k-')
ax2.plot(lgr2, rho2,'b--')
ax2.plot(lgr3, rho3,'r:')
#ax1.xlim(0,time.max())

ax3 = fig.add_subplot(3,1,3)
ax3.set_yscale('log')
ax3.plot(lgr1, sigma1,'k-')
ax3.plot(lgr2, sigma2,'b--')
ax3.plot(lgr3, sigma3,'r:')

#fig.savefig('./'+pname)
fig.savefig('./'+pname,dpi=80)
plt.close()
