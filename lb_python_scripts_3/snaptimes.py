import PyCos as pc
import numpy as np
import matplotlib.pyplot as plt

fp=open('/n/hernquistfs2/lblecha/illustris_data/Illustris-1/postprocess/output_times.txt','r')
ascale = np.loadtxt(fp,unpack=True,dtype=np.float)
fp.close()

z = 1.0/ascale - 1

cosmo = pc.Cosmology(0.27, 0.73, 0.0, -1.0, 0.704)

age = np.array([cosmo.age(elem) for elem in z])

dt = np.zeros((age.size))
dt[1:] = age[1:] - age[:-1]

for i in range(136):
    print "%d %g %g %g %g\n"%(i,ascale[i],z[i],age[i]/1.0e9,dt[i]/1.0e9)

a = 1/128.0
afof = np.array([a])
while a<=1.0:
    a = a*1.03
    afof = np.append(afof,a)
zfof = 1.0/afof - 1
agefof = np.array([cosmo.age(elem) for elem in zfof])
dtfof = np.zeros(agefof.size)
dtfof[1:] = agefof[1:] - agefof[:-1]
print dtfof/1.0e9

plt.clf()
plt.xlim(0,1)
ylim=(0,0.45)
plt.ylim(ylim)
#plt.plot(range(136),dt/1.0e9)
plt.plot(afof,dtfof/1.0e9,'bo-')
for tup in zip(afof,afof):
    plt.plot(tup,ylim,'b',linewidth=0.5)
plt.plot(ascale,dt/1.0e9,'ko-')
