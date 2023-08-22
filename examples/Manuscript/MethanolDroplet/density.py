import numpy as np

import matplotlib.pyplot as plt
from mscg import TrajReader



def plot(file, legend):
    reader = TrajReader(file)
    nread = 0
    data = []
    
    while reader.next_frame():

        com = reader.traj.x.sum(axis=0) / reader.traj.x.shape[0]
        #com = np.array([50., 50., 50.])
        com = np.tile(com, (reader.traj.x.shape[0], 1))
        dr = np.sqrt(np.square((reader.traj.x - com)).sum(axis=1))
        data.append(dr)
        nread += 1
        if nread % 100 == 0: print(nread)
  
    data = np.hstack(data)
    hist = np.histogram(data, bins=np.arange(15, 60.0, 0.5))

    x = (hist[1][:-1] + hist[1][1:]) * 0.5
    y = hist[0] / nread / np.square(x) * np.pi * 3.5
    plt.plot(x, y, label=legend)
    

plot('cg.lammpstrj', 'AA')
plot('mscg_sim.dcd', 'MSCG')
plot('ucg_sim.dcd', 'UCG')
plt.legend(loc='upper right')
plt.show()
