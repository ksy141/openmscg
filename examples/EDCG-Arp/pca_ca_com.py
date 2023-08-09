import MDAnalysis as mda
from MDAnalysis.analysis import pca,align
import pandas as pd
#from MDAnalysis.analysis.rms import RMSF
import matplotlib.pyplot as plt
plt.ioff()
fig = plt.figure()

u_traj = mda.Universe('data.data','GC.dcd',format="LAMMPS")
u_ref = mda.Universe('data.data','GC.dcd',format="LAMMPS")

u_traj.trajectory[-1]
u_ref.trajectory[0]

alignment = align.AlignTraj(u_traj, u_ref, select='all',filename='aligned_to_first_frame.dcd')
alignment.run()
aligned_u = mda.Universe('data.data','aligned_to_first_frame.dcd',format="LAMMPS")
pc = pca.PCA(aligned_u).run()

import numpy as np
np.save('pc.npy', pc.p_components)
np.save('ev.npy', pc.variance)
