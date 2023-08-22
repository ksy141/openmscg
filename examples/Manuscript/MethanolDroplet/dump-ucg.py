from mscg.cli import cgdump

import matplotlib.pyplot as plt
import numpy as np

cgdump.main(
    file = 'cg_models.p',
    dump = ['Pair_MeOH-MeOH,2.3,12.0,0.05']
)

cgdump.main(
    file = 'ucg_models.p',
    dump = ['Pair_Near-Near,2.3,12.0,0.05',
            'Pair_Near-Far,2.3,12.0,0.05',
            'Pair_Far-Far,2.3,12.0,0.05']
)

for tbl in ['Pair_Near-Near', 'Pair_Near-Far', 'Pair_Far-Far']:
    data = np.loadtxt(tbl + '.table', skiprows=5)
    plt.plot(data[:,1], data[:,2], label='UCG ' + tbl)

tbl = np.loadtxt('Pair_MeOH-MeOH.table', skiprows=5)
plt.plot(tbl[:,1], tbl[:,2], '--', label='MSCG Pair')

plt.legend(loc='upper right')
plt.xlim(3, 10)
plt.ylim(-3, 3)

plt.show()
