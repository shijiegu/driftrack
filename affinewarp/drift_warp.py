import numpy as np
import os
import sys
from .piecewisewarp import PiecewiseWarping
from .rank1_piecewisewarp import Rank1_PiecewiseWarping

def drift_warp(datapath):
    data = np.load(datapath).T
    data=np.log(1+data)
    data=np.expand_dims(data, axis=2)

    model = PiecewiseWarping(n_knots=4, warp_reg_scale=1e-9, smoothness_reg_scale=0)
    model.fit(data, iterations=20, warp_iterations=200)
    shifts=model.get_shift(data)*-1

    path, file = os.path.split(datapath)
    shift_datapath=os.path.join(path,'shift.npy')
    print('shift done...about to save')

    np.save(shift_datapath,shifts.squeeze())
    return shift_datapath

if __name__ == '__main__':
    datapath = sys.argv[1]
    sys.stdout.write(drift_warp(datapath))