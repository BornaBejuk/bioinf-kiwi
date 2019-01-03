import numpy as np
from utils import *

if __name__ == '__main__':
    path_c_r = '../data/EColi-synthetic/overlaps-c-r.paf'
    path_r_r = '../data/EColi-synthetic/overlaps-r-r.paf'

    overlaps_c_r = np.loadtxt(path_c_r, delimiter='\t', dtype='str')
    overlaps_r_r = np.loadtxt(path_r_r, delimiter='\t', dtype='str')

    print (get_avg_SI(overlaps_r_r))
    bla = overlaps_c_r[:,9:10].astype(float) / overlaps_c_r[:,10:11].astype(float) > 0.5
    print ((overlaps_c_r[ bla, :]).shape)
    print ((SI[SI[:,0:1] > 0.5]).shape)
    # print (overlaps_c_r[:][9]/overlaps_c_r[:][10])

    # print (overlaps_c_r[0][9])
    # print (overlaps_c_r[0][10])
