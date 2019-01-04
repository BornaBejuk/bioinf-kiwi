import numpy as np
from utils import *

if __name__ == '__main__':
    path_c_r = '../data/EColi-synthetic/overlaps-c-r.paf'
    path_r_r = '../data/EColi-synthetic/overlaps-r-r.paf'

    overlaps_c_r = np.loadtxt(path_c_r, delimiter='\t', dtype='str')
    # overlaps_r_r = np.loadtxt(path_r_r, delimiter='\t', dtype='str')

    # print (get_avg_SI(overlaps_r_r))
