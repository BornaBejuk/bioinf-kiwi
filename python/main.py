import numpy as np
from utils import *

if __name__ == '__main__':
    path_c_r = '../data/EColi-synthetic/overlaps-c-r.paf'
    path_r_r = '../data/EColi-synthetic/overlaps-r-r.paf'

    overlaps_c_r = load_data(path_c_r)

    # m = np.argmin(overlaps_r_r[:,10:11].astype(int))
    # print (m)
    # print (overlaps_r_r[m])
    # o = np.rec.fromarrays((overlaps_c_r), names=('a','b','c','d','e','f','g','h','u','r','w','q','z','i','o','p'), shape=overlaps_c_r.shape)
    # print (overlaps_c_r[0])
    # print (o[0])
    # print (get_OS(overlaps_c_r))
    # print (overlaps_c_r[0])
    # print (type(overlaps_c_r[0][0]))


    print (overlaps_c_r['query_end'])
