import numpy as np
from utils import *
from test_dfs import *
import pickle

if __name__ == '__main__':
    path_c_r = 'overlaps-c-r.paf'
    path_r_r = 'overlaps-r-r.paf'

    # overlaps_c_r = load_data(path_c_r)
    # overlaps_r_r = load_data(path_r_r)

    # m = np.argmin(overlaps_r_r[:,10:11].astype(int))
    # print (m)
    # print (overlaps_r_r[m])
    # o = np.rec.fromarrays((overlaps_c_r), names=('a','b','c','d','e','f','g','h','u','r','w','q','z','i','o','p'), shape=overlaps_c_r.shape)
    # print (overlaps_c_r[0])
    # print (o[0])
    # print (get_OS(overlaps_c_r))
    # print (overlaps_c_r[0])
    # print (type(overlaps_c_r[0][0]))


    # print (overlaps_c_r['query_end'])
    # print (get_ES(overlaps_r_r)[0].shape)


    # overlaps_c_r = append_scores(overlaps_c_r)
    # overlaps_r_r = append_scores(overlaps_r_r)
    # with open('overlaps_c_r_loaded', 'wb') as grouped_data_file:
    #     pickle.dump(overlaps_c_r, grouped_data_file)
    # with open('overlaps_c_r_loaded', 'rb') as grouped_data_file:
    #     overlaps_c_r = pickle.load(grouped_data_file)
    # overlaps_c_r = append_scores2(overlaps_c_r)
    # overlaps_r_r = append_scores2(overlaps_r_r)

    # with open('grouped_r_r_2', 'wb') as grouped_data_file:
    #     pickle.dump(overlaps_r_r, grouped_data_file)
    with open('grouped_r_r_2', 'rb') as grouped_data_file:
        overlaps_r_r_2 = pickle.load(grouped_data_file)
    print(overlaps_r_r_2[0])





    # print(overlaps_c_r[0])
    # print(np.array(len(overlaps_c_r[0])))
    # print(np.array(overlaps_c_r[0].shape))
    # print(get_EL(overlaps_c_r[0]))

    grouped_data = get_grouped_data(overlaps_r_r_2)
    # keys = grouped_data1.unique
    # with open('keys_r_r', 'wb') as grouped_data_file:
    #     pickle.dump(keys, grouped_data_file)
    # with open('grouped_data_c_r', 'wb') as grouped_data_file:
        # pickle.dump(grouped_data1, grouped_data_file)

    grouped_data = grouped_data.split(overlaps_r_r_2[['query_name','SI','OS','ES','extension_side', 'EL', 'OH']])
    with open('grouped_data_r_r_2', 'wb') as grouped_data_file:
        pickle.dump(grouped_data, grouped_data_file)
    # print(grouped_data1.unique)
    # print(grouped_data1[keys[0]])
    # print (keys)
    # print (grouped_data1[0])
