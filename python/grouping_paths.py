import pickle
import numpy as np

def divide_paths_into_anchor_groups(list_of_list_of_paths):
    map_of_paths = dict()
    for l in list_of_list_of_paths:
        for p in l:
            path = p[0]
            start_contig = path[0]
            end_contig = path[-1]
            contigs = sorted([start_contig, end_contig])
            key = contigs[0] + contigs[1]
            if key not in map_of_paths:
                map_of_paths[key] = []
                map_of_paths[key].append(path)
            else:
                map_of_paths[key].append(path)
    return map_of_paths


def get_all_lengths(map_of_paths, lengths_of_contigs, grouped_c_r, grouped_r_r, keys_c_r, keys_r_r):
    map_of_lengths = dict()
    for m in map_of_paths:
        for path in map_of_paths[m]:
            length_of_path = 0
            number_of_nodes = len(path)
            for i in range(number_of_nodes):
                if i < number_of_nodes - 1:
                    node_target = path[i-1]
                    node_cont = path[i]
                else:
                    node_target = path[i]
                    node_cont = path[i-1]
                if i == 0:
                    length_of_path += lengths_of_contigs[node_target]
                elif i == 1:
                    key = keys_c_r.tolist().index(node_target)
                    group = grouped_c_r[key]
                    read_cont = group[np.where(group['query_name'] == node_cont)]
                    length_of_path -= read_cont['OH'][0][0]
                    length_of_path += read_cont['EL'][0][1]
                elif i == number_of_nodes - 1:
                    key = keys_c_r.tolist().index(node_target)
                    group = grouped_c_r[key]
                    read_cont = group[np.where(group['query_name'] == node_cont)]
                    length_of_path -= read_cont['OH'][0][1]
                    length_of_path += read_cont['EL'][0][0]
                else:
                    key = keys_r_r.tolist().index(node_target)
                    group = grouped_r_r[key]
                    read_cont = group[np.where(group['query_name'] == node_cont)]
                    length_of_path -= read_cont['OH'][0][0]
                    length_of_path += read_cont['EL'][0][1]
            if m not in map_of_lengths:
                map_of_lengths[m] = []
                map_of_lengths[m].append((path, length_of_path))
            else:
                map_of_lengths[m].append((path, length_of_path))
    return map_of_lengths




if __name__ == '__main__':
    # all_paths = []

    # with open('mc_ctg1_left', 'rb') as mc_ctg1_right_file:
    #     mc_ctg1_left = pickle.load(mc_ctg1_right_file)
    # with open('mc_ctg1_right', 'rb') as mc_ctg1_right_file:
    #     mc_ctg1_right = pickle.load(mc_ctg1_right_file)
    # with open('mc_ctg1_left_2', 'rb') as mc_ctg1_right_file:
    #     mc_ctg1_left_2 = pickle.load(mc_ctg1_right_file)
    # with open('mc_ctg1_left_3', 'rb') as mc_ctg1_right_file:
    #     mc_ctg1_left_3 = pickle.load(mc_ctg1_right_file)
    # with open('mc_ctg2_left', 'rb') as mc_ctg1_right_file:
    #     mc_ctg2_left = pickle.load(mc_ctg1_right_file)
    # with open('mc_ctg2_right', 'rb') as mc_ctg1_right_file:
    #     mc_ctg2_right = pickle.load(mc_ctg1_right_file)
    # with open('mc_ctg2_right_2', 'rb') as mc_ctg1_right_file:
    #     mc_ctg2_right_2 = pickle.load(mc_ctg1_right_file)
    # with open('mc_ctg2_right_3', 'rb') as mc_ctg1_right_file:
    #     mc_ctg2_right_3 = pickle.load(mc_ctg1_right_file)
    # with open('mc_ctg3_left', 'rb') as mc_ctg1_right_file:
    #     mc_ctg3_left = pickle.load(mc_ctg1_right_file)
    # all_paths.append(mc_ctg1_left)
    # all_paths.append(mc_ctg1_right)
    # all_paths.append(mc_ctg1_left_2)
    # all_paths.append(mc_ctg1_left_3)
    # all_paths.append(mc_ctg2_left)
    # all_paths.append(mc_ctg2_right)
    # all_paths.append(mc_ctg2_right_2)
    # all_paths.append(mc_ctg2_right_3)
    # all_paths.append(mc_ctg3_left)

    # mapa = divide_paths_into_anchor_groups(all_paths)
    # with open('map_of_paths', 'wb') as mc_ctg1_right_file:
    #     pickle.dump(mapa, mc_ctg1_right_file)

    # treba na neki nacin dobit
    lengths_of_contigs = {'ctg1':1000000, 'ctg2':1800000, 'ctg3':1541652}

    with open('map_of_paths', 'rb') as mc_ctg1_right_file:
        map_of_paths = pickle.load(mc_ctg1_right_file)
    with open('grouped_data_c_r_2', 'rb') as mc_ctg1_right_file:
        grouped_c_r = pickle.load(mc_ctg1_right_file)
    # print(grouped_c_r[0])
    with open('keys_c_r', 'rb') as mc_ctg1_right_file:
        keys_c_r = pickle.load(mc_ctg1_right_file)
    with open('grouped_data_r_r_2', 'rb') as mc_ctg1_right_file:
        grouped_r_r = pickle.load(mc_ctg1_right_file)
    with open('keys_r_r', 'rb') as mc_ctg1_right_file:
        keys_r_r = pickle.load(mc_ctg1_right_file)
    # print(map_of_paths['ctg1ctg2'])
    # print(keys_c_r)
    map_of_lengths = get_all_lengths(map_of_paths, lengths_of_contigs, grouped_c_r, grouped_r_r, keys_c_r, keys_r_r)
    print(map_of_lengths['ctg1ctg2'])
    with open('map_of_lengths', 'wb') as grouped_data_file:
        pickle.dump(map_of_lengths, grouped_data_file)
    # print(map['ctg2ctg3'])
    # for key in map:
    #     print(key)
    # print(mapa['ctg1ctg3'])
    # lengths = []
    # for path in mapa['ctg1ctg3']:
    #     # print(len(path))
    #     lengths.append(len(path))
    # print(lengths)
    # ind = len(lengths)//2
    # print(ind)
    # print(mapa['ctg1ctg3'][ind])
    
    # print(map['Ctg1Ctg2'])
    # print(mc_ctg3_left)
    # print('#######################3')
    # print(mc_ctg2_left)

    # with open('mc_ctg3_left', 'rb') as mc_ctg_right_file:
    #     mc_ctg3_left = pickle.load(mc_ctg_right_file)
    # print(mc_ctg3_left[0][0][0])