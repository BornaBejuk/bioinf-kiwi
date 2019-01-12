import pickle
import numpy as np

def divide_paths_into_anchor_groups(list_of_list_of_paths):
    '''Takes as input an object that is in format:
        [ [ [path1], ..., [pathN] ],
          [ [path1]], ...,[pathN] ],
          ...
          [ [path1], ...,[pathN] ] ]
    where e.g. path1 == ['ctg1', 'read1', 'ctg2']
    Same as before, paths do not have to be inside list. i.e. input should probably look like
    [ [ path1, ..., pathN],
      [ path1, ..., pathN],
      ...
      [ path1, ..., pathN] ]
     That syntax was left over from first two approaches. 
     That can be changed if wanted when rewriting in cpp.

    Method outputs a map, e.g.
        {'ctg1ctg2': [[path1],..., [pathN]],
         'ctg1ctg3': [[path1],..., [pathN]],
         'ctg2ctg3': [[path1],..., [pathN]]}
    The goal of the method is to group all paths so that paths between same contigs are put together.
    This method could be modified so that the keys are tuples (ctg1, ctg2)
    instead of a string 'ctgctg2'.
    '''
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
    '''
    Takes as input a map of paths that the method 'divide_paths_into_anchor_groups' produced 
    and outputs a map in format:
    {'ctg1ctg2': [(path1, len_of_path1), ..., (pathN, len_of_pathN)],
     'ctg1ctg3': [(path1, len_of_path1), ..., (pathN, len_of_pathN)],
     'ctg2ctg3': [(path1, len_of_path1), ..., (pathN, len_of_pathN)]}
    '''
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


def divide_paths_into_groups(map_of_lengths, limit):
    '''
    Takes as input a map of lengths that the method 'get_all_lengths' produced 
    and outputs a map in format:
    {'ctg1ctg2': [group1, group2, ..., groupN],
     'ctg1ctg3': [group1, group2, ..., groupN],
     'ctg2ctg3': [group1, group2, ..., groupN]}
     where a group is a list of (path, len_of_path) tuples.
    '''
    map_of_groups = dict()
    for m in map_of_lengths:
        paths = map_of_lengths[m]
        path_lengths = [path[1] for path in paths]
        path_lengths_sorted = sorted(path_lengths)
        lowest = path_lengths_sorted[0]
        highest = path_lengths_sorted[-1]
        groups = [[] for i in range(len(paths))]
        group_window = (highest - lowest) / len(paths)
        if len(map_of_lengths[m]) < limit:
            first_group = map_of_lengths[m]
            map_of_groups[m] = [first_group]
            continue
        for path in paths:
            length = path[1]
            # da se ovaj najveci (highest) ne nadje uvijek sam u zadnjoj grupi
            if length == highest:
                length -= 1
            length = length - lowest
            position = int(length // group_window)
            groups[position].append(path)
        map_of_groups[m] = groups
    return map_of_groups


def choose_path_from_groups(map_of_groups):
    '''
    Takes as input a map of groups that the method 'divide_paths_into_groups' produced 
    and outputs a map in format:
    {'ctg1ctg2': chosen_path,
     'ctg1ctg3': chosen_path,
     'ctg2ctg3': chosen_path}
     where chosen path is a tuple (path, len_of_path)
    '''
    map_of_chosen_paths = dict()
    for m in map_of_groups:
        groups = map_of_groups[m]

        # uzmi grupu iz groups
        list_of_group_sizes = []
        for group in groups:
            list_of_group_sizes.append(len(group))

        # ako ih je vise koji su jednako veliki, argmax vraca prvi koji ima tu max velicinu
        # sto i zelim, pripazit na to
        max_indeks = np.argmax(list_of_group_sizes)
        paths = groups[max_indeks]

        length = len(paths)
        paths.sort(key = lambda x:x[1])
        index = length // 2
        
        map_of_chosen_paths[m] = paths[index]
    return map_of_chosen_paths

    


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
    # print(mc_ctg1_left_3)
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
    # lengths_of_contigs = {'ctg1':1000000, 'ctg2':1800000, 'ctg3':1541652}

    # with open('map_of_paths', 'rb') as mc_ctg1_right_file:
    #     map_of_paths = pickle.load(mc_ctg1_right_file)
    # with open('grouped_data_c_r_2', 'rb') as mc_ctg1_right_file:
    #     grouped_c_r = pickle.load(mc_ctg1_right_file)
    # # print(grouped_c_r[0])
    # with open('keys_c_r', 'rb') as mc_ctg1_right_file:
    #     keys_c_r = pickle.load(mc_ctg1_right_file)
    # with open('grouped_data_r_r_2', 'rb') as mc_ctg1_right_file:
    #     grouped_r_r = pickle.load(mc_ctg1_right_file)
    # with open('keys_r_r', 'rb') as mc_ctg1_right_file:
    #     keys_r_r = pickle.load(mc_ctg1_right_file)
    # # print(map_of_paths['ctg1ctg2'])
    # # print(keys_c_r)
    # print(map_of_paths)
    # map_of_lengths = get_all_lengths(map_of_paths, lengths_of_contigs, grouped_c_r, grouped_r_r, keys_c_r, keys_r_r)
    # print(map_of_lengths['ctg1ctg2'])
    # with open('map_of_lengths', 'wb') as grouped_data_file:
    #     pickle.dump(map_of_lengths, grouped_data_file)



    with open('map_of_lengths', 'rb') as grouped_data_file:
        map_of_lengths = pickle.load(grouped_data_file)
    map_of_groups = divide_paths_into_groups(map_of_lengths, 10)
    # print(map_of_groups['ctg1ctg3'])
    # for g in map_of_groups['ctg1ctg2']:
    #     print(g)
    #     print('\n')
    map_of_chosen_paths = choose_path_from_groups(map_of_groups)
    print(map_of_chosen_paths)
    # mapa = map_of_lengths['ctg1ctg2']
    # list_of_lengths = []
    # for path in mapa:
    #     print(path)
    #     list_of_lengths.append(path[1])
    # print(sorted(list_of_lengths))
    # print(len(list_of_lengths))
    # print(len(np.unique(list_of_lengths)))




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