from utils import *
import pickle
import numpy_indexed
import copy

def dfs(graph, start, goals):
    paths_to_goals = []
    stack = [(start, [start])]
    while stack:
        (vertex, path) = stack.pop()
        if vertex in graph:
            for next in graph[vertex] - set(path):
                if next in goals:
                    paths_to_goals.append(path + [next])
                else:
                    stack.append((next, path + [next]))
    return paths_to_goals


# version of dfs that builds the graph as it traverses it
def dfs_2(start, goals, side, grouped_c_r, keys_c_r, grouped_r_r, keys_r_r, max_depth = 12):
    graph = dict()

    # replace with a function that finds the connecting reads for this anchoring node
    graph[start] = set(get_n_best_connecting_reads_for_contig(start, grouped_c_r, keys_c_r, side, n = 10))

    paths_to_goals = []
    stack = [(start, [start])]
    while stack:
        (vertex, path) = stack.pop()
        if len(path) <= max_depth:
            if vertex in graph:
                for next in graph[vertex] - set(path):
                    if next in goals:
                        paths_to_goals.append(path + [next])
                        print(paths_to_goals)
                    else:
                        stack.append((next, path + [next]))
                        # find new connecting reads for this node if they exist
                        # if next in graph:
                        graph[next] = set(get_n_best_connecting_reads_for_read(next, side, grouped_r_r, keys_r_r, grouped_c_r, keys_c_r, n=2))
    return paths_to_goals


def get_n_best_connecting_reads_for_contig(contig, grouped, keys, side, n = 2, first_approach = True):
    group = grouped[keys.index(contig)]
    connecting_reads = group[np.where(group['extension_side'] == side)]
    if first_approach:
        ind = np.lexsort((connecting_reads['SI'], connecting_reads['OS']))
    else:
        ind = np.lexsort((connecting_reads['SI'], connecting_reads['ES']))
    return connecting_reads[ind][-n:]['query_name']


def get_n_best_connecting_reads_for_read(read, side, grouped_r_r, keys_r_r, grouped_c_r, keys_c_r, n = 2, first_approach = True):
    if side == 'right':
        other_side = 'left'
    else:
        other_side = 'right'
    
    second_group = []
    for i in range(len(keys_c_r)):
        group = grouped_c_r[i]
        for j in group:
            if j['query_name'] == read and j['extension_side'] == other_side:
                new_row = copy.deepcopy(j)
                new_row['query_name'] = keys_c_r[i]
                second_group.append(new_row)
    
    if read in keys_r_r:
        group = grouped_r_r[keys_r_r.index(read)]
        connecting_reads = group[np.where(group['extension_side'] == side)]
        final_group = connecting_reads
        if second_group != []:
            final_group = np.append(connecting_reads, second_group)
        if first_approach:
            ind = np.lexsort((final_group['SI'], final_group['OS']))
        else:
            ind = np.lexsort((final_group['SI'], final_group['ES']))
        return final_group[ind][-n:]['query_name']
    else:
        return []
    

if __name__ == '__main__':
    # data = {'A': set(['r1', 'r2', 'r3']),
    #      'r1': set(['r4', 'r5']),
    #      'r2': set(['r6', 'r7']),
    #      'r3': set(['r8', 'r9']),
    #      'r4': set(['A1', 'r10']),
    #      'r5': set(['r11', 'r12']),
    #      'r9': set(['r13', 'r14']),
    #      'r11': set(['A2', 'r16']),
    #      'r10': set(['r4']),
    #      'r13': set(['A3', 'r15']),
    #      'r15': set(['r17', 'r18'])}
    # start = 'A'
    # goals = ['A1', 'A2', 'A3']
    # paths = dfs_2(data, start, goals)
    # print(paths)
    path_c_r = 'overlaps-c-r.paf'
    path_r_r = 'overlaps-r-r.paf'

    # overlaps_c_r = load_data(path_c_r)
    # overlaps_r_r = load_data(path_r_r)
    # print(len(overlaps_c_r))
    # print(len(overlaps_r_r))

    with open('grouped_data_c_r', 'rb') as grouped_c_r_file:
        grouped_c_r = pickle.load(grouped_c_r_file)
    with open('grouped_data_r_r', 'rb') as grouped_r_r_file:
        grouped_r_r = pickle.load(grouped_r_r_file)
    with open('keys_c_r', 'rb') as keys_c_r_file:
        keys_c_r = pickle.load(keys_c_r_file).tolist()
    with open('keys_r_r', 'rb') as keys_r_r_file:
        keys_r_r = pickle.load(keys_r_r_file).tolist()

    # print(grouped_c_r[0][np.where(grouped_c_r[0]['extension_side'] == 'right')])

    # test
    # print(type(keys_c_r))
    start = keys_c_r[0]
    goals = keys_c_r
    side = 'left'
    paths = dfs_2(start, goals, side, grouped_c_r, keys_c_r, grouped_r_r, keys_r_r)
    with open('test_ctg1_paths_right_side', 'wb') as pahts_right_side_file:
        pickle.dump(paths, pahts_right_side_file)
    print(paths)
