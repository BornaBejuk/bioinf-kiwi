from utils import *
import pickle
import numpy_indexed
import copy

# version of dfs that builds the graph as it traverses it
def dfs_2(start, goals, side, grouped_c_r, keys_c_r, grouped_r_r, keys_r_r, max_depth = 11):
    '''
    Parameters:
    start : str
        name of the root node, e.g. 'ctg1'
    goals : list of str
        list of goal nodes, e.g. goals = ['ctg1', 'ctg2']
    side : str
        side to which we try to build path, right or left
        if side == 'right', continuing read is building path on the right side of target node
    grouped_c_r : numpy indexed object
        object that contains info about overlaps 
        e.g. grouped_c_r[0] gives info about continuing reads for contig1, grouped_c_r[1] for contig2, etc.
    grouped_r_r : numpy indexed object
        overlaps between reads, similair to grouped_c_r 
    keys_c_r : list of str
        a list of keys for grouped_r_r object, e.g. ['ctg1', 'ctg2']
        keys_r_r are used to get the index of contig in grouped_c_r object based on its name
    keys_r_r : list of str
        similair to keys_c_r

    Returns:
    paths_to_goals: list of lists of str
    e.g. [['ctg1', 'read1', 'ctg2'], ['ctg1', 'read2', 'ctg2']]

    a list of paths, where a path is a list of strings

    '''
    graph = dict()
    graph[start] = set(monte_carlo_extending_for_contig(start, side, grouped_c_r, keys_c_r))
    paths_to_goals = []
    stack = [(start, [start])]
    while stack:
        # print(stack)
        (vertex, path) = stack.pop()
        if len(path) <= max_depth:
            if vertex in graph:
                for next in graph[vertex] - set(path):
                    if next in goals:
                        paths_to_goals.append(path + [next])
                        # print(paths_to_goals)
                    else:
                        stack.append((next, path + [next]))
                        # find new connecting reads for this node if they exist
                        # if next in graph:
                        graph[next] = set(monte_carlo_extending_for_read(next, side, grouped_c_r, keys_c_r, grouped_r_r, keys_r_r))
    return paths_to_goals


# def get_n_best_connecting_reads_for_contig(contig, grouped, keys, side, n = 2, first_approach = True):
#     group = grouped[keys.index(contig)]
#     connecting_reads = group[np.where(group['extension_side'] == side)]
#     if first_approach:
#         ind = np.lexsort((connecting_reads['SI'], connecting_reads['OS']))
#     else:
#         ind = np.lexsort((connecting_reads['SI'], connecting_reads['ES']))
#     return connecting_reads[ind][-n:]['query_name']


# def get_n_best_connecting_reads_for_read(read, side, grouped_r_r, keys_r_r, grouped_c_r, keys_c_r, n = 2, first_approach = True):
#     if side == 'right':
#         other_side = 'left'
#     else:
#         other_side = 'right'
    
#     second_group = []
#     for i in range(len(keys_c_r)):
#         group = grouped_c_r[i]
#         for j in group:
#             if j['query_name'] == read and j['extension_side'] == other_side:
#                 new_row = copy.deepcopy(j)
#                 new_row['query_name'] = keys_c_r[i]
#                 second_group.append(new_row)
#     # print(type(second_group))
#     # print(second_group)
#     if read in keys_r_r:
#         group = grouped_r_r[keys_r_r.index(read)]
#         connecting_reads = group[np.where(group['extension_side'] == side)]
#         final_group = connecting_reads
#         # final_group = []
#         if second_group != []:
#             final_group = np.append(connecting_reads, second_group)
#         # if second_group != []:
#         #     final_group.append(second_group)
#         if first_approach:
#             ind = np.lexsort((final_group['SI'], final_group['OS']))
#         else:
#             ind = np.lexsort((final_group['SI'], final_group['ES']))
#         return final_group[ind][-n:]['query_name']
#     else:
#         return []


def monte_carlo_extending_for_contig(contig, side, grouped_c_r, keys_c_r,):
    '''
    Parameters:
    contig : str
        name of the starting contig, e.g. 'ctg1'
    side : str
        side to which we try to build path, right or left
        if side == 'right', continuing read is building path on the right side of target node
    grouped_c_r : numpy indexed object
        object that contains info about overlaps 
        e.g. grouped_c_r[0] gives info about continuing reads for contig1, grouped_c_r[1] for contig2, etc.
    keys_c_r : list of str
        a list of keys for grouped_r_r object, e.g. ['ctg1', 'ctg2']
        keys_r_r are used to get the index of contig in grouped_c_r object based on its name

    Returns:
    chosen_read : list of one str
    read that the monte carlo approach found as a continuing read for a given contig
    in a list because method dfs_2 expects a list of continuing reads. dfs_2 written that way
    because of first two approaches
    In case it doesnt find it, method returns empty list so that dfs_2 method doesnt crash

    subject to change ???!!!
    Karlo, feel free to change this so that it returns a string, and change dfs_2 so that is accepts a string
    if that suits you
    '''
    group = grouped_c_r[keys_c_r.index(contig)]
    group = group[np.where(group['extension_side'] == side)]
    group = group[np.where(group['ES'] >= 0)]
    if group.size != 0:
        reads_ES = group['ES']
        sum_ES = np.sum(reads_ES)
        probabilities = [x / sum_ES for x in reads_ES]
        chosen_read = np.random.choice(a = group, p = probabilities)
        chosen_read = [chosen_read['query_name']]
        return chosen_read
    else:
        return []
    

def monte_carlo_extending_for_read(read, side, grouped_c_r, keys_c_r, grouped_r_r, keys_r_r):
    '''
    Parameters:
    read : str
        name of the read for which we want to find continuing read, e.g. 'read00291'
    side : str
        side to which we try to build path, right or left
        if side == 'right', continuing read is building path on the right side of target node
    grouped_c_r/grouped_r_r/keys_c_r/keys_r_r
        look at previous comments in other methods

    Returns:
    chosen_read : list of one str
    read that the monte carlo approach found as a continuing read for a given read
    in a list because method dfs_2 expects a list of continuing reads. dfs_2 written that way
    because of first two approaches
    In case it doesnt find it, method returns empty list so that dfs_2 method doesnt crash

    subject to change ???
    Karlo, feel free to change this so that it returns a string, and change dfs_2 so that is accepts a string
    if that suits you
    '''
    if side == 'right':
        other_side = 'left'
    else:
        other_side = 'right'
    
    second_group = []
    for i in range(len(keys_c_r)):
        group = grouped_c_r[i]
        for j in group:
            if j['query_name'] == read and j['extension_side'] == other_side and j['ES'] > 0:
                new_row = copy.deepcopy(j)
                new_row['query_name'] = keys_c_r[i]
                second_group.append(new_row)
    if read in keys_r_r:
        group = grouped_r_r[keys_r_r.index(read)]
        group = group[np.where(group['extension_side'] == side)]
        group = group[np.where(group['ES'] >= 0)]
        if second_group != []:
            group = np.append(group, second_group)
        if group.size != 0:
            reads_ES = group['ES']
            sum_ES = np.sum(reads_ES)
            probabilities = [x / sum_ES for x in reads_ES]
            chosen_read = np.random.choice(a = group, p = probabilities)
            chosen_read = [chosen_read['query_name']]
            return chosen_read
        else:
            return []
    else:
        return []


def try_monte_carlo(start, goals, side, grouped_c_r, keys_c_r, grouped_r_r, keys_r_r, max_depth=50, n_times = 100):
    '''
    Parameters:
    start : str
        name of the root node, e.g. 'ctg1'
    goals : list of str
        list of goal nodes, e.g. goals = ['ctg1', 'ctg2']
    side : str
        side to which we try to build path, right or left
        if side == 'right', continuing read is building path on the right side of target node
    grouped_c_r : numpy indexed object
        object that contains info about overlaps 
        e.g. grouped_c_r[0] gives info about continuing reads for contig1, grouped_c_r[1] for contig2, etc.
    grouped_r_r : numpy indexed object
        overlaps between reads, similair to grouped_c_r 
    keys_c_r : list of str
        a list of keys for grouped_r_r object, e.g. ['ctg1', 'ctg2']
        keys_r_r are used to get the index of contig in grouped_c_r object based on its name
    keys_r_r : list of str
        similair to keys_c_r
    max_depth : int
        maximum number of nodes in a path
    n_times : int
        number of times that we try to get a path with monte carlo approach
        given a 1000 tries, we usually find around 5-10 paths

    Returns:
    all_paths: list of lists of paths
    e.g. [[['ctg1', 'read1', 'ctg2']], [['ctg1', 'read2', 'ctg2']]]
    probably one level of unnecessary complication, built that way because of first two approaches
    where we expected to find more than one path per dfs search, SUBJECT TO CHANGE, you could
    change this so that all_paths is a list of paths, not a list inside list
    '''
    all_paths = []
    for i in range(n_times):
        paths = dfs_2(start, goals, side, grouped_c_r, keys_c_r, grouped_r_r, keys_r_r, max_depth=max_depth)
        if paths != []:
            all_paths.append(paths)
    return all_paths


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
    
    # testing on fake data
    # grouped_c_r = [np.array([('read1', 1, 2, 3, 'right')], dtype=(np.record, [('query_name', '<U250'), ('SI', '<f8'), ('OS', '<f8'), ('ES', '<f8'), ('extension_side', '<U25')])), np.array([('read3', 1, 2, 3, 'left')], dtype=(np.record, [('query_name', '<U250'), ('SI', '<f8'), ('OS', '<f8'), ('ES', '<f8'), ('extension_side', '<U25')]))]
    # keys_c_r = ['ctg1', 'ctg2']
    # grouped_r_r = [np.array([('read2', 1, 2, 3, 'right')], dtype=(np.record, [('query_name', '<U250'), ('SI', '<f8'), ('OS', '<f8'), ('ES', '<f8'), ('extension_side', '<U25')])), np.array([('read3', 1, 2, 3, 'right')], dtype=(np.record, [('query_name', '<U250'), ('SI', '<f8'), ('OS', '<f8'), ('ES', '<f8'), ('extension_side', '<U25')])), np.array([], dtype=(np.record, [('query_name', '<U250'), ('SI', '<f8'), ('OS', '<f8'), ('ES', '<f8'), ('extension_side', '<U25')]))]
    # keys_r_r = ['read1', 'read2', 'read3']


    with open('grouped_data_c_r', 'rb') as grouped_c_r_file:
        grouped_c_r = pickle.load(grouped_c_r_file)
    with open('grouped_data_r_r', 'rb') as grouped_r_r_file:
        grouped_r_r = pickle.load(grouped_r_r_file)
    with open('keys_c_r', 'rb') as keys_c_r_file:
        keys_c_r = pickle.load(keys_c_r_file).tolist()
    with open('keys_r_r', 'rb') as keys_r_r_file:
        keys_r_r = pickle.load(keys_r_r_file).tolist()

    start = keys_c_r[1]
    goals = keys_c_r
    side = 'left'

    paths = try_monte_carlo(start, goals, side, grouped_c_r, keys_c_r, grouped_r_r, keys_r_r, max_depth=30, n_times = 500)
    # with open('mc_ctg3_right', 'wb') as paths_right_side_file:
    #     pickle.dump(paths, paths_right_side_file)
    print(paths)