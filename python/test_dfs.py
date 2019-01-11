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
def dfs_2(start, goals, side, grouped_c_r, keys_c_r, grouped_r_r, keys_r_r, max_depth = 11, first_approach=True, monte_carlo=False):
    graph = dict()

    # replace with a function that finds the connecting reads for this anchoring node
    if monte_carlo == False:
        graph[start] = set(get_n_best_connecting_reads_for_contig(start, grouped_c_r, keys_c_r, side, n = 5, first_approach=first_approach))
    else:
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
                        if monte_carlo == False:
                            graph[next] = set(get_n_best_connecting_reads_for_read(next, side, grouped_r_r, keys_r_r, grouped_c_r, keys_c_r, n=2, first_approach=first_approach))
                        else:
                            graph[next] = set(monte_carlo_extending_for_read(next, side, grouped_c_r, keys_c_r, grouped_r_r, keys_r_r))
    return paths_to_goals


def get_n_best_connecting_reads_for_contig(contig, grouped, keys, side, n = 2, first_approach = True):
    group = grouped[keys.index(contig)]
    connecting_reads = group[np.where(group['extension_side'] == side)]
    if first_approach:
        ind = np.lexsort((connecting_reads['SI'], connecting_reads['OS']))
    else:
        ind = np.lexsort((connecting_reads['SI'], connecting_reads['ES']))
    return connecting_reads[ind][-n:]['query_name']

# def get_n_best_connecting_reads_for_contig(contig, grouped, keys, side, n = 2, first_approach = True):
#     group = grouped[keys.index(contig)]
#     # connecting_reads = group[np.where(group['extension_side'] == side)]
#     connecting_reads = group[np.where(group[4] == side)]
#     if first_approach:
#         ind = np.lexsort((connecting_reads[1], connecting_reads[2]))
#         # ind = [connecting_reads[0]]
#     else:
#         ind = np.lexsort((connecting_reads[1], connecting_reads[3]))
#     return connecting_reads[ind][-n:][0]


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
    # print(type(second_group))
    # print(second_group)
    if read in keys_r_r:
        group = grouped_r_r[keys_r_r.index(read)]
        connecting_reads = group[np.where(group['extension_side'] == side)]
        final_group = connecting_reads
        # final_group = []
        if second_group != []:
            final_group = np.append(connecting_reads, second_group)
        # if second_group != []:
        #     final_group.append(second_group)
        if first_approach:
            ind = np.lexsort((final_group['SI'], final_group['OS']))
        else:
            ind = np.lexsort((final_group['SI'], final_group['ES']))
        return final_group[ind][-n:]['query_name']
    else:
        return []

def get_n_best_connecting_reads_for_read2(read, side, grouped_r_r, keys_r_r, grouped_c_r, keys_c_r, n = 2, first_approach = True):
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
    # print(type(second_group))
    # print(second_group)
    if read in keys_r_r:
        group = grouped_r_r[keys_r_r.index(read)]
        connecting_reads = group[np.where(group['extension_side'] == side)]
        final_group = connecting_reads
        # final_group = []
        # if second_group != []:
        #     final_group = np.append(connecting_reads, second_group)
        if second_group != []:
            final_group.append(second_group)
        
        # print(final_group)
        if first_approach:
            ind = np.lexsort((final_group[0][0][1], final_group[0][0][2]))
        else:
            ind = np.lexsort((final_group['SI'], final_group['ES']))
        a= final_group[ind]
        b = a[-n:]
        c = b[0]
        return [c[0]]
    else:
        return []


def monte_carlo_extending_for_contig(contig, side, grouped_c_r, keys_c_r,):
    group = grouped_c_r[keys_c_r.index(contig)]
    group = group[np.where(group['extension_side'] == side)]
    group = group[np.where(group['ES'] >= 0)]
    if group.size != 0:
        reads_ES = group['ES']
        sum_ES = np.sum(reads_ES)
        probabilities = [x / sum_ES for x in reads_ES]
        chosen_read = np.random.choice(a = group, p = probabilities)
        return [chosen_read['query_name']]
    else:
        return []
    

def monte_carlo_extending_for_read(read, side, grouped_c_r, keys_c_r, grouped_r_r, keys_r_r):
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
            return [chosen_read['query_name']]
        else:
            return []
    else:
        return []

        
def try_monte_carlo(start, goals, side, grouped_c_r, keys_c_r, grouped_r_r, keys_r_r, max_depth=50, n_times = 10):
    all_paths = []
    for i in range(n_times):
        paths = dfs_2(start, goals, side, grouped_c_r, keys_c_r, grouped_r_r, keys_r_r, max_depth=max_depth, monte_carlo=True)
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
    path_c_r = 'overlaps-c-r.paf'
    path_r_r = 'overlaps-r-r.paf'

    # overlaps_c_r = load_data(path_c_r)
    # overlaps_r_r = load_data(path_r_r)
    # print(len(overlaps_c_r))
    # print(len(overlaps_r_r))

    with open('grouped_data_c_r', 'rb') as grouped_c_r_file:
        grouped_c_r = pickle.load(grouped_c_r_file)
    # print(type(grouped_c_r[0]['ES']))
    # print((grouped_c_r[0]['ES']))
    # print(grouped_c_r[0]['extension_side'])
    # print(grouped_c_r)
    with open('grouped_data_r_r', 'rb') as grouped_r_r_file:
        grouped_r_r = pickle.load(grouped_r_r_file)
    with open('keys_c_r', 'rb') as keys_c_r_file:
        keys_c_r = pickle.load(keys_c_r_file).tolist()
    with open('keys_r_r', 'rb') as keys_r_r_file:
        keys_r_r = pickle.load(keys_r_r_file).tolist()
    # print(grouped_r_r[10])

    # read = monte_carlo_extending_for_read('read00246', 'left', grouped_c_r, keys_c_r, grouped_r_r, keys_r_r)
    # read = monte_carlo_extending_for_contig('ctg1', 'right', grouped_c_r, keys_c_r)
    # print(read)

    # print(grouped_c_r[0][np.where(grouped_c_r[0]['extension_side'] == 'right')])
    # print(grouped_c_r)
    # print(type(grouped_c_r[0]))
    # test
    # print(type(keys_c_r))
    # grouped_c_r = [np.array([('read1', 1, 2, 3, 'right')], dtype=(np.record, [('query_name', '<U250'), ('SI', '<f8'), ('OS', '<f8'), ('ES', '<f8'), ('extension_side', '<U25')])), np.array([('read3', 1, 2, 3, 'left')], dtype=(np.record, [('query_name', '<U250'), ('SI', '<f8'), ('OS', '<f8'), ('ES', '<f8'), ('extension_side', '<U25')]))]
    # keys_c_r = ['ctg1', 'ctg2']
    # grouped_r_r = [np.array([('read2', 1, 2, 3, 'right')], dtype=(np.record, [('query_name', '<U250'), ('SI', '<f8'), ('OS', '<f8'), ('ES', '<f8'), ('extension_side', '<U25')])), np.array([('read3', 1, 2, 3, 'right')], dtype=(np.record, [('query_name', '<U250'), ('SI', '<f8'), ('OS', '<f8'), ('ES', '<f8'), ('extension_side', '<U25')])), np.array([], dtype=(np.record, [('query_name', '<U250'), ('SI', '<f8'), ('OS', '<f8'), ('ES', '<f8'), ('extension_side', '<U25')]))]
    # keys_r_r = ['read1', 'read2', 'read3']
    start = keys_c_r[2]
    goals = keys_c_r
    side = 'right'

    # # # # print(type(grouped_c_r))
    # paths = dfs_2(start, goals, side, grouped_c_r, keys_c_r, grouped_r_r, keys_r_r, max_depth=50, monte_carlo=True)
    paths = try_monte_carlo(start, goals, side, grouped_c_r, keys_c_r, grouped_r_r, keys_r_r, max_depth=50, n_times = 1000)

    with open('mc_ctg3_right', 'wb') as pahts_right_side_file:
        pickle.dump(paths, pahts_right_side_file)
    print(paths)
