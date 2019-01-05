from utils import *

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
def dfs_2(data, start, goals):
    graph = dict()

    # replace with a function that finds the connecting reads for this anchoring node
    graph[start] = data[start]

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
                    # find new connecting reads for this node if they exist
                    if next in data:
                        graph[next] = data[next]
    return paths_to_goals


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

    overlaps_c_r = load_data(path_c_r)
    overlaps_r_r = load_data(path_r_r)
    print(len(overlaps_c_r))
    print(len(overlaps_r_r))
