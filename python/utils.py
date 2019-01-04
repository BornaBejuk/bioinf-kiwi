import numpy as np


def load_data(path):
    overlaps2 = np.loadtxt(path, delimiter='\t', dtype='str')

    overlaps = np.recarray(shape=(overlaps2.shape[0],),
                     dtype=[('query_name', 'U25'),
                            ('query_len', int),
                            ('query_start', float),
                             ('query_end', float),
                            ('strand', 'U25'),
                            ('target_name', 'U25'),
                            ('target_len', int),
                            ('target_start', float),
                            ('target_end', float),
                            ('res_matches', float),
                            ('block_len', float),
                            ('map_quality', int),
                            ('a', 'U25'),
                            ('b', 'U25'),
                            ('c', 'U25'),
                            ('d', 'U25'),
                            ])
    for i, row in zip(range(overlaps2.shape[0]),overlaps2):
        overlaps['query_name'][i] = row[0]
        overlaps['query_len'][i] = row[1]
        overlaps['query_start'][i] = row[2]
        overlaps['query_end'][i] = row[3]
        overlaps['strand'][i] = row[4]
        overlaps['target_name'][i] = row[5]
        overlaps['target_len'][i] = row[6]
        overlaps['target_start'][i] = row[7]
        overlaps['target_end'][i] = row[8]
        overlaps['res_matches'][i] = row[9]
        overlaps['block_len'][i] = row[10]
        overlaps['map_quality'][i] = row[11]
        overlaps['a'][i] = row[12]
        overlaps['b'][i] = row[13]
        overlaps['c'][i] = row[14]
        overlaps['d'][i] = row[15]

    return overlaps

def get_SI(overlaps):
    matches = overlaps[:,9:10].astype(float)
    block = overlaps[:,10:11].astype(float)
    SI = matches/block
    return SI


def get_avg_SI(overlaps):
    SI = get_SI(overlaps)
    return SI.mean()

def filter_SI(overlaps, SI_min=0.5):
    SI = get_SI(overlaps) > SI_min
    return overlaps[np.repeat(SI,overlaps.shape[1], axis=1)].reshape(-1,overlaps.shape[1])

def get_OL(overlaps):
    OL1 = overlaps[:,3:4] - overlaps[:,2:3]
    OL2 = overlaps[:,8:9] - overlaps[:,7:8]
    return OL1, OL2

def get_OH(overlaps):
    """
    a)  1 ---------|---|----
        2        --|---|-----------

    b)  1        ---|---|--------------
        2 ----------|---|----
    """
    OH1a = overlaps[:,1:2] - overlaps[:,3:4] - 1
    OH2a = overlaps[:,7:8]

    OH1b = overlaps[:,2:3]
    OH2b = overlaps[:,6:7] - overlaps[:,8:9] - 1

    if OH2a < OH2b:
        return OH1a, OH2a
    return OH1b, OH2b

def get_OS(overlaps):
    return get_OH(overlaps)
