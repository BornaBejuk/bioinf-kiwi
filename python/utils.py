import numpy as np
import numpy_indexed as npi

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
                            ('OS', float),
                            ('ES', float),
                            ('SI', float),
                            ('extension_side', 'U25'),
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
        # overlaps['a'][i] = row[12]
        # overlaps['b'][i] = row[13]
        # overlaps['c'][i] = row[14]
        # overlaps['d'][i] = row[15]

    overlaps = filter_contained(overlaps)
    overlaps = filter_SI(overlaps)
    return overlaps

def get_SI(overlaps):
    SI = overlaps['res_matches'] / overlaps['block_len']
    return SI

def get_avg_SI(overlaps):
    SI = get_SI(overlaps)
    return SI.mean()

def filter_SI(overlaps, SI_min=0.5):
    SI = get_SI(overlaps) > SI_min
    return overlaps[SI]

def get_OL(overlaps):
    OL1 = overlaps['query_end'] - overlaps['query_start']
    OL2 = overlaps['target_end'] - overlaps['target_start']
    return OL1, OL2

def extend_right(end1, len1, end2, len2):
    return (len2 - end2 - 1) < (len1 - end1 - 1)

def extend_left(start1, start2):
    return start2 < start1

def filter_contained(overlaps):
    to_delete = []
    for i in range(overlaps.shape[0]):
        if extend_right(overlaps['query_end'][i], overlaps['query_len'][i], overlaps['target_end'][i], overlaps['target_len'][i]):
            overlaps['extension_side'][i] = 'right'
            continue
        elif extend_left(overlaps['query_start'][i], overlaps['target_start'][i]):
            overlaps['extension_side'][i] = 'left'
            continue
        else:
            to_delete.append(i)

    return np.delete(overlaps, to_delete, None)

def get_OH(overlaps):
    """
    a)  1 ---------|---|+++
        2        ++|---|-----------

    b)  1        +++|---|--------------
        2 ----------|---|+++
    """
    OH1a = overlaps['query_len'] - overlaps['query_end'] - 1
    OH2a = overlaps['target_start']

    OH1b = overlaps['query_start']
    OH2b = overlaps['target_len'] - overlaps['target_end'] - 1

    final_OH = []
    for i in range(overlaps.shape[0]):
        if extend_right(overlaps['query_end'][i], overlaps['query_len'][i], overlaps['target_end'][i], overlaps['target_len'][i]):
            final_OH.append([OH1b[i], OH2b[i]])
        elif extend_left(overlaps['query_start'][i], overlaps['target_start'][i]):
            final_OH.append([OH1a[i], OH2a[i]])
        else:
            # TODO filtriraj one koji uopce ne produzuju nista ni sa jedne strane
            continue

    final_OH = np.array(final_OH).reshape(-1,2)

    return np.split(final_OH, 2, axis=1)

def get_EL(overlaps):
    """
    a)  1 +++++++++|---|---
        2        --|---|++++++++++

    b)  1        ---|---|+++++++++++++
        2 ++++++++++|---|---
    """
    EL1a = overlaps['query_start']
    EL2a = overlaps['target_len'] - overlaps['target_end'] - 1

    EL1b = overlaps['query_len'] - overlaps['query_end'] - 1
    EL2b = overlaps['target_start']

    final_EL = []
    for i in range(overlaps.shape[0]):
        if extend_right(overlaps['query_end'][i], overlaps['query_len'][i], overlaps['target_end'][i], overlaps['target_len'][i]):
            final_EL.append([EL1b[i], EL2b[i]])
        elif extend_left(overlaps['query_start'][i], overlaps['target_start'][i]):
            final_EL.append([EL1a[i], EL2a[i]])
        else:
            # TODO filtriraj one koji uopce ne produzuju nista ni sa jedne strane
            continue

    final_EL = np.array(final_EL).reshape(-1,2)

    return np.split(final_EL, 2, axis=1)


def get_OS(overlaps):
    OL1, OL2 = get_OL(overlaps)
    return 0.5 * (OL1 + OL2) * get_SI(overlaps)

def get_ES(overlaps):
    EL1, EL2 = get_EL(overlaps)
    OH1, OH2 = get_OH(overlaps)

    ES1 = get_OS(overlaps).reshape(-1,1) + .5 * EL2 - 0.5 * (OH1 + OH2)
    ES2 = get_OS(overlaps).reshape(-1,1) + .5 * EL1 - 0.5 * (OH1 + OH2)

    return ES1, ES2

def append_scores(overlaps):
    OS = get_OS(overlaps)
    ES1, ES2 = get_ES(overlaps)
    SI = get_SI(overlaps)
    for i in zip(range(overlaps.shape[0])):
        overlaps['OS'][i] = OS[i]
        overlaps['ES'][i] = ES2[i][0]
        overlaps['SI'][i] = SI[i]

    return overlaps

def get_grouped_data(overlaps):
    return npi.group_by(overlaps['target_name'])
