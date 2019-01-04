import numpy as np

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
