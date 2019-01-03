import numpy as np

def get_SI(overlaps):
    matches = overlaps[:,9:10].astype(float)
    block = overlaps[:,10:11].astype(float)
    SI = matches/block
    return SI


def get_avg_SI(overlaps):
    SI = get_SI(overlaps)
    return SI.mean()
