import numpy as np
from re import search
from os import system, path, makedirs

def recover_filename_from_pattern(input_pattern, t):
    match = search('\?+', input_pattern)
    digits = len(match.group())
    return input_pattern.replace('?'*digits, '%0{}d'.format(digits) % t)

def read_MDF_file(filename):
    with open(filename, "r") as fi:
        points = np.zeros(shape=(0, 3))
        for ln in fi:
            if ln.startswith("Point 1"):
                record = np.fromstring(ln[8:], dtype=np.float, sep=' ')
                points = np.vstack([points, record[[1, 0, 3]]])
    return(points)

def transform_points(input_folder_list):
    for n_file in range(0, len(input_folder_list)):
        segmentation_file = recover_filename_from_pattern(input_folder_list[n_file][0] + '/SPM00/TM??????/segmentation_curated.mdf', input_folder_list[n_file][1])
        points = read_MDF_file(segmentation_file)
        points[:, 2] = (points[:, 2] - 1) * 6
        np.savetxt(input_folder_list[n_file][0] + '/Warp/points_before.txt', points, fmt='%.8f')
        system('applyTrsfToPoints ' + input_folder_list[n_file][0] + '/Warp/points_before.txt ' + input_folder_list[n_file][0] + '/Warp/points_after.txt \
                -trsf ' + input_folder_list[n_file][0] + '/Warp/3_NL.inr')