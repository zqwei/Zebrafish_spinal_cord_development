
import numpy as np
from re import search
from os import system, path, makedirs
from IO import imread, imsave

def recover_filename_from_pattern(input_pattern, t):
    match = search('\?+', input_pattern)
    digits = len(match.group())
    return input_pattern.replace('?'*digits, '%0{}d'.format(digits) % t)

def pre_post_corr(input_folder_list):
    for n_file in range(0, len(input_folder_list)):
        fixed_file = recover_filename_from_pattern(input_folder_list[n_file][0] + '/SPM00/TM??????/SPM00_TM??????_CM01_CHN00.klb', input_folder_list[n_file][1])
        moving_file = recover_filename_from_pattern(input_folder_list[n_file][3] + '/SPM00/TM??????/SPM00_TM??????_CM01_CHN00.klb', input_folder_list[n_file][4])
        output_folder = input_folder_list[n_file][0] + '/Warp/'
        if not path.exists(output_folder):
            makedirs(output_folder)


        # im = imread(moving_file)
        # im.voxelsize = (1, 1, 6)
        # imsave(output_folder + 'A02.inr', im)
        #
        # im = imread(fixed_file)
        # im.voxelsize = (1, 1, 6)
        # imsave(output_folder + 'A01.inr', im)


        system('blockmatching -ref ' + fixed_file + ' -flo ' + moving_file + ' \
         -res ' + output_folder + 'after_ref_warped.klb -floating-voxel 1 1 6 -reference-voxel 1 1 6 \
         -res-trsf ' + output_folder + '3_NL.inr \
         -trsf-type vectorfield -estimator wlts \
         -lts-fraction 0.55 \
         -py-hl 6 -py-ll 1 \
         -py-gf -elastic-sigma 5.0 5.0 5.0 -fluid-sigma 5.0 5.0 5.0')

        # im = imread(output_folder + 'res.inr')
        # imsave(output_folder + 'after_ref_warped.klb', im)









