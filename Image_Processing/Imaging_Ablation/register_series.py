from drift_correction import drift_correction
from pre_post_corr import pre_post_corr
from transform_points import transform_points

input_folder_list = [
    ['/media/P/SV2/YW_18-02-03/fish1_before_20180203_191939.corrected/', 1225, 1225, '/media/P/SV2/YW_18-02-03/fish1_after2_20180203_200124.corrected/', 240, 3339],
    ['/media/P/SV2/YW_18-02-03/fish2_before_20180203_203907.corrected/', 1225, 1225, '/media/P/SV2/YW_18-02-03/fish2_after_20180203_205131.corrected/', 240, 3567],
    ['/media/P/SV2/YW_18-02-03/fish3_before_20180203_213738.corrected/', 1225, 1225, '/media/P/SV2/YW_18-02-03/fish3_after_20180203_215441.corrected/', 240, 7162],
    ['/media/P/SV2/YW_18-02-03/fish4_before_20180203_224436.corrected/', 1225, 1225, '/media/P/SV2/YW_18-02-03/fish4_after_20180203_225915.corrected/', 240, 7162],
    ['/media/P/SV2/YW_18-02-03/fish5_before_20180203_234605.corrected/', 1225, 1225, '/media/P/SV2/YW_18-02-03/fish5_after_20180203_235929.corrected/', 240, 2838],
    ['/media/P/SV2/YW_18-02-03/fish6_before_20180204_003041.corrected/', 1225, 1225, '/media/P/SV2/YW_18-02-03/fish6_after_20180204_004222.corrected/', 240, 3600]
]


# drift_correction(input_folder_list)
# pre_post_corr(input_folder_list)
transform_points(input_folder_list)