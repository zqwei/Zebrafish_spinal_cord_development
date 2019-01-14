%%
% *profile.mat* contains the following information:
% profile_all     is a matrix containing (*cells*) rows and (*time points*) columns of the graylevel intensity information. The temporal sampling is 4 Hz. 
% timepoints      correspondence of TM index to imaging data
% tracks_smoothed contains the xyz coordinates in imaging space of cells over time
% mnx             binary variable indicating the identity of mnx1 marker
% mnx_level_func  is the intensity level from the mnx channel during functional imaging
% side            is manually identified to indicate the cell on the left (1) or right (2)
% x, y, z         AP, ML and DV location of each cell on the spinal cord atlas
% lineage_id      correspondence to MaMuT file of curated cell tracking
% islet          (only available for IHC dataset) binary variable indicating the identity of islet1/2 marker
% birthtime      (only available for lineaging dataset) birth time information of individual cells

%% Name List of longitudinal functional imaging data
% functional imaging data
fileDirNames{1} = '../../Data/func_20150410';
fileDirNames{2} = '../../Data/func_20150417';
fileDirNames{3} = '../../Data/func_20161004';

% functional imaging data with IHC information
fileDirNames{4} = '../../Data/islet_20170202';
fileDirNames{5} = '../../Data/islet_20170216';

% functional imaging data with lineage information
fileDirNames{6} = '../../Data/lineage_20160328';
fileDirNames{7} = '../../Data/lineage_20170111';
fileDirNames{8} = '../../Data/lineage_20170925';


%% Name List of ablation data
% single ablation at different locations, before and after
fileDirNames{9} = '../../AblationData/single_cutA_before';
fileDirNames{10} = '../../AblationData/single_cutA_after';
fileDirNames{11} = '../../AblationData/single_cutM_before';
fileDirNames{12} = '../../AblationData/single_cutM_after';
fileDirNames{13} = '../../AblationData/single_cutP_before';
fileDirNames{14} = '../../AblationData/single_cutP_after';

% double ablation, before and after
fileDirNames{15} = '../../AblationData/double_fish_1_before';
fileDirNames{16} = '../../AblationData/double_fish_1_after';
fileDirNames{17} = '../../AblationData/double_fish_2_before';
fileDirNames{18} = '../../AblationData/double_fish_2_after';
fileDirNames{19} = '../../AblationData/double_fish_3_before';
fileDirNames{20} = '../../AblationData/double_fish_3_after';

