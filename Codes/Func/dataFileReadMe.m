%%
% *dff.mat* is a matrix containing 94 rows (*cells*) and 63768 columns
% (*time points*) of the dF/F information. The temporal sampling is 4 Hz. 
% *No guarantee* the identity_ of every single cell track >30 min.
% 
% *clustering.mat* is just a vector of tags that came out of our clustering 
% analysis to mark out the cells that finally evolved the patterns.
% 
% * 0 indicating a non-patterned neuron
% * 1 a left-patterned neuron
% * 2 a right-patterned neuron

%% Name List of data files
fileDirNames{1} = '../../Data/20140818'; %no atlas info
fileDirNames{2} = '../../Data/20141006'; %no atlas info
fileDirNames{3} = '../../Data/20150410';
fileDirNames{4} = '../../Data/20150417';
fileDirNames{5} = '../../Data/20160312';
fileDirNames{6} = '../../Data/20160328';
fileDirNames{7} = '../../Data/20161004';
fileDirNames{8} = '../../Data/20161026'; %no xyz data yet
fileDirNames{9} = '../../Data/20161027'; %no xyz data yet
fileDirNames{10} = '../../Data/20161028';
fileDirNames{11} = '../../Data/20170111';
fileDirNames{12} = '../../Data/20170112';
fileDirNames{13} = '../../Data/20170126';
fileDirNames{14} = '../../Data/20170201';
fileDirNames{15} = '../../Data/20170202';
fileDirNames{16} = '../../Data/20170216';

% data with MO_slc6a9 injection
fileDirNames{17} = '../../Data/20170315';
fileDirNames{18} = '../../Data/20170323';
fileDirNames{19} = '../../Data/20170412';
fileDirNames{20} = '../../Data/20170503';
fileDirNames{21} = '../../Data/20170517';
fileDirNames{22} = '../../Data/20170731';
fileDirNames{23} = '../../Data/20170828';
fileDirNames{24} = '../../Data/20170925';
fileDirNames{97} = '../../Data/20180420';

for nFish = 1:36
    fileDirNames{23 + nFish * 2} = ['../../AblationData/fish_' num2str(nFish) '_before'];
    fileDirNames{24 + nFish * 2} = ['../../AblationData/fish_' num2str(nFish) '_after'];
end

