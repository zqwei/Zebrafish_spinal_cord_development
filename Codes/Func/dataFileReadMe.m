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
% fileDirNames{1} = '../../Data/20140818';
% fileDirNames{2} = '../../Data/20141006';
% fileDirNames{3} = '../../Data/20150410';
fileDirNames{1} = '../../Data/20140818_curated';
fileDirNames{2} = '../../Data/20141006_curated';
fileDirNames{3} = '../../Data/20150410_curated';
fileDirNames{4} = '../../Data/20150417';
fileDirNames{5} = '../../Data/20160312';
fileDirNames{6} = '../../Data/20160328';
fileDirNames{7} = '../../Data/20161004';
fileDirNames{8} = '../../Data/20161028';
fileDirNames{9} = '../../Data/20170111';
fileDirNames{10} = '../../Data/20170112';