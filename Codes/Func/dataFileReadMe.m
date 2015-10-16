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
fileDirNames{1} = '../../Data/20140818';
fileDirNames{2} = '../../Data/20141006';
fileDirNames{3} = '../../Data/20150410';