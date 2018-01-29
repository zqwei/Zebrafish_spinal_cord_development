%
% Compute the number of more neurons included in fit as add the one with
% low correlations
%
%

% for nFile = 1:24
%     Neuron_selection_v3_1(nFile, thresTwichCor(nFile))
% end


addpath('../Func');
setDir;

numMat   = [];
groupInd = [];

for nFile = 1:24
    fileName      = fileNames{nFile}; 
    load([tempDatDir, 'numActiveNewNeuron', fileName, '.mat'], 'numActiveNewNeuron');
    numMat   = [numMat, numActiveNewNeuron];
    groupInd = [groupInd, ones(size(numActiveNewNeuron))*nFile];
end

figure
boxplot(numMat,groupInd)
xlabel('Fish Index')
ylabel('Number of more active neurons')
