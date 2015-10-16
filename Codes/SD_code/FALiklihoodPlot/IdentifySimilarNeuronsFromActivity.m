% Assumes dff variable is loaded
neuronNum = size(dff,1);
timePointNum = size(dff,2);
Z = linkage(dff,'single','correlation');
maximalAllowedCorrelation = 0.9;
neuronCorrelationGrouping = cluster(Z,'cutoff',1-maximalAllowedCorrelation);

newNeuronNum = length(unique(neuronCorrelationGrouping));
neuronInGroupInd = cell(newNeuronNum,1);
neuronInGroupNum = zeros(newNeuronNum,1);
newDff = zeros(newNeuronNum,timePointNum);

for ii=1:newNeuronNum
  neuronInGroupInd{ii} = find(neuronCorrelationGrouping==ii);
  neuronInGroupNum(ii) = length(neuronInGroupInd{ii});
  newDff(ii,:) = mean(dff(neuronCorrelationGrouping==ii,:),1);
end