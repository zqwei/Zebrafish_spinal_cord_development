%   Assumes useDFF is loaded

seed = 1;
rng(seed);

timePointNum = size(useDFF,1);
neuronNum = size(useDFF,2);

timePeriodLength = 5000;
windowShift = 5000;
timePointLim = 1:windowShift:(timePointNum-timePeriodLength-1);% Window style
timePeriodNum = length(timePointLim)-1;

finalTime = timePointLim(end-1):timePointLim(end);
Z = linkage(transpose(useDFF(finalTime,:)),'single','correlation');
[H,T,outperm] = dendrogram(Z,0);
useDFF = useDFF(:,outperm);

%% Clustering over time

cutoffCorrelation = 0.2;
minNeuronNumInCluster = 4;
clusterNumOverTime = zeros(1,length(timePointLim)-1);
neuronCorrelationGroupingMat = zeros(neuronNum,length(timePointLim)-1);
neuronsPerClusterMat = zeros(neuronNum,length(timePointLim)-1);
keepClustersInd = cell(length(timePointLim)-1,1);
keepClustersNeuronInd = cell(length(timePointLim)-1,1);

for ii=1:(length(timePointLim)-1)
  display(['Calculating models for time point index ' num2str(ii) ' of ' num2str(length(timePointLim)-1)])

  timeIndex = timePointLim(ii):(timePointLim(ii)+timePeriodLength-1);
  normalizedData = useDFF(timeIndex,:);
  normalizedData = normalizedData-repmat(mean(normalizedData),size(normalizedData,1),1);
  normalizedData = normalizedData./repmat(std(normalizedData),size(normalizedData,1),1);
  
  Z = linkage(normalizedData','single','correlation');
  
  neuronCorrelationGroupingMat(:,ii) = cluster(Z,'cutoff',1-cutoffCorrelation);
  clusterNumOverTime(ii) = length(unique(neuronCorrelationGroupingMat(:,ii)));
  
  for jj=1:clusterNumOverTime(ii)
    neuronsPerClusterMat(jj,ii) = sum(neuronCorrelationGroupingMat(:,ii)==jj);
  end
  [m,ix] = sort(neuronsPerClusterMat(:,ii),'descend');
  lastClustInd = find(m < minNeuronNumInCluster,1,'first');
  if lastClustInd>1; lastClustInd = lastClustInd-1; end
  keepClustersInd{ii} = ix(1:lastClustInd);
  
  keepClustersNeuronInd{ii} = cell(lastClustInd,1);
  for jj=1:lastClustInd
    keepClustersNeuronInd{ii}{jj} = find(neuronCorrelationGroupingMat(:,ii)==keepClustersInd{ii}(jj));
  end
end

fitCell = cell(length(timePointLim)-1,1);
overlapCell = cell(length(timePointLim)-1,1);

figure; hold on;

colorCell = {'r','b','g','c','y'};
baseArea = 10;
maxYSpan = 0.5;
colorVec = [0 0 0];
for ii=2:(length(timePointLim)-1)
  [overlapCell{ii},fitCell{ii}] = ...
    TempCompareTwoClusterings(keepClustersNeuronInd{ii},keepClustersNeuronInd{end});
  
  clusterNum = length(keepClustersNeuronInd{ii});
  for jj=1:clusterNum
  % Add text indicator of cluster
    colorVec = zeros(1,3)+(ones(1,3)*(jj-1)*(0.75/clusterNum));
    curUseNeuronInd = keepClustersNeuronInd{ii}{jj};
    scatter(curUseNeuronInd(1),ii+0.2,baseArea+80,colorVec,'filled','marker','pentagram');
    scatter(curUseNeuronInd,ones(length(curUseNeuronInd),1)*ii,baseArea+80,colorCell{fitCell{ii}(jj)},'filled');
  
    for kk=2:length(keepClustersNeuronInd{ii}{jj});
      DrawSimpleArc([keepClustersNeuronInd{ii}{jj}(kk-1) ii]',[keepClustersNeuronInd{ii}{jj}(kk) ii]',[(keepClustersNeuronInd{ii}{jj}(kk-1) + keepClustersNeuronInd{ii}{jj}(kk))/2 ii]',maxYSpan,colorVec);
    end
  end
  
%   for jj=1:length(keepClustersNeuronInd{end})
%     curUseNeuronInd = keepClustersNeuronInd{end}{jj};
%     scatter(curUseNeuronInd,ones(length(curUseNeuronInd),1)*ii,baseArea+40,colorCell{jj},'filled');
%   end
  
  useNeuronIndVec = 1:neuronNum;
  scatter(useNeuronIndVec,ones(length(useNeuronIndVec),1)*ii,baseArea,'k','filled');
  
  
end
% How about drawing arcs in between those that are of the same cluster,
% color coded according to their own cluster, perhaps just grey.

% This color within color visualization is no good.