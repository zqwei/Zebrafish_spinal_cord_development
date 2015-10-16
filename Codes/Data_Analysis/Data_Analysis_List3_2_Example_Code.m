addpath('../Func');
nFile               = 1;
setDir;    
fileDirName         = fileDirNames{nFile};
fileName            = fileNames{nFile};
load(['../../Data/' fileName, '/dff.mat'],'baseline')
load([tempDatDir, fileName, '.mat']);

baseline            = baseline(slicedIndex, :);
baseline            = baseline(leafOrder, :);


corrThres           = 0.4; 


nPlot           = 1; 
figure;
slicedDFF       = dff(:,timePoints(nPlot)+1:timePoints(nPlot+1));
corrDFF         = corr(slicedDFF'); 
% off-diagonal correlation matrix
corrDFF         = corrDFF - eye(size(corrDFF));        
% find node has strong correlations
maxCorr         = max(corrDFF, [], 2);
numNeuronPlot   = min(20, sum(maxCorr>corrThres));
[~, neuronIndex]= sort(maxCorr, 'descend');
xTrack          = squeeze(mean(tracks(:, timePoints(nPlot)+1:timePoints(nPlot+1), 1), 2));
yTrack          = squeeze(mean(tracks(:, timePoints(nPlot)+1:timePoints(nPlot+1), 2), 2)); 
zTrack          = squeeze(mean(tracks(:, timePoints(nPlot)+1:timePoints(nPlot+1), 3), 2)); 

distNeurons     = pdist(slicedDFF(neuronIndex(1:numNeuronPlot),:), 'correlation');
linkNeurons     = linkage(slicedDFF(neuronIndex(1:numNeuronPlot),:),'complete','correlation');
leafOrder       = optimalleaforder(linkNeurons, distNeurons);
neuronIndex(1:numNeuronPlot) = neuronIndex(leafOrder);
lineColors      = [jet(numNeuronPlot); 0.5 0.5 0.5];

figure;
hold on;
plot3(xTrack(neuronIndex(numNeuronPlot+1:end)), yTrack(neuronIndex(numNeuronPlot+1:end)), zTrack(neuronIndex(numNeuronPlot+1:end)), ...
    'o', 'markerFaceColor', lineColors(end,:), ...
    'markerEdgeColor', lineColors(end,:));
for nPoint      = 1:numNeuronPlot
    plot3(xTrack(neuronIndex(nPoint)), yTrack(neuronIndex(nPoint)), zTrack(neuronIndex(nPoint)), ...
        'o', 'markerFaceColor', lineColors(nPoint,:), ...
        'markerEdgeColor', lineColors(nPoint,:));
end
hold off;
box off
axis([0 1600 0 400 0 50]);



figure;
plot((timePoints(nPlot)+1:timePoints(nPlot+1)), ...
    slicedDFF(neuronIndex(3:6),:));
xlim([timePoints(nPlot)+1 timePoints(nPlot+1)]);
xlabel('Time (hr)');
ylabel('df/f');
box off


figure;
plot((timePoints(nPlot)+1:timePoints(nPlot+2))/3600/4, ...
    baseline(neuronIndex(3:6),timePoints(nPlot)+1:timePoints(nPlot+2)));
xlim([timePoints(nPlot)+1 timePoints(nPlot+1)]/3600/4);
xlabel('Time (hr)');
ylabel('df/f');
box off


figure;
plot((timePoints(nPlot)+1:timePoints(nPlot+2))/3600/4, ...
    (dff(neuronIndex(3:6),timePoints(nPlot)+1:timePoints(nPlot+2))+1).*baseline(neuronIndex(3:6),timePoints(nPlot)+1:timePoints(nPlot+2)));
xlim([timePoints(nPlot)+1 timePoints(nPlot+2)]/3600/4);
xlabel('Time (hr)');
ylabel('df/f');
box off

nPlot = 30;
figure;
plot((timePoints(nPlot)+1:timePoints(nPlot+2))/3600/4, ...
    dff(neuronIndex(3:6),timePoints(nPlot)+1:timePoints(nPlot+2)));
xlim([timePoints(nPlot)+1 timePoints(nPlot+2)]/3600/4);
xlabel('Time (hr)');
ylabel('df/f');
box off
