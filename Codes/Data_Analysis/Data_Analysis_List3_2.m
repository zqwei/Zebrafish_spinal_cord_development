%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.  Covariance analysis: cross-middle line activity neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 3.2 plot the neuron's raw df/f data with strong activity across them
%     max_example neurons (20 neurons)
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
%

function Data_Analysis_List3_2(nFile)
    % load data
    addpath('../Func');
    setDir;    
    fileDirName         = fileDirNames{nFile}; %#ok<USENS,NASGU>
    fileName            = fileNames{nFile}; %#ok<USENS>
    load([tempDatDir, fileName, '.mat']);

    corrThres           = 0.4; 

    if ~exist([plotDir, 'Strong_correlated_neurons_', fileName], 'dir')
        mkdir([plotDir, 'Strong_correlated_neurons_', fileName]);
    end
    
    
    
    
    
    for nPlot           = 1:20 
        figure;
        slicedDFF       = dff(:,timePoints(nPlot)+1:timePoints(nPlot+1));   %#ok<NODEF>
        corrDFF         = corr(slicedDFF'); 
        % off-diagonal correlation matrix
        corrDFF         = corrDFF - eye(size(corrDFF));        
        % find node has strong correlations
        maxCorr         = max(corrDFF, [], 2);
        numNeuronPlot   = min(20, sum(maxCorr>corrThres));
        [~, neuronIndex]= sort(maxCorr, 'descend');
        xTrack          = squeeze(mean(tracks(:, timePoints(nPlot)+1:timePoints(nPlot+1), 1), 2)); %#ok<NODEF>
        yTrack          = squeeze(mean(tracks(:, timePoints(nPlot)+1:timePoints(nPlot+1), 2), 2));        

        distNeurons     = pdist(slicedDFF(neuronIndex(1:numNeuronPlot),:), 'correlation');
        linkNeurons     = linkage(slicedDFF(neuronIndex(1:numNeuronPlot),:),'complete','correlation');
        leafOrder       = optimalleaforder(linkNeurons, distNeurons);
        neuronIndex(1:numNeuronPlot) = neuronIndex(leafOrder);
        lineColors      = [jet(numNeuronPlot); 0.5 0.5 0.5];

        subplot(4, 7, [1 2 8 9])
%         imagesc(corr(slicedDFF(neuronIndex(1:numNeuronPlot),:)'),[-1 1]);
        imagesc(corr(slicedDFF(neuronIndex,:)'),[-1 1]);
        colormap
        xlabel('Neuronal index')
        ylabel('Neuronal index')
        title(['Time from ' num2str(timePoints(nPlot)/3600/4) ' to ' num2str(timePoints(nPlot+1)/3600/4) 'hr']);
        box off
        
        subplot(4, 7, [15 16 22 23])
        hold on;
        plot(xTrack(neuronIndex(numNeuronPlot+1:end)), yTrack(neuronIndex(numNeuronPlot+1:end)), ...
            'o', 'markerFaceColor', lineColors(end,:), ...
            'markerEdgeColor', lineColors(end,:));
        for nPoint      = 1:numNeuronPlot
            plot(xTrack(neuronIndex(nPoint)), yTrack(neuronIndex(nPoint)), ...
                'o', 'markerFaceColor', lineColors(nPoint,:), ...
                'markerEdgeColor', lineColors(nPoint,:));
        end
        hold off;
        box off
        axis([0 1600 0 400]);

        
        for nPoint      = 1:numNeuronPlot
            subplot(4, 7, mod(nPoint-1,5)+3 + (nPoint-mod(nPoint-1,5)-1)/5*7)
            plot((timePoints(nPlot)+1:timePoints(nPlot+1))/3600/4, ...
                slicedDFF(neuronIndex(nPoint),:),'color',lineColors(nPoint,:));
            xlim([timePoints(nPlot)+1 timePoints(nPlot+1)]/3600/4);
            xlabel('Time (hr)');
            ylabel('df/f');
            box off
        end                
        
        setPrint(8*5, 6*4, [plotDir, 'Strong_correlated_neurons_', fileName, '/Time_at_' num2str(timePoints(nPlot),'%07d') ],'pdf')
        close all
    end

end