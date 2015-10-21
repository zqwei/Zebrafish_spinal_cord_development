% 
% Selecting neurons and storing their raw data, baseline, dff, tracks
% 
% -------------------------------------------------------------------------
%
% version 0.0
% preprocessing of orginal data using sgolayfilt
%
%
% --- Deleting the time points has unreasonbale colinearity
% --- Deleting the units without oscillatory activity (KS test)
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function Neuron_selection_v0(nFile)
    % load data
    addpath('../Func');
    setDir;
    
    fileDirName   = fileDirNames{nFile}; %#ok<USENS>
    fileName      = fileNames{nFile}; %#ok<USENS>
    
    lenWindow     = 511;
    numOrder      = 9;

    dirImageData  = [fileDirName '/'];

    load ([dirImageData 'profile.mat']) % profile_all, timepoints, nCells, tracks_smoothed, side

    baseline      = sgolayfilt(profile_all, numOrder, lenWindow, [], 2); % nDim = 2: neuron
    rawf          = profile_all;
    dff           = bsxfun(@rdivide,(profile_all - baseline), mean(baseline, 2));
    tracks        = tracks_smoothed;
    
    ksTestAllTime    = false(nCells, 1);
    
    for nNeuron      = 1:nCells
        dat          = dff(nNeuron, :);
        ksTestAllTime(nNeuron) = kstest2(-dat(dat<0), dat(dat>0), 'alpha', 0.01);
    end
    
    slicedIndex   = ksTestAllTime;
    timeStep      = 1200;
    numT          = size(dff, 2);
    timeEnd       = (floor(numT/timeStep) - 2)*timeStep;    

    baseline      = baseline(slicedIndex, timeStep:timeEnd+timeStep);
    rawf          = rawf(slicedIndex, timeStep:timeEnd+timeStep);
    dff           = dff(slicedIndex, timeStep:timeEnd+timeStep);
    tracks        = tracks(slicedIndex, timeStep:timeEnd+timeStep);
    side          = side(slicedIndex); %#ok<NODEF>
    
        
    timePoints    = 0:240:timeEnd-timeStep;     %#ok<NASGU>
%     logDetCorrAllDFF = zeros(length(timePoints),1);    
%     for nTime           = 1:length(timePoints)
%         slicedDFF       = dff(:, timePoints(nTime)+1:timePoints(nTime)+1200);
%         logDetCorrAllDFF(nTime) = log(det(corr(slicedDFF')));
%     end
%     figure;
%     plot(timePoints/4/60, logDetCorrAllDFF, '-')
%     xlim([0 timePoints(end)/4/60]);
%     xlabel('Time (min)')
%     ylabel('log|Corr|')
%     title('Colinearity index')
%     box off
%     setPrint(8, 6, [plotDir 'LogDetCorrIndex_' fileName], 'pdf')
    
    
    slicedDFF       = dff(:, end - timeStep*10:end);  
    distNeurons     = pdist(slicedDFF, 'correlation');
    linkNeurons     = linkage(slicedDFF,'single','correlation');
    leafOrder       = optimalleaforder(linkNeurons, distNeurons);   
    
%     m               = ceil(sqrt(length(timePoints)));

%     figure;
%     for nPlot           = 1:length(timePoints)
%         subplot(m, m, nPlot)
%         slicedDFF       = dff(:, timePoints(nPlot)+1:timePoints(nPlot)+1200);        
%         imagesc(corr(slicedDFF(leafOrder,:)'),[-1 1]);
%         xlabel('Neuronal index')
%         ylabel('Neuronal index')
%         title([num2str(nPlot) 'min']);
%         box off;
%     end
%     setPrint(m*8, m*6, [plotDir, 'CorrMat_', fileName], 'pdf');


%     figure;    
%     for nPlot           = 1:length(timePoints)
%         subplot(m, m, nPlot)
%         slicedDFF       = dff(:,timePoints(nPlot)+1:timePoints(nPlot)+1200);  
%         distNeurons     = pdist(slicedDFF, 'correlation');
%         linkNeurons     = linkage(slicedDFF,'single','correlation');
%         leafOrderLocal  = optimalleaforder(linkNeurons, distNeurons);       
%         imagesc(corr(slicedDFF(leafOrderLocal,:)'),[-1 1]);
%         xlabel('Neuronal index')
%         ylabel('Neuronal index')
%         title([num2str(nPlot) 'min']);
%         box off;
%     end
%     setPrint(m*8, m*6, [plotDir, 'CorrMatLocal_', fileName], 'pdf');
        
    baseline      = baseline(leafOrder, :); %#ok<NASGU>
    rawf          = rawf(leafOrder, :); %#ok<NASGU>
    dff           = dff(leafOrder, :); %#ok<NASGU>
    tracks        = tracks(leafOrder, :); %#ok<NASGU>
    side          = side(leafOrder, :); %#ok<NASGU>

    save([tempDatDir, fileName, '.mat'], 'dff', 'tracks', 'leafOrder', 'slicedIndex','baseline', 'rawf', 'side', 'timePoints');
    if exist('mnx', 'var')
        mnx       = mnx(slicedIndex);
        mnx       = mnx(leafOrder);
        save([tempDatDir, fileName, '.mat'], 'mnx', '-append')
    end
end