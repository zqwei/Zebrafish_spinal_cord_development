%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.  Covariance analysis: factor analysis
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.2 Contunity check of FA
%     See if the FA is changed hugely at the adjacent time points
%     Based on result 2.1
%     FA dimesion should be between 3 - 10
%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 

function Data_Analysis_List2_2(nFile)

    if nargin<1
        nFile = 1;
    end
        
    addpath('../Func');
    setDir;    
    fileDirName       = fileDirNames{nFile}; %#ok<NASGU,USENS>
    fileName          = fileNames{nFile}; %#ok<USENS>
    
    maxNumFactor      = 10; % This should be changed on some points...
    minNumFactor      = 2;
    
    load([tempDatDir, fileName, '.mat']);
        
    numPlot       = length(timePoints)-1;
    matEV         = nan(numPlot, numPlot, length(minNumFactor:maxNumFactor));
    numUnit       = size(dff,1); %#ok<NODEF>
    matEVSingleUnit = nan(numPlot, numPlot, length(minNumFactor:maxNumFactor), numUnit);

    for nPlot1       = numPlot:-1:1
        refDFF       = dff(:,timePoints(nPlot1)+1:timePoints(nPlot1+1));
        refDFF       = bsxfun(@minus, refDFF, mean(refDFF,2));
        refDFF       = bsxfun(@rdivide, refDFF, std(refDFF,[],2))';
        for nFactor         = minNumFactor:maxNumFactor
            % get lambda and psi from refDFF
            [lambda,psi]    = factoran(refDFF, nFactor, 'scores','regression', 'rotate', 'none');
            for nPlot2      = numPlot:-1:1
                % compute EV of slicedDFF based on lambda and psi from refDFF
                slicedDFF   = dff(:,timePoints(nPlot2)+1:timePoints(nPlot2+1));
                slicedDFF   = bsxfun(@minus, slicedDFF, mean(slicedDFF,2));
                slicedDFF   = bsxfun(@rdivide, slicedDFF, std(slicedDFF,[],2))';                
                matEV(nPlot1, nPlot2, nFactor-minNumFactor+1) = LONOFA(slicedDFF, lambda, psi);
                matEVSingleUnit(nPlot1, nPlot2, nFactor-minNumFactor+1, :) = LONOFASingleUnitEV(slicedDFF, lambda, psi);
                display([nPlot1, nPlot2, nFactor])
            end
        end
    end
    save([tempDatDir, fileName, '_FAEVCrossTime.mat'], 'matEV', 'matEVSingleUnit'); 
    figure;
    imagesc(squeeze(max(matEV,[],3)));
end
