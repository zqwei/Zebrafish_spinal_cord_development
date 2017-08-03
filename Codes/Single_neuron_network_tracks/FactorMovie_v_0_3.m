%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Factor movie
% 
% based on code -- 
% FactorMovie_v_0_2 no movie version
% 
% -------------------------------------------------------------------------
% 
% Ziqiang Wei
% weiz@janelia.hhmi.org
%
%


function FactorMovie_v_0_3(nFile)            
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile};   %#ok<*USENS>
    load([tempDatDir, fileName, '.mat'], 'dff', 'sideSplitter', 'side', 'tracks', 'timePoints', 'new_x', 'new_y', 'new_z'); 
    load([tempDatNetDir, 'LONOLoading_' fileName, '_v_0_1.mat'], 'new_activeNeuronMat')
    load([tempDatNetDir, 'LONOLoading_' fileName, '_v_0_2.mat'], 'preLMat', 'preLMatIndex', 'preLMatTime');
    
    if ~exist('new_x', 'var'); return; end
    
    x                 = new_x;
    y                 = new_y;
    z                 = new_z;
    numTime           = size(new_activeNeuronMat, 2);
    numNeuron         = length(side);
    halfActTime       = nan(numNeuron, 1);
    timeBin           = 15;
    activeThres       = 0.65;
    
    [~, neworder]     = sort(x);
    neworder = [neworder(side(neworder)==1); neworder(side(neworder)==2)];

    mColor = cbrewer('qual', 'Dark2',  max(preLMatIndex), 'cubic');
    linew = 1.25;
    z = z/max(z) * 1.8;
    y = y/2;

    figure;
    hold on

    for period = 1:numel(timePoints)
        timeRange = timePoints(period)+1:timePoints(period)+1200;
        radius = 0.05 + period*0.0005;
        LMat        = preLMat(:, preLMatTime == period);
        factorIndex = preLMatIndex(preLMatTime == period);
        activeTag   = new_activeNeuronMat(:, period);
        hold on
        
        if size(LMat,2) >= 1
            for nFactor = 1:size(LMat, 2)
                neuronFactor = LMat(:, nFactor)>0;
                dominateSide = ceil(mean(side(neuronFactor)) - 0.5);
                otherSide    = 3 - dominateSide;
                dominateSideNeuron = neuronFactor & side == dominateSide;
                otherSideNeuron    = neuronFactor & side == otherSide;
                if sum(dominateSideNeuron) < 4
                    CHPoints = smoothedBoundary(x(dominateSideNeuron), y(dominateSideNeuron), radius);
                    plot(CHPoints(:,1), CHPoints(:,2), '-', 'color', mColor(factorIndex(nFactor), :));
                    if sum(otherSideNeuron)>0
                        if sum(otherSideNeuron) > 0
                            CHPoints = smoothedBoundary(x(otherSideNeuron), y(otherSideNeuron), radius);
                            plot(CHPoints(:,1), CHPoints(:,2), '-', 'color', mColor(factorIndex(nFactor), :));
                        end
                    end
                end
            end
        end
    end
    plot(x, y, 'ok', 'MarkerFaceColor','k');
    colormap(mColor)
    colorbar
    xlim([0 ceil(max(x))+1]);
    ylim([-1 1]);
    gridxy(1:ceil(max(x)), 0, 'color', 'k', 'linestyle', '--')   
    box off
    setPrint(20, 12, [plotNetDir 'FactorEvolutionSpaceEarlyTime_' fileName], 'pdf')
end