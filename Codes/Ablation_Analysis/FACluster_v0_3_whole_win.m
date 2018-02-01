%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.  Covariance analysis: Factor analysis -- number of island using non
% LONO methods with active neurons
%
%
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
%


function FACluster_v0_3_whole_win(nFile)
    addpath('../Func');
    setDir;
    fileName          = fileNames{nFile}; %#ok<USENS>
    load([tempDatDir, fileName, '.mat'], 'dff','timePoints','activeNeuronMat', 'timeStep', 'new_x', 'new_y', 'new_z'); 
    load([tempDatDir, 'FALONO_', fileName, '.mat'], 'uncorrectedLONOM', 'numFactors');

    LONOM         = mode(uncorrectedLONOM);
    activeNeuron  = sum(activeNeuronMat, 2) >0;
%     activeNeuron  = mean(activeNeuronMat, 2) >0.5;
    LMat          = nan(size(activeNeuronMat, 1), LONOM);
    PsiMat        = nan(size(activeNeuronMat, 1), 1);

%     slicedDFF     = dff(activeNeuron,timePoints(1)+1:timePoints(end)+timeStep);
    slicedDFF     = dff(activeNeuron,timePoints(1)+1:timePoints(end));
    slicedDFF     = bsxfun(@minus, slicedDFF, mean(slicedDFF,2));
    slicedDFF     = bsxfun(@rdivide, slicedDFF, std(slicedDFF,[],2))';
    [Lambda, Psi] = factoran(slicedDFF, LONOM,'maxit',10000);

    LMat(activeNeuron, :)   = Lambda;
    PsiMat(activeNeuron) = Psi;
    
    new_y = new_y/2;
    new_z = new_z/max(new_z) * 1.2;
    figure;
    LMat(LMat<0) = 0;
    LMat(isnan(LMat)) = 0;
    thresL       = 0.3;
    subplot(3, 2, [1 3 5])
    imagesc(LMat)
    subplot(3, 2, 2)
    hold on
    plot(new_x(~activeNeuron), new_y(~activeNeuron), 'sk')
    plot(new_x(activeNeuron), new_y(activeNeuron), 'ok')
    for nFactor = 1:LONOM
        scatter(new_x(LMat(:,nFactor)>thresL), new_y(LMat(:,nFactor)>thresL), 'filled')
    end
    subplot(3, 2, 4)
    hold on
    plot(new_x(~activeNeuron & new_y<0), new_z(~activeNeuron & new_y<0), 'sk')
    plot(new_x(activeNeuron & new_y<0), new_z(activeNeuron & new_y<0), 'ok')
    for nFactor = 1:LONOM
        scatter(new_x(LMat(:,nFactor)>thresL & new_y<0), new_z(LMat(:,nFactor)>thresL & new_y<0), 'filled')
    end
    subplot(3, 2, 6)
    hold on
    plot(new_x(~activeNeuron & new_y>0), new_z(~activeNeuron & new_y>0), 'sk')
    plot(new_x(activeNeuron & new_y>0), new_z(activeNeuron & new_y>0), 'ok')
    for nFactor = 1:LONOM
        scatter(new_x(LMat(:,nFactor)>thresL & new_y>0), new_z(LMat(:,nFactor)>thresL & new_y>0), 'filled')
    end
    
    setPrint(8*2, 6, [plotDir, 'AverageFactorAblation_', fileName], 'pdf')
end
