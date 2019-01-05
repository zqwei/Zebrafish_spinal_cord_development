%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.  Covariance analysis: Factor analysis -- loadin matrix evolution --
% reordering
% 
% center of the factor
% connectivity strength to each node
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function FACluster_v0_8(nFile, nTime) 
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>    
    load([tempDatDir, fileName, '.mat'], 'sideSplitter', 'side', 'tracks', 'timePoints', 'activeNeuronMat', 'dff'); 
    load([tempDatDir, 'EvoLoading_' fileName, '.mat'], 'networkMat', 'neuronXLoc', 'neuronYLoc', 'neuronZLoc'); 
    xtracks           = neuronXLoc(:, nTime);
    ytracks           = neuronYLoc(:, nTime);
    ztracks           = neuronZLoc(:, nTime);
    factorSet         = networkMat{nTime};
    numNeuron         = length(xtracks);
    
    slicedDFF         = dff(:,timePoints(nTime)+1:timePoints(nTime)+1200); %#ok<NODEF>
    slicedDFF         = bsxfun(@minus, slicedDFF, mean(slicedDFF,2));
    slicedDFF         = bsxfun(@rdivide, slicedDFF, std(slicedDFF,[],2))';
    
    mColor            = cbrewer('qual', 'Dark2',  16, 'cubic');
    numDelay          = 200;
    maxDelay          = 25;
    dColor            = cbrewer('div', 'Spectral',  numDelay, 'cubic');
    
    figure;
    subplot(1, 2, 1)
    hold on;
    for nFactor              = 1:length(factorSet)
        if nFactor == 1
            neuronIndex = factorSet{nFactor}.neuronIndex;
            x           = xtracks(neuronIndex);
            y           = ytracks(neuronIndex);
            z           = ztracks(neuronIndex);
            plot3(x, y, z, 'o', 'color', [0.3 0.3 0.3]);
        elseif nFactor<length(factorSet)
            neuronIndex = factorSet{nFactor}.neuronIndex;
            if length(neuronIndex) == 1 % single neuron FA
                x      = factorSet{nFactor}.x;
                y      = factorSet{nFactor}.y;
                z      = factorSet{nFactor}.z;
                plot3(x, y, z, 'o', 'color', mColor(nFactor, :), 'MarkerFaceColor', mColor(nFactor, :), 'MarkerSize', 10);
            else % multi-neuron FA
                x      = factorSet{nFactor}.x;
                y      = factorSet{nFactor}.y;
                z      = factorSet{nFactor}.z;
                plot3(x, y, z, '^', 'color', mColor(nFactor, :), 'MarkerFaceColor', mColor(nFactor, :), 'MarkerSize', 10);
                plot3(xtracks(neuronIndex), ytracks(neuronIndex), ztracks(neuronIndex), 'ok', 'MarkerFaceColor', 'k');
                neuronCCMat = factorSet{nFactor}.neuronCCMat;
                for nNeuron = 1:numNeuron
                    if ~isnan(neuronCCMat(nNeuron, 1)) && neuronCCMat(nNeuron, 2)>0
%                         colorIndex = ceil(neuronCCMat(nNeuron, 1)/10*100)+100;
                        colorIndex = ceil(abs(neuronCCMat(nNeuron, 1))/maxDelay*(numDelay-1))+1;
                        xs     = [x, xtracks(nNeuron)];
                        ys     = [y, ytracks(nNeuron)];
                        zs     = [z, ztracks(nNeuron)];
                        plot3(xs, ys, zs, '-', 'linewid', neuronCCMat(nNeuron, 2)*4, 'color', dColor(colorIndex,:))
%                         drawArrow(xs, ys, zs, {'-', 'linewid', neuronCCMat(nNeuron, 2)*4, 'color', dColor(colorIndex,:)})
                    end
                end
            end
        else % FA connection
%             for pFactor          = 1:length(factorSet)-2
%                 for qFactor      = pFactor+1:length(factorSet)-2
%                     if ~isnan(factorSet{nFactor}.delayMat(pFactor, qFactor)) && factorSet{nFactor}.corrMat(pFactor, qFactor)>0
% %                         colorIndex = ceil(neuronCCMat(nNeuron, 1)/10*100)+100;
%                         colorIndex = ceil(abs(neuronCCMat(nNeuron, 1))/10*199)+1;
%                         xs         = [factorSet{pFactor+1}.x, factorSet{qFactor+1}.x];
%                         ys         = [factorSet{pFactor+1}.y, factorSet{qFactor+1}.y];
%                         zs         = [factorSet{pFactor+1}.z, factorSet{qFactor+1}.z];
%                         plot3(xs, ys, zs, '-', 'linewid', neuronCCMat(nNeuron, 2)*4, 'color', dColor(colorIndex,:))
% %                         drawArrow(xs, ys, zs, {'-', 'linewid', factorSet{nFactor}.corrMat(pFactor, qFactor)*4, 'color', dColor(colorIndex,:)})
%                     end
%                 end
%             end
        end
    end
    
    colormap(dColor)
    hcb = colorbar;
    set(hcb, 'YTick', 0:0.2:1, 'YTickLabel', num2str((0:0.2:1)'*maxDelay))
    ylabel(hcb, 'Delay (sec)');
    xlabel('x-axis')
    ylabel('y-axis')
    zlabel('z-axis')
    box off
    
    subplot(1, 2, 2)
    hold on;
    for nFactor              = 1:length(factorSet)
        if nFactor == 1
            neuronIndex = factorSet{nFactor}.neuronIndex;
            x           = xtracks(neuronIndex);
            y           = ytracks(neuronIndex);
            z           = ztracks(neuronIndex);
            plot3(x, y, z, 'o', 'color', [0.3 0.3 0.3]);
        elseif nFactor<length(factorSet)
            neuronIndex = factorSet{nFactor}.neuronIndex;
            if length(neuronIndex) == 1 % single neuron FA
                x      = factorSet{nFactor}.x;
                y      = factorSet{nFactor}.y;
                z      = factorSet{nFactor}.z;
                plot3(x, y, z, 'o', 'color', mColor(nFactor, :), 'MarkerFaceColor', mColor(nFactor, :), 'MarkerSize', 10);
            else % multi-neuron FA
                x      = factorSet{nFactor}.x;
                y      = factorSet{nFactor}.y;
                z      = factorSet{nFactor}.z;
                plot3(x, y, z, '^', 'color', mColor(nFactor, :), 'MarkerFaceColor', mColor(nFactor, :), 'MarkerSize', 10);
                plot3(xtracks(neuronIndex), ytracks(neuronIndex), ztracks(neuronIndex), 'ok', 'MarkerFaceColor', 'k');
%                 neuronCCMat = factorSet{nFactor}.neuronCCMat;
%                 for nNeuron = 1:numNeuron
%                     if ~isnan(neuronCCMat(nNeuron, 1)) && neuronCCMat(nNeuron, 2)>0
% %                         colorIndex = ceil(neuronCCMat(nNeuron, 1)/10*100)+100;
%                         colorIndex = ceil(abs(neuronCCMat(nNeuron, 1))/10*199)+1;
%                         xs     = [x, xtracks(nNeuron)];
%                         ys     = [y, ytracks(nNeuron)];
%                         zs     = [z, ztracks(nNeuron)];
%                         plot3(xs, ys, zs, '-', 'linewid', neuronCCMat(nNeuron, 2)*4, 'color', dColor(colorIndex,:))
% %                         drawArrow(xs, ys, zs, {'-', 'linewid', neuronCCMat(nNeuron, 2)*4, 'color', dColor(colorIndex,:)})
%                     end
%                 end
            end
        else % FA connection
            for pFactor          = 1:length(factorSet)-2
                for qFactor      = pFactor+1:length(factorSet)-2
                    if ~isnan(factorSet{nFactor}.delayMat(pFactor, qFactor)) && factorSet{nFactor}.corrMat(pFactor, qFactor)>0
%                         colorIndex = ceil(neuronCCMat(nNeuron, 1)/10*100)+100;
                        colorIndex = ceil(abs(factorSet{nFactor}.delayMat(pFactor, qFactor))/maxDelay*(numDelay-1))+1;
                        xs         = [factorSet{pFactor+1}.x, factorSet{qFactor+1}.x];
                        ys         = [factorSet{pFactor+1}.y, factorSet{qFactor+1}.y];
                        zs         = [factorSet{pFactor+1}.z, factorSet{qFactor+1}.z];
                        plot3(xs, ys, zs, '-', 'linewid', factorSet{nFactor}.corrMat(pFactor, qFactor)*4, 'color', dColor(colorIndex,:))
%                         drawArrow(xs, ys, zs, {'-', 'linewid', factorSet{nFactor}.corrMat(pFactor, qFactor)*4, 'color', dColor(colorIndex,:)})
                    end
                end
            end
        end
    end
    
    colormap(dColor)
    hcb = colorbar;
    set(hcb, 'YTick', 0:0.2:1, 'YTickLabel', num2str((0:0.2:1)'*maxDelay))
    ylabel(hcb, 'Delay (sec)');
    xlabel('x-axis')
    ylabel('y-axis')
    zlabel('z-axis')
    box off
    setPrint(32, 12, [plotDir, 'CorrSpace_', fileName, '_' num2str(nTime, '%03d')], 'pdf');

end