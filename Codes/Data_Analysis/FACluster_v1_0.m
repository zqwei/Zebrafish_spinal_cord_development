%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.  Covariance analysis: Factor analysis -- loadin matrix evolution --
% reordering
% 
% size of factor against time
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function FACluster_v1_0(nFile)            
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>    
    mnx = [];
    load([tempDatDir, fileName, '.mat'], 'sideSplitter', 'side', 'tracks', 'timePoints', 'activeNeuronMat','mnx'); 
    load([tempDatDir, 'LONOLoading_' fileName, '.mat'], 'CorrectedLMat'); 
    numTime           = length(CorrectedLMat); %#ok<USENS>
    numNeuron         = length(side);
    radius            = 10;
    
    xlimMin           = min(min(tracks(:,:,1))); %#ok<NODEF>
    xlimMax           = max(max(tracks(:,:,1)));
    ylimMin           = min(min(tracks(:,:,2)));
    ylimMax           = max(max(tracks(:,:,2)));    
    
    mColor            = [     0    0.4470    0.7410
                         0.8500    0.3250    0.0980];
                     
    mColor            = [mColor; cbrewer('qual', 'Dark2',  128, 'cubic')];
    figure;
    
    video             = VideoWriter([plotDir, 'FactorEvolutionV3_', fileName, '.mp4'],'MPEG-4');%
    open(video)
    
    preLMat           = nan(numNeuron, 1);
    
    for nTime         = 1:numTime
        LMat          = CorrectedLMat{nTime}; % threshold by 0.3
        LMat(isnan(LMat)) = 0;  
        LMat(:, sum(LMat, 1)==0) = [];
        cla reset
        hold on
        xtracks   = squeeze(mean(tracks(:, timePoints(nTime)+(1:1200), 1), 2)); 
        ytracks   = squeeze(mean(tracks(:, timePoints(nTime)+(1:1200), 2), 2));            
        plot(xtracks(activeNeuronMat(:, nTime)), ytracks(activeNeuronMat(:, nTime)), 'ok', 'MarkerFaceColor','k') %#ok<NODEF>
        plot(xtracks(~activeNeuronMat(:, nTime)), ytracks(~activeNeuronMat(:, nTime)), 'ok')
%         plot(xtracks(mnx==0), ytracks(mnx==0), 'or')

        if size(LMat,2) > 1
            if sum(~isnan(preLMat(:))) == 0
                factorIndex  = 1:size(LMat, 2);
                preLMat      = LMat;
            else
%                 sideRemoveList = [];
%                 for nFactor = 1:size(LMat, 2)
%                     neuronFactor = LMat(:, nFactor)>0; 
%                     if length(unique(side(neuronFactor)))>1
%                         sideRemoveList   = [sideRemoveList; nFactor]; %#ok<AGROW>
%                     end
%                 end

                [~, maxFactorPerNeuronIndex] = max(LMat(sum(LMat, 2)>0, :), [], 2);
                sideRemoveList  = histc(maxFactorPerNeuronIndex, 1:size(LMat, 2)) <2; % remove the factor has no dominate factors
                
                LMat(:, sideRemoveList) = [];
                LMat            = LMat > 0;               
                sizeLMat        = sum(LMat, 1);
                [~, indexLMat]  = sort(sizeLMat, 'descend');
                LMat            = LMat(:, indexLMat);  
                factorIndex     = zeros(size(LMat, 2), 1);                

                % compute similarity matrix
                similarityScore = zeros(size(LMat, 2), size(preLMat, 2));
                for nFactor     = 1:size(LMat, 2)
                    if sum(isnan(LMat(:)))>0; keyboard();end
                    similarityScore(nFactor, :) = sum(bsxfun(@and, LMat(:, nFactor), preLMat));                    
                end
                
                % check if any prefactor has no connection with new factors
                % decide which factor is not included in preLMatIndex
                % maxIndex is the factor with the maximum coverage with the prefactors 
                [~, maxIndex]   = max(similarityScore, [], 1);
                                
                % check if any prefactor is merged (factor has maximum coverages with more than one prefactors, pick the larger one as its index)
                for nFactor     = 1:size(LMat, 2)
                    nFacotrNumPreFactor = sum(maxIndex == nFactor);
                    switch nFacotrNumPreFactor
                        case 0
                            preLMat = [preLMat, LMat(:, nFactor)]; %#ok<AGROW>
                            factorIndex(nFactor) = size(preLMat, 2);
                        case 1
                            if similarityScore(nFactor, maxIndex == nFactor) == 0
                                preLMat = [preLMat, LMat(:, nFactor)]; %#ok<AGROW>
                                factorIndex(nFactor) = size(preLMat, 2);
                            else
                                preLMatIndex         = find(maxIndex == nFactor);
                                factorIndex(nFactor) = preLMatIndex;
                                preLMat(:, preLMatIndex) = preLMat(:, preLMatIndex) | LMat(:, nFactor);
                            end
                        otherwise
                            preLMatIndex         = find(maxIndex == nFactor);
                            [~, nFactorMaxIndex] = max(similarityScore(nFactor, preLMatIndex));
                            factorIndex(nFactor) = preLMatIndex(nFactorMaxIndex);
                            preLMat(:, preLMatIndex(nFactorMaxIndex)) = preLMat(:, preLMatIndex(nFactorMaxIndex)) | LMat(:, nFactor);
                    end
                end
            end
                            
            for nFactor = 1:size(LMat, 2)
                neuronFactor = LMat(:, nFactor)>0;
                if length(unique(side(neuronFactor)))==1
                    CHPoints = smoothedBoundary(xtracks(neuronFactor), ytracks(neuronFactor), radius);
                    patch(CHPoints(:,1), CHPoints(:,2), mColor(factorIndex(nFactor), :), 'facealpha', 0.6, 'edgecolor', 'none');
                else
                    if sum(side(neuronFactor)==1) == 1 || sum(side(neuronFactor)==2) == 1 
                        plot(xtracks(neuronFactor), ytracks(neuronFactor), 'o', 'Color', mColor(factorIndex(nFactor), :), 'MarkerFaceColor', mColor(factorIndex(nFactor), :))
                    else
                        plot(xtracks(neuronFactor), ytracks(neuronFactor), 's', 'Color', mColor(factorIndex(nFactor), :), 'MarkerFaceColor', mColor(factorIndex(nFactor), :), 'MarkerSize', 10)
                    end
                end
            end       
        end
        text((xlimMin+xlimMax)/2, ylimMax-10, [num2str(nTime) ' min'],'fontsize', 24)
        xlim([xlimMin-5 xlimMax+5])
        ylim([ylimMin-5 ylimMax+5])
        axis off
        hold off;
        set(gcf,'color','w');
        pause(0.1);
        frame        = getframe;
        writeVideo(video, frame);
    end
    
    close(video)
    close all
    
end