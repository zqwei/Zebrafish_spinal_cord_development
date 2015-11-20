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


function FACluster_v0_6(nFile)            
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>    
    load([tempDatDir, fileName, '.mat'], 'sideSplitter', 'side', 'tracks', 'timePoints', 'activeNeuronMat'); 
    load([tempDatDir, 'LONOLoading_' fileName, '.mat'], 'CorrectedLMat') 
%     CorrectedLMat     = LMat;
    numTime           = length(CorrectedLMat);
    numNeuron         = length(side);
    
    xlimMin           = min(min(tracks(:,:,1))); %#ok<NODEF>
    xlimMax           = max(max(tracks(:,:,1)));
    ylimMin           = min(min(tracks(:,:,2)));
    ylimMax           = max(max(tracks(:,:,2)));    
    
    mColor            = [     0    0.4470    0.7410
                         0.8500    0.3250    0.0980];
                     
    mColor            = [mColor; cbrewer('qual', 'Dark2',  16, 'cubic')];
    figure;
    
    mkdir([plotDir, 'Video_', fileName])
    
    for nTime         = 1:numTime
        LMat          = CorrectedLMat{nTime};
        LMat(isnan(LMat)) = 0;  
        LMat(:, sum(LMat, 1)==0) = [];
        
        cla reset
        hold on
        xtracks   = squeeze(mean(tracks(:, timePoints(nTime)+(1:1200), 1), 2)); 
        ytracks   = squeeze(mean(tracks(:, timePoints(nTime)+(1:1200), 2), 2));            
        plot(xtracks(activeNeuronMat(:, nTime)), ytracks(activeNeuronMat(:, nTime)), 'ok', 'MarkerFaceColor','k') %#ok<NODEF>
        plot(xtracks(~activeNeuronMat(:, nTime)), ytracks(~activeNeuronMat(:, nTime)), 'ok')    
        text((xlimMin+xlimMax)/2, ylimMax-10, [num2str(nTime) ' min'],'fontsize', 24)
        
        for nFactor = 1:size(LMat, 2)
            if sum(LMat(:,nFactor)~=0)>1
                cX      = xtracks'*abs(LMat(:,nFactor))/sum(abs(LMat(:,nFactor)));
                cY      = ytracks'*abs(LMat(:,nFactor))/sum(abs(LMat(:,nFactor)));            
                for nNeuron = 1:numNeuron
                    if LMat(nNeuron, nFactor)>0
                        plot([cX xtracks(nNeuron)], [cY ytracks(nNeuron)], '-', 'linewid', 4*LMat(nNeuron, nFactor), 'color', mColor(nFactor, :));
                    elseif LMat(nNeuron, nFactor)<0
                        disp(LMat(nNeuron, nFactor))
                        plot([cX xtracks(nNeuron)], [cY ytracks(nNeuron)], '--', 'linewid', -4*LMat(nNeuron, nFactor), 'color', mColor(nFactor, :));
                    end
                end
                
                plot(cX, cY, '^', 'color', mColor(nFactor, :), 'MarkerFaceColor', mColor(nFactor, :));

%                 if sum(side==1 & LMat(:,nFactor)~=0)>1 && sum(side==2 & LMat(:,nFactor)~=0)>1
%                     disp(nTime)
%                     plot(cX, cY, 's', 'color', mColor(nFactor, :), 'MarkerFaceColor', mColor(nFactor, :));
%                 else
%                     plot(cX, cY, '^', 'color', mColor(nFactor, :), 'MarkerFaceColor', mColor(nFactor, :));
%                 end
            end
        end
        

        xlim([xlimMin-5 xlimMax+5])
        ylim([ylimMin-5 ylimMax+5])
        axis off
        hold off;
        
        setPrint(16, 12, [plotDir, 'Video_', fileName, '/' num2str(nTime, '%03d')], 'pdf');
    end
    close all
end