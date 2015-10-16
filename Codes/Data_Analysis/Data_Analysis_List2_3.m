%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.  Covariance analysis: factor analysis
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.2 Contunity check of FA (Plots based on 2.2)
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

function Data_Analysis_List2_3(nFile)

    if nargin<1
        nFile = 1;
    end
        
    addpath('../Func');
    setDir;    
    fileDirName       = fileDirNames{nFile}; %#ok<NASGU,USENS>
    fileName          = fileNames{nFile}; %#ok<USENS>
    maxNumFactor      = 10;
    minNumFactor      = 2;
    load([tempDatDir, fileName, '.mat']);
    load([tempDatDir, fileName, '_FAEVCrossTime.mat'], 'matEV', 'matEVSingleUnit'); 
%     load([tempDatDir, fileName, '_FAEV.mat'], 'opt1Dim');
    load([tempDatDir, fileName, '_FAEV.mat'], 'opt2Dim');
    if exist('opt2Dim','var') && ~exist('opt1Dim','var'); opt1Dim = opt2Dim; end
    figure;
    subplot(2,1,1)
%     imagesc(squeeze(max(matEV,[],3)));
    plot_EV_Time_Time(timePoints, matEV)
    
%     figure;
%     plot_EVDim_Time_Time(timePoints, matEV, minNumFactor)
    
    opt1Dim = min(opt1Dim,maxNumFactor);
    subplot(2,1,2)
    plot_EVopt1Dim_Time_Time(timePoints, matEV, opt1Dim, minNumFactor)
end

function plot_EV_Time_Time(timePoints, matEV)
    imagesc(timePoints(1:end-1)/4/3600, timePoints(1:end-1)/4/3600, max(matEV,[],3)');
    xlim([0 timePoints(end)'/4/3600])
    ylim([0 timePoints(end)'/4/3600])
    box off
    title ('Max Total EV')
    xlabel('Time (hour)') % reference
    ylabel('Time (hour)') % test
    colorbar
    caxis([-0.3 0.7])
end


function plot_EVDim_Time_Time(timePoints, matEV, minNumFactor)
    [~, EVDim] = max(matEV,[],3);
    EVDim      = EVDim'  + minNumFactor -1;
    imagesc(timePoints(1:end-1)/4/3600, timePoints(1:end-1)/4/3600, EVDim);
    xlim([0 timePoints(end)'/4/3600])
    ylim([0 timePoints(end)'/4/3600])
    box off
    title ('Max EV Dim')
    xlabel('Time (hour)') % reference
    ylabel('Time (hour)') % test
end


function plot_EVopt1Dim_Time_Time(timePoints, matEV, opt1Dim, minNumFactor)
    maxEV      = zeros(size(matEV,1), size(matEV,2));
    
    for  nPlot = 1:length(timePoints)-1
        maxEV(:,nPlot) = squeeze(matEV(nPlot, :, opt1Dim(nPlot)-minNumFactor+1));
    end
    
    imagesc(timePoints(1:end-1)/4/3600, timePoints(1:end-1)/4/3600, maxEV);
    xlim([0 timePoints(end)'/4/3600])
    ylim([0 timePoints(end)'/4/3600])
    box off
    title ('Optimal Total EV')
    xlabel('Time (hour)') % reference
    ylabel('Time (hour)') % test
    colorbar
    caxis([-0.3 0.7])
end