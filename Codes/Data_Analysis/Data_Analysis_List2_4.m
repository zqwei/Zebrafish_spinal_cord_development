%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.  Covariance analysis: factor analysis
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.3 EV vs location
%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function Data_Analysis_List2_4(nFile)

    if nargin<1
        nFile = 1;
    end
        
    addpath('../Func');
    setDir;    
    fileDirName       = fileDirNames{nFile}; %#ok<NASGU,USENS>
    fileName          = fileNames{nFile}; %#ok<USENS>
    load([tempDatDir, fileName, '.mat']);
    load([tempDatDir, fileName, '_FAEVCrossTime.mat'], 'matEV', 'matEVSingleUnit'); 
    load([tempDatDir, fileName, '_FAEV.mat'], 'opt2Dim');
    
    % Plot EVs as a function of x-y location
    numPlot          = length(timePoints)-1;
    mCol             = 8;
    mRow             = ceil(numPlot/mCol);
    minNumFactor     = 2;
    opt1Dim          = opt2Dim-minNumFactor+1;
    opt1Dim          = min(opt1Dim, 8);
        
    figure;
    h = suptitle(fileName);
    set(h,'Interpreter','none');
    
    
    
    color_max        = 0.7;
    color_min        = -0.3;
    
    color_nTime      = zeros(size(matEVSingleUnit,4),3); 
    
    for nTime           = 1:length(timePoints)-1
        ev_units_nTime  = squeeze(matEVSingleUnit(nTime,nTime,opt1Dim(nTime),:));
        color_nTime(:,1)= (ev_units_nTime - color_min)/(color_max - color_min);
        color_nTime(:,2)= 1 - (ev_units_nTime - color_min)/(color_max - color_min);
        subplot(mRow, mCol, nTime)
        xTrack          = squeeze(mean(tracks(:, timePoints(nTime)+1:timePoints(nTime+1), 1), 2)); %#ok<NODEF>
        yTrack          = squeeze(mean(tracks(:, timePoints(nTime)+1:timePoints(nTime+1), 2), 2));
        scatter(xTrack, yTrack, [], color_nTime)
        title(['Time from: ' num2str(timePoints(nTime)/4/3600) ' to ' num2str(timePoints(nTime+1)/4/3600) 'hr'])
    end   

end