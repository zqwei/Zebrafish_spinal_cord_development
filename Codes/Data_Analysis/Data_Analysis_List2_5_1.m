%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.  Covariance analysis: factor analysis
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.4 FA vs islands
%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org


function Data_Analysis_List2_5_1(nFile)

    if nargin<1
        nFile = 1;
    end
        
    addpath('../Func');
    setDir;    
    fileDirName       = fileDirNames{nFile}; %#ok<NASGU,USENS>
    fileName          = fileNames{nFile}; %#ok<USENS>
    load([tempDatDir, fileName, '.mat'],'timePoints');
    load([tempDatDir, fileName, '_LMatTime.mat'], 'LMatTime')
    load([tempDatDir, fileName, '_FAEV.mat'], 'opt2Dim')
    coeffThres        = 0.3;
    minSizeCluster    = 3;
        
    figure;
    hold on;
    for nTime           = 1:length(timePoints)-1        
        sizeCluster     = sum(abs(LMatTime{nTime})>coeffThres, 1); %#ok<USENS>
        sizeCluster     = sizeCluster(sizeCluster>minSizeCluster);
        plot(timePoints(nTime)/4/3600*ones(length(sizeCluster),1), sizeCluster+rand(size(sizeCluster)), '+b')
    end  
    xlabel('Time (hour)')
    ylabel('# neurons loading mat.')
    xlim([0 timePoints(end)/4/3600]); %#ok<COLND>
    setPrint(6, 4.5, [plotDir, 'FALMatNumNeuron_', fileName], 'pdf')
end

