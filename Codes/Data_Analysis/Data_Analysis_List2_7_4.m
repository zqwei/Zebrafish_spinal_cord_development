%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.  Covariance analysis: factor analysis
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.4 FA vs islands -- clustering -- Computing unrotated Loading matrix,
% Psi, correlation matrix etc. -- Lambda as function of time
%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function Data_Analysis_List2_7_4(nFile)

    if nargin<1
        nFile = 1;
    end
        
    addpath('../Func');
    setDir;    
    fileDirName         = fileDirNames{nFile}; %#ok<NASGU,USENS>
    fileName            = fileNames{nFile}; %#ok<USENS>
    load([tempDatDir, fileName, '.mat'], 'timePoints');
    load([tempDatDir, fileName, '_LMatTime.mat'], 'LMatTime');
    
    figure;
    hold on
    
    numTimes             = length(LMatTime); %#ok<USENS>
    
    colorTimes           = linspecer(numTimes);
    
    for nTime            = 1:numTimes
        LMat             = LMatTime{nTime};
        plot(sort(eig(LMat' * LMat)/size(LMat, 1),'descend'), '-o', 'color', colorTimes(nTime, :))
    end
    
    colormap(linspecer); hcb = colorbar; set(hcb,'YTick',[0, 1], 'YTickLabel', {num2str(timePoints(1)/3600/4,'%.2f'), num2str(timePoints(end)/3600/4,'%.2f')})
    
    ylim([0 0.4])
    ylabel('EV Loading Mat.')
    xlabel('# Factors')
    setPrint(8, 6, [plotDir, 'FADecompositionEVLMatCrossTime', fileName], 'pdf');
        
end