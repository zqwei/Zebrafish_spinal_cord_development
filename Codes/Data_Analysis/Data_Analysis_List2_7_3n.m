%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.  Covariance analysis: factor analysis
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.4 FA vs islands -- clustering -- Computing unrotated Loading matrix,
% Psi, correlation matrix etc. -- Psi as function of time
%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function Data_Analysis_List2_7_3n(nFile)

    if nargin<1
        nFile = 1;
    end
        
    addpath('../Func');
    setDir;    
    fileDirName         = fileDirNames{nFile}; %#ok<NASGU,USENS>
    fileName            = fileNames{nFile}; %#ok<USENS>
    load([tempDatDir, fileName, '.mat'], 'timePoints');
    load([tempDatDir, fileName, '_LMatTime.mat'], 'PsiMat');
    
    figure;
    
    plot(timePoints(1:end-1)/4/3600, mean(PsiMat, 2),'o','linewid',2)
%     ylim([0 1])
    xlim([0 timePoints(end)'/4/3600])
    ylabel('EV uncorrelated noise')
    xlabel('Time (hour)')
    setFigureSize(8, 6);
        
end