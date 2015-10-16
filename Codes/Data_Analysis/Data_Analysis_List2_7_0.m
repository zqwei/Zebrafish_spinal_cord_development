%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.  Covariance analysis: factor analysis
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.4 FA vs islands -- clustering -- Computing unrotated Loading matrix,
% Psi, correlation matrix etc.
%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function Data_Analysis_List2_7_0(nFile)

    if nargin<1
        nFile = 1;
    end
        
    addpath('../Func');
    setDir;    
    fileDirName         = fileDirNames{nFile}; %#ok<NASGU,USENS>
    fileName            = fileNames{nFile}; %#ok<USENS>
    load([tempDatDir, fileName, '.mat']);
%     load([tempDatDir, fileName, '_FAEV.mat'], 'opt2Dim');
    load([tempDatDir, fileName, '_PSDPeakTimeFAEV.mat'], 'opt2Dim');
    load([tempDatDir, fileName, '_PSDPeakTime.mat'], 'PSDPeakTime')
    
    
    LMatTime            = cell(length(timePoints)-1, 1);
    numNeuron           = size(dff, 1); %#ok<NODEF>
    PsiMat              = ones(length(timePoints)-1, numNeuron);
    corrMat             = zeros(length(timePoints)-1, numNeuron, numNeuron);
    
    for nTime           = 1: length(timePoints)-1     
        activeIndex     = PSDPeakTime(:, nTime)>0 & ~isnan(PSDPeakTime(:, nTime)); %#ok<NODEF>
        allDFF          = dff(:,timePoints(nTime)+1:timePoints(nTime+1));
        slicedDFF       = dff(activeIndex,timePoints(nTime)+1:timePoints(nTime+1));
        slicedDFF       = bsxfun(@minus, slicedDFF, mean(slicedDFF,2));
        slicedDFF       = bsxfun(@rdivide, slicedDFF, std(slicedDFF,[],2))';  
        if opt2Dim(nTime) >= 1 && ~isnan(opt2Dim(nTime))
            LambdaAll       = zeros(numNeuron, opt2Dim(nTime));
            [Lambda, Psi]   = factoran(slicedDFF, opt2Dim(nTime),'maxit',10000);
            PsiMat(nTime, activeIndex) = Psi;
            corrMat(nTime, :, :) = corr(allDFF');
            LambdaAll(activeIndex, :) = Lambda;
            LMatTime{nTime} = LambdaAll;
        else
            LMatTime{nTime} = zeros(numNeuron, 1);
        end
    end   
    
    save([tempDatDir, fileName, '_LMatTime.mat'], 'LMatTime', 'PsiMat', 'corrMat');
    
end