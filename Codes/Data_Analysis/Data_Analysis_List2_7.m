%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.  Covariance analysis: factor analysis
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.4 FA vs islands -- clustering
%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function Data_Analysis_List2_7(nFile)

    if nargin<1
        nFile = 1;
    end
        
    addpath('../Func');
    setDir;    
    fileDirName         = fileDirNames{nFile}; %#ok<NASGU,USENS>
    fileName            = fileNames{nFile}; %#ok<USENS>
    load([tempDatDir, fileName, '.mat']);
    load([tempDatDir, fileName, '_FAEV.mat'], 'opt2Dim');
    
    mCol                = 8;
    mRow                = ceil((length(timePoints)-1)/mCol);
        
    Thres               = 0.3;
    
    LMatTime            = cell(length(timePoints)-1,1);
    
    for nTime           = 1: length(timePoints)-1     
        slicedDFF       = dff(:,timePoints(nTime)+1:timePoints(nTime+1)); %#ok<NODEF>
        slicedDFF       = bsxfun(@minus, slicedDFF, mean(slicedDFF,2));
        slicedDFF       = bsxfun(@rdivide, slicedDFF, std(slicedDFF,[],2))';                
        LMatTime{nTime} = factoran(slicedDFF, opt2Dim(nTime),'maxit',10000);
    end   
    
    save([tempDatDir, fileName, '_LMatTime.mat'], 'LMatTime')
    
end