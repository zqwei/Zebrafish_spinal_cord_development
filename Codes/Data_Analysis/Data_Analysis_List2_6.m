%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.  Covariance analysis: factor analysis
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.4 FA vs islands -- Tree plot
%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function Data_Analysis_List2_6(nFile)

    if nargin<1
        nFile = 1;
    end
        
    addpath('../Func');
    setDir;    
    fileDirName         = fileDirNames{nFile}; %#ok<NASGU,USENS>
    fileName            = fileNames{nFile}; %#ok<USENS>
    load([tempDatDir, fileName, '.mat']);
    load([tempDatDir, fileName, '_FAEV.mat'], 'opt2Dim');
    
    numPlot             = length(timePoints)-1;
    mCol                = 8;
    mRow                = ceil(numPlot/mCol);
    
    nodeRoot            = 1;
    nodes               = nan(1,1000);
    nodes(1)            = 0;
    nNodes              = 1;
    Thres               = 0.3;
    
    for nTime           = length(timePoints)-1:-1:1      
        slicedDFF       = dff(:,timePoints(nTime)+1:timePoints(nTime+1)); %#ok<NODEF>
        slicedDFF       = bsxfun(@minus, slicedDFF, mean(slicedDFF,2));
        slicedDFF       = bsxfun(@rdivide, slicedDFF, std(slicedDFF,[],2))';                
        lMatTime        = factoran(slicedDFF, opt2Dim(nTime),'maxit',10000)>Thres;
        if nTime == length(timePoints)-1
            parentNodes = nan(size(lMatTime,2),1);
            for iNode   = 1:size(lMatTime,2)
                if sum(lMatTime(:,iNode))>0
                    nNodes             = nNodes + 1;
                    parentNodes(iNode) = nNodes;
                    nodes(nNodes)      = nodeRoot;
                end
            end
            parentNodes = parentNodes(~isnan(parentNodes));
        else
            parentNodes = nan(size(lMatTime,2),1);
            for iNode   = 1:size(lMatTime,2)
                if sum(lMatTime(:,iNode))>0
                    nNodes             = nNodes + 1;
                    parentNodes(iNode) = nNodes;
                    [~, pIndex]        = max(double(lMatTime(:,iNode))'*double(lMatPreTime));
                    nodes(nNodes)      = parentPreNodes(pIndex);
                end
            end
            parentNodes = parentNodes(~isnan(parentNodes)');
        end
        lMatPreTime      = lMatTime;  
        parentPreNodes   = parentNodes;
    end   
    
    treeplot(nodes(~isnan(nodes)))
end