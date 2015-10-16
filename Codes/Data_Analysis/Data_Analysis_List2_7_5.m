%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.  Covariance analysis: factor analysis
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.4 FA vs islands -- clustering -- Computing unrotated Loading matrix,
% Psi, correlation matrix etc. -- Raw data by FA communities
%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function Data_Analysis_List2_7_5(nFile, nTime)

    if nargin<1
        nFile = 1;
    end
        
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>
    
    load([tempDatDir, fileName, '.mat']);
    load([tempDatDir, fileName, '_FAEV.mat'], 'opt2Dim');
    
    totVar              = size(dff, 1); %#ok<NODEF>    
    
    slicedDFF           = dff(:,timePoints(nTime)+1:timePoints(nTime+1));
    distNeurons         = pdist(slicedDFF, 'correlation');
    linkNeurons         = linkage(slicedDFF,'single','correlation');
    leafOrder           = optimalleaforder(linkNeurons, distNeurons); 
    slicedDFF           = slicedDFF(leafOrder,:)';        
%     slicedDFF           = bsxfun(@minus, slicedDFF, mean(slicedDFF,2));
%     slicedDFF           = bsxfun(@rdivide, slicedDFF, std(slicedDFF,[],2))';
%     numDim              = opt2Dim(nTime);
%     Lambda              = factoran(slicedDFF, numDim, 'maxit', 10000);
    numDim              = opt2Dim(nTime);
    Lambda              = factoran(slicedDFF, numDim, 'maxit', 10000);
    [newLambda, ~]      = nnmf(Lambda, numDim);

    colorTimes          = linspecer(numDim+1);
    
    figure;   
    
    hold on;
    
    for nNeuron = 1:totVar
        nDim    = find(newLambda(nNeuron, :)>0.2, 1, 'first');  
        if isempty(nDim)
            nDim = 1;
        else
            nDim = nDim + 1;
        end
        plot(slicedDFF(:, nNeuron) + nNeuron, 'Color', colorTimes(nDim, :));
    end
    
    colormap(colorTimes)
    
end