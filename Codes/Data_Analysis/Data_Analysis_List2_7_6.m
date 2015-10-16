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


function Data_Analysis_List2_7_6(nFile, nTime)

    if nargin<1
        nFile = 1;
    end
        
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>
    
    load([tempDatDir, fileName, '.mat']);
    load([tempDatDir, fileName, '_FAEV.mat'], 'opt2Dim');
    
%     totVar              = size(dff, 1); %#ok<NODEF>    
        
    dff_threshould      = 0.03;
    numDim              = opt2Dim(nTime);
    slicedDFF           = dff(:,timePoints(nTime)+1:timePoints(nTime+1));
    activeIndex         = std(slicedDFF, [], 2) > dff_threshould;
    df                  = .5*((sum(activeIndex)-numDim)^2 - (sum(activeIndex)+numDim));
    if (sum(activeIndex) <= numDim) || (df<0); return; end
    slicedDFF           = slicedDFF(activeIndex, :);
    distNeurons         = pdist(slicedDFF, 'correlation');
    linkNeurons         = linkage(slicedDFF,'single','correlation');
    leafOrder           = optimalleaforder(linkNeurons, distNeurons); 
    slicedDFF           = slicedDFF(leafOrder,:)';        
%     slicedDFF           = bsxfun(@minus, slicedDFF, mean(slicedDFF,2));
%     slicedDFF           = bsxfun(@rdivide, slicedDFF, std(slicedDFF,[],2))';
%     numDim              = opt2Dim(nTime);
%     Lambda              = factoran(slicedDFF, numDim, 'maxit', 10000);
    
    Lambda              = factoran(slicedDFF, numDim, 'maxit', 10000);
    [newLambda, ~]      = nnmf(Lambda, numDim);

    colorTimes          = linspecer(numDim);
    
    figure;   
    
    hold on;
    
    for nNeuron = 1:sum(activeIndex)
        nDim    = find(newLambda(nNeuron, :)==max(newLambda(nNeuron, :)), 1, 'first');  
        plot(slicedDFF(:, nNeuron) + nNeuron, 'Color', colorTimes(nDim, :));
    end
    
    
    figure;
    hold on;
    
    slicedDFF           = bsxfun(@minus, slicedDFF, mean(slicedDFF,1));
    slicedDFF           = bsxfun(@rdivide, slicedDFF, std(slicedDFF,[],1));
    
    for nNeuron = 1:sum(activeIndex)
        nDim    = find(newLambda(nNeuron, :)==max(newLambda(nNeuron, :)), 1, 'first');  
        plot(slicedDFF(:, nNeuron) + nNeuron, 'Color', colorTimes(nDim, :));
    end
    
    colormap(colorTimes)
    
end