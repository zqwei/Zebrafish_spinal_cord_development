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


function Data_Analysis_List2_7_1n(nFile, nTimePoints)

    if nargin<1
        nFile = 1;
    end
        
    addpath('../Func');
    setDir;    
    fileDirName         = fileDirNames{nFile}; %#ok<NASGU,USENS>
    fileName            = fileNames{nFile}; %#ok<USENS>
    load([tempDatDir, fileName, '.mat']);
    load([tempDatDir, fileName, '_FAEV.mat'], 'opt2Dim');
%     load([tempDatDir, fileName, '_LMatTime.mat'], 'LMatTime', 'PsiMat', 'corrMat');
    
%     peakTime            = find(opt2Dim == max(opt2Dim), 1);    
%     nTimePoints         = [3, peakTime, length(opt2Dim) - 3];
    totVar              = size(dff, 1); %#ok<NODEF>
    
    nCol                = max(opt2Dim) + 3;
    nRow                = length(nTimePoints);
    
    figure;
    
    for nPlot           = 1:nRow
        nTime           = nTimePoints(nPlot);
        slicedDFF       = dff(:,timePoints(nTime)+1:timePoints(nTime+1));
        distNeurons     = pdist(slicedDFF, 'correlation');
        linkNeurons     = linkage(slicedDFF,'single','correlation');
        leafOrder       = optimalleaforder(linkNeurons, distNeurons); 
        slicedDFF       = slicedDFF(leafOrder,:);        
        slicedDFF       = bsxfun(@minus, slicedDFF, mean(slicedDFF,2));
        slicedDFF       = bsxfun(@rdivide, slicedDFF, std(slicedDFF,[],2))';
        
        [Lambda, Psi]   = factoran(slicedDFF, opt2Dim(nTime),'maxit',10000);
        Psi             = diag(Psi);
        corrMat         = corr(slicedDFF);
        newLambda       = pcaRotationFactors(Lambda);
        resMat          = corrMat - (newLambda * newLambda' + Psi);
        
        % Correlation matrix plot
        subplot(length(nTimePoints), nCol, (nPlot-1)*nCol + 1)
        imagesc(corrMat); caxis([-1 1]); title('Correlation Matrix')
        
        for nDim       = 1:opt2Dim(nTime)
            subplot(length(nTimePoints), nCol, (nPlot-1)*nCol + 1 + nDim)
            corrFactorMat = newLambda(:,nDim) * newLambda(:,nDim)';
            imagesc(corrFactorMat); caxis([-1 1]); 
            title(['Factor #' num2str(nDim) ': %EV ' ...
                num2str(newLambda(:,nDim)' * newLambda(:,nDim)/totVar*100,'%.2f')])
        end
        
        subplot(length(nTimePoints), nCol, (nPlot-1)*nCol + 1 + nDim + 1)
        imagesc(Psi); caxis([-1 1]); 
        title(['Psi: %EV ' num2str(trace(Psi)/totVar*100,'%.2f')])
            
        subplot(length(nTimePoints), nCol, (nPlot-1)*nCol + 1 + nDim + 2)
        imagesc(resMat); caxis([-1 1]); 
        title(['Residual: %EV ' num2str(trace(resMat)/totVar*100,'%.2f')])                
    end   
    
    setFigureSize(8*nCol, 6*nRow);
        
end