%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.  Covariance analysis: factor analysis
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.4 FA vs islands -- clustering -- Computing unrotated Loading matrix,
% Psi, correlation matrix etc. -- Rotaion
%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function Data_Analysis_List2_7_2(nFile)

    if nargin<1
        nFile = 1;
    end
        
    addpath('../Func');
    setDir;    
    fileDirName         = fileDirNames{nFile}; %#ok<NASGU,USENS>
    fileName            = fileNames{nFile}; %#ok<USENS>
    load([tempDatDir, fileName, '.mat']);
    load([tempDatDir, fileName, '_FAEV.mat'], 'opt2Dim');
    
%     peakTime            = find(opt2Dim == max(opt2Dim), 1);    
    totVar              = size(dff, 1); %#ok<NODEF>
    
    nTime               = length(opt2Dim) - 3;
    
    slicedDFF           = dff(:,timePoints(nTime)+1:timePoints(nTime+1));
    distNeurons         = pdist(slicedDFF, 'correlation');
    linkNeurons         = linkage(slicedDFF,'single','correlation');
    leafOrder           = optimalleaforder(linkNeurons, distNeurons); 
    slicedDFF           = slicedDFF(leafOrder,:);        
    slicedDFF           = bsxfun(@minus, slicedDFF, mean(slicedDFF,2));
    slicedDFF           = bsxfun(@rdivide, slicedDFF, std(slicedDFF,[],2))';
    numDim              = opt2Dim(nTime);
    Lambda              = factoran(slicedDFF, numDim, 'maxit', 10000);
    
    nRow                = 7;
    nCol                = numDim + 3;
    
    figure;
    
    % varimax
    nPlot               = 2;
    newLambda           = Lambda;
    subplotDecompostionLambaMat(newLambda, nPlot, nCol, nRow, totVar, 'VARIMAX');
    
    % PCA: EV-max
    nPlot               = 1;
    newLambda           = pcaRotationFactors(Lambda);
    subplotDecompostionLambaMat(newLambda, nPlot, nCol, nRow, totVar, 'EV-MAX');
    
    % quartimax
    nPlot               = 3;
    newLambda           = rotatefactors(Lambda, 'Method', 'quartimax', 'Maxit', 10000);
    subplotDecompostionLambaMat(newLambda, nPlot, nCol, nRow, totVar, 'QUARTIMAX');
    
    % equamax
    nPlot               = 4;
    newLambda           = rotatefactors(Lambda, 'Method', 'equamax', 'Maxit', 10000);
    subplotDecompostionLambaMat(newLambda, nPlot, nCol, nRow, totVar, 'EQUAMAX');
    
    % parsimax
    nPlot               = 5;
    newLambda           = rotatefactors(Lambda, 'Method', 'parsimax', 'Maxit', 10000);
    subplotDecompostionLambaMat(newLambda, nPlot, nCol, nRow, totVar, 'PARSIMAX');
    
    % oblique
    nPlot               = 6;
    newLambda           = rotatefactors(Lambda, 'Method', 'promax');
    H                   = (newLambda' * newLambda) \ (newLambda' * Lambda);
    subplotDecompostionLambaMat(newLambda, nPlot, nCol, nRow, totVar, 'OBLIQUE-PROMAX', H);
    
    % NNMF
    nPlot               = 7;
    [newLambda, ~]      = nnmf(Lambda, numDim);
    H                   = (newLambda' * newLambda) \ (newLambda' * Lambda);
    subplotDecompostionLambaMat(newLambda, nPlot, nCol, nRow, totVar, 'NNMF', H);
    
    setPrint(8*nCol, 6*nRow, [plotDir, 'FADecompositionRotation', fileName], 'pdf');
        
end

function subplotDecompostionLambaMat(newLambda, nPlot, nCol, nRow, totVar, titleName, H)
    
    subplot(nRow, nCol, (nPlot-1)*nCol + 1)
    imagesc(newLambda * newLambda'); caxis([-1 1]); title('Lambda Corr. Mat.')
    for nDim       = 1:size(newLambda, 2)
        subplot(nRow, nCol, (nPlot-1)*nCol + 1 + nDim)
        corrFactorMat = newLambda(:,nDim) * newLambda(:,nDim)';
        imagesc(corrFactorMat); caxis([-1 1]); 
        title(['Factor #' num2str(nDim) ': %EV ' ...
            num2str(newLambda(:,nDim)' * newLambda(:,nDim)/totVar*100,'%.2f')])
    end
    
    subplot(nRow, nCol, (nPlot-1)*nCol + 1 + nDim + 1)
    imagesc(newLambda); caxis([-1 1]); title([titleName ' Lambda'])
    
    if nargin > 6
        subplot(nRow, nCol, (nPlot-1)*nCol + 1 + nDim + 2)
        imagesc(H*H'); caxis([-1 1]); title([titleName ' H'])
    end
end