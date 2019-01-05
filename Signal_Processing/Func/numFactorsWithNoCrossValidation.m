% 
% compute number of factors using different criteria
%
% Input:
% 
% X             ---  N x p
% maxM          ---  max number of factors in search
%
% Output:
%
% eigenvalue based result:
% kgM           --- # of factors using Kaiser-Guttman criterion
% varM          --- # of factors using total variance (95%)
% paM           --- # of factors using parallel analysis
% vMapM         --- # of factors using Velicer's MAP (4th order)
%
% likelihood based result:
% pChisqTest    --- p-values from chi-square test
% AIC
% BIC
% CAIC
% 
% estimation error results:
% SRMR          --- signficance of residuals
% RMSEA         --- root mean square error of approximation
% AGFI          --- adjusted goodness of fit index (dropped; not stable)
% CFI           --- comparative fit index
%
% plots:
% Scree Plot
% \chi^2 Plot
% Likelihood Plot
% Estimation Error Plot
%
% -------------------------------------------------------------------------
% % Example
% % generate artificial data
% 
% maxM                    = 20;
% numFactor               = 4;
% numTrials               = 100;
% realLoadingMat          = tril(randn(maxM, numFactor));
% commonFactors           = randn(numTrials, numFactor);
% noise                   = randn(numTrials, maxM) * 0.6;
% 
% X                       = commonFactors * realLoadingMat' + noise;
% 
% isScreePlot             = true;
% numFactors              = numFactorsWithNoCrossValidation(X, 10, isScreePlot);
% 
% disp('Based on eigenvalue analysis')
% disp(['# of factors using Kaiser-Guttman criterion: ', num2str(numFactors.kgM)]);
% disp(['# of factors using total variance (95%): ', num2str(numFactors.varM)]);
% disp(['# of factors using parallel analysis: ', num2str(numFactors.paM)]);
% disp(['# of factors using Velicer MAP (4th order): ', num2str(numFactors.vMapM)]);
% 
% disp('Based on likelihood analysis')
% disp(['# of factors using chi square test: ', num2str(numFactors.pChisqTestM)]);
% disp(['# of factors using AIC: ', num2str(numFactors.AICM)]);
% disp(['# of factors using BIC: ', num2str(numFactors.BICM)]);
% disp(['# of factors using CAIC: ', num2str(numFactors.CAICM)]);
% 
% disp('Based on estimation error analysis')
% disp(['# of factors using signficance of residuals: ', num2str(numFactors.SRMRM)]);
% disp(['# of factors using root mean square error of approximation: ', num2str(numFactors.RMSEAM)]);
% disp(['# of factors using comparative fit index: ', num2str(numFactors.CFIM)]);
%
% -------------------------------------------------------------------------
% version 1.0
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org


function numFactors = numFactorsWithNoCrossValidation(X, maxM, isScreePlot)

    % Potential warning comes from factoran when Psi has close to zeros variance
    warning('off','all'); 

    if nargin < 3; isScreePlot = false; end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%% eigenvalue based result:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [ntrial, nvars] = size(X);
    r           = corr(X);
    [eigvect, eigval] = eig(r);
    eigval      = diag(eigval);
    [eigval, k] = sort(eigval,'descend'); % sort the eigenvalues & get the indicies
    eigvect     = eigvect(:,k); % sort the eigenvectors based on the indicies
    
    % maxM from K-G criterion
    kgM         = sum(eigval > 1.0); 
    numFactors.kgM = kgM;
    
    % 0.95 variance criterion
    perVar      = 0.95;
    varM        = sum(cumsum(eigval) < nvars * perVar) + 1; 
    numFactors.varM = varM;
    
    
    % Velicer's MAP test -- 4th order
    % code from O'Connor, B. P. (2000)
    %
    loadings    = eigvect * sqrt(diag(eigval,0));
    fm4         = [(1:nvars); (1:nvars)]';
    fm4(1,2)    = (sum(sum(r.^4))-nvars)/(nvars*(nvars-1));
    for m       = 1:nvars - 1
        biga    = loadings(:,1:m);
        partcov = r - (biga * biga');
        d       = diag((1 ./ sqrt(diag(partcov))), 0);
        pr      = d * partcov * d;
        fm4(m+1,2) = (sum(sum(pr.^4))-nvars)/(nvars*(nvars-1));
    end
    %
    % identifying the smallest fm value & its location
    %
    vMapM       = find(fm4(:,2) == min(fm4(:,2)), 1, 'first') - 1;
    numFactors.vMapM = vMapM;
    
    % parallel analysis
    kind        = 2; % FA analysis
    randtype    = 2; % permutation of data
    percent     = perVar * 100;
    ndatsets    = 100;
    paM         = parallelAnalysis(X, kind, randtype, percent, ndatsets);
    numFactors.paM = paM;
    
    % summary of number of FA from eigenvalue analysis
    eigMaxM     = max([kgM, vMapM, paM]);    
    if maxM < eigMaxM
%         disp(['max number of factor is reset to ' num2str(eigMaxM) ', using eigvalue criterion.'])
        maxM = eigMaxM;
    end
    
    if maxM > nvars + 0.5 -sqrt(2*nvars + 0.25) || maxM > nvars        
        maxM = floor(min(nvars-1, nvars + 0.5 -sqrt(2*nvars + 0.25)));
%         disp(['max number of factor is reset to ' num2str(maxM) ', using freedom criterion.'])
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%% likelihood based results:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
    pChisqTest   = zeros(maxM, 1);
    valueAIC     = zeros(maxM, 1);
    valueBIC     = zeros(maxM, 1);
    valueCAIC    = zeros(maxM, 1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%% estimation error results:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    SRMR         = zeros(maxM, 1); % signficance of residuals
    RMSEA        = zeros(maxM, 1); % root mean square error of approximation
%     AGFI         = zeros(maxM, 1); % 
    CFI          = zeros(maxM, 1); % 
    
    numParams0   = 0.5 * nvars * (nvars+1);
    numParams00  = 0.5 * nvars * (nvars-1);
    
    for m = 1:maxM
        [Lambda, Psi] = factoran(X, m, 'scores','regression', 'rotate', 'none');
        S                       = diag(Psi) + Lambda*Lambda';
        logDetS                 = log(det(S));
        logDetR                 = log(det(r));
        traceSR                 = trace(S\r);
        logLik                  = -ntrial * (nvars*log(2*pi) + logDetS + traceSR);
        numParams               = (m+1)*nvars - 0.5*m*(m-1);
        
        nll                     = logDetS - logDetR;
        logLikRatio             = -nll;
        chisq                   = -(ntrial - (2*nvars + 11)/6 - 2*m/3) * logLikRatio;
        dfe                     = .5*((nvars-m)^2 - (nvars+m));
        pChisqTest(m)           = chi2pval(chisq, dfe);
        valueAIC(m)             = -2*logLik + 2* numParams;
        valueBIC(m)             = -2*logLik + log(ntrial)* numParams;
        valueCAIC(m)            = -2*logLik + (log(ntrial)+1)* numParams;
        
        fHat                    = logDetS - logDetR + traceSR - nvars;
        fHat0                   = sum(log(diag(r))) - logDetR;
        
        SRMR(m)                 = sqrt(sum(sum((tril(r-S)).^2))/numParams0);
        RMSEA(m)                = sqrt(max(fHat/numParams - 1/ntrial, 0));
%         GFI                     = 1 - sum(sum((S\(r-S)).^2))/sum(sum((S\r).^2));
%         AGFI(m)                 = 1 - (1-GFI)*nvars*(nvars+1)/(nvars*(nvars+1) - 2*numParams);
        CFI(m)                  = 1 - (fHat - numParams/ntrial)/(fHat0 - numParams00/ntrial);
    end
    
%     numFactors.pChisqTest   = pChisqTest;
%     numFactors.AIC          = valueAIC;
%     numFactors.BIC          = valueBIC;
%     numFactors.CAIC         = valueCAIC;
    
    pThres                  = 0.05;
    
    numFactors.pChisqTestM  = find(pChisqTest > pThres, 1, 'first');
    numFactors.AICM         = find(valueAIC == min(valueAIC), 1, 'first');
    numFactors.BICM         = find(valueBIC == min(valueBIC), 1, 'first');
    numFactors.BICM         = find(valueBIC == min(valueBIC), 1, 'first');
    numFactors.CAICM        = find(valueCAIC == min(valueCAIC), 1, 'first');
    
    
%     numFactors.SRMR         = SRMR;
%     numFactors.RMSEA        = RMSEA;
% %     numFactors.AGFI         = AGFI;
%     numFactors.CFI          = CFI;
    
    SRMRThres               = 0.05;
    RMSEAThres              = 0.08;
    CFIThres                = 0.90;

    numFactors.SRMRM        = find(SRMR < SRMRThres, 1, 'first');
    numFactors.RMSEAM       = find(RMSEA < RMSEAThres, 1, 'first');
%     numFactors.AGFIM        = find(AGFI == max(AGFI), 1, 'first');
    numFactors.CFIM         = find(CFI > CFIThres, 1, 'first');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%% Plots:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    
    if isScreePlot
        figure;
        %%% Scree plot
        subplot(2, 2, 1)
        plot(1: maxM, eigval(1: maxM), '-ob', 'linewidth',2);
        hold on;
        gridxy([], 1.0, 'Color','k','Linestyle','--', 'linewidth',2) ;
        h(1) = gridxy(kgM,'Color','c','Linestyle','--', 'linewidth',2) ;
        h(2) = gridxy(varM,'Color','r','Linestyle','--', 'linewidth',2) ;
        h(3) = gridxy(vMapM,'Color','g','Linestyle','--', 'linewidth',2) ;
        h(4) = gridxy(paM,'Color','m','Linestyle','--', 'linewidth',2) ;
        legend(h, {'Kaiser-Guttman', 'Variance', 'MAP', 'PA'})
        xlim([0.5 maxM+0.5])
        title('Scree Plot')        
        
        %%% chi-square plot
        subplot(2, 2, 2)
        plot(1: maxM, pChisqTest, '-ob', 'linewidth',2);
        hold on;
        gridxy([], pThres, 'Color','k','Linestyle','--', 'linewidth',2) ;
        xlim([0.5 maxM+0.5])
        title('\chi^2 Plot')        

        %%% chi-square plot
        subplot(2, 2, 3)
        plot(1: maxM, valueAIC, '-o', 'linewidth',2);
        hold on;
        plot(1: maxM, valueBIC, '-o', 'linewidth',2);
        plot(1: maxM, valueCAIC, '-o', 'linewidth',2);
        legend({'AIC', 'BIC', 'CAIC'})
        xlim([0.5 maxM+0.5])
        title('Likelihood Plot')  
        
        %%% chi-square plot
        subplot(2, 2, 4)
        plot(1: maxM, SRMR, '-o', 'linewidth',2);
        hold on;
        plot(1: maxM, RMSEA, '-o', 'linewidth',2);
        plot(1: maxM, CFI, '-o', 'linewidth',2);
        legend({'SRMR', 'RMSEA', 'CFI'})
        gridxy([], [SRMRThres, RMSEAThres, CFIThres], 'Color','k','Linestyle','--', 'linewidth',2) ;
        xlim([0.5 maxM+0.5])
        title('Estimation Error Plot')  
    end

    