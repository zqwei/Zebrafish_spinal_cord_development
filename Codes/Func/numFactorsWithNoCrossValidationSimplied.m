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
% 
% estimation error results:
% SRMR          --- signficance of residuals
% CFI           --- comparative fit index
%
%
% -------------------------------------------------------------------------
% version 1.0
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org


function numFactors = numFactorsWithNoCrossValidationSimplied(X, maxM)

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
        
    % summary of number of FA from eigenvalue analysis
    eigMaxM     = kgM;    
    if maxM < eigMaxM
%         disp(['max number of factor is reset to ' num2str(eigMaxM) ', using eigvalue criterion.'])
        maxM = eigMaxM;
    end
    
    if maxM > nvars + 0.5 -sqrt(2*nvars + 0.25) || maxM > nvars        
        maxM = floor(min(nvars-1, nvars + 0.5 -sqrt(2*nvars + 0.25)));
%         disp(['max number of factor is reset to ' num2str(maxM) ', using freedom criterion.'])
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%% estimation error results:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    SRMR         = zeros(maxM, 1); % signficance of residuals
    CFI          = zeros(maxM, 1); % 
    
    numParams0   = 0.5 * nvars * (nvars+1);
    numParams00  = 0.5 * nvars * (nvars-1);
    
    for m = 1:maxM
        [Lambda, Psi] = factoran(X, m, 'scores','regression', 'rotate', 'none');
        S                       = diag(Psi) + Lambda*Lambda';
        logDetS                 = log(det(S));
        logDetR                 = log(det(r));
        traceSR                 = trace(S\r);
        numParams               = (m+1)*nvars - 0.5*m*(m-1);
        
        fHat                    = logDetS - logDetR + traceSR - nvars;
        fHat0                   = sum(log(diag(r))) - logDetR;
        
        SRMR(m)                 = sqrt(sum(sum((tril(r-S)).^2))/numParams0);
        CFI(m)                  = 1 - (fHat - numParams/ntrial)/(fHat0 - numParams00/ntrial);
    end
        
    SRMRThres               = 0.05;
    RMSEAThres              = 0.08;
    CFIThres                = 0.90;

    numFactors.SRMRM        = find(SRMR < SRMRThres, 1, 'first');
    numFactors.CFIM         = find(CFI > CFIThres, 1, 'first');
    