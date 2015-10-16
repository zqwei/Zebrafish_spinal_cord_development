% 
% parallelAnalysis
% 
% code adopted from O'Connor, B. P. (2000)
% [realeval, means, percentiles] = parallelAnalysis(rawData, kind, randtype, percent, ndatsets)
% -------------------------------------------------------------------------
% version 1.0
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org

function numFactor = parallelAnalysis(rawData, kind, randtype, percent, ndatsets)

    % Specify the desired kind of parellel analysis, where:
    % 1 = principal components analysis
    % 2 = principal axis / common factor analysis
    if nargin < 2; kind     = 2; end
    
    % Enter either
    %  1 for normally distributed random data generation parallel analysis, or
    %  2 for permutations of the raw data set (more time consuming).
    if nargin < 3; randtype = 2; end
    
     % Enter the desired percentile here
    if nargin < 4; percent  = 95; end
    
    % Enter the desired number of parallel data sets here
    if nargin < 5; ndatsets  = 1000; end    
    
    [ncases, nvars]  = size(rawData);
    r                = corrcoef(rawData);
    realeval         = flipud(sort(eig(r)));
    evals            = zeros(nvars, ndatsets);
    
    % principal components analysis & random normal data generation
    if  kind == 1 && randtype == 1
%         realeval = flipud(sort(eig(corrcoef(rawData))));
        for nds  = 1:ndatsets; evals(:,nds) = eig(corrcoef(randn(ncases,nvars)));end
    end

    % principal components analysis & raw data permutation
    if  kind == 1 && randtype == 2
%         realeval = flipud(sort(eig(corrcoef(rawData))));
        for nds = 1:ndatsets 
            x = rawData;
            for lupec = 1:nvars
                for luper = 1:(ncases -1)
                    k = fix( (ncases - luper + 1) * rand(1) + 1 )  + luper - 1;
                    d = x(luper,lupec);
                    x(luper,lupec) = x(k,lupec);
                    x(k,lupec) = d;
                end
            end
            evals(:,nds) = eig(corrcoef(x));
        end
    end

    % PAF/common factor analysis & random normal data generation
    if  kind == 2 && randtype == 1
        smc = 1 - (1 ./ diag(inv(r)));
        for ii=1:size(r,1); r(ii,ii) = smc(ii,1); end;
        realeval = flipud(sort(eig(r)));
        for nds = 1:ndatsets 
            r = corrcoef(randn(ncases,nvars));
            smc = 1 - (1 ./ diag(inv(r)));
            for ii=1:size(r,1); r(ii,ii) = smc(ii,1); end;
            evals(:,nds) = eig(r);
        end
    end

    % PAF/common factor analysis & raw data permutation
    if  kind == 2 && randtype == 2
        smc = 1 - (1 ./ diag(inv(r)));
        for ii=1:size(r,1); r(ii,ii) = smc(ii,1); end;
        realeval = flipud(sort(eig(r)));
        for nds = 1:ndatsets; 
            x = rawData;
            for lupec = 1:nvars;
                for luper = 1:(ncases -1);
                    k = fix( (ncases - luper + 1) * rand(1) + 1 )  + luper - 1;
                    d = x(luper,lupec);
                    x(luper,lupec) = x(k,lupec);
                    x(k,lupec) = d;
                end
            end
            r = corrcoef(x);
            smc = 1 - (1 ./ diag(inv(r)));
            for ii=1:size(r,1); r(ii,ii) = smc(ii,1); end;
            evals(:,nds) = eig(r);
        end
    end


    evals = flipud(sort(evals, 1));
%     means = (mean(evals,2));   % mean eigenvalues for each position.
    evals = sort(evals,2);     % sorting the eigenvalues for each position.
    percentiles = (evals(:,round((percent*ndatsets)/100)));  % percentiles.
    
    numFactor  = sum(realeval > percentiles & realeval>0);
    
    
    