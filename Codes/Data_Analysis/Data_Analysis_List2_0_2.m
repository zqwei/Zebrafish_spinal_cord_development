%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.  Covariance analysis: factor analysis
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.0 Figure summary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 

function Data_Analysis_List2_0_2(nFile)
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>
    
    load([tempDatDir, fileName, '.mat']);
    minNumFactor      = 2;
    load([tempDatDir, fileName, '_FAEV.mat'], 'opt2Dim');
    load([tempDatDir, fileName, '_FAEVCrossTime.mat'], 'matEV', 'matEVSingleUnit');
    
    opt1Dim          = opt2Dim-minNumFactor+1;
    opt1Dim          = min(opt1Dim, 8);    
    
    max_ev_units        = nan(size(matEVSingleUnit,4), length(timePoints)-1); %#ok<NODEF>
    for nTime           = 1:length(timePoints)-1
        max_ev_units(:,nTime) = squeeze(matEVSingleUnit(nTime,nTime,opt1Dim(nTime),:));
    end 
    numPlot          = size(max_ev_units, 1);
        
    RSquare          = zeros(numPlot, 1);
    halfTime         = zeros(numPlot, 1);
    for nUnit        = 1:numPlot
        [fitParams, fitResult] = sigm_fit((timePoints(1:end-1)/4/3600), max_ev_units(nUnit,:));
        RSquare(nUnit)         = 1 - mean((max_ev_units(nUnit,:)' - fitResult.ypred).^2)./var(max_ev_units(nUnit,:));
        halfTime(nUnit)        = fitParams(3);
    end
    
    save([tempDatDir, fileName, '_HalfTime.mat'], 'RSquare', 'halfTime');
end