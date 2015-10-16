%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.  Covariance analysis: factor analysis
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.4 FA vs islands -- clustering -- overlay to 
%     following 2.7
%     For this plot, I ignored the points not in a cluster not a
%     combination of right or left groups
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function Data_Analysis_List2_8_1(nFile)

    if nargin<1
        nFile = 1;
    end
        
    addpath('../Func');
    setDir;    
    fileName            = fileNames{nFile}; %#ok<USENS>
    load([tempDatDir, fileName, '.mat']);
    load([tempDatDir, fileName, '_FAEV.mat'], 'opt2Dim');
    load([tempDatDir, fileName, '_LMatTime.mat'], 'LMatTime')
    load(['/Volumes/Zebrafish_Imaging_Data_by_Yinan/' fileName '/periods.mat']);
    flist               = dir(['/Volumes/Zebrafish_Imaging_Data_by_Yinan/' fileName '/image_data/*.tif']);

    Thres               = 0.3;
    FinalGroup          = false(size(LMatTime{1},1), 2);                                %#ok<USENS>
    [~, TimeG1]         = max(cell2mat(arrayfun(@(x) ...
                          max(sum(LMatTime{x}>Thres,1)), 1:length(LMatTime), ...
                          'UniformOutput', false))); 
    [~, FAIndexG1]      = max(sum(LMatTime{TimeG1}>Thres,1));
    FinalGroup(:,1)     = LMatTime{TimeG1}(:,FAIndexG1)>Thres;                      
    [~, TimeG2]         = max(cell2mat(arrayfun(@(x) ...
                          max(sum(LMatTime{x}(~FinalGroup(:,1),:)>Thres,1)), ...
                          1:length(LMatTime), ...
                          'UniformOutput', false)));
    [~, FAIndexG2]      = max(sum(LMatTime{TimeG2}(~FinalGroup(:,1),:)>Thres,1));
    FinalGroup(:,2)     = LMatTime{TimeG1}(:,FAIndexG2)>Thres;  
    
    maxColorValue       = 0.9;
                
    for nTime           = 1:length(timePoints)-1
        xyzTrack        = squeeze(mean(tracks(:, timePoints(nTime)+1:timePoints(nTime+1), :), 2)); %#ok<NODEF>
        LMatnTime       = LMatTime{nTime};
        coordinates     = [];
        colorTable      = [];
        radiusTable     = [];
        
        for nFactor     = size(LMatnTime, 2):-1:1
            LMatTimeFactorPos = LMatnTime(:,nFactor)>Thres;
            nNeuronPosInG1    = sum(LMatTimeFactorPos & FinalGroup(:,1))/sum(LMatTimeFactorPos);
            nNeuronPosInG2    = sum(LMatTimeFactorPos & FinalGroup(:,2))/sum(LMatTimeFactorPos);
            LMatTimeFactorNeg = LMatnTime(:,nFactor)<-Thres;
            nNeuronNegInG1    = sum(LMatTimeFactorNeg & FinalGroup(:,1))/sum(LMatTimeFactorNeg);
            nNeuronNegInG2    = sum(LMatTimeFactorNeg & FinalGroup(:,2))/sum(LMatTimeFactorNeg);
            
            if nNeuronPosInG1+nNeuronPosInG2 > 0
                hValue        = (nNeuronPosInG1*0+nNeuronPosInG2*1)*maxColorValue;
                sizeValue     = LMatnTime(LMatTimeFactorPos,nFactor);
                coordinates   = [coordinates; xyzTrack(LMatTimeFactorPos,:)]; %#ok<AGROW>
                colorTable    = [colorTable; hsv2rgb(hValue*ones(sum(LMatTimeFactorPos),1),...
                                                     ones(sum(LMatTimeFactorPos),1),...
                                                     ones(sum(LMatTimeFactorPos),1))];%#ok<AGROW>
                radiusTable   = [radiusTable; sizeValue]; %#ok<AGROW>
            end
            
            if nNeuronNegInG1+nNeuronNegInG2 > 0
                hValue        = (nNeuronNegInG1*0+nNeuronNegInG2*1)*maxColorValue; % 1 - (nNeuronNegInG1*0+nNeuronNegInG2*1)*0.5;
                sizeValue     = -LMatnTime(LMatTimeFactorNeg,nFactor);
                coordinates   = [coordinates; xyzTrack(LMatTimeFactorNeg,:)]; %#ok<AGROW>
                colorTable    = [colorTable; hsv2rgb(hValue*ones(sum(LMatTimeFactorNeg),1),...
                                                     ones(sum(LMatTimeFactorNeg),1),...
                                                     ones(sum(LMatTimeFactorNeg),1))];                 %#ok<AGROW>
                radiusTable   = [radiusTable; sizeValue]; %#ok<AGROW>
            end
            
        end
        nTimePoint            = (timePoints(nTime)+timePoints(nTime+1))/2;
        nImage                = sum(periodCenters<=nTimePoint);
        exportMaskRadius(['/Volumes/Zebrafish_Imaging_Data_by_Yinan/' fileName ...
            '/image_data/' flist(nImage).name], ['/Volumes/Zebrafish_Imaging_Data_by_Yinan/' fileName ...
            '/Overlay_image_data/' flist(nImage).name], coordinates, colorTable, radiusTable)
    end
end