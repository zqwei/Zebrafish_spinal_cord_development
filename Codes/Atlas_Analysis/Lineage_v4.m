%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lineage analysis 4: export features to csv files, for use with Tulip
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%


function Lineage_v4(nFile)
addpath('../Func');
setDir;
fileName       = fileNames{nFile};
fileDirName    = fileDirNames{nFile};

load([fileDirName '/', 'dev_data.mat'], 'leafID');
load([tempDatDir, fileName, '.mat'], 'slicedIndex', 'leafOrder', 'new_x', 'new_y', 'new_z', 'side', 'mnx');
load([tempDatDir, 'Leader_' fileName, '.mat'], 'activeTime', 'patternTime', 'birthtime', 'factorSize', 'neuralPlateLoc', 'birthLoc', 'divAngle');


leafID = leafID(slicedIndex);
leafID = leafID(leafOrder);

% name of csv, variable, default/NaN value
LUT = { 1:numel(leafID), 'cell_id', -1; ...
        new_x,  'AP_location',  -1; ...
        new_y,  'LR_location',  -9999; ...
        new_z,  'DV_location',  -1; ...
        side,   'side',         -1; ...
        mnx,    'mnx',          -1; ...
        birthtime, 'birthtime', -1; ...
        factorSize, 'initial_fs', -1; ...
        activeTime, 'active_time', -1; ...
        patternTime, 'pattern_time', -1; ...
        neuralPlateLoc(:, 1), 'neural_plate_AP', -9999; ...
        neuralPlateLoc(:, 2), 'neural_plate_LR', -9999; ...
        birthLoc(:, 1), 'birth_loc_AP', -1; ...
        birthLoc(:, 2), 'birth_loc_LR', -1; ...
        divAngle(:, 1), 'div_angle_theta', -1;
        divAngle(:, 2), 'div_angle_phi', -1; ...
       };
   
csvFolder = [tempDatDir '/csv_' fileName];
if ~exist(csvFolder, 'dir')
    mkdir(csvFolder);
end
for i = 1:size(LUT, 1)
    csvFile = fopen([csvFolder '/' LUT{i, 2} '.csv'], 'w');
    fprintf(csvFile, ['default, ' num2str(LUT{i, 3}) ',\n']);
    metric = LUT{i, 1};
    pID    = leafID(~isnan(metric));
    metric = metric(~isnan(metric));
    for nLine = 1:numel(pID)
        fprintf(csvFile, [num2str(pID(nLine)), ',', num2str(metric(nLine)), ',\n']);
    end
    fclose(csvFile);
end


