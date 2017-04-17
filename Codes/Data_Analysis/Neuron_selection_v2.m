% 
% Selecting neurons and storing their raw data, baseline, dff, tracks
% 
% -------------------------------------------------------------------------
%
% add new x, y, z (Yinan's new coordinates) to TempDat File
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function Neuron_selection_v2(nFile)
    % load data
    addpath('../Func');
    setDir;
    fileDirName   = fileDirNames{nFile}; %#ok<USENS>
    fileName      = fileNames{nFile}; %#ok<USENS>
    dirImageData  = [fileDirName '/'];
    load ([dirImageData 'profile.mat'], 'x', 'y', 'z')
    
    if ~exist('x', 'var')
        return;
    end    
    
    if ~exist([tempDatDir, fileName, '.mat'], 'file')
        return;
    end
    
    load([tempDatDir, fileName, '.mat'], 'leafOrder', 'slicedIndex');
    
    new_x         = x(slicedIndex);
    new_y         = y(slicedIndex);
    new_z         = z(slicedIndex);
    
    new_x         = new_x(leafOrder);
    new_y         = new_y(leafOrder);
    new_z         = new_z(leafOrder);
    
    save([tempDatDir, fileName, '.mat'], 'new_x', 'new_y', 'new_z', '-append')
    
end