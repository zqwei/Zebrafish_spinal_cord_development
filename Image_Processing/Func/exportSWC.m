%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% exportSWC.m
% return tracks in swc format: 
% id, type, x, y, z, radius, pid, time, confidence, lineage_id
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%
function matrix = exportSWC(tracks, timepoints, varargin)
flagKill = 0;
if ~isempty(varargin)
    flagKill = 1;
    terminationTimepoint = varargin{1};
end

nTimepoints = size(tracks, 1);
nCells = size(tracks{1}, 1);


matrix = zeros(nCells*nTimepoints, 10);
matrix(:, 6) = 1;


currentIds = 1:nCells;
matrixLabel = false(nCells*nTimepoints, 1);
for i = 1:nTimepoints
    matrix(currentIds, 1) = currentIds;
    matrix(currentIds, 3:5) = tracks{i};
    matrix(currentIds, 5) = matrix(currentIds, 5);
    if i==1
        matrix(currentIds, 7) = -1;
    else
        matrix(currentIds, 7) = currentIds - nCells;
    end
    matrix(currentIds, 8) = timepoints(i);
    matrix(currentIds, 10) = 1:nCells;
    if flagKill
        matrixLabel(currentIds(terminationTimepoint<=i)) = 1;
    end
    currentIds = currentIds + nCells;
end
matrix(matrixLabel, :) = [];

