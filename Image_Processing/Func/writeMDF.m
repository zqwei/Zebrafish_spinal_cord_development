%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% writeMDF.m
% 
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%
function writeMDF(filename, pointList)
nPoints = size(pointList, 1);
fid = fopen(filename, 'w');
fprintf(fid, 'MTrackJ 1.5.0 Data File\n');
fprintf(fid, 'Displaying true true true 1 2 0 4 100 10 0 0 0 2 1 12 0 true true true false\n');
fprintf(fid, 'Assembly 1 FF0000\n');
fprintf(fid, 'Cluster 1 FF0000\n');
for i = 1:nPoints
    fprintf(fid, ['Track ' num2str(i) ' FF0000 true\n']);
    fprintf(fid, ['Point 1 ' num2str(pointList(i, 1), '%.1f') ' ' num2str(pointList(i, 2), '%.1f') ' 1.0 ' num2str(pointList(i, 3), '%.1f') ' 1.0\n']);
end

fprintf(fid, 'End of MTrackJ Data File\n');
fclose(fid);