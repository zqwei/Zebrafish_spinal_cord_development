%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% readMDF.m
% 
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%
function [points, varargout] = readMDF(inputFile)

fid = fopen(inputFile);

lineIdx = 0;
headerIdx = 0;
colorHeader = cell(lineIdx, 1);
fileHeader = cell(headerIdx, 1);
points = zeros(lineIdx, 3);
tline = fgetl(fid);
while ischar(tline)
    words = strsplit(tline, ' ');
    if (strcmp(words{1}, 'Point'))
        points(lineIdx, 1) = str2double(words{3});
        points(lineIdx, 2) = str2double(words{4});
        points(lineIdx, 3) = str2double(words{6});
    elseif (strcmp(words{1}, 'Track'))
        lineIdx = lineIdx + 1;
        colorHeader{lineIdx} = tline;
    else
        headerIdx = headerIdx + 1;
        fileHeader{headerIdx} = tline;
    end
    tline = fgetl(fid);
end

fclose(fid);

if (nargout>1)
    varargout{1} = fileHeader;
    varargout{2} = colorHeader;
end