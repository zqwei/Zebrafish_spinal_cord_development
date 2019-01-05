function imageStack = readImage(filename, nThreads)

if exist(filename, 'file') == 0
    error 'File does not exist.'
end;

% determine file extension to select appropriate reader module
[~, ~, ext] = fileparts(filename);

switch(ext)
    case '.jp2'
        if nargin < 2
            nThreads = min(feature('numcores'), 8);
        end;
        imageStack = readJP2stack(filename, nThreads);
    case {'.tif', '.tiff'}
        imageStack = readTIFstack(filename);
    case '.klb'
        if nargin < 2
            nThreads = feature('numcores');
        end;
        imageStack = readKLBstack(filename, nThreads);
    otherwise
        error 'Format not recognized.'
end;