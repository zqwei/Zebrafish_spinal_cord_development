function writeImage(imageStack, filename, varargin)

% determine file extension to select appropriate write module
[~, ~, extension] = fileparts(filename);

switch(extension)
    case '.jp2'
        
        % set all options to default settings
        nThreads = min(feature('numcores'), 8);
        compression = 0;
        
        % read optional input parameters
        if ~isempty(varargin)
            for c = 1:2:length(varargin)
                switch lower(varargin{c})
                    case {'nthreads'}
                        nThreads = varargin{c + 1};
                    case {'compression'}
                        compression = varargin{c + 1};
                    otherwise
                        error(['Invalid optional argument for JP2 writer: ' varargin{c}]);
                end; % switch
            end; % for
        end; % if
        
        writeJP2stack(imageStack, filename, nThreads, compression);
        
    case {'.tif', '.tiff'}
        
        writeTIFstack(imageStack, filename, 2^31);
        
    case '.klb'
        
       % set all options to default settings
        nThreads = feature('numcores');
        pixelSize = [];
        blockSize = [];
        compressionType = [];
        metadata = [];
        
        %read optional input parameters
        if ~isempty(varargin)
            for c = 1:2:length(varargin)
                switch lower(varargin{c})
                    case {'numthreads'}
                        nThreads = varargin{c + 1};
                    case {'blocksize'}
                        blockSize = varargin{c + 1};
                    case {'pixelsize'}
                        pixelSize = varargin{c + 1};
                    case {'compressiontype'}
                        compressionType = varargin{c + 1};
                    case {'metadata'}
                        metadata = varargin{c + 1};
                    otherwise
                        error(['Invalid optional argument for KLB writer: ' varargin{c}]);
                end; % switch
            end; % for
        end; % if
        
        writeKLBstack(imageStack, filename, nThreads, pixelSize, blockSize, compressionType, metadata);
        
    otherwise
        error 'Format not recognized.'
end;