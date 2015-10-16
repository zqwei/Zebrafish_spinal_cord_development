% 
% Data Analysis List
% 
% -------------------------------------------------------------------------
% version 1.0
% some missing data existing for the initially generated dataset.
%
% version 1.1
% using only the orginal dataset.
%
% Analysis list
%
% 1.  Covariance analysis: factor analysis
% 2.  Time series analysis: LDS analysis
% 3.  Cluster in the LDS connectivity matrix
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.  Covariance analysis: clustering analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data_Analysis_List1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.  Covariance analysis: factor analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data_Analysis_List2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.  Time series analysis: LDS analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data_Analysis_List3
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4.  Cluster in the LDS connectivity matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data_Analysis_List4

function Data_Analysis_List0_1(nFile)

    addpath('../Func');

    % load data
    setDir;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 0.1 Loading data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % dirResult         = [fileDirName, '/Fitted_result/'];
    fileDirName       = fileDirNames{nFile};
    fileName          = fileNames{nFile}; %#ok<USENS>

    dirImageData      = [fileDirName '/']; % [fileDirName, '/Data/'];
    load ([dirImageData 'dff.mat']);
    load ([dirImageData 'clustering.mat'])
    
    dffWithNoDuplicates       = dff;  %#ok<NASGU>
    tracksWithNoDuplicates    = tracks_smoothed;  %#ok<NASGU>
    idxWithNoDuplicates       = idx; % make idx = [1 2 3] % 1: left, 2: right, 3: non-selective
    unique_idx                = unique(idxWithNoDuplicates);
    if isequal(unique_idx, [0 1 2]')
        idxWithNoDuplicates (idxWithNoDuplicates==0) = 3;    %#ok<NASGU>
    end

    save([tempDatDir, fileName, '.mat'], 'dffWithNoDuplicates','tracksWithNoDuplicates','idxWithNoDuplicates');
end