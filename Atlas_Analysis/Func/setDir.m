%% Data Folders
% data with functional imaging + atlas annotation
DirNames{3} = '../Data/20150410';
DirNames{4} = '../Data/20150417';
DirNames{5} = '../Data/20160312';
DirNames{7} = '../Data/20161004';
DirNames{8} = '../Data/20161026';

DirNames{12} = '../Data/20170112';


% data with functional imaging + dev info (birth time) "birthtime.mat"
DirNames{6} = '../Data/20160328';
DirNames{9} = '../Data/20161027'; % missing birthtime
DirNames{11} = '../Data/20170111'; % missing birthtime
DirNames{14} = '../Data/20170201'; % missing birthtime

% data with functional imaging + IHC info "islet.mat"
DirNames{10} = '../Data/20161028';
DirNames{13} = '../Data/20170126';
DirNames{15} = '../Data/20170202';
DirNames{16} = '../Data/20170216';


% data with functional imaging but no atlas information
DirNames{1} = '../Data/20140818';
DirNames{2} = '../Data/20141006';

% data with MO_slc6a9 injection
DirNames{17} = '../Data/20170315';
DirNames{18} = '../Data/20170323';
DirNames{19} = '../Data/20170412';

%% Plot Folders
PlotDir = '../Plot';
if ~exist(PlotDir, 'dir')
    mkdir(PlotDir);
end

%% TempData Folder
TempDataDir = '../TempData';
if ~exist(TempDataDir, 'dir')
    mkdir(TempDataDir);
end
%% Dataset Name
dataset = DirNames;
for i = 1:length(DirNames)
    dataset{i} = DirNames{i}(9:end);
end
