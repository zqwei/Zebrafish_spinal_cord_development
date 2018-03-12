%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot active time, pattern time with cell type info
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
% 

function Leader_v3(nFile)
addpath('../Func');
setDir;
fileName          = fileNames{nFile};

load([tempDatDir, fileName, '.mat'], 'mnx');
load([tempDatDir, 'Leader_', fileName, '.mat'], 'activeTime', 'patternTime');
figure, hold on 
plot(activeTime(mnx>0), patternTime(mnx>0), 'ok');
plot(activeTime(mnx==0), patternTime(mnx==0), 'sr');
plot([0, nanmax([patternTime, activeTime])], [0, nanmax([patternTime, activeTime])], 'k--');
xlabel('Active time (hour)')
ylabel('Pattern time (hour)');
setPrint(8, 6, [plotDir,  'MNXActivePatternTime_' fileName], 'pdf');
close
