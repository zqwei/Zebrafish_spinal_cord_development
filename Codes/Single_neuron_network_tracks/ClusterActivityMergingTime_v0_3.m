%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LMatCorrection_v_0_1
% PreLMatTracker_v_0_1
% Spectrogram_v0_3
% Spectrogram_v0_4
%
% Based on the merging cluster identified by ClusterActivityMergingTime_v0_0
% summary of
% ClusterActivityMergingTime_v0_1
% ClusterActivityMergingTime_v0_2
%
% -------------------------------------------------------------------------
% 
% Ziqiang Wei
% weiz@janelia.hhmi.org
%
%

addpath('../Func');
setDir;    

fileToAnalysis = [3 7 10 16]; % 

for nFile = fileToAnalysis
    ClusterActivityMergingTime_v0_1(nFile)
    ClusterActivityMergingTime_v0_2(nFile)
end

markerStyle = {'', '', 'o', '', '', '', 's', '', '', '^', '', '', '', '', '', '>'};
mColor         = cbrewer('qual', 'Dark2',  32, 'cubic');

figure

for nFile = fileToAnalysis

    fileName   = fileNames{nFile};   %#ok<*USENS>
    load([tempDatNetDir, 'LocalCommunityMerging_' fileName, '_v_0_1.mat'], 'corr_prepost_mat', 'corr_mat')
    subplot(1, 2, 1)
    hold on
    scatter(corr_mat(:,1), corr_mat(:,3), corr_mat(:,2)*200, corr_mat(:,4), markerStyle{nFile}, 'filled', 'MarkerFaceAlpha', 0.7)
    
    subplot(1, 2, 2)
    hold on
    scatter(corr_prepost_mat(:,1), corr_prepost_mat(:,2), [], corr_prepost_mat(:,3), markerStyle{nFile}, 'filled', 'MarkerFaceAlpha', 0.7)
    
end

subplot(1, 2, 1)
set(gca, 'TickDir', 'out')
caxis([0 1])
xlim([0 1])
ylim([0 1])
refline(1)
xlabel('Local community corr.')
ylabel('Global community corr.')
box off

subplot(1, 2, 2)
caxis([0 0.2])
xlim([0 0.4])
ylim([0 0.4])
refline(1)
xlabel('Local community xcorr.')
ylabel('Global community xcorr.')
box off

setPrint(8*2, 6, [plotNetDir 'LocalCommunityMergingCorr_summary'], 'pdf')