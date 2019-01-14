addpath('../Data_Analysis/')
lineage_datasets = [6, 11, 24];
control_datasets = [3, 4, 7, 12, 10, 15, 16];
islet_datasets = [10 13 15 16];

%lineage datasets
for nFile = lineage_datasets
    disp(['processing dataset #' num2str(nFile)]);
    FAEV_v0_1(nFile); % estimate firstActiveTime
    Leader_v0(nFile, 10); % calculate activeTime, patternTime, birthtime
    Leader_v1_2(nFile); % plot leader metrics on polar plot (dep)
    Leader_v4_1(nFile);  % initial pair composition analysis
    Leader_v4_2(nFile); % neuron-wise preActLevel and prePowerLevel
    Lineage_v1(nFile); % plot simple birth-func matrix
    Lineage_v3(nFile); % calculate lineage tree info from trackingM
    Lineage_v3_1_1(nFile); % calculate dev statistics for each cell
    Lineage_v1_2_1(nFile); % plot extended birth-func matrix
    Lineage_v3_1(nFile); % calculate dev statistics for each cell
    Lineage_v1_2(nFile); % plot extended birth-func matrix
    Lineage_v3_2(nFile); % plot lineage tree
end
% % 
Lineage_v1_3(lineage_datasets) % conserved correlation matrix
Lineage_v2_2(lineage_datasets); %group scattered plot
Lineage_v5(lineage_datasets); % sibling analysis
Lineage_v2_1(lineage_datasets); %segmental birth-act, summary

% visualize lineage trees
% visOptions = {'ActMat', 'EVMat', 'patternTime'};
visOptions = {'EVMat'};
for nFile = lineage_datasets
    for nVis = 1:numel(visOptions)
        Lineage_v3_3(nFile, visOptions{nVis});
    end
end

% leader analysis
for nFile = control_datasets
    disp(['processing dataset #' num2str(nFile)]);
    Leader_v4_1(nFile); % initial pair composition analysis
    Leader_v4_2_1_2(nFile); % leader pair IEI distribution
    Leader_v4_2_3(nFile); % leader pair power spectrum analysis
end
preActThres = Leader_v4_3(control_datasets);% determine preActLevel threshold

preActThres = 0.6;
Leader_v5(control_datasets, preActThres); % summarize initail pair L-N composition stacked bar plot
Leader_v4_5_1(control_datasets);% summarize IEI similarity
Leader_v4_5_3(control_datasets); % summarize spectral correlation
Leader_v5_2(control_datasets, preActThres); % initial pair L-N composition errorbar plot summarize fish stats

Leader_v2_1(control_datasets, preActThres);% summarize leader seg location on preActLevel

% islet analysis
for nFile = islet_datasets
    Leader_v0(nFile, 20);
end
Leader_v2_2(islet_datasets, preActThres); % bar plot leader identity with islet and mnx
Leader_v2_3(islet_datasets); % violin plot leader identity with islet and mnx