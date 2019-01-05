%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lineage analysis 5: sibling statistics
%
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%

function Lineage_v5(datasets)
addpath('../Func');
setDir;

metricList = {'birthtime/(1.5*60)', 0,  'birth time';... %1
    'abs(neuralPlateLoc(:, 2))', 2, 'neural plate ML';...%2
    'neuralPlateLoc(:, 1)', 2, 'neural plate AP';...%3
    'divAngle(:, 2)', 2, 'division angle AP'; ...%4
    'birthLoc(:, 1)', 2, 'birth location AP'; ...%5
    'birthLoc(:, 2)', 2, 'birth location LR'; ...%6
    'siblingCloneSize', 0, 'sibling clone size'; ...%7
    'activeTime', 1, 'activation time'; ...%8
    'patternTime', 1, 'patterned time'; ...%9
    'preActLevel', 1, 'pre-pattern activity'; ...%10
    'new_x', 2, 'AP location'; ...%11
    'new_y', 2, 'LR location'; ...%12
    'new_z', 2, 'DV location'; ...%13
    'mnx', 3 , 'mnx'; ...%14
    };
nMetrics = size(metricList, 1);
metricsAll = [];
EVsimAll = [];  % 7 columns: 
                % 1 - similarity pair, 
                % 2 - similarity ctrl_mean
                % 3 - similarity ctrl_std
                % 4 - p-value for ranksign
                % 5 - rank percentile in control, 
                % 6 - number of pairs in control
                % 7 - ipsilateral(1) or contralateral(0) pair
EVctrlAll = cell(0, 1);
vecInd = 1;
plotID = 0;
figure,
for nFish = 1:numel(datasets)
    nFile = datasets(nFish);
    fileName          = fileNames{nFile};
    fileDirName       = fileDirNames{nFile};
    
    % load and calculate all metrics
    load([fileDirName, '/profile.mat'], 'siblings');
    load([tempDatDir, fileName, '.mat'], 'new_x', 'new_y', 'new_z', 'mnx', 'leafOrder', 'slicedIndex', 'side');
    load([tempDatDir, 'Leader_' fileName, '.mat']);
    load([tempDatDir, 'EV_' fileName, '.mat'], 'EVMat');
    
    EVMat(isnan(EVMat) | EVMat<0.05) = 0;
    ori_id = 1:numel(slicedIndex);
    ori_id = ori_id(slicedIndex);
    ori_id = ori_id(leafOrder);
    EVsimMat = corr(EVMat');
    EVsim    = nan(size(siblings, 1), 6);
    EVprctle = nan(size(siblings, 1), 6);
    
    for i = 1:size(siblings, 1)
        siblings(i, 1) = find(siblings(i, 1)==ori_id);
        siblings(i, 2) = find(siblings(i, 2)==ori_id);
    end
    
    metrics = nan(size(siblings, 1), 2, nMetrics+1);
    for i = 1:size(metricList, 1)
        metric = eval(metricList{i, 1});
        metrics(:, :, i) = metric(siblings);
    end
    %     for i = 1:size(siblings, 1)
    %         control = abs(birthtime - birthtime(siblings(i, 1))) < 20/1.5;
    %         EVcontrol = tril(EVsimMat(control, control), -1);
    %         idxTrue = tril(true(size(EVcontrol)), -1);
    %         EVcontrolVector = EVcontrol(idxTrue(:));
    %         EVsim(i, 1) = EVsimMat(siblings(i, 1), siblings(i, 2));
    %         EVcontrolVector(EVcontrolVector == EVsim(i, 1)) = [];
    %
    %         EVsim(i, 2) = nanmean(EVcontrolVector);
    %         EVsim(i, 3) = nanstd(EVcontrolVector);
    %         if ~isnan(EVsim(i, 1))
    %             EVsim(i, 4) = signrank(EVcontrolVector, EVsim(i, 1), 'tail', 'left');
    %         end
    %         EVsim(i, 5) = side(siblings(i, 1));
    %         EVsim(i, 6) = side(siblings(i, 2));
    %         EVsim(i, 7) = sum(control);
    %     end
    
    for i = 1:size(siblings, 1)
        control = abs(birthtime - birthtime(siblings(i, 1))) < 60/1.5;
        EVcontrol = tril(EVsimMat(control, control), -1);
        idxTrue = false(size(EVcontrol));
        
        idxTrue(side(control)==side(siblings(i, 1)), side(control)==side(siblings(i, 2))) = true;
        idxTrue(side(control)==side(siblings(i, 2)), side(control)==side(siblings(i, 1))) = true;
        %         idxTrue(side==side(sibling(i, 2)), side==side(sibling(i, 1))) = true;
        idxTrue = tril(idxTrue, -1);
        EVcontrolVector = EVcontrol(idxTrue(:));
        EVsim(i, 1) = EVsimMat(siblings(i, 1), siblings(i, 2)); %similarity index
        EVcontrolVector(EVcontrolVector == EVsim(i, 1) | isnan(EVcontrolVector)) = [];
        
        EVsim(i, 2) = nanmean(EVcontrolVector); % average similarity of control
        EVsim(i, 3) = nanstd(EVcontrolVector);  % std similarity of control
        if ~isnan(EVsim(i, 1)) && ~isempty(EVcontrolVector)
            EVsim(i, 4) = signrank(EVcontrolVector, EVsim(i, 1)); % p value of sib vs. control
            EVsim(i, 5) = sum(EVsim(i, 1)>EVcontrolVector)/numel(EVcontrolVector);
        end
        EVsim(i, 6) = numel(EVcontrolVector); % number of pairs to control
        EVsim(i, 7) = side(siblings(i, 1)) == side(siblings(i, 2)); %ipsilateral pair(1), or contraleteral pair(0)
        EVctrlAll{vecInd, 1} = EVcontrolVector;
        vecInd = vecInd+1;
        
        if EVsim(i, 7)==1
            plotID = plotID+1;
            t = (1:size(EVMat, 2))/60;
            subplot(3, 3, plotID)
            hold on
            currCtrl =  abs(birthtime - birthtime(siblings(i, 1))) < 60/1.5;
            currCtrl(siblings(i, 1)) = false;
            currCtrl(siblings(i, 2)) = false;
            h3 = plot(t, EVMat(currCtrl, :)', 'color', [.8, .8, .8]);
            h1 = plot(t, EVMat(siblings(i, 1), :)', 'r');
            h2 = plot(t, EVMat(siblings(i, 2), :)', 'g');
            xlim([0, t(end)]);
            ylim([0, 1]);
            xlabel('Time (hour)');
            ylabel('EV');
%             legend([h1; h2; h3], {'sibling 1', 'sibling 2', 'control'});
            title(['dataset #' num2str(nFish) ' sibling #' num2str(i)]);
            hold off
        end
    end
    
    metrics(:, :, end) = nFish;
    metricsAll = [metricsAll; metrics];
    EVsimAll = [EVsimAll; EVsim];
    
end

setPrint(8*3, 6*3, [plotDir '/sibling_EV_curves'], 'pdf');

figure,
for i = 8:12
    subplot(1, 5, i-7);
    gscatter(metricsAll(:, 1, i), metricsAll(:, 2, i), metricsAll(:, 1, end));
    xlabel('sibling 1');
    ylabel('sibling 2');
    legend off
    title(metricList{i, 3});
    hold on,
    axLimU = nanmax(nanmax(metricsAll(:, :, i)));
    axLimL = min(0, nanmin(nanmin(metricsAll(:, :, i))));
    plot([axLimL, axLimU], [axLimL, axLimU], 'k--');
    xlim([axLimL, axLimU]);
    ylim([axLimL, axLimU]);
    box off
end
setPrint(40, 6, [plotDir,  'SiblingMetricGroupScatterPlot'], 'pdf');

ipSibList = EVsimAll(EVsimAll(:, 7)==1 & ~isnan(EVsimAll(:, 1)), :);
ipCtrlVec = EVctrlAll(EVsimAll(:, 7)==1 & ~isnan(EVsimAll(:, 1)));
figure,     hold on
for i  = 1:size(ipSibList, 1)
    currCtrl = ipCtrlVec{i};
    scatter(i+(rand(numel(currCtrl), 1)-0.5)*0.5, currCtrl, [], [.8, .8, .8], 'filled');
    scatter(i, ipSibList(i, 1), 'r', 'filled');
    if ipSibList(i, 4)<0.05 && ipSibList(i, 1)>ipSibList(i, 2)
        text(i, 1.01, '*');
    end
end
set(gca, 'xTick', 1:size(ipSibList, 1));
ylim([0, 1.1]);
xlabel('Sibling #');
ylabel('EV history similarity score');
setPrint(8, 6, [plotDir,  'ipsiSiblingSimMetric'], 'pdf');


