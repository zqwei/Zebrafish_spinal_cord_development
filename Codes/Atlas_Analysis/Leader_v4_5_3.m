%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial factor correlation analysis - Batch processing
% Determine dominance of LPA-SPA in merging events, use average trace of
% post merging
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%
function Leader_v4_5_3(fileList)
% % calculate similarity statistics and accuracy
addpath('../Func');
setDir;
leaderSpecCorr = cell(numel(fileList), 1);

tag = '_avgPost';

for i = 1:numel(fileList)
    nFile = fileList(i);
    fileName = fileNames{nFile};
    load([tempDatDir, 'LeaderPairs_PSD' tag  fileName, '.mat'], 'tagList');
    leaderSpecCorr{i, 1} = tagList;
end

accuracy = zeros(3, 2); %columns: LL, LS/SL, SS; rows: miss, total
LSAccuracy = zeros(4, 1);
pairwiseSim = []; %LPre-LPost, SPre-LPost, LPre-SPost, SPre-SPost
fishNumber = [];
for i = 1:numel(fileList)
    c = leaderSpecCorr{i, 1};
    for nFactor = 1:numel(c)
        switch sum(c(nFactor).preActLevel)
            case 1
                accuracy(2, 1) = accuracy(2, 1) + any(c(nFactor).preActLevel ~= c(nFactor).spectralCorr);
                accuracy(2, 2) = accuracy(2, 2)+1;
                LPA = find(c(nFactor).preActLevel);
                SPA = 3 - LPA;
                M = c(nFactor).distM;
                pairwiseSim = [pairwiseSim; [M(LPA, 3), M(SPA, 3)]];
                fishNumber  = [fishNumber; i];
                if any(isnan(c(nFactor).spectralCorr))
                    LSAccuracy(4) = LSAccuracy(4) + 1;
                elseif (sum(c(nFactor).preActLevel ~= c(nFactor).spectralCorr)==2)
                    LSAccuracy(2) = LSAccuracy(2)+1;
                    disp(['Disagree: dataset #' num2str(fileList(i)) ', Factor #' num2str(nFactor)]);
                elseif (sum(c(nFactor).preActLevel ~= c(nFactor).spectralCorr)==1)
                    LSAccuracy(3) = LSAccuracy(3)+1;
                    disp(['No leader: dataset #' num2str(fileList(i)) ', Factor #' num2str(nFactor)]);
                else
                    LSAccuracy(1) = LSAccuracy(1)+1;
                    %                     disp(['Correct: dataset #' num2str(fileList(i)) ', Factor #' num2str(nFactor)]);
                end
            case 0
                accuracy(3, 1) = accuracy(3, 1) + any(c(nFactor).preActLevel ~= c(nFactor).spectralCorr);
                accuracy(3, 2) = accuracy(3, 2)+1;
            case 2
                accuracy(1, 1) = accuracy(1, 1) + any(c(nFactor).preActLevel ~= c(nFactor).spectralCorr);
                accuracy(1, 2) = accuracy(1, 2)+1;
        end
    end
end

%% plot example pair similarity
example = [1, 5]; % example correlation - dataset#3, factor#5
fileName = fileNames{fileList(example(1))};
load([tempDatDir '/' fileName '.mat'],'dff', 'activeNeuronMat', 'new_x', 'timePoints');
load([tempDatDir, 'Leader_', fileName, '.mat'], 'leaderPairMetrics', 'activeTime');
leaderPairInfo = leaderPairMetrics(example(2));
leaderPair      = leaderPairInfo.neuronIndex;
appearTime      = leaderPairInfo.appearTime;
firstActiveTime = leaderPairInfo.firstActiveTime;
preActLevel     = leaderPairInfo.preActLevel;
tStart = timePoints(firstActiveTime)+1; %timePoints(firstActiveTime+1);
tMerge = timePoints(appearTime+1)+1;
tEnd   = min(tMerge + tMerge-tStart, timePoints(end));
L = find(preActLevel>=0.6);
S = find(preActLevel<0.6);
LPAPre = zscore(dff(leaderPair(L), tStart:tMerge));
SPAPre = zscore(dff(leaderPair(S), tStart:tMerge));
post = zscore(mean(dff(leaderPair, tMerge+1:tEnd), 1));
w = floor(min(numel(LPAPre), numel(post))/2)*2;
fThres = [0, 1]; % signal with frequency above fThres considered as noise
nw = 101;
[Y1pre, f] = pmtm(LPAPre, nw, w, 4);
Y2pre      = pmtm(SPAPre, nw, w, 4);
Ypost     = pmtm(post, nw, w, 4);

fRange     = f<fThres(2) & f>fThres(1);
spec = [Y1pre.*f, Y2pre.*f, Ypost.*f];
spec = spec(fRange, :)';
f = f(fRange);
spec = spec./repmat(sum(spec, 2), 1, sum(fRange));
spec = cumsum(spec, 2);
% nColors = 256;
% endColor = [76/256, 0/256, 33/256];
% midColor = [1, 1, 1];
% colorMap = zeros(nColors, 3);
% gamma = 0.6;
% for r = 1:3
%     colorMap(:, r) = linspace(midColor(r), endColor(r), nColors)'.^gamma;
% end

c = leaderSpecCorr{example(1), 1};
C = c(example(2)).distM;
if ~isempty(tag)
    C(isnan(C)) = 0;
    C = C+C';
end
C = C([L, S, 3], [L, S, 3]);
I = ones(3, 3);

mColor = lines(7);
figure,
subplot(1, 2, 1)
hold on
plot(f, spec(1, :), '--', 'color', mColor(1, :));
plot(f, spec(2, :), '--', 'color', mColor(2, :));
plot(f, spec(3, :), 'color', mColor(1, :));
% plot([0, 1], [0, 1], 'k');
xlabel('Frequency (Hz)')
ylabel('Cum. Boosted PMTM');
legend({'LPA-Pre', 'SPA-Pre', 'Post'}, 'Location', 'southeast');
ylim([0 1]);
set(gca, 'ytick', 0:0.2:1);

subplot(1, 2, 2)
imagesc(C, 'alphadata', triu(I, 1), [-0.6, 1]);
colorbar;
% colormap(colorMap)
set(gca, 'XTick', 1:4, 'XTickLabel', {'L-pre', 'S-pre', 'Post'});
set(gca, 'YTick', 1:4, 'YTickLabel', {'L-pre', 'S-pre', 'LPost'});
setPrint(32, 12, [plotDir  'InitialLeaderPairs_ExampleFactor_SpecSim' tag], 'pdf');

% pie plot for agreement level
figure, pie(LSAccuracy);
setPrint(8, 6, [plotDir  'InitialLeaderPairs_Accuracy_PSD' tag], 'pdf');


%% Pairwise distance
pd = pairwiseSim(~isnan(pairwiseSim(:, 1)), :);
figure, hold on
plot([2, 1], pd(:, 1:2), 'color', [.8, .8, .8]);
plot([2, 1], median(pd(:, 1:2)), 'or-', 'linew', 3);
set(gca, 'XTick', 1:4, 'XTickLabel', {'SPre-Post', 'LPre-Post'});
ylabel('PSD similarity');
if isempty(tag) || strcmp(tag, 'fftcorr') || strcmp(tag, '_avgPost')
    p1 = signrank(pd(:, 1), pd(:, 2), 'tail', 'right');
else
    p1 = signrank(pd(:, 1), pd(:, 2), 'tail', 'left');
end
xlim([0.5, 2.5])
title(['p=' num2str(p1)]);
set(gca, 'ytick', -0.8:0.2:1);
ylim([-1 1.2])
setPrint(16, 12, [plotDir  'PSD_Similarity' tag], 'pdf');

%% Pairwise distance - with fish color annotation and HB correction
pd = pairwiseSim(~isnan(pairwiseSim(:, 1)), :);
figure, hold on
for i = 1:size(pd, 1)
    plot([2, 1], pd(i, 1:2), 'Color', mColor(fishNumber(i), :));
end
pFish = numel(fileList, 1);
for i = 1:numel(fileList)
    pFish(i) = signrank(pd(fishNumber==i, 1), pd(fishNumber==i, 2), 'tail', 'right');
end
pFish = sort(pFish, 'descend');
pMax = 0;
for i = 1:numel(fileList)
    pMax = max(pMax, pFish(i)*i);
end

plot([2, 1], median(pd(:, 1:2)), 'ok-', 'linew', 3);
set(gca, 'XTick', 1:4, 'XTickLabel', {'SPre-Post', 'LPre-Post'});
ylabel('PSD similarity');
set(gca, 'ytick', -0.8:0.2:1);
ylim([-1 1.2])
xlim([0.5, 2.5])
title(['p1=' num2str(pMax)]);
setPrint(16, 12, [plotDir  'PSD_Similarity_fishn' tag], 'pdf');


