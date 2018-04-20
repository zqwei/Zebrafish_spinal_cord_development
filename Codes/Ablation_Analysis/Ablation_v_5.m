%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5  (double cut) evaluate spatio-temporal distribution of stable factors
%
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%

function Ablation_v_5(fishList)
addpath('../Func');
setDir;

secDistriList = zeros(numel(fishList), 3); % each side of each fish, relative length of factors
for i = 1:numel(fishList);
    nFile = 24+fishList(i)*2;
    fileDirName  = fileDirNames{nFile};
    fileName          = fileNames{nFile}; %#ok<USENS>

    load([fileDirName '/', 'profile.mat'], 'segAblation_A', 'segAblation_P');
    load([tempDatDir, fileName, '.mat'], 'new_x', 'new_y'); 
    load([tempDatDir, 'LONOLoading_' fileName, '_v2.mat'], 'factorComp', 'factorSizes');
    domSecVote = zeros(3, 1); %AP-M 
%     domSecVote = zeros(3, 2); %AP-M x LR
    for nFactor = 1:size(factorSizes, 2)
        currentFactor = factorComp(:, nFactor);
        for nTime = 1:numel(currentFactor)
            if isempty(currentFactor{nTime})
                continue;
            end
            xFac = new_x(currentFactor{nTime});
            sideFac = unique(new_y(currentFactor{nTime})>0);
            voteAMP = histcounts(xFac, [-Inf, segAblation_A(2), segAblation_P(1), Inf]);
            [~, domSec] = max(voteAMP);
            if numel(domSec)==1
                domSecVote(domSec) = domSecVote(domSec)+1;
%                 domSecVote(domSec, sideFac+1) = domSecVote(domSec, sideFac+1)+1;
%                 domSecVote(mod(domSec, 2)+1, sideFac+1) = domSecVote(mod(domSec, 2)+1, sideFac+1)+1; %dom in M->1, dom in AP->2
            end
        end
    end
    secDistriList(i, :) = domSecVote';
%     secDistriList(2*i-1:2*i, :) = domSecVote';
end
lifeTime_M = secDistriList(:, 2);
lifeTime_AP = max(secDistriList(:, [1, 3]), [], 2);
figure, scatter(lifeTime_M, lifeTime_AP, 'filled');
hold on
plot([0, 45], [0, 45], '--k');
xlabel('Life time of communities in M (min)');
ylabel('Life time of communities in A/P (min)');
xlim([0 45]);
ylim([0 45]);
box off
setPrint(8, 6, [plotDir '/DoubleCutFactorLife'], 'pdf');