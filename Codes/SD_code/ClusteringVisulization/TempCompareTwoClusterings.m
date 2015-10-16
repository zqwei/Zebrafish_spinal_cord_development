function [bestOverlapForOne,bestFitForOne] = ...
  TempCompareTwoClusterings(clusteringOneIndCell,clusteringTwoIndCell)
%   Outputs should be:
%   (adjusted) Rand index
%   Which of the clusters in the first clustering is most similar to a
%   given cluster in the second

clusteringOneNum = length(clusteringOneIndCell);
clusteringTwoNum = length(clusteringTwoIndCell);


% I am assuming for now we have the same amount of clusters
bestFitForOne = zeros(clusteringTwoNum,1);
bestOverlapForOne = zeros(clusteringTwoNum,1);
clusterSize = zeros(clusteringOneNum,1);

for ii=1:clusteringOneNum
  curClustOneInd = clusteringOneIndCell{ii};
  clusterSize(ii) = length(curClustOneInd);
  matchVec = zeros(clusteringTwoNum,1);
  for jj=1:clusteringTwoNum
    commonInds = intersect(clusteringOneIndCell{ii},clusteringTwoIndCell{jj});
    inBothNum = length(commonInds)./clusterSize(ii);
    missingNum = length(setdiff(clusteringOneIndCell{ii},commonInds))./clusterSize(ii);
    
    matchVec(jj) = inBothNum - 0.1*missingNum;
  end
  [bestOverlapForOne(ii),bestFitForOne(ii)] = max(matchVec);
end

end

% bestFitForTwo = zeros(clusterNum,1);
% bestOverlapForTwo = zeros(clusterNum,1);
% 
% for ii=1:clusterNum
%   matchVec = zeros(clusterNum,1);
%   curClustTwoLogical = clusteringTwoIndCell==ii;
% 
%   for jj=1:clusterNum
%     inBothNum = sum(curClustTwoLogical.*(clusteringOneIndCell == jj))./sum(curClustTwoLogical);
%     missing = sum((clusteringOneIndCell == jj).*(clusteringTwoIndCell ~= ii))./sum(clusteringOneIndCell==jj);
%     matchVec(jj) = inBothNum - missing;
%   end
%   [bestOverlapForTwo(ii),bestFitForTwo(ii)] = max(matchVec);
% end



% compareVec = zeros(unitNum,1);
% for ii=1:clusterNum
%   compareVec(clusterTwoInd==ii) = bestFitForOne(ii); 
% end
% 
% clusterColors = {'k','r','c','g','b','y'};
% 
% [sortedClusterSize,sortedClusterInd] = sort(clusterSize,'descend');
% 
% 
% figure; hold on;
% indVec = [];
% ind = 0;
% 
% firstLineHeight = 1;
% 
% for ii=1:clusterNum
%   indVec = [indVec sort(find(clusterOneInd==sortedClusterInd(ii)))'];
%   line([ind ind+sortedClusterSize(ii)],[firstLineHeight firstLineHeight],...
%     'color',clusterColors{ii},'lineWidth',3);
%   ind = ind+sortedClusterSize(ii);
%   bestFitForTwo(ii) = find(sortedClusterInd == bestFitForTwo(ii));
% end
% 
% textHeight = 0;
% secondLineHeight = 2;
% textShift = 0.75;
% for ii=1:unitNum
%   equivalentCluster = bestFitForTwo(clusterTwoInd(indVec(ii)));
% %   text(ii-textShift,textHeight,num2str(indVec(ii)))
%   text(ii-textShift,textHeight,lineNames{indVec(ii)},'rotation',90,'interpreter','none');
% 
%   line([ii-1 ii],[secondLineHeight secondLineHeight],...
%     'color',clusterColors{equivalentCluster},'lineWidth',3);
%   line([ii ii],[textHeight secondLineHeight],'linestyle','--','color','k','linewidth',0.25)
% end
% 
% axis([0 unitNum+1 textHeight-1 secondLineHeight+1])
% 
% axis off;
% end