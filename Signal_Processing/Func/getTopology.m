function [lineage, annotation] = getTopology(trackingM, priorityChild)
% return lineage tree structures with only root, branch and leaf
% input: 
% trackingM     - swc format of 4D lineages
% priorityChild - Px1 vector if Node ID in this list, enforce it to be left child
% output: 
% lineage: N x 1 struct, each containing a separate lineage tree with fields: 
%     lineage_ID - lineage ID from trackingM(:, 10)
%     tree       - Mx5 table with: parent, child_left, child_right, leafOrder, time
%                  leafOrder: leaf level 1~N-1; parent = mean(L,R)
%                  time: time point info from trackingM(:, 8)
%     MaMuT_ID   - MaMuT ID from trackingM(:, 1), corresponding to rows topoM
% annotation: P x 2 matrix, each correspond to an entry in priorityChild
%             columns are: lineage number, leaf number (entry in lineage)
%
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%
lineageList = unique(trackingM(:, 10));
lineageList(lineageList<0) = [];
lineage = [];
annotation = nan(numel(priorityChild), 2); %lineage number, node number

for nLineage = 1:numel(lineageList)
    swc = trackingM(trackingM(:, 10)==lineageList(nLineage), :);    
    rootNode = swc(swc(:, 7)<1, 1);
    if isempty(rootNode)
        continue;
    end
%     if numel(rootNode)>1
%         disp(['lineage #' num2str(lineageList(nLineage)) ' has more than one root: ' num2str(rootNode')]);
%         continue
%     end
    numChild = hist(categorical(swc(:, 7)), categorical(swc(:, 1)));
    leafNodes = swc(numChild == 0, 1);
    branchNodes = swc(numChild == 2, 1);
    branchNodes(branchNodes==rootNode) = [];
    
    id = [leafNodes; branchNodes; rootNode];
    
    if sum(ismember(leafNodes, priorityChild))==0
        warning(['Lineage #' num2str(lineageList(nLineage)) ' does not have an anotated node']);
    end
    
    topoM = -ones(numel(id), 5);
    priorityScore = zeros(numel(id), 1);
    childList     = cell(numel(id), 1);
    for i = 1:numel(leafNodes)
        % travel from each leaf to node, add up levelScore and priority
        currNode_ID = leafNodes(i);
        priorityTag = 0;
        if ismember(currNode_ID, priorityChild)
            priorityTag = 1;
            priorityScore(i) = priorityScore(i) + 1;
            annotation(currNode_ID==priorityChild, 1) = nLineage;
            annotation(currNode_ID==priorityChild, 2) = i;
        end
        lastNode_ID = currNode_ID;
        topoM(i, 5) = swc(swc(:, 1)==currNode_ID, 8);
        try
        while currNode_ID>0
            currNode_ID = swc(swc(:, 1)==currNode_ID, 7);
            if ismember(currNode_ID, branchNodes) || currNode_ID==rootNode
                currNode = find(id==currNode_ID); % parent
                lastNode = find(id==lastNode_ID); % child
                topoM(lastNode, 1) = currNode;
                topoM(currNode, 5) = swc(swc(:, 1)==currNode_ID, 8);
                childList{currNode} =  [childList{currNode}, lastNode];
                priorityScore(currNode) = priorityScore(currNode) + priorityTag;
                lastNode_ID = currNode_ID;
            end
        end
        catch
            disp('');
        end
    end
    for i = 1:numel(id)
        children = unique(childList{i});
        if numel(children) == 2
            if priorityScore(children(1))>priorityScore(children(2))
                topoM(i, 2) = children(1); % left child has higher priority
                topoM(i, 3) = children(2);
            else
                topoM(i, 2) = children(2);
                topoM(i, 3) = children(1);
            end
        elseif numel(children) == 1
            topoM(i, 2) = children;
        elseif numel(children) > 2
            disp('node has more than 2 children');
        end
    end
    % calculate leafOrder, replace child in L->R from root
    nodeList = numel(id);
    leafSeq = numel(id);
    while ~isempty(nodeList)
        currNode = nodeList(1);
        nodeList = nodeList(2:end);
        currLoc = find(leafSeq==currNode);
        children = topoM(currNode, [2, 3]);
        children(children<0) = [];
        if ~isempty(children)
            leafSeq = [leafSeq(1:currLoc-1), children, leafSeq(currLoc+1:end)];
            nodeList = [nodeList, children];
        end
    end
    topoM(leafSeq, 4) = 1:numel(leafNodes);
    for i = numel(leafNodes)+1:numel(id)
        topoM(i, 4) = getLeafOrder(topoM, i);
    end
    currLineage.lineageID = lineageList(nLineage);
    currLineage.MaMuT_ID = id;
    currLineage.tree = topoM;
    lineage = [lineage; currLineage];
end

end

function pos = getLeafOrder(M, i)
if M(i, 2)<0 % leaf node
    pos = M(i, 4);
elseif M(i, 3)<0 % a single child
    pos = getLeafOrder(M, M(i, 2));
else % both children
    pos = (getLeafOrder(M, M(i, 2))+ getLeafOrder(M, M(i, 3)))/2;
end
end