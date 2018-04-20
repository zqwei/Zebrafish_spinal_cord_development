function lineage = getSibling(trackingM, leafNode, birth)
% return sibling sub-tree of leafNode and division triplets ID in trackingM
% siblingTrackingM: sub-tree of sibling in swc format;
% divisionTriplet: [Mother, daughter1(gives birth to leafNode), daughter2(gives birth to sibling)]
% siblingCloneSize: size of the sibling lineage tree

% find birthNode from leafNode bottom-up
currentNode = leafNode;
% childNode   = leafNode;
currentTime = trackingM(trackingM(:, 1)==leafNode, 8);
sibling     = NaN;
while currentNode>0
    parentNode = trackingM(trackingM(:, 1)== currentNode, 7);
    currentTime = trackingM(trackingM(:, 1)== parentNode, 8);
    childList   = trackingM(trackingM(:, 7)==parentNode, 1);
    if numel(childList)>1
        sibling = childList(childList~=currentNode);
        break;
    end
%     childNode   = currentNode;
    currentNode = parentNode;
%     disp(['parentNode=' num2str(parentNode) ', currentTime=' num2str(currentTime)]);
end
if (currentTime ~= birth)
    disp(['Detected last division at ' num2str(currentTime) ', actual birth time at ' num2str(birth)]);
end
lineage.divisionTriplet = [parentNode, currentNode, sibling];

% find all sub-trees of sibling
siblingLeafList = [];
currentNodeList = sibling;
visitTag    = false(size(trackingM, 1), 1);
while ~isempty(currentNodeList)
    currentNode = currentNodeList(1);
    visitTag(trackingM(:, 1)==currentNode) = 1;
    currentNodeList(1) = [];
    childList = trackingM(trackingM(:, 7)==currentNode, 1);
    if isempty(childList)
        siblingLeafList = [siblingLeafList, currentNode];
%         disp(['leaf node in sibling tree: ' num2str(currentNode) ', timePoint= ' num2str(trackingM(trackingM(:, 1)==currentNode, 8))]);
    end
    currentNodeList = [currentNodeList; childList(:)];
end
siblingTrackingM = trackingM(visitTag, :);
siblingTrackingM(siblingTrackingM(:, 1)==sibling, 7) = -1;
lineage.siblingTrackingM = siblingTrackingM;
lineage.siblingLeafList = siblingLeafList;

end