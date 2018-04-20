function trackingM_splitted = splitTracksbyID(trackingM, leaf_id, lineage_id)
% split the lineages into single tracks starting from leaf_id backwards
% new Track_ID will become leaf_id, and sister cells will share the same
% trajectory
% lienage_id: final ID of the nodes corresponding to each leaf_id
trackingM_splitted = [];
for i = 1:numel(leaf_id)
    currentLeaf = leaf_id(i);
    if ~isnan(lineage_id(i))
        while currentLeaf > 0
            currentNode = trackingM(trackingM(:,1) == currentLeaf, :);
            currentNode(10) = lineage_id(i);
            trackingM_splitted = [trackingM_splitted;  currentNode];
            currentLeaf = trackingM(trackingM(:,1) == currentLeaf, 7);
        end
    end
end