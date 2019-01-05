% plot lineage tree graph from left to right, with annotation node
% highlighted by colors of "colorGroup" and "fillGroup"
%
% topoM: 1.parent, 2.left_child, 3.right_child, 4.leafOrder, 5.time 
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%
function plotLineageTree(trees, annotation, colorGroup, fillGroup)
treeSpacing = 3;
nodeSpacing = 2;

xOffset = 0;
figure, hold on
set(gca, 'YDir', 'Reverse');
nodes2Draw = nan(size(annotation, 1), 2);%x, y
for nTree = 1:numel(trees)
    topoM = trees(nTree).tree;
    for i = 1:size(topoM, 1)-1
        % draw vertical line to parent
        nodeX = xOffset + (topoM(i, 4)-1)*nodeSpacing;
        plot([nodeX, nodeX], [topoM(i, 5), topoM(topoM(i, 1), 5)], 'k-');      
        % if not leaf, draw horizontal lines
        if topoM(i, 3)>0
            leftX  = xOffset + (topoM(topoM(i, 2), 4)-1)*nodeSpacing;
            rightX = xOffset + (topoM(topoM(i, 3), 4)-1)*nodeSpacing;
            plot([leftX, rightX], [topoM(i, 5), topoM(i, 5)], 'k');
        end
        % if i is annotation point, add it to list
        noteEntry = find(annotation(:, 1)==nTree & annotation(:, 2)==i, 1);
        if ~isempty(noteEntry)
            nodes2Draw(noteEntry, 1:2) = [nodeX, topoM(i, 5)];
        end     
    end
    xOffset = xOffset + max(topoM(:, 4)-1)*nodeSpacing + treeSpacing;
end
xlim([-treeSpacing, xOffset+treeSpacing]);
axis off

scatter(nodes2Draw(fillGroup>0, 1), nodes2Draw(fillGroup>0, 2), [], colorGroup(fillGroup>0), 's', 'filled', 'MarkerEdgeColor', 'flat');
scatter(nodes2Draw(fillGroup==0, 1), nodes2Draw(fillGroup==0, 2), [], colorGroup(fillGroup==0), 's');
colormap jet