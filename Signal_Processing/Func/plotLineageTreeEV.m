% plot lineage tree graph from left to right, with annotation node
% highlighted by colors of "colorGroup" and "fillGroup"
% at the bottom of the trees, plot activation history in bars
%
% topoM: 1.parent, 2.left_child, 3.right_child, 4.leafOrder, 5.time 
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%
function plotLineageTreeEV(trees, annotation, colorGroup, fillGroup, EVMat)
treeSpacing = 3;
nodeSpacing = 2;

xOffset = 0;
figure, hold on
set(gca, 'YDir', 'Reverse');
nodes2Fill = nan(size(annotation, 1), 2);%x, y
nodes2Draw = [];%x, y
endTime = trees(1).tree;
endTime = max(endTime(:, 5));


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
            nodes2Fill(noteEntry, 1:2) = [nodeX, topoM(i, 5)];
        elseif topoM(i, 5)==endTime && topoM(i, 3)<=0
            nodes2Draw = [nodes2Draw; nodeX, topoM(i, 5)];
        end     
    end
    xOffset = xOffset + max(topoM(:, 4)-1)*nodeSpacing + treeSpacing;
end
xlim([-treeSpacing, xOffset+treeSpacing]);
axis off
for i = 1:size(nodes2Fill, 1)
    xNode = nodes2Fill(i, 1);
    xSquare = [xNode- nodeSpacing/4, xNode - nodeSpacing/4, xNode+ nodeSpacing/4, xNode+ nodeSpacing/4];
    plot(xSquare, [endTime, endTime+size(EVMat, 2), endTime+size(EVMat, 2), endTime], 'k');
    for j = 1:size(EVMat, 2)
        fill(xSquare, [endTime+j-1, endTime+j, endTime+j, endTime+j-1], 1-EVMat(i, j), 'LineStyle', 'None');
    end
    if colorGroup(i) == 1
        text(xNode, endTime+size(EVMat, 2)+50, 'R', 'FontSize', 20, 'HorizontalAlignment', 'center');
    else
        text(xNode, endTime+size(EVMat, 2)+50, 'L', 'FontSize', 20, 'HorizontalAlignment', 'center');
    end
    if fillGroup(i) == 0
        text(xNode, endTime+size(EVMat, 2)+100, '*', 'FontSize', 40, 'HorizontalAlignment', 'center');
    end
end
if ~isempty(nodes2Draw)
    scatter(nodes2Draw(:, 1), nodes2Draw(:, 2),18, [.5, .5, .5], 'filled');
end
colormap hot
caxis([0 1])


