%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lineage analysis 3.1: developmental trajectory analysis
%
% Plot trajectory of progenitors in neural plate space
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%
function Lineage_v3_1(nFile)
addpath('../Func');
setDir;
fileName       = fileNames{nFile};

load([tempDatDir, 'DevelopmentStat_' fileName, '.mat'], 'cellTr', 'birthStats', 'trackingM');
load([tempDatDir, fileName, '.mat'], 'slicedIndex', 'leafOrder', 'new_x', 'new_y', 'new_z', 'side', 'mnx');


cellTr = cellTr(slicedIndex, :, :);
cellTr = cellTr(leafOrder, :, :);

% fit body axes in developmental space
nCells = sum(slicedIndex);
lastPos = squeeze(cellTr(:, end, :));
initialPos = squeeze(cellTr(:, 1, :));
[AP, LR, ori] = fitBodyAxes(lastPos, new_x, new_y);
DV = cross(AP, LR);

% measure neural plate coordinates
neuralPlateLoc = (initialPos - repmat(ori, nCells, 1)) * [AP', LR'];
neuralPlateLocTime = cellTr(:,:, 1:2);
for i = 1:size(cellTr, 2)
    currPos = squeeze(cellTr(:, i, :));
    neuralPlateLocTime(:, i, :) = (currPos - repmat(ori, nCells, 1)) * [AP', LR'];
end

% [coeff, score, ~] = pca(initialPos - repmat(ori, nCells, 1));
% [~, dimAP] = max(abs(AP*coeff));
% [~, dimLR] = max(abs(LR*coeff));
% neuralPlateLoc = score(:, [dimAP, dimLR]);

% clone size, birthlocation and division angle (mapped in neural plate space)
siblingCloneSize = nan(sum(slicedIndex), 1);
divAngle         = nan(sum(slicedIndex), 2);
birthLoc         = nan(sum(slicedIndex), 2);
birthStats = birthStats(slicedIndex);
birthStats = birthStats(leafOrder);
for i = 1:sum(slicedIndex)
    if ~isempty(birthStats(i).siblingLeafList)
%         siblingCloneSize(i) = NaN;
%         divAngle(i, :)      = NaN;
%     else
        siblingCloneSize(i) = numel(birthStats(i).siblingLeafList);
        PMother             = trackingM(trackingM(:, 1)==birthStats(i).divisionTriplet(1), 3:5);
        PSibl1              = trackingM(trackingM(:, 1)==birthStats(i).divisionTriplet(2), 3:5);
        PSibl2              = trackingM(trackingM(:, 1)==birthStats(i).divisionTriplet(3), 3:5);
        divV                = PSibl2-PSibl1;
%         divAngle(i, 1)      = acosd(abs(dot(PSibl1-PMother, AP)/norm(PSibl1-PMother)));
%         divAngle(i, 2)      = acosd(abs(dot(PSibl2-PMother, AP)/norm(PSibl2-PMother)));
        divAngle(i, 1)      = 90 - acosd(abs(dot(divV, DV)/norm(divV)));
        divVPlate           = project2Plane(divV, DV);
        divAngle(i, 2)      = acosd(abs(dot(divVPlate, AP)/norm(divVPlate)));
        birthLoc(i, :)      = (birthStats(i).birthLocation - ori) * [AP', LR'];
    end
end


%% plotting
figure, h1=plot(squeeze(neuralPlateLocTime(side==1, :, 1))', 'r');
hold on, h2=plot(squeeze(neuralPlateLocTime(side==2, :, 1))', 'g');
xlabel('time')
ylabel('AP location neural plate')
legend([h1(1), h2(1)], {'left', 'right'});
setPrint(8, 6, [plotDir, 'NeuralPlateAP_side_', fileName], 'pdf');
close;

mColor = cool(256);
figure,hold on
for i = 1:nCells
    ci = round((new_x(i)-min(new_x)+1)/(max(new_x)-min(new_x)+1)*255);
    plot(squeeze(neuralPlateLocTime(i, :, 1))', 'Color', mColor(ci, :));
end
xlabel('time')
ylabel('AP location neural plate')
setPrint(8, 6, [plotDir, 'NeuralPlateAP_APsegments_', fileName], 'pdf');
close


figure, h1=plot(squeeze(neuralPlateLocTime(side==1, :, 2))', 'r');
hold on, h2=plot(squeeze(neuralPlateLocTime(side==2, :, 2))', 'g');
xlabel('time')
ylabel('LR location neural plate')
legend([h1(1), h2(1)], {'left', 'right'});
setPrint(8, 6, [plotDir, 'NeuralPlateLR_side_', fileName], 'pdf');
close

figure, h1=plot(squeeze(neuralPlateLocTime(mnx==0, :, 2))', 'r');
hold on, h2=plot(squeeze(neuralPlateLocTime(mnx==1, :, 2))', 'g');
xlabel('time')
ylabel('LR location neural plate')
legend([h1(1), h2(1)], {'mnx-', 'mnx+'});
setPrint(8, 6, [plotDir, 'NeuralPlateLR_mnx_', fileName], 'pdf');
close

mColor = winter(256);
figure,hold on
for i = 1:nCells
    ci = round((new_z(i)-min(new_z)+1)/(max(new_z)-min(new_z)+1)*255);
    plot(squeeze(neuralPlateLocTime(i, :, 2))', 'Color', mColor(ci, :));
end
xlabel('time')
ylabel('LR location neural plate')
setPrint(8, 6, [plotDir, 'NeuralPlateLR_DV_', fileName], 'pdf');
close

%% save metrics
save([tempDatDir, 'Leader_' fileName, '.mat'], 'neuralPlateLoc', 'siblingCloneSize', 'birthLoc', 'divAngle', '-append');
save([tempDatDir, 'DevelopmentStat_' fileName, '.mat'], 'AP', 'LR', 'ori', '-append');
end

