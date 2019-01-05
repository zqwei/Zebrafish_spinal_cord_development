%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lineage analysis 3.1.1: developmental trajectory analysis
%
% Plot trajectory of progenitors in neural plate space
% estimate neural plate location based on polar plot using centroid of yolk
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%
function Lineage_v3_1_1(nFile)
addpath('../Func');
setDir;
fileName       = fileNames{nFile};
fileDirName    = fileDirNames{nFile};
load([fileDirName '/', 'dev_data.mat'], 'centroid');
load([tempDatDir, 'DevelopmentStat_' fileName, '.mat'], 'cellTr', 'birthStats', 'trackingM');
load([tempDatDir, fileName, '.mat'], 'slicedIndex', 'leafOrder', 'new_x', 'new_y', 'new_z', 'side', 'mnx');

cellTr = cellTr(slicedIndex, :, :);
cellTr = cellTr(leafOrder, :, :);

% fit plane saparating L and R in developmental space
nCells = sum(slicedIndex);
lastPos = squeeze(cellTr(:, end, :));
MdlLinear = fitcdiscr(lastPos,side);
K = MdlLinear.Coeffs(1,2).Const;
L = MdlLinear.Coeffs(1,2).Linear;
% figure,
% hold on
% scatter3(lastPos(side==1, 1), lastPos(side==1, 2), lastPos(side==1, 3), 'r', 'filled');
% scatter3(lastPos(side==2, 1), lastPos(side==2, 2), lastPos(side==2, 3), 'g', 'filled');
% f = @(x1,x2) (-K - L(1)*x1 - L(2)*x2)/L(3);
% ezsurf(f ,[min(lastPos(:, 1)), max(lastPos(:, 1)), min(lastPos(:, 2)), max(lastPos(:, 2))]);
LR_centroid = L'/norm(L);

% define AP
[AP_centroid, ~, ~] = fitBodyAxes(lastPos, new_x, new_y);
AP_centroid = project2Plane(AP_centroid, LR_centroid);
DV_centroid = cross(AP_centroid, LR_centroid);

% measure neural plate coordinates: phi(AP),theta(LR), r(DV) to centroid

neuralPlateLocTime_centroid = cellTr;
for i = 1:size(cellTr, 2)
    currPos = squeeze(cellTr(:, i, :)) - repmat(centroid, nCells, 1);
    neuralPlateLocTime_centroid(:, i, 3) = sqrt(sum(currPos.^2, 2));
    currPos = currPos ./ repmat(sqrt(sum(currPos.^2, 2)), 1, 3);
    neuralPlateLocTime_centroid(:, i, 1) = 180-acosd(project2Plane(currPos, LR_centroid) * AP_centroid');
    neuralPlateLocTime_centroid(:, i, 2) = 90-acosd(currPos * LR_centroid');
end
neuralPlateLoc_centroid = squeeze(neuralPlateLocTime_centroid(:, 1, :));
% clone size, birthlocation and division angle (mapped in neural plate space)
siblingCloneSize = nan(sum(slicedIndex), 1);
divAngle_centroid         = nan(sum(slicedIndex), 2); 
birthLoc_centroid         = nan(sum(slicedIndex), 3);
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
        
        % method 1: use radial direction and centroid
%         absV                = PMother - centroid;
%         absV                = absV/norm(absV);
%         divAngle(i, 1)      = acosd(abs(dot(divV, absV)/norm(divV))); % radial angle
%         divVPlate           = project2Plane(divV, absV);
%         divAngle(i, 2)      = acosd(abs(dot(divVPlate, AP)/norm(divVPlate))); % AP angle of projection on tangential plane

        % method2: use AP/DV axis
        divAngle_centroid(i, 1)      = 90 - acosd(abs(dot(divV, DV_centroid)/norm(divV)));
        divVPlate           = project2Plane(divV, DV_centroid);
        divAngle_centroid(i, 2)      = acosd(abs(dot(divVPlate, AP_centroid)/norm(divVPlate)));
        
        currPos = birthStats(i).birthLocation - centroid;        
        birthLoc_centroid(i, 3)      = norm(currPos);
        currPos             = currPos/norm(currPos);
        birthLoc_centroid(i, 1)      = 180-acosd(project2Plane(currPos, LR_centroid) * AP_centroid');
        birthLoc_centroid(i, 2)      = 90- acosd(currPos * LR_centroid');
    end
end


%% plotting
figure, h1=plot(squeeze(neuralPlateLocTime_centroid(side==1, :, 1))', 'r');
hold on, h2=plot(squeeze(neuralPlateLocTime_centroid(side==2, :, 1))', 'g');
xlabel('time')
ylabel('AP location neural plate')
legend([h1(1), h2(1)], {'left', 'right'});
setPrint(8, 6, [plotDir, 'NeuralPlateAP_side_', fileName '_centroid'], 'pdf');
close;

mColor = cool(256);
figure,hold on
for i = 1:nCells
    ci = round((new_x(i)-min(new_x)+1)/(max(new_x)-min(new_x)+1)*255);
    plot(squeeze(neuralPlateLocTime_centroid(i, :, 1))', 'Color', mColor(ci, :));
end
xlabel('time')
ylabel('AP location neural plate')
setPrint(8, 6, [plotDir, 'NeuralPlateAP_APsegments_', fileName '_centroid'], 'pdf');
close


figure, h1=plot(squeeze(neuralPlateLocTime_centroid(side==1, :, 2))', 'r');
hold on, h2=plot(squeeze(neuralPlateLocTime_centroid(side==2, :, 2))', 'g');
xlabel('time')
ylabel('LR location neural plate')
legend([h1(1), h2(1)], {'left', 'right'});
setPrint(8, 6, [plotDir, 'NeuralPlateLR_side_', fileName '_centroid'], 'pdf');
close

figure, h1=plot(squeeze(neuralPlateLocTime_centroid(mnx==0, :, 2))', 'r');
hold on, h2=plot(squeeze(neuralPlateLocTime_centroid(mnx==1, :, 2))', 'g');
xlabel('time')
ylabel('LR location neural plate')
legend([h1(1), h2(1)], {'mnx-', 'mnx+'});
setPrint(8, 6, [plotDir, 'NeuralPlateLR_mnx_', fileName '_centroid'], 'pdf');
close

mColor = winter(256);
figure,hold on
for i = 1:nCells
    ci = round((new_z(i)-min(new_z)+1)/(max(new_z)-min(new_z)+1)*255);
    plot(squeeze(neuralPlateLocTime_centroid(i, :, 3))', 'Color', mColor(ci, :));
end
xlabel('time')
ylabel('DV location neural plate')
disp('');
setPrint(8, 6, [plotDir, 'NeuralPlateDV_DV_', fileName '_centroid'], 'pdf');
close

%% save metrics
save([tempDatDir, 'Leader_' fileName, '.mat'], 'neuralPlateLoc_centroid', 'siblingCloneSize', 'birthLoc_centroid', 'divAngle_centroid', '-append');
save([tempDatDir, 'DevelopmentStat_' fileName, '.mat'], 'AP_centroid', 'LR_centroid', 'centroid', '-append');
end

