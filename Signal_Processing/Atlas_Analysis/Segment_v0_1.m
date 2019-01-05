%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Segmental analysis: 
% 0. export the location of segmental leaders ranked by activeTime
% 
% 
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
% 

function Segment_v0_1(nFile)
addpath('../Func');
setDir;
fileName          = fileNames{nFile};

load([tempDatDir, 'Leader_', fileName, '.mat'], 'activeTime', 'patternTime');
load([tempDatDir, fileName, '.mat'], 'timePoints', 'tracks', 'new_x', 'new_y');

x = new_x;
y = new_y;
refTP = timePoints(find(mod(timePoints, 1200)==0, 1, 'last'));
points = squeeze(tracks(:, refTP, :));
leaderPoints_rank1 = nan(0, 3);
leaderPoints_rank2 = nan(0, 3);
points(:, 3) = points(:, 3) * 6;
for seg = floor(min(x-0.5)):ceil(max(x-0.5))
    currentSegLeft = find(y<0 & x>=seg+0.5 & x<seg+1.5);
    [~, ~, ic] = unique(activeTime(currentSegLeft));
    leaderPoints_rank1 = [leaderPoints_rank1; points(currentSegLeft(ic==1), :)];
    leaderPoints_rank2 = [leaderPoints_rank2; points(currentSegLeft(ic==2), :)];
    currentSegRight = find(y>0 & x>=seg+0.5 & x<seg+1.5);
    [~, ~, ic] = unique(activeTime(currentSegRight));
    leaderPoints_rank1 = [leaderPoints_rank1; points(currentSegRight(ic==1), :)];
    leaderPoints_rank2 = [leaderPoints_rank2; points(currentSegRight(ic==2), :)];
end

leaderPoints_rank1(:, 3) = leaderPoints_rank1(:, 3);
leaderPoints_rank2(:, 3) = leaderPoints_rank2(:, 3);
writeMarkerFile(leaderPoints_rank1, [plotDir, 'segmentleaders_rank1_' fileName '_TM' num2str(refTP, '%6d') '.marker']);
writeMarkerFile(leaderPoints_rank2, [plotDir, 'segmentleaders_rank2_' fileName '_TM' num2str(refTP, '%6d') '.marker']);
writeMarkerFile(points, [plotDir 'allcells_' fileName '_TM' num2str(refTP, '%6d') '.marker']);