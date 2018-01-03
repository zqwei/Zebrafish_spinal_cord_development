function Export_Segment_Leaders(nFile)
% export marker files containing locations of the first active cell in each
% segment (segmental boundary here defined by 0.5 in x)
addpath('../Func');
setDir;

load([TempDataDir '/tmp_' dataset{nFile} '.mat']);
load([DirNames{nFile} '\data.mat'], 'tracks', 'timePoints');

refTP = timePoints(find(mod(timePoints, 1200)==0, 1, 'last'));
points = squeeze(tracks(:, refTP, :));
leaderPoints_rank1 = nan(0, 3);
leaderPoints_rank2 = nan(0, 3);
points(:, 3) = points(:, 3) * 6;
for seg = floor(min(x-0.5)):ceil(max(x-0.5))
    currentSegLeft = find(y<0 & x>=seg+0.5 & x<seg+1.5);
    [~, ~, ic] = unique(halfActTime(currentSegLeft));
    leaderPoints_rank1 = [leaderPoints_rank1; points(currentSegLeft(ic==1), :)];
    leaderPoints_rank2 = [leaderPoints_rank2; points(currentSegLeft(ic==2), :)];
    currentSegRight = find(y>0 & x>=seg+0.5 & x<seg+1.5);
    [~, ~, ic] = unique(halfActTime(currentSegRight));
    leaderPoints_rank1 = [leaderPoints_rank1; points(currentSegRight(ic==1), :)];
    leaderPoints_rank2 = [leaderPoints_rank2; points(currentSegRight(ic==2), :)];
end

leaderPoints_rank1(:, 3) = leaderPoints_rank1(:, 3);
leaderPoints_rank2(:, 3) = leaderPoints_rank2(:, 3);
writeMarkerFile(leaderPoints_rank1, [PlotDir '/segmentleaders_rank1_' dataset{nFile} '_TM' num2str(refTP, '%6d') '.marker']);
writeMarkerFile(leaderPoints_rank2, [PlotDir '/segmentleaders_rank2_' dataset{nFile} '_TM' num2str(refTP, '%6d') '.marker']);
writeMarkerFile(points, [PlotDir '/segmentall_' dataset{nFile} '_TM' num2str(refTP, '%6d') '.marker']);