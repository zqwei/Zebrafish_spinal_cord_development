% This analysis is for the temporal correlation of the single neurons at 
% each stage in Yinan's dataset

setDir;
load ([dirImageData, flistImageData(1).name]);
load ([dirImageData, 'clustering.mat']);
load ([dirImageData, 'data.mat']);
timeNow      = 56000;
Y            = dff(:, timeNow:timeNow+5000); % Suggested by Yinan

Z            = linkage(Y, 'ward', 'correlation');
C            = cluster(Z, 'maxclust', 10);

figure;

time_Points  = timeNow + 5000;
scatter3(tracks_smoothed(:,time_Points,1),...
         tracks_smoothed(:,time_Points,2),...
         tracks_smoothed(:,time_Points,3),...
         [],C);
% colormap(cool)
     

load(['/Volumes/Zebrafish_Imaging_Data_by_Yinan/Data_Dre_E1_BTXinjHuCH2BGCaMP6f_TL_20140818_045650_corrected_signal/EV_Results/result_',...
      num2str(timeNow),'.mat']);
figure;
scatter3(tracks_smoothed(:,time_Points,1),...
         tracks_smoothed(:,time_Points,2),...
         tracks_smoothed(:,time_Points,3),...
         [],ev_Y);
% colormap(cool)

idxLR        = (idx == 1) | (idx ==2);
Motor_Y      = Y(idxLR,:);
figure;
imagesc(corr(Motor_Y'))

figure;
imagesc(corr(Y'))