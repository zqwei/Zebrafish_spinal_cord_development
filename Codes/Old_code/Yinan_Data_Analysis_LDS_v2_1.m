% This analysis is for the temporal correlation of the single neurons at 
% each stage in Yinan's dataset

setDir;
load ([dirImageData, flistImageData(1).name]);
load ([dirImageData, 'clustering.mat']);

rng('shuffle');

for nFit   = 10101 : 12020
    nStart_point = floor(nFit/10000);
    nxDim        = floor((nFit - nStart_point*10000)/100);
    nTrial       = nFit - nStart_point*10000 - nxDim * 100;
    
    if nxDim <=20 && nxDim >=1 && nTrial>=1 && nTrial<=20 ...
            && ~exist(['U:/Yinan_Fitted_result/Latent_system_',...
                       num2str((nStart_point+2)*5000+1000,'%05d'),...
                       '_',num2str(nxDim,'%02d'),'_',...
                       num2str(nTrial,'%02d'),'.mat'],'file')
                   
        Yinan_Data_Analysis_LDS_v2((nStart_point+2)*5000+1000, nxDim, dff, idx, nTrial);
    end    
end