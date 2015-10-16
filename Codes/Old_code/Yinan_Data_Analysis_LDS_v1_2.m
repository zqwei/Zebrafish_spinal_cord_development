% This analysis is for the temporal correlation of the single neurons at 
% each stage in Yinan's dataset
% This code deal with distribution of pink and purple cells in each fit

setDir;
load ([dirImageData, flistImageData(1).name]);
load ([dirImageData, 'clustering.mat']);
idxLR        = (idx == 1) | (idx ==2);
Y            = dff(idxLR, :); % Suggested by Yinan

% load the file with maximum EV in each dimensionality
LT           = 5001;




Start_point        = 16000:5000:56000;
xDim               = 1:20;
numTrial           = 20;
Best_fit           = nan(length(Start_point),length(xDim));
sum_EV             = nan(length(Start_point),length(xDim));

figure;
for nStart_point   = 1:length(Start_point)
    for nxDim      = 1:length(xDim)
        value_matrix  = nan(numTrial,1);
        for nTrial = 1:numTrial
            if exist(['Fitted_result/Latent_system_',...
                       num2str(Start_point(nStart_point),'%05d'),...
                       '_',num2str(xDim(nxDim),'%02d'),'_',...
                       num2str(nTrial,'%02d'),'.mat'],'file')
               load(['Fitted_result/Latent_system_',...
                       num2str(Start_point(nStart_point),'%05d'),...
                       '_',num2str(xDim(nxDim),'%02d'),'_',...
                       num2str(nTrial,'%02d'),'.mat']);
               value_matrix(nTrial) = correct;
            end
        end
        
        if ~isempty(find(value_matrix == nanmax(value_matrix)))
            Best_fit(nStart_point, nxDim) = find(value_matrix == nanmax(value_matrix));
        end
        
        if ~isnan(Best_fit(nStart_point, nxDim))
            
            load(['Fitted_result/Latent_system_',...
                       num2str(Start_point(nStart_point),'%05d'),...
                       '_',num2str(xDim(nxDim),'%02d'),'_',...
                       num2str(Best_fit(nStart_point, nxDim),'%02d'),'.mat']);
            
            subplot(length(Start_point), length(xDim), (nStart_point-1)*length(xDim) + nxDim)
            Y_now                     = Y(:,Start_point(nStart_point):Start_point(nStart_point)+5000);
            mean_Y                    = mean(Y_now, 2);
            var_Y                     = sum(bsxfun(@minus, Y_now, mean_Y).^2, 2);
            [~, y_est, ~]             = loo (Y_now, Ph, [0, LT]);
            err_Y                     = sum((Y_now-y_est).^2, 2);
            ev_Y                      = 1 - err_Y./var_Y;           
            hist(ev_Y,-0.4:0.1:1);    
            xlim([-0.4 1.0])
            ylim([0 30])
            sum_EV(nStart_point, nxDim) = sum(ev_Y);
        end
        
    end
end

figure;
imagesc(sum_EV)