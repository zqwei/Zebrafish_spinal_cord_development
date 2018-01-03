numFiles = 21;
for nFile = fileList
    disp(['processing dataset nFile = ' num2str(nFile)]);
%     % step 1: plot calcium traces and atlas
%     disp('plotting calcium traces and atlas');
%     Plot_Traces_Map(nFile);
    % step 2: plot metrics onto atlas and save the metrics
    % half-Act, half-EV, tEV-tAct, mnx level, factorSize, birthtime, islet
    disp('plotting metrics onto atlas');
    Plot_Metrics_Atlas(nFile);
%     % step 3: plot segmental model
%     disp('plotting metrics onto segmental model');
%     Plot_Segmental_Model(nFile);
end


