for nFile = 1:6
    addpath('Func');
    disp(['Processing fish #' num2str(nFile)]);
    estimateDrift(nFile);
    generateCoordinates(nFile);
    dbList = getSignal(nFile);
    
    for i = 1:2
        collectResult(dbList{i});
    end
    buildAtlas(nFile);
    calculateDFF(nFile);
%     ablationAnalysis(nFile);
end