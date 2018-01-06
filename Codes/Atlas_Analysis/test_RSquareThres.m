RThresList = 0:0.05:0.9;
outlierList = zeros(numel(RThresList), 24);
fitList = zeros(numel(RThresList), 24);
for i = 1:numel(RThresList)
    for nFile = 1:24
        [outlierList(i, nFile), fitList(i, nFile)] = Leader_v0_1(nFile, RThresList(i));
    end
end
    