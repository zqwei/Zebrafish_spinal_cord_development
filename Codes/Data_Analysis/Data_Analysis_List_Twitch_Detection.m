
thresTwichCor = 0.35:0.05:0.5;

performanceMat = nan(length(thresTwichCor), 24);

for nFile = 1:24
    for nThres = 1:length(thresTwichCor)
        performanceMat(nThres, nFile) = Neuron_selection_v4(nFile, thresTwichCor(nThres));
    end
end