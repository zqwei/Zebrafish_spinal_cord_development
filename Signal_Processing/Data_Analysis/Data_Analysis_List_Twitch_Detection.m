
thresTwichCor = 0.35:0.05:0.6;

performanceMat = nan(length(thresTwichCor), 24);

for nFile = 1:24
    for nThres = 5:6%1:length(thresTwichCor)
        performanceMat(nThres, nFile) = Neuron_selection_v4(nFile, thresTwichCor(nThres));
        disp(performanceMat(nThres, nFile))
    end
end

% result

% performanceMat = [333,119,9,220,203,4,19,80,0,0,4,0,1,0,2,3,0,0,0,1,1,1,0,7;...
%                   207,49,2,123,70,0,11,43,0,0,2,0,1,0,2,3,0,0,0,1,0,1,0,6;...
%                   80,7,1,37,11,0,6,13,0,0,0,0,1,0,0,2,0,0,0,0,0,1,0,3;...
%                   36,2,1,14,9,0,6,3,0,0,0,0,1,0,0,1,0,0,0,0,0,1,0,2;...
%                   17,1,1,6,4,0,4,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,2;...
%                   4,0,1,5,1,0,3,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,2];