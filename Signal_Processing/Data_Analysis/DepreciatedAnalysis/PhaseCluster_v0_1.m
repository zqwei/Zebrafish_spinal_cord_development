nFile = 1;
nTime = 60;

addpath('../Func');
setDir;    
fileName          = fileNames{nFile};
load([tempDatDir, fileName, '.mat'], 'dff', 'side', 'timePoints', 'activeNeuronMat'); 

slicedDFF     = dff(activeNeuronMat(:, nTime), timePoints(nTime)+1:timePoints(nTime)+1200); 
slicedDFF     = bsxfun(@minus, slicedDFF, mean(slicedDFF,2));
slicedDFF     = bsxfun(@rdivide, slicedDFF, std(slicedDFF,[],2))';

numNeuron     = sum(activeNeuronMat(:, nTime));

peakCorrMat   = zeros(numNeuron);
peakLagMat    = nan(numNeuron);

maxLags       = 40;
numStd        = 3;
fs            = 4;

for nNeuron   = 1:numNeuron
    for mNeuron = nNeuron+1:numNeuron
        [xcf, lags, bounds] = crosscorr(slicedDFF(:,nNeuron), slicedDFF(:,mNeuron), maxLags, numStd);
        [pks, locs]         = max(xcf);
        locs                = lags(locs);
        if locs>=0
            [corrValue, corrSig]        = corr(slicedDFF(1:end-locs,nNeuron), slicedDFF(1+locs:end,mNeuron));
        else
            [corrValue, corrSig]        = corr(slicedDFF(1-locs:end,nNeuron), slicedDFF(1:end+locs,mNeuron));
        end
        if corrSig<0.01
            peakLagMat(nNeuron, mNeuron)  = locs/fs;
            peakCorrMat(nNeuron, mNeuron) = pks;
        end
    end
end

