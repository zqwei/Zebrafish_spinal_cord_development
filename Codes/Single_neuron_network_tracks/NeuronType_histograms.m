neuronType_all = [];
neuronActTime_all = [];
neuronOscTime_all = [];
control_datasets = [3, 4, 7, 12, 10, 11, 13, 15, 16];
addpath('../Func');
setDir;
for nFile = control_datasets
    fileName          = fileNames{nFile};   %#ok<*USENS>
    load([tempDatNetDir, 'LONOLoading_' fileName, '_v_0_1.mat'], 'neuronOscTime', 'neuronType', 'neuronActTime');
    neuronType_all = [neuronType_all; neuronType];
    neuronActTime_all = [neuronActTime_all; neuronActTime];
    neuronOscTime_all = [neuronOscTime_all; neuronOscTime];
end

figure, hold on
for nType = 1:3
    [fout, xout] = ksdensity(neuronActTime_all(neuronType_all == nType), 0:0.02:1, 'Bandwidth', 0.02);
    stairs(xout, fout/sum(fout), 'linewid', 2)
    xlim([0 1])
end
legend({'type 1', 'type 2', 'type 3'});
title('neuron activity level');
setPrint(8, 6, [plotNetDir 'SingleNeuronActLevelFactorTimeCellType_Summary.pdf'])


figure, hold on
for nType = 1:3
    [fout, xout] = ksdensity(neuronOscTime_all(neuronType_all == nType), 0:0.02:1, 'Bandwidth', 0.02);
    stairs(xout, fout/sum(fout), 'linewid', 2)
    xlim([0 1])
end
legend({'type 1', 'type 2', 'type 3'});
title('neuron activity level');
setPrint(8, 6, [plotNetDir 'SingleNeuronOscLevelFactorTimeCellType_Summary.pdf'])


