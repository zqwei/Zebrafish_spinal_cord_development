addpath('../Func');
setDir;    

fileToAnalysis = [3, 4, 7, 12, 10, 11, 13, 15, 16];

neuronActTimeAll = [];
neuronFASizeAll  = [];

for nFile = fileToAnalysis
    fileName   = fileNames{nFile};   %#ok<*USENS>
    load([tempDatNetDir, 'NeuronActTime', fileName, '.mat'], 'neuronActTime', 'neuronTime', 'neuronFASize');
    
    neuronActTimeAll = [neuronActTimeAll; neuronActTime];
    neuronFASizeAll  = [neuronFASizeAll; neuronFASize];
    
end


neuronActTimeThres = 0.6;
bins  = 1:10;
figure;
subplot(1, 3, 1)
[fout, xout] = hist(neuronActTimeAll(~isnan(neuronActTimeAll)), 0:0.05:1);
stairs(xout, fout, 'linewid', 2)
xlim([0 1.1])
xlabel('Peak fraction of activation period', 'fontsize', 24)
ylabel('Number Neuron', 'fontsize', 24)
box off

subplot(1, 3, 2)
[fout, xout] = hist(neuronFASizeAll(neuronActTimeAll>neuronActTimeThres), bins);
stairs(xout, fout, 'linewid', 2)
xlabel('Factor size', 'fontsize', 24)
ylabel('Number Neuron', 'fontsize', 24)

box off

subplot(1, 3, 3)
[fout, xout] = hist(neuronFASizeAll(neuronActTimeAll<=neuronActTimeThres), bins);
stairs(xout, fout, 'linewid', 2)
xlabel('Factor size', 'fontsize', 24)
ylabel('Number Neuron', 'fontsize', 24)
box off