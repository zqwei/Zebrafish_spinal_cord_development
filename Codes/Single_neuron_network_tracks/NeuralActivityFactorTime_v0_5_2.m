addpath('../Func');
setDir;    

fileToAnalysis = [3, 4, 7, 12, 10, 15, 16];

neuronActTimeAll = [];
neuronFASizeAll  = [];
locationAll      = [];
group            = [];
neuronTimeAll    = [];

for nFile = fileToAnalysis
    fileName   = fileNames{nFile};   %#ok<*USENS>
    load([tempDatDir, fileName, '.mat'], 'new_x', 'new_y', 'new_z');
    load([tempDatDir, 'EV_' fileName, '.mat'], 'halfActTime');
    load([tempDatNetDir, 'NeuronActTime', fileName, '.mat'], 'neuronActTime', 'neuronTime', 'neuronFASize');    
    neuronActTimeAll = [neuronActTimeAll; neuronActTime];
    neuronFASizeAll  = [neuronFASizeAll; neuronFASize];   
    locationAll      = [locationAll; new_x];
    group            = [group; nFile * ones(size(new_x))]; 
    neuronTimeAll    = [neuronTimeAll; halfActTime];
end


neuronActTimeThres = 0.6;
% bins  = 1:10;
% figure;
% subplot(1, 3, 1)
% [fout, xout] = hist(neuronActTimeAll(~isnan(neuronActTimeAll)), 0:0.05:1);
% stairs(xout, fout, 'linewid', 2)
% xlim([0 1.1])
% xlabel('Peak fraction of activation period', 'fontsize', 24)
% ylabel('Number Neuron', 'fontsize', 24)
% box off
% 
% subplot(1, 3, 2)
% % [fout, xout] = hist(neuronFASizeAll(neuronActTimeAll>neuronActTimeThres & neuronTimeAll<1), bins);
% [fout, xout] = hist(neuronFASizeAll(neuronActTimeAll>neuronActTimeThres), bins);
% stairs(xout, fout, 'linewid', 2)
% xlabel('Factor size', 'fontsize', 24)
% ylabel('Number Neuron', 'fontsize', 24)
% 
% box off
% 
% subplot(1, 3, 3)
% [fout, xout] = hist(neuronFASizeAll(neuronActTimeAll<=neuronActTimeThres), bins);
% stairs(xout, fout, 'linewid', 2)
% xlabel('Factor size', 'fontsize', 24)
% ylabel('Number Neuron', 'fontsize', 24)
% box off

figure;
scatter(locationAll(neuronActTimeAll>neuronActTimeThres), neuronFASizeAll(neuronActTimeAll>neuronActTimeThres), [], group(neuronActTimeAll>neuronActTimeThres), 'filled');
xlabel('location')
ylabel('Factor Size')


figure;
scatter(neuronTimeAll(neuronActTimeAll>neuronActTimeThres), neuronFASizeAll(neuronActTimeAll>neuronActTimeThres), [], group(neuronActTimeAll>neuronActTimeThres), 'filled');
xlabel('first factor time')
ylabel('Factor Size')
% xlim([0 2])