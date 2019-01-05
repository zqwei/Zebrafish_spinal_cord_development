%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.  Covariance analysis: Factor analysis -- number of island using non
% LONO methods with active neurons
%
%
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
%


function FACluster_v0_3_whole_win(nFile, thresL, thresPsi)
addpath('../Func');
setDir;
fileName          = fileNames{nFile}; %#ok<USENS>
load([tempDatDir, fileName, '.mat'], 'dff','timePoints','activeNeuronMat', 'timeStep', 'new_x', 'new_y', 'new_z');
load([tempDatDir, 'FALONO_', fileName, '.mat'], 'uncorrectedLONOM', 'numFactors');

if mod(nFile, 2) == 1
    LONOM = 2;
else
    LONOM         = mode(uncorrectedLONOM);
end
activeNeuron  = sum(activeNeuronMat, 2) >0;
%     activeNeuron  = mean(activeNeuronMat, 2) >0.5;
LMat          = nan(size(activeNeuronMat, 1), LONOM);
PsiMat        = ones(size(activeNeuronMat, 1), 1);

slicedDFF     = dff(activeNeuron,timePoints(1)+1:timePoints(end)+timeStep);
%     slicedDFF     = dff(activeNeuron,timePoints(1)+1:timePoints(end));
slicedDFF     = bsxfun(@minus, slicedDFF, mean(slicedDFF,2));
slicedDFF     = bsxfun(@rdivide, slicedDFF, std(slicedDFF,[],2))';
[Lambda, Psi] = factoran(slicedDFF, LONOM,'maxit',10000, 'rotate', 'varimax');

LMat(activeNeuron, :)   = Lambda;
PsiMat(activeNeuron) = Psi;

uncorrectedLMat = LMat;


% based on SNR and LMat
LMat(LMat<0) = 0;
LMat(isnan(LMat)) = 0;

LMat(PsiMat>thresPsi, :) = 0;
LMat(LMat < thresL)   = 0;

% merge overlap
LMat(LMat<repmat(max(LMat, [], 2), 1, size(LMat, 2))) = 0;

% side correction: for bilateral factors with a single neuron on one
% side, remove that neuron
for nFactor = 1:size(LMat, 2)
    fNeuronLeft = find(new_y>0 & LMat(:, nFactor)>0);
    if sum(new_y>0 & LMat(:, nFactor)>0) ==1
        LMat(new_y>0 & LMat(:, nFactor)>0, nFactor) = 0;
    end
    if sum(new_y<0 & LMat(:, nFactor)>0) ==1
        LMat(new_y<0 & LMat(:, nFactor)>0, nFactor) = 0;
    end
end

% remove single neuron factor
LMat(:, sum(LMat, 1)<=1) = 0;


% remove empty factors
LMat(:, sum(LMat)==0) = [];
LONOM = size(LMat, 2);


new_y = new_y/2;
new_z = new_z/max(new_z) * 1.2;
figure;
subplot(3, 2, [1 3 5])
imagesc(LMat)
subplot(3, 2, 2)
hold on
plot(new_x(~activeNeuron), new_y(~activeNeuron), 'sk')
plot(new_x(activeNeuron), new_y(activeNeuron), 'ok')
for nFactor = 1:LONOM
    scatter(new_x(LMat(:,nFactor)>0), new_y(LMat(:,nFactor)>0), 'filled')
end
subplot(3, 2, 4)
hold on
plot(new_x(~activeNeuron & new_y<0), new_z(~activeNeuron & new_y<0), 'sk')
plot(new_x(activeNeuron & new_y<0), new_z(activeNeuron & new_y<0), 'ok')
for nFactor = 1:LONOM
    scatter(new_x(LMat(:,nFactor)>0 & new_y<0), new_z(LMat(:,nFactor)>0 & new_y<0), 'filled')
end
subplot(3, 2, 6)
hold on
plot(new_x(~activeNeuron & new_y>0), new_z(~activeNeuron & new_y>0), 'sk')
plot(new_x(activeNeuron & new_y>0), new_z(activeNeuron & new_y>0), 'ok')
for nFactor = 1:LONOM
    scatter(new_x(LMat(:,nFactor)>0 & new_y>0), new_z(LMat(:,nFactor)>0 & new_y>0), 'filled')
end

save([tempDatDir, 'FALONO_Average_', fileName, '.mat'], 'uncorrectedLMat', 'LMat', 'PsiMat');
setPrint(8*2, 6, [plotDir, 'AverageFactorAblation_', fileName], 'pdf')
end
