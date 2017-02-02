%% Plots -- L*T*L+Psi

nFile = 2;
nTime = 301;

addpath('../Func/')
setDir;

fileName = fileNames{nFile};
load([tempDatDir, fileName, '.mat'], 'sideSplitter', 'side', 'tracks', 'timePoints', 'activeNeuronMat','mnx'); 
load([tempDatDir, 'LONOLoading_' fileName, '.mat'], 'LMat', 'PsiMat', 'CorrectedLMat');

lambda = CorrectedLMat{nTime};
lambdaT = LMat{nTime};
lambda(isnan(lambda)) = 0;
lambdaT(isnan(lambdaT)) = 0;

figure;
imagesc(lambda,[-1 1])
ylabel('Neuron index')
xlabel('Func. Com.')
% setPrint(4, 6, 'corrected_lambda', 'pdf')

T = lambdaT'/lambda';
imagesc(T'*T,[-1 1]), axis xy
ylabel('Func. Com.')
xlabel('Func. Com.')
% setPrint(4, 3, 'corrected_T', 'pdf')

imagesc(diag(PsiMat{nTime}),[-1 1]), axis xy
ylabel('Neuron index')
xlabel('Neuron index')
addpath('../Func/')
% setPrint(8, 6, 'corrected_Psi', 'pdf')

%% Plots -- estDFF
nNeuron = 40;
psi     = PsiMat{nTime};
load([tempDatDir, fileName, '.mat'], 'dff', 'timePoints','sideSplitter');
slicedDFF    = dff(:,timePoints(nTime)+1:timePoints(nTime)+1200); %#ok<NODEF>
slicedDFF    = bsxfun(@minus, slicedDFF, mean(slicedDFF,2));
slicedDFF    = bsxfun(@rdivide, slicedDFF, std(slicedDFF,[],2))';

DFFEst  = LONOFASingleUnit (slicedDFF, lambdaT, psi, nNeuron);

figure;
plot((1:1200)/4, slicedDFF(:, nNeuron), '-k', (1:1200)/4, DFFEst, '-r')
ylabel('z-DF/F')
xlabel('Time (s)')
xlim([0 300])
box off
% setPrint(8, 6, 'DFF_EstDFF_1', 'pdf')

nNeuron = 25;
DFFEst  = LONOFASingleUnit (slicedDFF, lambdaT, psi, nNeuron);

figure;
plot((1:1200)/4, slicedDFF(:, nNeuron), '-k', (1:1200)/4, DFFEst, '-r')
ylabel('z-DF/F')
xlabel('Time (s)')
xlim([0 300])
box off
% setPrint(8, 6, 'DFF_EstDFF_2', 'pdf')

%% LMat plots
m = ceil(sqrt(length(CorrectedLMat)/10));

figure;
for nTime = 1:floor(length(CorrectedLMat)/10)
    lambda = CorrectedLMat{nTime*10-9};
    lambda(isnan(lambda)) = 0;
    lambda(sum(lambda,1)==0) = [];
    subplot(m, m, nTime)
    hold on
    imagesc(lambda,[-1 1])
    axis xy
    plot([0.5 7.5], [sideSplitter sideSplitter], '--w')
    ylabel('Neuron index')
    xlabel('Func. Com.')
    xlim([0.5 7.5])
    ylim([0.5 length(lambda)+0.5])
    title([num2str(nTime*10-9) 'min'])
    box off
end
setPrint(8*m, 6*m, 'LMat_Time', 'pdf')


%% LMat plot in space
load([tempDatDir, fileName, '.mat'], 'sideSplitter', 'side', 'tracks', 'timePoints', 'activeNeuronMat');
numNeuron         = length(side);
nTime             = 301;
radius            = 10;
xlimMin           = min(min(tracks(:,:,1)));
xlimMax           = max(max(tracks(:,:,1)));
ylimMin           = min(min(tracks(:,:,2)));
ylimMax           = max(max(tracks(:,:,2)));    
mColor            = [     0    0.4470    0.7410
                     0.8500    0.3250    0.0980];
mColor            = [mColor; cbrewer('qual', 'Dark2',  128, 'cubic')];
figure;
LMat          = CorrectedLMat{nTime}; % threshold by 0.3
LMat(isnan(LMat)) = 0;  
LMat(:, sum(LMat, 1)==0) = [];
cla reset
hold on
xtracks   = squeeze(mean(tracks(:, timePoints(nTime)+(1:1200), 1), 2)); 
ytracks   = squeeze(mean(tracks(:, timePoints(nTime)+(1:1200), 2), 2));            
plot(xtracks(activeNeuronMat(:, nTime)), ytracks(activeNeuronMat(:, nTime)), 'ok', 'MarkerFaceColor','k') %#ok<NODEF>
plot(xtracks(~activeNeuronMat(:, nTime)), ytracks(~activeNeuronMat(:, nTime)), 'ok')

[~, maxFactorPerNeuronIndex] = max(LMat(sum(LMat, 2)>0, :), [], 2);
sideRemoveList  = histc(maxFactorPerNeuronIndex, 1:size(LMat, 2)) <2; % remove the factor has no dominate factors

LMat(:, sideRemoveList) = [];
LMat            = LMat > 0;               
sizeLMat        = sum(LMat, 1);
[~, indexLMat]  = sort(sizeLMat, 'descend');
LMat            = LMat(:, indexLMat);  
factorIndex     = 1:size(LMat, 2);

for nFactor = 1:size(LMat, 2)
    neuronFactor = LMat(:, nFactor)>0;
    if length(unique(side(neuronFactor)))==1
        CHPoints = smoothedBoundary(xtracks(neuronFactor), ytracks(neuronFactor), radius);
        patch(CHPoints(:,1), CHPoints(:,2), mColor(factorIndex(nFactor), :), 'facealpha', 0.6, 'edgecolor', 'none');
    else
        if sum(side(neuronFactor)==1) == 1 || sum(side(neuronFactor)==2) == 1 
            plot(xtracks(neuronFactor), ytracks(neuronFactor), 'o', 'Color', mColor(factorIndex(nFactor), :), 'MarkerFaceColor', mColor(factorIndex(nFactor), :))
        else
            plot(xtracks(neuronFactor), ytracks(neuronFactor), 's', 'Color', mColor(factorIndex(nFactor), :), 'MarkerFaceColor', mColor(factorIndex(nFactor), :), 'MarkerSize', 10)
        end
    end
end       
text((xlimMin+xlimMax)/2, ylimMax-10, [num2str(nTime) ' min'],'fontsize', 24)
xlim([xlimMin-5 xlimMax+5])
ylim([ylimMin-5 ylimMax+5])
box off
hold off;
xlabel('x-axis')
ylabel('y-axis')
% setPrint(16, 12, 'ExampleFASpace', 'pdf')

%% Plot z-DFF
nTime  = 301;
lambda = CorrectedLMat{nTime};
lambda(isnan(lambda)) = 0;
psi     = PsiMat{nTime};
slicedDFF    = dff(:,timePoints(nTime)+1:timePoints(nTime)+1200); %#ok<NODEF>
slicedDFF    = bsxfun(@minus, slicedDFF, mean(slicedDFF,2));
slicedDFF    = bsxfun(@rdivide, slicedDFF, std(slicedDFF,[],2))';
figure;
hold on
plot((1:1200)/4, bsxfun(@plus, slicedDFF(:, sum(lambda,2)==0), find(sum(lambda,2)==0)'*5), '-', 'color',[0.5 0.5 0.5])
plot((1:1200)/4, bsxfun(@plus, slicedDFF(:, lambda(:,1)~=0), find(lambda(:,1)~=0)'*5), '-b')
plot((1:1200)/4, bsxfun(@plus, slicedDFF(:, lambda(:,2)~=0), find(lambda(:,2)~=0)'*5), '-r')
% ylabel('z-DF/F')
% xlabel('Time (s)')
axis off
box off
setPrint(16, 20, 'Trace', 'png')

DFF    = estFactor(lambda, psi, slicedDFF);
figure;
hold on
plot((1:1200)/4, bsxfun(@plus, DFF(1, :), 1*5), '-b')
plot((1:1200)/4, bsxfun(@plus, DFF(2, :), 2*5), '-r')
axis off
hold off
setPrint(16, 4, 'DFFTrace', 'png')

%%
close all
