function Plot_Metrics_Atlas(nFile)
% factorSize taken from mean of a timeWindow=20 with outliers (>3sigma
% away) removed
addpath('../Func');
setDir;
load([DirNames{nFile} '\EV.mat'], 'halfActTime', 'halfEVTime', 'RSquare', 'validFitIndex');
load([DirNames{nFile} '\data.mat'], 'new_x', 'new_y', 'new_z', 'side', 'timePoints', 'activeNeuronMat', 'slicedIndex', 'leafOrder', 'mnx');
load([DirNames{nFile} '\LONOLoading.mat'], 'CorrectedLMat');
load([DirNames{nFile} '\profile_mnx_r8.mat'], 'profile_all');

if (exist([DirNames{nFile} '\birthtime.mat'], 'file'))
    load([DirNames{nFile} '\birthtime.mat']);
end

if (exist([DirNames{nFile} '\islet.mat'], 'file'))
    load([DirNames{nFile} '\islet.mat']);
end



% throw away bad fitted tact and tEV

halfActTime(~validFitIndex) = NaN;
halfActTime(halfActTime<0 | halfActTime>timePoints(end)/(4*3600)) = NaN;
halfEVTime(~validFitIndex) = NaN;
halfEVTime(halfEVTime<0 | halfEVTime>timePoints(end)/(4*3600)) = NaN;

nNeurons = numel(new_x);

% calculate factor size with window=20
windowSize = 20;
fs = nan(nNeurons, numel(timePoints));

for period = 1:numel(timePoints)
    LMat = CorrectedLMat{period};
    for i = 1:nNeurons
        if activeNeuronMat(i, period)
            % determine factor belongings - a single active neuron
            if any(isnan(LMat(i, :))) || ~any(LMat(i, :)>0)
                fs(i, period) = 1;
            else
                % factor belongings - determine dominate factor
                fs(i, period) = sum(LMat(:, LMat(i, :)==max(LMat(i, :)))>0);
            end
        end
    end
end

% for debugging
plotFactorSize(fs);

factorSize = nan(nNeurons, 1);
for i = 1:nNeurons
    ft = fs(i, ~isnan(fs(i, :)));
    if ~isempty(ft) && numel(ft)>windowSize
        f = ft(1:windowSize);
        f(abs(f-mean(f))>3*std(f)) = [];
        factorSize(i) = mean(f);
    end
end
% factorSize with invalid fitting curve will also no count
% factorSize(~validFitIndex) = NaN;
% factorSize(isnan(halfActTime) | isnan(halfEVTime)) = NaN;
% factorSize(RSquare<0) = NaN;

% calculate other metrics
mnx_level = profile_all(slicedIndex, :);
mnx_level = mnx_level(leafOrder, 1);
diffHalfTime = halfEVTime - halfActTime;

% plot atlas with different metrics
x = new_x;
y = new_y;
z = new_z;
metrics = {'halfActTime', 'halfEVTime', 'tEV-tAct', 'mnx level', 'factorSize', 'birth time', 'islet'};
for it = 1:7
    validMetrics = 1;
    switch it
        case 1
            me = halfActTime;
        case 2
            me = halfEVTime;
        case 3
            me = diffHalfTime;
        case 4
            
            me = mnx_level;
        case 5
            me = factorSize;
        case 6
            if exist('birthtime', 'var');
                birthtime = birthtime(slicedIndex);
                birthtime = birthtime(leafOrder, :);
                me = birthtime;
            else
                validMetrics = 0;
            end
        case 7
            if exist('islet', 'var');
                islet = islet(slicedIndex);
                islet = islet(leafOrder, :);
                me = islet;
            else
                validMetrics = 0;
            end
    end
    if ~validMetrics
        continue;
    end
    show_atlas(x, y, z, 1:9, me, mnx>0);
    export_fig([PlotDir '/' metrics{it} '_' dataset{nFile} '.pdf'], '-nocrop');
    close
    
    save([TempDataDir '/tmp_' dataset{nFile} '.mat'], 'halfActTime', 'halfEVTime', 'diffHalfTime', 'mnx_level', 'factorSize', 'mnx', 'x', 'y', 'z');
    if exist('birthtime', 'var')
        save([TempDataDir '/tmp_' dataset{nFile} '.mat'], 'birthtime', '-append');
    end
    if exist('islet', 'var')
        save([TempDataDir '/tmp_' dataset{nFile} '.mat'], 'islet', '-append');
    end
end