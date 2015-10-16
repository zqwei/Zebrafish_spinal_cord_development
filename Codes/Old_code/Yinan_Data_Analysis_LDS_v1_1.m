% This analysis is for the temporal correlation of the single neurons at 
% each stage in Yinan's dataset

setDir;
load ([dirImageData, flistImageData(1).name]);
load ([dirImageData, 'clustering.mat']);


numTrial      = 20;
Start_point   = 16000:5000:56000;
xDim          = 1:20;

Best_fit      = nan(length(Start_point),length(xDim));

for nStart_point   = 1:length(Start_point)
    for nxDim      = 1:length(xDim)
        value_matrix  = nan(numTrial,1);
        for nTrial = 1:numTrial
            if exist(['Fitted_result/Latent_system_',...
                       num2str(Start_point(nStart_point),'%05d'),...
                       '_',num2str(xDim(nxDim),'%02d'),'_',...
                       num2str(nTrial,'%02d'),'.mat'],'file')
               load(['Fitted_result/Latent_system_',...
                       num2str(Start_point(nStart_point),'%05d'),...
                       '_',num2str(xDim(nxDim),'%02d'),'_',...
                       num2str(nTrial,'%02d'),'.mat']);
               value_matrix(nTrial) = correct;
            end
        end
        
        if ~isempty(find(value_matrix == nanmax(value_matrix)))
            Best_fit(nStart_point, nxDim) = find(value_matrix == nanmax(value_matrix));
        end
        
    end
end



Y            = dff(:, 56000:61000); % Suggested by Yinan
idxLR        = (idx == 1) | (idx ==2);
Y            = Y(idxLR, :); % Suggested by Yinan

nStart_point = 9;
nxDim        = 10;              
load(['Fitted_result/Latent_system_',...
       num2str(Start_point(nStart_point),'%05d'),...
       '_',num2str(xDim(nxDim),'%02d'),'_',...
       num2str(Best_fit(nStart_point, nxDim),'%02d'),'.mat']);
figure;
m           = ceil(sqrt(size(Y, 1)));

idx_LR      = idx(idxLR);

for nNeuron = 1:size(Y,1)
    subplot(m, m, nNeuron);
    hold on
    switch idx_LR(nNeuron)
        case 1
            plot(Y(nNeuron,:),'-k');
        case 2
            plot(Y(nNeuron,:),'-b');
        case 3
            plot(Y(nNeuron,:),'-g');
    end
    plot(y_est(nNeuron,:),'-r');
    xlim([0 5000])
    hold off
end

setPrintPdf(8*m, 6*m, 'Fig_xDim_10_LR_uints.pdf');

figure;
m           = ceil(sqrt(size(Y, 1)));

for nNeuron = 1:size(Y,1)
    subplot(m, m, nNeuron);
    hold on
    switch idx_LR(nNeuron)
        case 1
            plot(400:900, Y(nNeuron,400:900),'-k');
        case 2
            plot(400:900, Y(nNeuron,400:900),'-b');
        case 3
            plot(400:900, Y(nNeuron,400:900),'-g');
    end
    plot(400:900, y_est(nNeuron,400:900),'-r');
    xlim([400 900])
    hold off
end

setPrintPdf(8*m, 6*m, 'Fig_xDim_10_LR_uints_short.pdf');


nStart_point = 9;
nxDim        = 2;              
load(['Fitted_result/Latent_system_',...
       num2str(Start_point(nStart_point),'%05d'),...
       '_',num2str(xDim(nxDim),'%02d'),'_',...
       num2str(Best_fit(nStart_point, nxDim),'%02d'),'.mat']);
figure;
m           = ceil(sqrt(size(Y, 1)));

idx_LR      = idx(idxLR);

for nNeuron = 1:size(Y,1)
    subplot(m, m, nNeuron);
    hold on
    switch idx_LR(nNeuron)
        case 1
            plot(Y(nNeuron,:),'-k');
        case 2
            plot(Y(nNeuron,:),'-b');
        case 3
            plot(Y(nNeuron,:),'-g');
    end
    plot(y_est(nNeuron,:),'-r');
    xlim([0 5000])
    hold off
end

setPrintPdf(8*m, 6*m, 'Fig_xDim_02_LR_uints.pdf');

figure;
m           = ceil(sqrt(size(Y, 1)));

for nNeuron = 1:size(Y,1)
    subplot(m, m, nNeuron);
    hold on
    switch idx_LR(nNeuron)
        case 1
            plot(400:900, Y(nNeuron,400:900),'-k');
        case 2
            plot(400:900, Y(nNeuron,400:900),'-b');
        case 3
            plot(400:900, Y(nNeuron,400:900),'-g');
    end
    plot(400:900, y_est(nNeuron,400:900),'-r');
    xlim([400 900])
    hold off
end

setPrintPdf(8*m, 6*m, 'Fig_xDim_02_LR_uints_short.pdf');


CMatrix    = zeros(size(Y,1)+nxDim, size(Y,1)+nxDim);
CMatrix(size(Y,1)+1:end, :) = [Ph.C', Ph.A];

ids        = cell(size(Y,1)+nxDim, 1);

for nCell  = 1:size(Y,1)
    ids(nCell)  = {['Cell',num2str(nCell,'%02d')]};
end

for nDim   = 1:nxDim
    ids(nCell + nDim) = {['Latent',num2str(nDim,'%02d')]};
end

nodeColors     = [repmat([1 0 0], size(Y,1), 1); repmat([0 1 0], nxDim, 1)];
graphViz4Matlab('-adjMat',CMatrix,'-nodeLabels',ids,'-layout',Treelayout,'-nodeColors',nodeColors)

