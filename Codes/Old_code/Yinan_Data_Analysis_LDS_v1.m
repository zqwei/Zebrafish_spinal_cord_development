% This analysis is for the temporal correlation of the single neurons at 
% each stage in Yinan's dataset

setDir;
timeNow      = 56000;
Y            = dff(:, timeNow:timeNow+5000); % Suggested by Yinan
nLength      = 2000;
nTrial       = 1000;
rand_Y       = randomCropSingleTrial(Y, nLength, nTrial);

mean_type    = 'Constant_mean';
tol          = 1e-5;
cyc          = 10000;

xDim         = 10;
[yDim, T, K] = size(rand_Y);


is_suc = false;

while ~is_suc
    try
        Ph      = lds_uni_latent(rand_Y, xDim, 'mean_type',mean_type,'tol',tol,'cyc',cyc);
        is_suc  = true;
    catch
        is_suc  = false;
    end
end

LT                   = size(Y,2);

[err, y_est, rand_y] = loo (Y, Ph, [0, LT]);


figure;
m           = ceil(sqrt(size(Y, 1)));

for nNeuron = 1:size(Y,1)
    subplot(m, m, nNeuron);
    hold on
    switch idx(nNeuron)
        case 1
            plot(Y(nNeuron,:),'-k');
        case 2
            plot(Y(nNeuron,:),'-b');
        case 3
            plot(Y(nNeuron,:),'-g');
    end
    plot(y_est(nNeuron,:),'-r');
    hold off
    xlim([0 5000])
end


setPrintPdf(8*m, 6*m, ['Fig_xDim_10_all_units_',num2str(timeNow),'.pdf']);


figure;
m           = ceil(sqrt(size(Y, 1)));

for nNeuron = 1:size(Y,1)
    subplot(m, m, nNeuron);
    hold on
    switch idx(nNeuron)
        case 1
            plot(400:700,Y(nNeuron,400:700),'-k');
        case 2
            plot(400:700,Y(nNeuron,400:700),'-b');
        case 3
            plot(400:700,Y(nNeuron,400:700),'-g');
    end
    plot(400:700, y_est(nNeuron,400:700),'-r');
    hold off
    xlim([400 700])
end


setPrintPdf(8*m, 6*m, ['Fig_xDim_10_all_units_short_',num2str(timeNow),'.pdf']);

% figure;
% mean_Y       = mean(Y, 2);
% var_Y        = sum(bsxfun(@minus, Y, mean_Y).^2, 2);
% err_Y        = sum((Y-y_est).^2, 2);
% ev_Y         = 1 - err_Y./var_Y;
% 
% histg(ev_Y, idx)
% box off
% xlabel('Explained Variance','fontsize',12);
% ylabel('Frequency', 'fontsize', 12)
% xlim([-0.2 1.0]);
% ylim([0 0.2]);
% setPrintPdf(8, 6, ['Fig_xDim_10_all_units_EV',num2str(timeNow),'.pdf']);
% 
% figure;
% scatter(log10(var_Y), ev_Y, [], idx, 'filled')
% box off
% xlabel('Total Variance','fontsize',12);
% ylabel('Explained Variance', 'fontsize', 12)
% set(gca,'XTick',0:1:3,'XTickLabel',{'1','','100',''})
% xlim([0 3]);
% ylim([0 1.0]);
% setPrintPdf(8, 6, ['Fig_xDim_10_all_units_EV_scatter_',num2str(timeNow),'.pdf']);

figure;
hold on;
scatter(log10(var_Y(idx==1)), ev_Y(idx==1), [], ev_Y(idx==1), 'o','filled')
scatter(log10(var_Y(idx==2)), ev_Y(idx==2), [], ev_Y(idx==2), 's','filled')
scatter(log10(var_Y(idx==3)), ev_Y(idx==3), [], ev_Y(idx==3), 'v','filled')
hold off
box off
xlabel('Total Variance','fontsize',12);
ylabel('Explained Variance', 'fontsize', 12)
set(gca,'XTick',-1:1:2,'XTickLabel',{'0.1','1','10','100'})
xlim([-1 3]);
ylim([-0.4 1.0]);
colormap(cool)
caxis([-0.4 1.0]); colorbar
setPrintPdf(8, 6, ['Fig_xDim_10_all_units_EV_scatter_',num2str(timeNow),'.pdf']);


% % % idxLR        = (idx == 1) | (idx ==2);
% % % Y            = dff(idxLR, 15000:20000); % Suggested by Yinan
% % % nLength      = 2000;
% % % nTrial       = 1000;
% % % rand_Y       = randomCropSingleTrial(Y, nLength, nTrial);
% % % mean_type    = 'Constant_mean';
% % % tol          = 1e-5;
% % % cyc          = 10000;
% % % xDim         = 10;
% % % [yDim, T, K] = size(rand_Y);
% % % is_suc = false;
% % % 
% % % while ~is_suc
% % %     try
% % %         Ph      = lds_uni_latent(rand_Y, xDim, 'mean_type',mean_type,'tol',tol,'cyc',cyc);
% % %         is_suc  = true;
% % %     catch
% % %         is_suc  = false;
% % %     end
% % % end
% % % 
% % % LT                   = size(Y,2);
% % % [err, y_est, rand_y] = loo (Y, Ph, [0, LT]);
% % % 
% % % 1 - err
% % % 
% % % figure;
% % % m           = ceil(sqrt(size(Y, 1)));
% % % 
% % % idx_LR      = idx(idxLR);
% % % 
% % % for nNeuron = 1:size(Y,1)
% % %     subplot(m, m, nNeuron);
% % %     hold on
% % %     switch idx_LR(nNeuron)
% % %         case 1
% % %             plot(Y(nNeuron,:),'-k');
% % %         case 2
% % %             plot(Y(nNeuron,:),'-b');
% % %         case 3
% % %             plot(Y(nNeuron,:),'-g');
% % %     end
% % %     plot(y_est(nNeuron,:),'-r');
% % %     xlim([0 5000])
% % %     hold off
% % % end
% % %           
% % % % setPrintPdf(8*m, 6*m, 'Fig4.pdf');