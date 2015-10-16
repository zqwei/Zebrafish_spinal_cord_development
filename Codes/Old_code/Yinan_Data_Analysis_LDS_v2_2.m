% This analysis is for the temporal correlation of the single neurons at 
% each stage in Yinan's dataset

setDir;

numTrial      = 20;

Start_point   = 16000:5000:56000;
xDim          = 1:16;

Exp_Var       = nan(length(Start_point)+1,length(xDim)+1);

for nStart_point   = 1:length(Start_point)
    for nxDim      = 1:length(xDim)
        value_matrix  = nan(numTrial,1);
        for nTrial = 1:numTrial
            if exist([dirResult, 'Latent_system_',...
                       num2str(Start_point(nStart_point),'%05d'),...
                       '_',num2str(xDim(nxDim),'%02d'),'_',...
                       num2str(nTrial,'%02d'),'.mat'],'file')
               load([dirResult, 'Latent_system_',...
                       num2str(Start_point(nStart_point),'%05d'),...
                       '_',num2str(xDim(nxDim),'%02d'),'_',...
                       num2str(nTrial,'%02d'),'.mat']);
               value_matrix(nTrial) = correct;
            end
        end
        
        Exp_Var(nStart_point, nxDim) = nanmax(value_matrix);
        
    end
end

[X, Y] = meshgrid([Start_point/1000 max(Start_point)/1000+5], [xDim max(xDim)+1]);
contourfWithNoLine(X, Y, Exp_Var', 128);
caxis([0 1])
colormap(jet(128))
colorbar
box off
xlabel('Starting points of the trial (k)','fontsize',12)
ylabel('# of latent dimenion','fontsize',12)
title('Amount of explained variance','fontsize',12)
setPrintPdf(8,6,'Fig1.pdf')

figure;
pcolorWithNoLine(X, Y, Exp_Var');
caxis([0 1])
colormap(jet(128))
colorbar
box off
xlim([16 61])
xlabel('Starting points of the trial (k)','fontsize',12)
ylabel('# of latent dimenion','fontsize',12)
title('Amount of explained variance','fontsize',12)
setPrintPdf(8,6,'Fig2.pdf')

Best_dim = [6 5 5 4 2 2 2 2 2];
plot(Start_point/1000, Best_dim,'-o','linewid',2);
box off
xlim([16 61])
xlabel('Starting points of the trial (k)','fontsize',12)
ylabel('Optimal # of latent dimenion','fontsize',12)
title('Amount of explained variance','fontsize',12)
setPrintPdf(8,6,'Fig3.pdf')