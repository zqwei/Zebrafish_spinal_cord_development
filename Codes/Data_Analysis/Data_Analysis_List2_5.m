%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.  Covariance analysis: factor analysis
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.4 FA vs islands
%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


% function Data_Analysis_List2_5(nFile)
% 
%     if nargin<1
%         nFile = 1;
%     end
%         
%     addpath('../Func');
%     setDir;    
%     fileDirName       = fileDirNames{nFile}; %#ok<NASGU,USENS>
%     fileName          = fileNames{nFile}; %#ok<USENS>
%     load([tempDatDir, fileName, '.mat']);
%     load([tempDatDir, fileName, '_FAEV.mat'], 'opt2Dim');
%     
%     % Plot EVs as a function of x-y location
%     numPlot          = length(timePoints)-1;
%     mCol             = 8;
%     mRow             = ceil(numPlot/mCol);
% %     minNumFactor     = 3;
% %     opt1Dim          = opt2Dim-minNumFactor+1;
% %     opt1Dim          = min(opt1Dim, 8);
%         
%     figure;
%     h = suptitle(fileName);
%     set(h,'Interpreter','none');
%     coeffThres       = 0.5;
%     
% %     for nTime           = 1:length(timePoints)-1        
% %         slicedDFF       = dff(:,timePoints(nTime)+1:timePoints(nTime+1));
% %         slicedDFF       = bsxfun(@minus, slicedDFF, mean(slicedDFF,2));
% %         slicedDFF       = bsxfun(@rdivide, slicedDFF, std(slicedDFF,[],2))';                
% %         lMatTime        = factoran(slicedDFF, opt2Dim(nTime),'maxit',10000);
% %         subplot(mRow, mCol, nTime)
% %         imagesc(lMatTime)
% %         caxis([-1 1])
% %         title(['Time from: ' num2str(timePoints(nTime)/4/3600,'%.2f') ' to ' num2str(timePoints(nTime+1)/4/3600,'%.2f') ' hr'])
% %         xlim([1 opt2Dim(nTime)]);
% %         ylim([1 size(slicedDFF,2)])        
% %     end
% %     setPrint(mCol*6, mRow*4.5, [plotDir, 'FALMat_', fileName], 'pdf')
%     
%     for nTime           = 1:length(timePoints)-1        
%         slicedDFF       = dff(:,timePoints(nTime)+1:timePoints(nTime+1)); %#ok<NODEF>
%         slicedDFF       = bsxfun(@minus, slicedDFF, mean(slicedDFF,2));
%         slicedDFF       = bsxfun(@rdivide, slicedDFF, std(slicedDFF,[],2))';                
%         lMatTime        = factoran(slicedDFF, opt2Dim(nTime),'maxit',10000);
%         subplot(mRow, mCol, nTime)
% %         imagesc([lMatTime>0.5,sum(lMatTime>0.5,2)>1])
%         imagesc(lMatTime>coeffThres)
%         colormap(gray(2))
%         title(['Time from: ' num2str(timePoints(nTime)/4/3600,'%.2f') ' to ' num2str(timePoints(nTime+1)/4/3600,'%.2f') ' hr'])
%         xlim([0.5 opt2Dim(nTime)+0.5]);
%         set(gca, 'XTick', 1:opt2Dim(nTime));
%         ylim([1 size(slicedDFF,2)])
%     end   
%     setPrint(mCol*6, mRow*4.5, [plotDir, 'FALMatTres05_', fileName], 'pdf')
% 
% end


function Data_Analysis_List2_5(nFile)

    if nargin<1
        nFile = 1;
    end
        
    addpath('../Func');
    setDir;    
    fileDirName       = fileDirNames{nFile}; %#ok<NASGU,USENS>
    fileName          = fileNames{nFile}; %#ok<USENS>
    load([tempDatDir, fileName, '.mat'],'timePoints');
    load([tempDatDir, fileName, '_LMatTime.mat'], 'LMatTime')
%     load([tempDatDir, fileName, '_PSDPeakTimeFAEV.mat'], 'opt2Dim');
    
    numPlot          = length(timePoints)-1;
    mCol             = 8;
    mRow             = ceil(numPlot/mCol);
        
%     figure;
%     h = suptitle(fileName);
%     set(h,'Interpreter','none');
%     coeffThres       = 0.3;
%     
%     for nTime           = 1:length(timePoints)-1        
%         subplot(mRow, mCol, nTime)
% %         imagesc([lMatTime>0.5,sum(lMatTime>0.5,2)>1])
%         imagesc(abs(LMatTime{nTime})>coeffThres) %#ok<USENS>
%         colormap(gray(2))
%         title(['Time from: ' num2str(timePoints(nTime)/4/3600,'%.2f') ' to ' num2str(timePoints(nTime+1)/4/3600,'%.2f') ' hr'])
%         xlim([0.5 opt2Dim(nTime)+0.5]);
%         set(gca, 'XTick', 1:opt2Dim(nTime));
%         ylim([1 size(LMatTime{nTime},1)])
%     end   
%     setPrint(mCol*6, mRow*4.5, [plotDir, 'FALMatTres05_', fileName], 'pdf')
    
    figure;
    h = suptitle(fileName);
    set(h,'Interpreter','none');
    
    for nTime           = 1:length(timePoints)-1        
        subplot(mRow, mCol, nTime)
%         imagesc([lMatTime>0.5,sum(lMatTime>0.5,2)>1])
        LMat            = LMatTime{nTime};
        if size(LMatTime{nTime}, 2) > 1
            LMat        = nnmf(LMat, size(LMatTime{nTime}, 2));
        end
        imagesc(LMat); 
        colormap(parula)
        caxis([-1 1])
        title(['Time from: ' num2str(timePoints(nTime)/4/3600,'%.2f') ' to ' num2str(timePoints(nTime+1)/4/3600,'%.2f') ' hr'])
        xlim([0.5 size(LMatTime{nTime}, 2)+0.5]);
        set(gca, 'XTick', 1:size(LMatTime{nTime}, 2));
        ylim([1 size(LMatTime{nTime},1)])
    end   
    setPrint(mCol*6, mRow*4.5, [plotDir, 'FALMat_', fileName], 'pdf')

end

