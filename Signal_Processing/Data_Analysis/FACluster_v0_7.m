%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.  Covariance analysis: Factor analysis -- loadin matrix evolution --
% reordering
%
% center of the factor
% connectivity strength to each node
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
%


function FACluster_v0_7(nFile)
    addpath('../Func');
    setDir;
    fileName          = fileNames{nFile}; %#ok<USENS>
    load([tempDatDir, fileName, '.mat'], 'sideSplitter', 'side', 'tracks', 'timePoints', 'activeNeuronMat', 'dff', 'timeStep'); 
    load([tempDatDir, 'LONOLoading_' fileName, '.mat'], 'CorrectedLMat', 'PsiMat');
    numTime           = length(CorrectedLMat); %#ok<USENS>
    numNeuron         = length(side);
    networkMat        = cell(numTime, 1);

    neuronXLoc        = zeros(numNeuron, numTime);
    neuronYLoc        = zeros(numNeuron, numTime);
    neuronZLoc        = zeros(numNeuron, numTime);


    for nTime         = 1:numTime
        LMat                     = CorrectedLMat{nTime};
        LMat(isnan(LMat))        = 0;
        Ph                       = PsiMat{nTime}; %#ok<USENS>
        LMat(:, sum(LMat, 1)==0) = [];
        xtracks   = squeeze(mean(tracks(:, timePoints(nTime)+(1:timeStep), 1), 2));  %#ok<NODEF>
        ytracks   = squeeze(mean(tracks(:, timePoints(nTime)+(1:timeStep), 2), 2));
        ztracks   = squeeze(mean(tracks(:, timePoints(nTime)+(1:timeStep), 3), 2));
        nAct      = activeNeuronMat(:, nTime); %#ok<NODEF>
        slicedDFF    = dff(:,timePoints(nTime)+1:timePoints(nTime)+timeStep); %#ok<NODEF>
        slicedDFF    = bsxfun(@minus, slicedDFF, mean(slicedDFF,2));
        slicedDFF    = bsxfun(@rdivide, slicedDFF, std(slicedDFF,[],2))';
        neuronXLoc(:, nTime) = xtracks;
        neuronYLoc(:, nTime) = ytracks;
        neuronZLoc(:, nTime) = ztracks;
        % non-active neurons
        nFactor.x            = nan;
        nFactor.y            = nan;
        nFactor.z            = nan;
        nFactor.DFF          = nan(timeStep, 1);
        nFactor.neuronIndex  = find(~nAct);
        nFactor.neuronCCMat  = nan(length(nFactor.neuronIndex), 2); % cross correlation mat % first col: corr; second col: time

        factorSet            = {nFactor};

        % single neuron factor
        singleNeuronSet      = [];
        for nNeuron          = find(nAct)'
            % if a neuron is a single neuron on one side of the fish
            idVec            = LMat(nNeuron, :)~=0;
            idSide           = side(nNeuron);
            nFactor.x            = xtracks(nNeuron);
            nFactor.y            = ytracks(nNeuron);
            nFactor.z            = ztracks(nNeuron);
            nFactor.DFF          = slicedDFF(:, nNeuron);
            nFactor.neuronIndex  = nNeuron;
            nFactor.neuronCCMat  = nan(length(nFactor.neuronIndex), 2); % cross correlation mat % first col: corr; second col: time
            if sum(side(sum(LMat(:, idVec)~=0, 2)>0) == idSide)<=1
                factorSet        = [factorSet, nFactor]; %#ok<*AGROW>
                singleNeuronSet  = [singleNeuronSet, nNeuron];
            end
        end

        LMat(singleNeuronSet, :) = 0;
        LMat(:, sum(LMat, 1)==0) = [];
        Ph(sum(LMat, 2)==0)  = 1; % put the noise of kicking-out neurons to ones
        DFF                  = estFactor(LMat, Ph, slicedDFF);

        fs                   = 4;
        maxLags              = 25 * fs;
        numStd               = 3;


        % nodes inside factors
        for mFactor          = 1:size(LMat, 2)
            lVec             = LMat(:, mFactor);
            numSide1         = sum(side(lVec~=0)==1);
            numSide2         = sum(side(lVec~=0)==2);
            sideFactor       = 0;
            if numSide2<2 && numSide1>=2
                sideFactor   = 1;
            elseif numSide1<2 && numSide2>=2
                sideFactor   = 2;
            end

            if sideFactor>0
                nFactor.x            = xtracks(side==sideFactor)'*abs(lVec(side==sideFactor))/sum(abs(lVec(side==sideFactor)));
                nFactor.y            = ytracks(side==sideFactor)'*abs(lVec(side==sideFactor))/sum(abs(lVec(side==sideFactor)));
                nFactor.z            = ztracks(side==sideFactor)'*abs(lVec(side==sideFactor))/sum(abs(lVec(side==sideFactor)));
                mDFF                 = DFF(mFactor, :)';
                nFactor.DFF          = mDFF;
                nFactor.neuronIndex  = find(lVec~=0);
                neuronCCMat          = nan(numNeuron, 2);
                for nNeuron          = find(lVec~=0)'
                    [xcf, lags, bound] = crosscorr(mDFF, slicedDFF(:,nNeuron), maxLags, numStd);
                    [pks, locs]      = findpeaks(xcf, lags,'MinPeakDistance', 5,'MinPeakHeight', bound(1), 'MinPeakProminence', bound(1));
%                     [pks, locs]         = max(xcf);
                    pklocs           = abs(locs) == min(abs(locs));
                    locs             = locs(pklocs);
                    pks              = pks(pklocs);
                    [~,  pklocs]     = max(pks);
                    locs             = locs(pklocs);
                    pks              = pks(pklocs);
                    if ~isempty(locs)
                        if locs>=0
                            [~, corrSig]  = corr(mDFF(1:end-locs), slicedDFF(1+locs:end,nNeuron));
    %                         if abs(corrValue - pks)>abs(pks*0.01); keyboard();end
                        else
                            [~, corrSig]  = corr(mDFF(1-locs:end), slicedDFF(1:end+locs,nNeuron));
    %                         if abs(corrValue - pks)>abs(pks*0.01); keyboard();end
                        end
                        if corrSig<0.01
                            neuronCCMat(nNeuron, 1)     = locs/fs;
                            neuronCCMat(nNeuron, 2)     = pks;
                        end
                    end
                end
                nFactor.neuronCCMat  = neuronCCMat;
                factorSet        = [factorSet, nFactor];
            else
                for nSide  = 1:2
                    nFactor.x            = xtracks(side==nSide)'*abs(lVec(side==nSide))/sum(abs(lVec(side==nSide)));
                    nFactor.y            = ytracks(side==nSide)'*abs(lVec(side==nSide))/sum(abs(lVec(side==nSide)));
                    nFactor.z            = ztracks(side==nSide)'*abs(lVec(side==nSide))/sum(abs(lVec(side==nSide)));
                    mDFF                 = DFF(mFactor, :)';
                    nFactor.DFF          = mDFF;
                    nFactor.neuronIndex  = find(lVec~=0 & side==nSide);
                    neuronCCMat          = nan(numNeuron, 2);
                    for nNeuron          = find(lVec~=0 & side==nSide)'
                        [xcf, lags, bound]   = crosscorr(mDFF, slicedDFF(:,nNeuron), maxLags, numStd);
                        [pks, locs]      = findpeaks(xcf, lags,'MinPeakDistance', 5,'MinPeakHeight', bound(1), 'MinPeakProminence', bound(1));
%                         [pks, locs]         = max(xcf);
                        pklocs           = abs(locs) == min(abs(locs));
                        locs             = locs(pklocs);
                        pks              = pks(pklocs);
                        [~,  pklocs]     = max(pks);
                        locs             = locs(pklocs);
                        pks              = pks(pklocs);
                        if ~isempty(locs)
                            if locs>=0
                                [~, corrSig]        = corr(mDFF(1:end-locs), slicedDFF(1+locs:end,nNeuron));
                            else
                                [~, corrSig]        = corr(mDFF(1-locs:end), slicedDFF(1:end+locs,nNeuron));
                            end
                            if corrSig<0.01
                                neuronCCMat(nNeuron, 1)     = locs/fs;
                                neuronCCMat(nNeuron, 2)     = pks;
                            end
                        end
                    end
                    nFactor.neuronCCMat  = neuronCCMat;
                    factorSet        = [factorSet, nFactor];
                end
            end
        end

        % among factors
        if length(factorSet)>1
            factorsMat.corrMat   = nan(length(factorSet)-1);
            factorsMat.delayMat  = nan(length(factorSet)-1);
            for pFactor          = 1:length(factorSet)-1
                for qFactor      = pFactor+1:length(factorSet)-1
                    pFactorDFF   = factorSet{pFactor+1}.DFF;
                    qFactorDFF   = factorSet{qFactor+1}.DFF;
                    if ~isequal(pFactorDFF, qFactorDFF)
                        [xcf, lags, bound] = crosscorr(pFactorDFF, qFactorDFF, maxLags, numStd);
                        [pks, locs]      = findpeaks(xcf, lags,'MinPeakDistance', 5,'MinPeakHeight', bound(1), 'MinPeakProminence', bound(1));
    %                     [pks, locs]         = max(xcf);
                        pklocs           = abs(locs) == min(abs(locs));
                        locs             = locs(pklocs);
                        pks              = pks(pklocs);
                        [~,  pklocs]     = max(pks);
                        locs             = locs(pklocs);
                        pks              = pks(pklocs);
                        if ~isempty(locs)
                            if locs>=0
                                [~, corrSig]        = corr(pFactorDFF(1:end-locs), qFactorDFF(1+locs:end));
                            else
                                [~, corrSig]        = corr(pFactorDFF(1-locs:end), qFactorDFF(1:end+locs));
                            end
                            if corrSig<0.01
                                factorsMat.delayMat(pFactor, qFactor)    = locs/fs;
                                factorsMat.corrMat(pFactor, qFactor)     = pks;
                            end
                        end
                    end
                end
            end
            factorSet        = [factorSet, factorsMat];
        end

        networkMat{nTime}    = factorSet;
    end

    save([tempDatDir, 'EvoLoading_' fileName, '.mat'], 'networkMat', 'neuronXLoc', 'neuronYLoc', 'neuronZLoc');
end
