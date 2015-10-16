%% Loading data and preprocessing
load ../Data/Data_Dre_E1_BTXinjHuCH2BGCaMP6f_TL_20140818_045650_corrected_signal/Data/dff.mat

IdentifySimilarNeuronsFromActivity;
useDFF = newDff';
clear newDff;

%%  Run analysis
seed = 1;
rng(seed);

timePointNum = size(useDFF,1);
neuronNum = size(useDFF,2);

timePeriodLength = 5000;
timePointLim = 1:timePeriodLength:timePointNum;
% timePointLim = flipdim(timePointLim,2);
timePeriodNum = length(timePointLim)-1;

finalTime = timePointLim(end-1):timePointLim(end);
figure;
Z = linkage(transpose(useDFF(finalTime,:)),'single','correlation');
[H,T,outperm] = dendrogram(Z,0);
useDFF = useDFF(:,outperm);

%% Two factor model over time

factorNum = 2;
lambdaArray = zeros(neuronNum,factorNum,timePeriodNum);
psiArray = zeros(neuronNum,timePeriodNum);
sampleCovArray = zeros(neuronNum,neuronNum,timePeriodNum);
communalityArray = zeros(neuronNum,timePeriodNum);

trainingFract = 0.5;
trainingPointNum = ceil(timePeriodLength/2);

trainingLogLikeVec = zeros(1,timePeriodNum);
matlabFunctionTrainingLogLikeVec = zeros(1,timePeriodNum);
testingLogLikeVec = zeros(1,timePeriodNum);

for ii=1:(length(timePointLim)-1)
  display(['Calculating models for time point index ' num2str(ii)])

  timeIndex = timePointLim(ii):timePointLim(ii+1);
  normalizedData = useDFF(timeIndex,:);
  normalizedData = normalizedData-repmat(mean(normalizedData),size(normalizedData,1),1);
  normalizedData = normalizedData./repmat(std(normalizedData),size(normalizedData,1),1);
  
  [lambdaArray(:,:,ii),psiArray(:,ii),T,stats,F] = factoran(normalizedData,factorNum);
  
  sampleCovArray(:,:,ii) = cov(normalizedData);
  r = randperm(timePeriodLength);
  
  communalityArray(:,ii) = sum(squeeze(lambdaArray(:,:,ii)).^2,2);

  trainingData = normalizedData(r(1:trainingPointNum),:);
  trainingData = trainingData-repmat(mean(trainingData),size(trainingData,1),1);
  trainingData = trainingData./repmat(std(trainingData),size(trainingData,1),1);
  
  testData = normalizedData(r((trainingPointNum+1):end),:);
  testData = testData-repmat(mean(testData),size(testData,1),1);
  testData = testData./repmat(std(testData),size(testData,1),1);
  
  [trainLambda,trainPsi,T,stats,F] = factoran(trainingData,factorNum);
  matlabFunctionTrainingLogLikeVec(ii) = stats.loglike;
  
  sigma = trainLambda*trainLambda'+diag(trainPsi);
  STrain = cov(trainingData);
  STest = cov(testData);
  
  trainingLogLikeVec(ii) = -(log(det(sigma)) + trace(STrain/sigma) - log(det(STrain)) - neuronNum);
  testingLogLikeVec(ii) = -(log(det(sigma)) + trace(STest/sigma) - log(det(STest)) - neuronNum);
end

overTimeLogLikeMat = zeros(timePeriodNum,timePeriodNum);
for ii=1:timePeriodNum
  sigma = squeeze(lambdaArray(:,:,ii))*squeeze(lambdaArray(:,:,ii))'+diag(psiArray(:,ii));
  for jj=1:timePeriodNum
%     FALogLike = log(det(sigma)) + trace(squeeze(sampleCovArray(:,:,jj))*inv(sigma)) - log(det(squeeze(sampleCovArray(:,:,jj)))) - neuronNum;
    FALogLike = log(det(sigma)) + trace(squeeze(sampleCovArray(:,:,jj))/sigma) - log(det(squeeze(sampleCovArray(:,:,jj)))) - neuronNum;

    overTimeLogLikeMat(ii,jj) = -FALogLike;
  end
end

% ll = flip(flip(overTimeLogLikeMat,1),2);
% ll = flipdim(flipdim(overTimeLogLikeMat,1),2);
ll = overTimeLogLikeMat;

colorVec = [1 0 0];
colorDiff = 1/(timePeriodNum-1);
colorMat = zeros(timePeriodNum,3);
for ii=1:timePeriodNum
  colorMat(ii,:) = colorVec+[-1 0 1]*(ii-1)*colorDiff;
end

figure('position',[350 700 850 650]); hold on;
for ii=1:timePeriodNum
  plot(ll(ii,:),'color',colorMat(ii,:),'linewidth',2);
end
scatter(1:timePeriodNum,diag(ll),ones(timePeriodNum,1)*60,colorMat,'fill');

xt = get(gca,'xtick');
if sum(ismember(xt,0));
  xt(1) = 1;
end
set(gca,'xtick',xt);
set(gca,'xticklabel',timePointLim(xt));
xlabel(['Model time point. Time length ' num2str(timePeriodLength)]);
ylabel('Log likelihood')
title(['Across time factor analysis. Factor number is ' num2str(factorNum)]);

figure, imagesc(communalityArray,[0 1])
set(gca,'xticklabel',timePointLim(get(gca,'xtick')));
title('Across time fraction variance explained by other neurons');
xlabel(['Model time point. Time length ' num2str(timePeriodLength)]);
ylabel('Neuron index')

%% Factor number analysis
% I still don't have a good way to do this

% factorNumVec = 1:15;
% modelNum = length(factorNumVec);
% dfMat = zeros(timePeriodNum,modelNum);
% llMat = zeros(timePeriodNum,modelNum);
% chisqMat = zeros(timePeriodNum,modelNum);
% for ii=11:(length(timePointLim))
%   display(['Calculating models for time point index ' num2str(ii)])
%   timeIndex = timePointLim(ii):timePointLim(ii+1);
%   normalizedData = useDFF(timeIndex,1:42);
%   normalizedData = normalizedData-repmat(mean(normalizedData),size(normalizedData,1),1);
%   normalizedData = normalizedData./repmat(std(normalizedData),size(normalizedData,1),1);
%   for jj=1:length(factorNumVec)
%     currentModelDF = ceil(0.5*((neuronNum-factorNumVec(jj)).^2-(neuronNum+factorNumVec(jj))));
%     [lambda,psi,T,stats,F] = factoran(normalizedData,factorNumVec(jj),'maxit',1000,'rotate','promax');
%     dfMat(ii,jj) = stats.dfe;
%     llMat(ii,jj) = stats.loglike;
%     if ~isfield(stats,'chisq')
%       chisqMat(ii,jj) = -llMat(ii,jj)*(timePeriodLength-1);
%     else
%       chisqMat(ii,jj) = stats.chisq;
%     end    
%   end
% end