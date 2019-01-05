%
%
% Explained variance from leave-one neuron-out using Factor analysis for
% each unit
% 
%
% Ziqiang Wei
% version 1.0
%
%
% Input
%
% nFoldDFFTest  --- n x yDim
%


function EVSingleUnit = LONOFASingleUnitSingleFactorEV(nFoldDFFTest, lambda, psi)
    
    numUnit           = size(nFoldDFFTest, 2);
    EVSingleUnit      = zeros(numUnit, size(lambda,2)); % numUnit, numFactor
    numTime           = size(nFoldDFFTest,1);
    
    for nUnit         = 1:numUnit
         DFFEstPerFactor = LONOFASingleUnitSingleFactor (nFoldDFFTest, lambda, psi, nUnit);
         nDFF            = nFoldDFFTest(:, nUnit);
         EVSingleUnit(nUnit, :) = 1 - sum((bsxfun(@minus, nDFF, DFFEstPerFactor)).^2, 1)/(var(nDFF)*(numTime-1));
    end
    
%     EVSingleUnit = 1 - sum((nFoldDFFTest-nFoldDFFEst).^2)./(var(nFoldDFFTest,[],1)*size(nFoldDFFTest,1));
    
end


% function DFFEst = LONOFASinglueUnit (nFoldDFFTest, lambda, psi, nUnit)
% 
% 
%     LONODFF      = nFoldDFFTest(:, [1:nUnit-1 nUnit+1:end]);
%     L            = lambda([1:nUnit-1 nUnit+1:end], :); % yDim -1 x xDim
%     Ph           = psi([1:nUnit-1 nUnit+1:end]);
%     
%     Currlambda   = lambda(nUnit, :);
%     
%     estX         = L'/(L*L'+diag(Ph)) * LONODFF';
%     % xDim x n
%     % This method of estX is equivalent to F from following code
%     % [lambda,psi, ~, ~, F] = factoran(data, m ,'scores','regression');
%     
%     DFFEst       = (Currlambda * estX)'; 
%     
% 
% end