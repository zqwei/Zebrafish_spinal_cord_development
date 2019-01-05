%
%
% Explained variance from leave-one neuron-out using Factor analysis
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


function EV = LONOFA(nFoldDFFTest, lambda, psi)
    
    numUnit     = size(nFoldDFFTest, 2);
    nFoldDFFEst = zeros(size(nFoldDFFTest));
    
    for nUnit = 1:numUnit
         nFoldDFFEst(:,nUnit) = LONOFASingleUnit (nFoldDFFTest, lambda, psi, nUnit);       
    end
    
    EV        = 1 - sum(sum((nFoldDFFTest-nFoldDFFEst).^2))/sum(var(nFoldDFFTest,[],1)*size(nFoldDFFTest,1));

end