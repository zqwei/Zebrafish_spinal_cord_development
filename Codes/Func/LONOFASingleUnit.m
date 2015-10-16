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






function DFFEst = LONOFASingleUnit (nFoldDFFTest, lambda, psi, nUnit)


    LONODFF      = nFoldDFFTest(:, [1:nUnit-1 nUnit+1:end]);
    L            = lambda([1:nUnit-1 nUnit+1:end], :); % yDim -1 x xDim
    Ph           = psi([1:nUnit-1 nUnit+1:end]);
    
    Currlambda   = lambda(nUnit, :);
    
    estX         = L'/(L*L'+diag(Ph)) * LONODFF';
    % xDim x n
    % This method of estX is equivalent to F from following code
    % [lambda,psi, ~, ~, F] = factoran(data, m ,'scores','regression');
    
    DFFEst       = (Currlambda * estX)'; 
    

end