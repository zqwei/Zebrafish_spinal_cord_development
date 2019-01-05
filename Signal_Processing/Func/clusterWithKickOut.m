% 
% Grouping data with kicking out low correlated units
%
% Input:
% 
% X             ---  N x m
% nGroup        ---  number of groups
% minSizeGroup  ---  a minimal size of group (a group having units less than this would discard)
%
% Output:
% groupInd      --- value from 0 to nGroup, where 0 means kicking-out units
% 
% -------------------------------------------------------------------------
% version 1.0
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 

function groupInd  = clusterWithKickOut(X, nGroup, minSizeGroup)
    
    numUnit        = size(X, 1);
    groupInd       = nan(numUnit, 1);
    validSet       = find(isnan(groupInd));
    
%     imagesc(corr(X(validSet,:)'),[-1 1]) 
%     pause()
    
    while ~isempty(validSet)    
        linkNeurons    = linkage(X(validSet,:),'single','correlation');
        tGroupInd      = cluster(linkNeurons,'maxclust',nGroup);
        lengthGroup    = cell2mat(arrayfun(@(tGroup) sum(tGroupInd==tGroup),...
                                            1:nGroup,'UniformOutput',false));
%         imagesc(corr(X(validSet,:)'),[-1 1])
%         pause()
        if sum(lengthGroup < minSizeGroup)
            kickoutInd                 = find(lengthGroup < minSizeGroup);
            for nKickOut               = 1:length(kickoutInd)
                invalidSet             = validSet(tGroupInd == kickoutInd(nKickOut));
                groupInd(invalidSet)   = 0;
            end
        else
            groupInd(validSet)         = tGroupInd;
        end
        
        validSet                       = find(isnan(groupInd));        
    end

end