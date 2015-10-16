% 
% List group index from 1 to nGroup
%
% 
% -------------------------------------------------------------------------
% version 1.0
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 

function regroupInd = reGroupIndex(groupInd, nGroup)
    
    regroupInd      = [];

    for tGroup      = 1:nGroup
        regroupInd  = [regroupInd; find(groupInd==tGroup)]; %#ok<AGROW>
    end
    

end