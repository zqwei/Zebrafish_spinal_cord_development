% 
% Data Analysis List
% 
% -------------------------------------------------------------------------
% version 1.0
%
% Analysis list
%
% 1.  Covariance analysis: factor analysis
% 2.  Time series analysis: LDS analysis
% 3.  Cluster in the LDS connectivity matrix
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.  Covariance analysis: clustering analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data_Analysis_List1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.  Covariance analysis: factor analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data_Analysis_List2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.  Time series analysis: LDS analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data_Analysis_List3
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4.  Cluster in the LDS connectivity matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data_Analysis_List4

function Data_Analysis_List0(nFile)

    addpath('../Func');

    % load data
    setDir;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 0.1 Loading data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % dirResult         = [fileDirName, '/Fitted_result/'];
    
    fileDirName       = fileDirNames{nFile};
    fileName          = fileNames{nFile};

    dirImageData      = [fileDirName '/']; % [fileDirName, '/Data/'];
    % dff
    % timepoints
    % tracks_smoothed
    load ([dirImageData 'dff.mat']);
    load ([dirImageData 'clustering.mat'])
    % idx (from PCA result)
    % load ([dirImageData 'clustering.mat']);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 0.2 Preprocessing the data and ruling out the duplicated units
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % neuronNum = size(dff,1);
    % timePointNum = size(dff,2);
    distNeurons               = linkage(dff,'single','correlation');
    % 
    % Below this value, we consider that a neuron could be a duplicated unit in the
    % population,where these neurons have a strong correlation averaged across
    % time.
    minDist                   = 0.15; 
    % We use the average value of these duplicated unit as one unit
    %
    % 1. Grouping the duplicated units into subgroups
    duplicatedUnitsSubGroup   = cluster(distNeurons, 'cutoff', minDist);
    % Index of subgroups
    % 
    numGroups                = length(unique(duplicatedUnitsSubGroup));
    % Check if merging units are reasonable in space.
    % -- Visually
    % figure;
    % hold on;
    % cJet                      = jet();
    % nMerge                    = 1;
    % for nNeuron               = 1:numGroups
    %     if sum(duplicatedUnitsSubGroup == nNeuron)>1
    %         plot3Traj(squeeze(mean(tracks_smoothed(duplicatedUnitsSubGroup == nNeuron, 50000:end, :), 2)),'.','color',cJet(nMerge*3,:),'markersize',8);
    %         nMerge            = nMerge + 1;
    %     end
    % end
    % hold off
    % view(3)
    %
    % -- Computing the radius of each subgroub
    % radiusSubGroup            = nan(numNeurons, 1);
    % for nNeuron               = 1:numNeurons
    %     if sum(duplicatedUnitsSubGroup == nNeuron)>1
    %         radiusSubGroup(nNeuron)   = getRadius(squeeze(mean(tracks_smoothed(duplicatedUnitsSubGroup == nNeuron, 50000:end, :), 2)));
    %     end
    % end
    % -- Computing the average distance between two random points in the whole
    % group
    %
    numStim                   = 10000;
    % Perform numStim times computing...
    meanDist                  = distRandomTwoPoint(squeeze(mean(tracks_smoothed(:, 50000:end, :), 2)), numStim);
    % Merging distance threshould (This threshould should be tuned)
    mergingDist               = 0.1 * meanDist;

    % New df/f data
    dffWithNoDuplicates       = nan(size(dff));
    % Index corresponding to the old data set
    indWithNoDuplicates       = cell(size(dff,1),1);
    % new tracks_smoothed
    tracksWithNoDuplicates    = nan(size(tracks_smoothed));
    % new idx
    idxWithNoDuplicates       = nan(size(idx));

    nNeuron                   = 1;
    for nGroup                = 1:numGroups
        if sum(duplicatedUnitsSubGroup == nNeuron)>1 && ...
                getRadius(squeeze(mean(tracks_smoothed(duplicatedUnitsSubGroup == nNeuron, 50000:end, :), 2))) < mergingDist
            dffWithNoDuplicates(nNeuron,:)                            = mean(dff(duplicatedUnitsSubGroup == nNeuron, :),1);
            tracksWithNoDuplicates(nNeuron,:,:)                       = mean(tracks_smoothed(duplicatedUnitsSubGroup == nNeuron,:,:),1);
            indWithNoDuplicates(nNeuron)                              = {find(duplicatedUnitsSubGroup == nNeuron)};
            idxWithNoDuplicates(nNeuron)                              = max(idx(duplicatedUnitsSubGroup == nNeuron));
            nNeuron           = nNeuron + 1;
        else
            numNeuronInGroup  = sum(duplicatedUnitsSubGroup == nNeuron);
            dffWithNoDuplicates(nNeuron:nNeuron+numNeuronInGroup-1,:)      = dff(duplicatedUnitsSubGroup == nNeuron, :);
            tracksWithNoDuplicates(nNeuron:nNeuron+numNeuronInGroup-1,:,:) = tracks_smoothed(duplicatedUnitsSubGroup == nNeuron, :, :);
            idxWithNoDuplicates(nNeuron:nNeuron+numNeuronInGroup-1)        = idx(duplicatedUnitsSubGroup == nNeuron);
            indWithNoDuplicates(nNeuron:nNeuron+numNeuronInGroup-1)        = num2cell(find(duplicatedUnitsSubGroup == nNeuron));        
    %         if numNeuronInGroup > 1
    %             keyboard();
    %         end
            nNeuron           = nNeuron + numNeuronInGroup;
        end
    end

    dffWithNoDuplicates       = dffWithNoDuplicates(1:nNeuron-1,:);
    indWithNoDuplicates       = indWithNoDuplicates(1:nNeuron-1);
    tracksWithNoDuplicates    = tracksWithNoDuplicates(1:nNeuron-1,:,:);
    idxWithNoDuplicates       = idxWithNoDuplicates(1:nNeuron-1);

    % Check again if merging units are reasonable in space.
    % -- Visually
    % figure;
    % cColor                    = ones(1, 3);
    % numMerge                  = sum(cell2mat(...
    %                             arrayfun(@(tIndex) length(indWithNoDuplicates{tIndex})>1,...
    %                             1:length(indWithNoDuplicates), 'UniformOutput', false)));
    % m                         = ceil(sqrt(numMerge));
    % 
    % nMerge                    = 1;
    % for nNeuron               = 1:length(indWithNoDuplicates)
    %     if length(indWithNoDuplicates{nNeuron})>1
    %         subplot(m, m, nMerge);
    %         hold on;
    %         totUnits          = indWithNoDuplicates{nNeuron};
    %         for nUnit         = 1 : length(totUnits)
    %             plot3Traj(squeeze(tracks_smoothed(totUnits(nUnit), :, :)),'-','linewid',1.0,'color',cColor/(length(totUnits)+1)*nUnit);
    %         end
    %         hold off;
    %         view(3) 
    %         xlabel('x');
    %         ylabel('y');
    %         zlabel('z');
    %         nMerge            = nMerge + 1;
    %     end
    % end
    % hold off

    save([tempDatDir, fileName, '.mat'], 'dffWithNoDuplicates', 'indWithNoDuplicates','tracksWithNoDuplicates','idxWithNoDuplicates');
    % setPrint(m*8, m*6, [plotDir, 'Duplicated_Units_', fileName], 'pdf');
end