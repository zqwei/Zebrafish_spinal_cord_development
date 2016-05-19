%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5.  MNX neurons -- active neurons
% 
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function MNX_v0_1(nFile)            
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>    
    load([tempDatDir, fileName, '.mat'], 'mnx', 'side','tracks'); 
    load([tempDatDir, 'EV_' fileName, '.mat'], 'halfEVTime', 'RSquare', 'halfActTime', 'validFitIndex')
    
    RSquareThres      = 0.7;
    meantracks        = squeeze(mean(tracks, 2));
    xTrack            = meantracks(RSquare>RSquareThres, 1);
    yTrack            = meantracks(RSquare>RSquareThres, 2);
    zTrack            = meantracks(RSquare>RSquareThres, 3);
    mnxSig            = mnx(RSquare>RSquareThres);
    halfEVTimeSig     = halfEVTime(RSquare>RSquareThres);
    
%     [~, ~, ~, ~, coeff] = classify([xTrack, yTrack, zTrack, halfEVTimeSig],[xTrack, yTrack, zTrack, halfEVTimeSig],mnxSig);
%     
%     coeff             = coeff(1,2).linear;
%     locs              = [xTrack, yTrack, zTrack] * coeff(1:3);
%     scatter(locs, halfEVTimeSig, [], mnxSig, 'filled');
    figure;
    plot(xTrack(mnxSig==0), halfEVTimeSig(mnxSig==0),'+r', xTrack(mnxSig==1), halfEVTimeSig(mnxSig==1),'ok');
    legend({'MNX-', 'MNX+'})
    legend('location', 'best')
%     colormap(jet)
%     xlabel('Location')
    xlabel('x-axis')
    ylabel('Half EV time (hour)');
    setPrint(8, 6, [plotDir, 'MNXNumNeuronsClassifier_', fileName], 'pdf');
    close all
end