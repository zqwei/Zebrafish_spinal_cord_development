% simple algorithm transfer ca++ to spike events
%
% Ziqiang Wei


function dffs_marks = ca2spike(dffs, activeNeuronTime)
    dffs_mark = ones(240, 1) * activeNeuronTime;
    dffs_mark = dffs_mark(:);
    [peaks, locs] = findpeaks(dffs, 'minPeakHeight', 0.0, 'minPeakDistance', 8);
    dffs_mark_locs = dffs_mark(locs);
    locs = locs(peaks > max([peaks(dffs_mark_locs==0), 0]) + 0.01);
    dffs_marks = false(size(dffs));
    dffs_marks(locs) = true;
    dffs_marks = double(dffs_marks);
end