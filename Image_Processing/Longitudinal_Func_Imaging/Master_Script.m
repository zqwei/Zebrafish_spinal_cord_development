%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Master Script for processing longitudinal functional imaging data
% Step-by-step guide can be found on wiki
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
inputFolder = 'Q:\SV2\YW_17-09-25\Dre_E1_HuCH2BGC6f-mnxTagRFP-H2BhaloTagJF635-0-1-2_func_20170926_033034.corrected';
timepoints = 0:64000; %range of time points to be analyzed

backgroundThres = 110; % background intensity for image segmentation
aspectRatio = 6; % aspect ratio of image stacks
preselectTimeRange = 48000; % for pre-selection: time window from the end to extract signal
FARange = 30001:30001+1200; % for pre-selection: time window to run factor analysis

dsFactor = 200; % factor for resampling in time
FAThres = 0.2; % for pre-selection: threshold of factor loading
rlowessWindow = 5; % number of time points to smooth curated tracks


inputFilePattern = [inputFolder '\SPM00\TM??????\SPM00_TM??????_CM00_CHN00.klb'];
resultFolder = [inputFolder '.simpleTrack\'];

%% Part 1: Pre-processing
Pre_Processing_s1(inputFilePattern, resultFolder, timepoints, dsFactor); % Downsampling data in time
Pre_Processing_s2(resultFolder, backgroundThres); % Automated image segmentation
% ** Manual Step **
%  Manual curation of automated segmentation result using MTrackJ, save result as Period_??????_curated.mdf in [resultFolder '\signal']

%% Part 2: Cell tracking and signal extraction
Cell_Tracking_s1(resultFolder, aspectRatio); % Automated cell tracking
Signal_Extraction_s1(inputFolder, aspectRatio, preselectTimeRange); % Extraction of automated cell tracking signal for pre-selectoin
Signal_Extraction_s2(inputFolder, FARange, FAThres) % Pre-selection based on activation test and factor analysis
Signal_Extraction_s3(inputFolder, aspectRatio); % Generate MaMuT Files of pre-selected tracks
% ** Manual Step **
% Generate data.xml using MaMuT in [resultFolder '\signal']
% manual curation of selected tracks using MaMuT, save as data-mamut.xml in [resultFolder '\signal']

Signal_Extraction_s4(inputFolder, aspectRatio, rlowessWindow); % Interpolation and smooth curated MaMuT tracks
Signal_Extraction_s5(inputFolder, aspectRatio); % Extraction of intensity from curated MaMut tracks

%% Part 3: Atlas construction and annotation
% ** Manual Step **
% 1. Mark motor nerve roots on image and store them as Atlas.tif_resampled.marker (Vaa3D) in [resultFolder '\signal']
% 2. Annotate mnx+/- cells based on channel 01, append 'mnx' variable to profile.mat in [resultFolder '\profile']
% 3. (optional) If islet staining data is available, append 'islet' variable to profile.mat in [resultFolder '\profile']
% 4. (optional) If birth time annotation is avialble, append 'birthtime' variable to profile.mat in [resultFolder '\profile']
Atlas_Annotation_s1(resultFolder);% Calculate AP, LR, DV coordiantes on atlas
