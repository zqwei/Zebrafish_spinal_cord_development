%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Signal Extraction Step 3: Generating MaMuT XML for manual curation
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%
function Signal_Extraction_s3(inputFolder, ratio3D)
addpath('../Func');
trackFile = [inputFolder '.simpleTrack\tracking\tracks.mat'];
selectFile = [inputFolder '.simpleTrack\active_neuron_extraction\pre_selection.mat'];
dataXML = [inputFolder '.simpleTrack\signal\data.xml'];

load(trackFile);
load(selectFile);
tracks_selected = tracks;
for i = 1:numel(tracks)
    tracks_selected{i} = tracks{i}([cells_common; cells_active; cells_pattern], :);
end
terminationTimepoint = terminationTimepoint([cells_common; cells_active; cells_pattern]);
swc = exportSWC(tracks_selected, timepoints, terminationTimepoint);

swc(:, [3 4 5]) = swc(:, [4 3 5]) - 1;

writeMamutXML(swc, [dataXML '_autotracks.xml'], dataXML, [1 1 ratio3D]);

