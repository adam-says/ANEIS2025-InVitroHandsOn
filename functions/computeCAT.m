%% Center of Activity Trajectory (CAT) analysis

% INPUTS:
%   t1 - in ms, the beginning of the analysis window
%   t2 - in ms, the ending of the analysis window
%   allspks - matrix with [spiketimes(ms) electrode_nr]
%   binningWindow - Binning window (in ms) for the analysis
%   doPlot - set to 1 if you want the plot, set to 0 otherwise

% Apr 2025, Adam Armada-Moreira

function CAT = computeCAT(t1, t2, allspks, metafile, binningWindow)
% check if time points are well defined
if t2 < t1
    warning('The ending timepoint is smaller than the starting!');
    CAT = [];
    return
end
% Make the STH according to the defined bin size
time_bins = t1:binningWindow:t2;
activeElectrodes = unique(abs(allspks(:,2)));
nrElectrodes = numel(activeElectrodes);
% TEMPORARY HARD-CODED DIMENSION - works with 120ch MCS MEA
electrodeSpace = NaN(12,12);
for ch=1:nrElectrodes
    [elec_x, elec_y] = get_electrode_xy(activeElectrodes(ch),metafile);%TODO: Dictionary with xy location for each electrode,
    %      in integers (1,1), (1,2),...
    mask = abs(allspks(:,2)) == activeElectrodes(ch);
    local_spk_times = allspks(mask,1);
    electrodeSpace(elec_y,elec_x) = activeElectrodes(ch);
    spatialSTH(elec_y,elec_x,:) = histcounts(local_spk_times,time_bins); % output matrix: y * x * time
end %end STH loop

x_dim = size(spatialSTH,2);
y_dim = size(spatialSTH,1);
[x, y] = meshgrid(1:x_dim,1:y_dim);

% TEMPORARY HARD-CODED VARIABLE
ELECTRODE_DISTANCE = 0.1;
x = x * ELECTRODE_DISTANCE; 
y = y * ELECTRODE_DISTANCE; 


for t = 1:size(spatialSTH,3)
    frame = spatialSTH(:,:,t);
    weightedx = x .* frame;
    weightedy = y .* frame;
    CAx(t) = sum(weightedx(:)) / sum(frame(:));
    CAy(t) = sum(weightedy(:)) / sum(frame(:));
end

CAT.x = CAx;
CAT.y = CAy;
CAT.t = time_bins(1:end-1);
CAT.electrodeSpace = electrodeSpace;
CAT.spatialSTH = spatialSTH;

end %end function