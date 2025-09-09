%% Extract the ISI 

% INPUT:
%   allspks - matrix with [spiketimes(ms) electrode_nr]

% OUTPUT:
%   out - struct containing:
%       time_ms: bin centers, in msec
%       isi: all inter-spike intervals computed, in msec
%       isiDist: The ISI distribution per bin, as probability

% Aug 2025, Alessio Di Clemente & Adam Armada-Moreira

% TODO: Add "method" to consider linear and log binning scales

function out = computeISI(allspks)

isi = [];
electrodes = unique(abs(allspks(:,2)));

for el = 1:numel(electrodes)
    spksPerElectrode = allspks(abs(allspks(:,2)) == electrodes(el), 1);
    el_isi = diff(spksPerElectrode); % Calculate ISI for the current electrode
    isi = vertcat(isi, el_isi);
end % for el = 1:numel(electrodes)


% Define the bin edges for the histcount
toplimit = ceil(max(isi));
btmlimit = floor(min(isi));

% Get order of magnitude of limit
topOrderOfMagnitude = floor(log10(toplimit));
btmOrderOfMagnitude = floor(log10(btmlimit));

% log edges
log_edges = logspace(btmOrderOfMagnitude, topOrderOfMagnitude, 20*(topOrderOfMagnitude-btmOrderOfMagnitude));

% Create histograms for logarithmic edges
logHist = histcounts(isi, log_edges, 'Normalization','probability');

binCenters = [];
for i = 1:numel(log_edges)-1
    binCenters(i) = (log_edges(i+1)+log_edges(i)) / 2;
end % for i = 1:numel(log_edges)-1

% output data
out.isi = isi;
out.isiDist = logHist;
out.time_ms = binCenters;


end % function out = computeISI(allspks, threshold, scale)