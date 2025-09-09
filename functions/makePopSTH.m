%% Compute the spiketime histogram

% INPUTS:
%   allspks - matrix with [spiketimes(ms) electrode_nr]
%   meta - metadata file created by opto.loadDataset
%   bin - binning window for histogram generation, in msec
%   method - normalization method for the histogram
%       'count': simple count of spikes in each bin
%       'fr': count/bin, general firing rate in each bin
%       'normfr': firing rate normalized by the number of overall active
%       electrodes

% OUTPUT:
%   STH - a struct containing:
%       bin: binning window used
%       fs: sampling rate for the histogram (Hz)
%       time: time vector for the histogram (msec)
%       method: normalization method used
%       network: population-wide spiketime histogram

% Aug 2025, Adam Armada-Moreira

function STH = makePopSTH(allspks,meta,bin,method)

% Make spike time histogram
STH = [];
STH.bin = bin; % in msec, default 1
STH.fs = 1/(STH.bin/1000);
STH.time = 0:STH.bin:(meta.duration_ms);
STH.method = method;


% Exclude inactive electrodes (those with a firing rate below 0.02 Hz)
allElectrodes = unique(abs(allspks(:,2)));
spkPerElectrode = histcounts(abs(allspks(:,2)),(0.5:1:120.5));
threshold = 0.02;
activeElectrodes = allElectrodes(spkPerElectrode > (meta.duration_s*threshold));
NactiveElectrodes = numel(activeElectrodes);

STH.network = histcounts(allspks(:,1),STH.time);
STH.time = STH.time(1:end-1);

switch method
    case 'count'
        STH.network = STH.network;
    case 'fr'
        STH.network = STH.network ./ (0.001*bin); % Get firing rate in Hz
    case 'normfr'
        STH.network = STH.network ./ (0.001*bin*NactiveElectrodes);
    otherwise
        error("Calculation method not recognized! Use 'count', 'fr', or 'normfr'.")
end % switch method


% TODO: calculate the number of active electrodes for each bin
end % function STH = makePopSTH(in,meta,bin,method)
