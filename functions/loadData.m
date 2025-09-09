%% Aux function to load SpiQ spiketime data into a MATLAB database

% INPUT:
%       datadir - Directory ending in ".dat" outputted from SpiQ

% OUTPUT:
%       allspks - Matrix with spiketimes and electrodes
%       metadata - Struct containing:
%           chan_names: labeled names of each electrode
%           fs: sampling rate
%           nsamples: total number of acquired samples
%           duration_s: recording duration in s
%           duration_ms: recording duration in ms

% Aug 2025, Adam Armada-Moreira

function [allspks metadata] = loadData(datadir)

% Initialize the output variables
allspks = [];
metadata = struct();

if ~isempty(datadir)
    if endsWith(datadir, '.dat')
        allspks = readmatrix(fullfile(datadir,"spk.txt"));
        allspks(:,1) = allspks(:,1) * 1000; % SpiQ output is in s, this converts it to ms

        metafile = fullfile(datadir,'meta.toml');
        toml_tmp = toml.read(metafile);
        tmp2 = toml.map_to_struct(toml_tmp);
        metadata.chan_names = tmp2.channels_names;
        metadata.fs = double(tmp2.info.Nsamples) / (double(tmp2.info.Data_Recording_0_Duration) / 1000000);
        metadata.nsamples = double(tmp2.info.Nsamples);
        metadata.duration_s = metadata.nsamples / metadata.fs;
        metadata.duration_ms = metadata.duration_s * 1000;
    else
        return
    end % if endsWith(datadir, '.dat')
end % if ~isempty(datadir)
end % function