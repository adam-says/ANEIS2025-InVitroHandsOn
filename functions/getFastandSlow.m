% Separate fast and slow components of the burst decay via a montecarlo-like process
%
% INPUT:
% 1) bursts -> table with info on all detected bursts. in particular needs
% to have variables 'burst_peak_ms','burst_end_ms','burst_durations_ms',
% 'TimeToPeak_ms'
% 2) spiketrain -> vector with all spike times from the entire recording
% 3) binsize -> width of the temporal bin to use for burst profiling
%
% OUTPUT:
% A matrix containing only fast oscillations and another one with slow
% oscillations


function [fastOscillations,slowOscillations] = getFastandSlow(bursts,spiketrain,binsize)

    N_bursts = height(bursts);
    spiketrain = spiketrain(:,1);

    burst_starts = bursts.burst_start_ms;
    burst_ends = bursts.burst_end_ms;

    plot_ext_ms = 150;
    burst_max_ext = max(bursts.burst_durations_ms)/binsize;
    burst_max_ext = round(burst_max_ext+2*(plot_ext_ms/binsize));

    rep = 150;

    allmeanfakebursts = zeros(N_bursts,burst_max_ext);

    burst_profiles = zeros(N_bursts,burst_max_ext);

    for n = 1:height(bursts)

        if burst_starts(n) - plot_ext_ms <= 0
            ext_start = burst_starts(n);
            ext_end = burst_ends(n)+plot_ext_ms;

        elseif burst_ends(n) + plot_ext_ms >= max(spiketrain)
            ext_start = burst_starts(n)-plot_ext_ms;
            ext_end = burst_ends(n);

        else
            ext_start = burst_starts(n)-plot_ext_ms;
            ext_end = burst_ends(n)+plot_ext_ms;

        end

        burst_spikes = spiketrain(spiketrain >= ext_start & spiketrain <= ext_end);
        edges = ext_start:binsize:ext_end;
        burst_profile = histcounts(burst_spikes,edges);

        burst_spikes_repmat = repmat(burst_spikes,1,rep);
        shift_mat = normrnd(0,15,size(burst_spikes_repmat));
        fake_bursts = burst_spikes_repmat+shift_mat;

        fakeBurst_profiles = zeros(rep,length(burst_profile));

        for f = 1:rep               
            fakeBurst_profiles(f,:) = histcounts(fake_bursts(:,f),edges);     
        end

        meanfakeburst = mean(fakeBurst_profiles,1);

        meanfakeburst(end+1:burst_max_ext) = 0;
        burst_profile(end+1:burst_max_ext) = 0;

        allmeanfakebursts(n,:) = meanfakeburst;
        burst_profiles(n,:) = burst_profile;

    end

    fastOscillations = burst_profiles-allmeanfakebursts;
    slowOscillations = allmeanfakebursts;

end