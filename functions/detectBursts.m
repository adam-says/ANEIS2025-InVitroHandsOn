% Burst detection starting from spike time histogram
% 
% INPUT:
% 1) STH -> struct file with fields 'network'= spike times histogram and
% 'bin' = size of the temporal bin used to produce the histogram
% 2) peakThr -> scalar number, minimum burst peak amplitude expressed as fraction of the
% active electrodes
% 3) detectThr -> scalar number, minimum threshold to consider a bin as
% a candidate member of a burst expressed as fraction of active channels
% 4) wait_time -> scalar number, expressed in ms and used to determine
% burst ends before after the peak. A burst ends when the STH stays below
% detectThr for at least wait_time
% 5) minDuration -> scalar number, minimum burst duration, expressed in ms
% 6) N_activeElec -> number of active electrodes in the recording
%
% OUTPUT:
% table with: burst number, index of start, peak and end of each burst in
% the STH, total number of spikes in the burst, burst duration in ms, burst
% start and end in ms, Time to peak, half decay time and full width at half
% maximum in ms




function bursts = detectBursts(STH,peakThr,detectThr,wait_time,minDuration,N_activeElec)

    binsize = STH.bin;
    spikeHistogram = STH.network;
    min_burst_dist = wait_time/binsize;       
    min_burst_dur = minDuration/binsize;      
    min_burst_peak = N_activeElec*peakThr;    
    detectionThreshold = N_activeElec*detectThr;

    % -firt rough identification of bursts as bins > detectionThreshold
    bursts_bins = find(spikeHistogram>detectionThreshold)';   %--create a variable with all the bins crossing the treshold (expressed as idx relative to bincounts)

    if isempty(bursts_bins) == 1
        disp('no burst detected')
        bursts = 'no burst detected';
        return
    end

    bursts_gap = find(diff(bursts_bins)>min_burst_dist);    %--creates a variable with the burst_bins indexes of all the gaps between a burst and the next one (that is the ends of the bursts)
    idxbursts_start = [1; bursts_gap+1];     %--indexes of the start of the burst from burst_gap
    idxbursts_end = [bursts_gap; length(bursts_bins)];   %--as above but with the ends (add length(burst_bins) or the last one is missing)
    
    %-Remove oscillations where more than 50% of the bins are below
    %threshold
    isABurst = false(length(idxbursts_start),1);
    
    for j = 1:length(idxbursts_start)
    
        burstCandidate = bursts_bins(idxbursts_start(j):idxbursts_end(j));
        binsAbove = length(burstCandidate);
        totBins = (burstCandidate(end)-burstCandidate(1))+1;
        isABurst(j) = binsAbove/totBins >= 0.5;
    
    end
    
    idxbursts_start = idxbursts_start(isABurst);

    if isempty(idxbursts_start) == 1
        disp('no burst detected')
        bursts = 'no burst detected';
        return
    end
    
    idxbursts_end = idxbursts_end(isABurst);
    burst_start_idx = bursts_bins(idxbursts_start);
    burst_end_idx = bursts_bins(idxbursts_end);
    
    % -Remove bursts not reaching "minimum burst duration"
    bursts_duration = (burst_end_idx-burst_start_idx);
    burst_idx = find(bursts_duration>min_burst_dur);

    if isempty(burst_idx) == 1
        disp('no burst detected')
        bursts = 'no burst detected';
        return
    end

    burst_start_idx = burst_start_idx(burst_idx);
    burst_end_idx = burst_end_idx(burst_idx);
    
    
    % -Remove bursts not reaching "minimum burst peak"
    burst_peaks = zeros(length(burst_start_idx),1);

    for p = 1:length(burst_start_idx)
        burst_peaks(p,1) = max(spikeHistogram(burst_start_idx(p):burst_end_idx(p)));
    end

    burst_idx = find(burst_peaks>min_burst_peak);

    if isempty(burst_idx) == 1
        disp('no burst detected')
        bursts = 'no burst detected';
        return
    end

    burst_start_idx = burst_start_idx(burst_idx);
    burst_end_idx = burst_end_idx(burst_idx);

    %- sum up and put data into table
    burst_N = 1:1:length(burst_start_idx);
    burst_N = burst_N';

    amplitude_Nspikes = zeros(length(burst_N),1);
    burst_peak_idx = zeros(length(burst_N),1);
    TimeToPeak_ms = zeros(length(burst_N),1);
    N_spikes = zeros(length(burst_N),1);
    halfDecay_ms = zeros(length(burst_N),1);
    FWHM_ms = zeros(length(burst_N),1);

    for n = 1:length(burst_N)

        N_spikes(n) = sum(spikeHistogram(burst_start_idx(n):burst_end_idx(n)));

        [maxSpikes, maxSpikes_idx] = max(spikeHistogram(burst_start_idx(n):burst_end_idx(n)));
        amplitude_Nspikes(n) = maxSpikes;
        TimeToPeak_ms(n) = maxSpikes_idx*binsize;
        burst_peak_idx(n) = burst_start_idx(n)+maxSpikes_idx-1;

        halfAmplitude = amplitude_Nspikes(n)/2;
        burstDecay = spikeHistogram((burst_peak_idx(n)):burst_end_idx(n));

        if isempty(burstDecay) == 1
            halfDecay_ms(n) = nan;
        else
            [~, half_dec_idx] = min(abs(burstDecay - halfAmplitude));
            halfDecay_ms(n) = half_dec_idx*binsize;
        end

        burstOnset = spikeHistogram(burst_start_idx(n):burst_peak_idx(n));

        if isempty(burstOnset) == 1
            FWHM_ms(n) = nan;
        else
            [~, halfonset_idx] = min(abs(burstOnset - halfAmplitude));
            FWHM_ms(n) = ((length(burstOnset)+half_dec_idx)-halfonset_idx)*binsize;
        end

    end

    burst_start_ms = (burst_start_idx-1)*binsize; % -1 because bin 1 goes from 0 to binsize ms
    burst_end_ms = (burst_end_idx-1)*binsize;
    burst_peak_ms = (burst_peak_idx-1)*binsize;
    burst_durations_ms = burst_end_ms - burst_start_ms;
    idx_StartPeakEnd = [burst_start_idx burst_peak_idx burst_end_idx];


    bursts = table(burst_N, idx_StartPeakEnd, N_spikes, amplitude_Nspikes, burst_durations_ms, burst_start_ms, ...
        burst_peak_ms, burst_end_ms, TimeToPeak_ms, halfDecay_ms, FWHM_ms);


end
