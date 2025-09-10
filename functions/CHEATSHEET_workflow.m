%% ANEIS 2025
% Hands-on topic

clc; clear; close all;
%% =±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±
% Session 1
%  =±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±

%% Load data

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


% Suggestion: use uigetdir so you don't have to type the whole path
datadir = uigetdir; % DIY

[allspks metadata] = loadData(datadir);

%% Do it yourself: Raster plot 
scatter(allspks(:,1) ./ 1000,abs(allspks(:,2)),5,'k','filled')
xlabel('Time (s)')
ylabel('Electrode number')

%% ===========================================
% Identifying and excluding inactive electrodes
%  ===========================================

%% Do it yourself: Identify active electrodes based on the spiketimes
allElectrodes = unique(abs(allspks(:,2)));

spkPerElectrode = histcounts(abs(allspks(:,2)),(0.5:1:120.5));

% Exclude inactive electrodes (those with a firing rate below 0.02 Hz)
threshold = 0.02;
activeElectrodes = allElectrodes(spkPerElectrode > (metadata.duration_s*threshold));
NactiveElectrodes = numel(activeElectrodes);


%% Spiketime histogram

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

%% Do it yourself: use the different methods for the STH and plot them all
% Just the spike count
STH = makePopSTH(allspks, metadata, 5, 'count');

% Firing rate
STH_fr = makePopSTH(allspks, metadata, 5,'fr'); % DIY

% Normalized firing rate
STH_normfr = makePopSTH(allspks, metadata, 5, 'normfr'); % DIY

% Make a plot with the three different STH
tiledlayout(3,1)

nexttile
plot(STH.time ./ 1000, STH.network) %DIY

nexttile
plot(STH_fr.time ./ 1000, STH_fr.network) %DIY

nexttile
plot(STH_normfr.time ./ 1000, STH_normfr.network) %DIY


%% (optional) Do it yourself: Raster plot with STH
tiledlayout(2,1,'TileSpacing','tight')
nexttile
plot(STH_normfr.time ./ 1000, STH_normfr.network,'r')

set(gca,'XTick',[],'YTick',[],'XColor','none','YColor','none','Color','none')

nexttile
scatter(allspks(:,1) ./ 1000,abs(allspks(:,2)),3,'k','filled')
xlabel('Time (s)')
set(gca,'YTick',[],'YColor','none','Color','none')


%% ===========================================
% Network parameters
%  ===========================================

%% Inter-spike interval (ISI)
% Option 1: run this function to get the ISI
% INPUT:
%   allspks - matrix with [spiketimes(ms) electrode_nr]

% OUTPUT:
%   out - struct containing:
%       time_ms: bin centers, in msec
%       isi: all inter-spike intervals computed, in msec
%       isiDist: The ISI distribution per bin, as probability
ISI = computeISI(allspks);

% Option 2: write some code to calculate the ISI

%% Plot ISI
% Since the time intervals given by the previous function are in logspace,
% try to set xscale('log')
figure
plot(ISI.time_ms, ISI.isiDist)
xscale('log')
xlim([0 10000])
set(gca,'XTick', [10 100 1000 10000],'XTickLabel',[10 100 1000 10000])

%% Do it yourself: Extract the main stats from ISI (mean and standard deviation)
meanISI = mean(ISI.isi); %DIY
sdISI = std(ISI.isi); % DIY

%% Based on the plotted distribution, we can identify two predominant ISI intervals 
% To make this easier, we can first smooth our distribution
smoothISIdist = movmean(ISI.isiDist, 6);

% We will find the main peaks in a subset of the ISI distribution, 
% excluding the first elements, where ISI is < 4 ms
[pks,locs] = findpeaks(smoothISIdist(13:end), ISI.time_ms(13:end),...
    'NPeaks', 2,...
    'MinPeakDistance',100,...
    'WidthReference','halfheight');

ISI_pks = pks;
ISI_interval = locs;

%% Now we can determine the area under the curve to compute the cumulative
% probability of finding ISI within the main ranges
% For that, we will use the half-height as our boundaries
for i = 1:numel(ISI_interval)
    pk_idx = find(ISI.time_ms == ISI_interval(i));
    halfheight = ISI_pks(i) / 2;
    range_start = find(ISI.isiDist(1:pk_idx) < halfheight,1,'last');
    range_end = find(ISI.isiDist(pk_idx:end) < halfheight,1,'first');
    range_end = range_end + pk_idx - 1;
    ISI_rangeStart(i) = ISI.time_ms(range_start);
    ISI_rangeEnd(i) = ISI.time_ms(range_end);
    prob = sum(ISI.isiDist(range_start:range_end));
    ISI_probability(i) = prob;
end
%
plot(ISI.time_ms, smoothISIdist)
hold on
scatter(ISI_interval, ISI_pks,'filled')
for i = 1:2
    xline(ISI_rangeStart(i))
    xline(ISI_rangeEnd(i))
    text(ISI_interval(i), ISI_pks(i)/2,...
        sprintf('Peak %i\nISI = %d ms\nPercentage = %d %%',i,round(ISI_interval(i)), round(ISI_probability(i)*100)))
end
xscale('log')
xlim([0 10000])
set(gca,'XTick', [10 100 1000 10000],'XTickLabel',[10 100 1000 10000])


%% Do it yourself: Firing rate (FR)
% Directly
FR_pop = numel(allspks(:,1)) / metadata.duration_s;
FR_pop_norm = FR_pop / NactiveElectrodes;

% As a mean of the FR of each electrode, the result should be similar but
% this strategy allows us to see the FR distribution across the MEA
FRperElectrode = [];
for i=1:numel(activeElectrodes)
    FRperElectrode(i) = numel(allspks(abs(allspks(:,2)) == activeElectrodes(i),1)) / metadata.duration_s;
end
FR_elec_mean = mean(FRperElectrode);
FR_elec_median = median(FRperElectrode);

% Using the STH
FR_sth = makePopSTH(allspks, metadata, 100, 'normfr');
FR_sth_mean = mean(FR_sth.network);
FR_sth_median = median(FR_sth.network);
% Why do we see such a striking difference in the median here, depending on the binning window?
% How do we normalize in relation to the number of active electrodes?

% FR as the inverse of the ISI
FR_isi = (1 / mean(ISI.isi)) * 1000;
% Does this last strategy require normalization?
% Why is the value calculated here slightly different?

% Using the ISI analysis, we can even estimate the mean bursting and
% non-bursting firing rate
FR_burst = 1 / ISI_interval(1) * 1000;
FR_nonburst = 1 / ISI_interval(2) * 1000;


%% =±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±
% Session 2
%  =±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±

%% Burst detection: use the provided function to detect bursts from the STH of the recording

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

detectThr = 0.01;
peakThr = 0.15;
minDuration = 100;
wait_time = 50;

BURST = detectBursts(STH,peakThr,detectThr,wait_time,minDuration,NactiveElectrodes);

%% Do it yourself: check burst detection, extract stats and plot one or more burst parameters

%-plot and check what you detected using BURST variables
burstCheck = figure;
plot(STH.time,STH.network,'b-')
hold on

for p = 1:height(BURST)

    plot(BURST.burst_peak_ms(p),BURST.amplitude_Nspikes(p),'*r')
    plot([BURST.burst_start_ms(p) BURST.burst_start_ms(p)],[0 120],'-g')
    plot([BURST.burst_end_ms(p) BURST.burst_end_ms(p)],[0 120],'-r')
    plot([0 max(STH.time)],[detectThr*NactiveElectrodes detectThr*NactiveElectrodes],'g--')
    plot([0 max(STH.time)],[peakThr*NactiveElectrodes peakThr*NactiveElectrodes],'r--')

end

%-extract statistics and plot variables
meanDuration = mean(BURST.TimeToPeak_ms);
medianDuration = median(BURST.TimeToPeak_ms);
sdDuration = std(BURST.TimeToPeak_ms);
iqrDuration = iqr(BURST.TimeToPeak_ms);

meanAmpl = mean(BURST.halfDecay_ms);
medianAmpl = median(BURST.halfDecay_ms);
sdAmpl = std(BURST.halfDecay_ms);
iqrAmpl =  iqr(BURST.halfDecay_ms);

violins = figure;
violinplot(BURST(:,["TimeToPeak_ms","halfDecay_ms"]))

ax = gca;
ax.TickDir = "out";
ax.FontWeight = "bold";
ax.FontSize = 12;
set(gca,'TickLabelInterpreter', 'none');
ylabel('Time (ms)','FontSize',14)

%% (optional) do it yourself: perform burst detection on the other dataset 
% and investigate the difference in one or more burst parameters
% easiest strategy: run the code twice and save burst table with different
% names
% better strategy: do it programmatically in loop


%% =±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±
% Session 3
%  =±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±

%% BURST Profiling: use the provided function to generate an array with all the bursts profiles 
% (spike time histograms) from a single recording alignining them at their start or at their peak.
% 
% INPUT
% 1) allspks -> column vector where the first column are the spiketimes and
% the second one the ID of the electrode
% 2) bursts -> table with start and end of all the bursts in ms
% (burst_start_ms, burst_end_ms)
% 3) binsize -> size of the temporal bin to use generate burst profiles
% 4) alignment -> string scalar or char vector defining the alignment
% method for the burst profiles: "start" 0 aligned at their start, "peak" =
% aligned at their peak
% 5) show -> if show == 0 don't show figures, if show == 1 show figure and
% pause
%
% OUTPUT:
% 1) an array with all the bursts profiles (spike times histogram
% and their plot
% 2) (optional) a table with the parameters used to generate the profiles
% 3) (optional) an array with the average burst calculated as point-bypoint average
% of all the bursts aligned at their start. Their relative plot will show
% the average one in red

[burst_profiles,profiling_parameters,mean_burst] = burstProfiling(allspks,BURST,STH.bin,"peak",0);

%% Do it yourself see how the average burst and the alignment change with the two methods
%try to put show == 1

%% Do it yourself: smooth burst profiles
% suggestion: use smoothdata function
% try different strategies to smooth burst profiles and plot to see the
% result

% option 1: moving mean or median
smoothedbursts = smoothdata(burst_profiles,2,'movmean',21);
% option 2: smoothing by fitting data with lowess/loess or other models
smoothedbursts = smoothdata(burst_profiles,2,'lowess',21);
% option 3: smoothig with savitzky-golay filter
smoothedbursts = sgolayfilt(burst_profiles',5,21)';

%-plot bursts and smoothed bursts to check, pausing on each one
for b = 1:height(BURST)

    time_ms = (0:STH.bin:(size(burst_profiles,2)*STH.bin)-1);
    plot(time_ms, burst_profiles(b,:),'k-')
    hold on
    plot(time_ms, smoothedbursts(b,:),'r-')
    pause
    close

end

%% Burst slope: use the provided function to calculate burst onset slope.
% do it with bursts normalized to their peak (relative slope) or not
% normalized (absolute slope)
% INPUT:
% 1) smoothedBursts -> a cell array of burst profiles smoothed to avoid
% noisy oscillations
% 2) binsize -> size of the temporal bins used to generate the burst
% profiles
% 3) normalize -> 'yes' or 'no' depending if you want to calculate the
% rising slopes on burst normalized by their amplitude or not
%
% OUTPUT:
% a cell array with the rising slopes (i.e. first order derivatives) of each burst (relative or absolute)
% over time. Each cell hosts the rising slope of one burst from burst profiles


%option 1: without normalization
RiseSlopes = getRiseSlopes(smoothedbursts,STH.bin,'No');

%option 2: normalizing the bursts at their peak
NormRiseSlopes = RiseSlopes(smoothedBursts,2,'yes');

%% Do it yourself: take the maximum of the derivatives and store them in an array

MaxRiseSlopes = cellfun(@max,RiseSlopes);


%% (optional) Do It Yourself: plot burst, smoothed burst, derivative and max slope to check

for j = 1:height(BURST)

    [~,peakIdx] = max(smoothedbursts(j,:));
    slope = RiseSlopes{j};
    [~,maxslopeIdx] = max(slope);

    tiledlayout(2,1,'TileSpacing','tight')
    
    nexttile
    plot(burst_profiles(j,1:peakIdx),'k-')
    hold on
    plot(smoothedbursts(j,1:peakIdx),'r-','Linewidth',1.5)
    plot([maxslopeIdx maxslopeIdx], [0 max(burst_profiles(j,1:peakIdx))],'g--','LineWidth',1.5)
    ax = gca;
    ax.TickDir = "out";
    ax.FontWeight = "bold";
    ax.FontSize = 12;
    set(gca,'TickLabelInterpreter', 'none');
    legend('burst','smoothed burst','max rise slope','Location','northwest')
    xlim([0 length(smoothedbursts(1:peakIdx))])

    nexttile
    plot(slope,'r-','LineWidth',1.5)
    hold on
    xlim([0 length(smoothedbursts(1:peakIdx))])
    pause
    
end

%% Do it yourself: calculate rising slopes in the two dataset and check if they are different
% do it with and without normalization
% again 2 strategies:
% easier: do it 2 times assigning different names to the riseslope
% variables
% better: do it in a loop

%% Burst oscillations: use the provided function to obtain fast and slow oscillations from the bursts
% Separate fast and slow components of the burst decay via surrogate
% spike-train jittering approach
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

[fastOscillations, slowOscialltions] = getFastandSlow(BURST,allspks,STH.bin);

%% Do it Yourself: plot original burst (fast+slow oscillations), fast and slow oscillations

for j = 1:height(BURST)

    fast = fastOscillations(j,:);
    slow = slowOscialltions(j,:);
    burst = slow+fast;
    time_ms = 0:STH.bin:(length(fast)*binsize)-1;

    tiledlayout(2,1,'TileSpacing','tight')

    nexttile
    plot(time_ms,burst,'k-')
    hold on
    plot(time_ms,slow,'r--','LineWidth',1.5)
    xlabel('time (ms)')
    ylabel('# spikes')
    ax = gca;
    ax.TickDir = "out";
    ax.FontWeight = "bold";
    ax.FontSize = 12;
    set(gca,'TickLabelInterpreter', 'none');
 

    nexttile
    plot(time_ms,fast,'b-')
    xlabel('time (ms)')
    ylabel('# spikes')
    ax = gca;
    ax.TickDir = "out";
    ax.FontWeight = "bold";
    ax.FontSize = 12;
    set(gca,'TickLabelInterpreter', 'none');

    pause
 
end

%% Spectrograms and power spectrum density

% example of spectrogram function usage
Fs = STH.fs;
windowLength = 20;
overlap = 10;
nfft = 256;

fast = fastOscillations(b,:);
fast_time_s = (0:binsize:length(fast)*binsize-1)/1000;
[~,Fr,Time,power] = spectrogram(fast,windowLength,overlap,nfft,Fs,'yaxis');
PSD = mean(power,2);
[~,dominantFrequency_idx] = max(PSD);
dominantFrequency = Fr(dominantFrequency_idx);

%--plot
tiledlayout(3,1,'TileSpacing','tight')

nexttile
imagesc(Time, Fr, power);
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title(sprintf('Spectrogram (power) — burst %d', b));
colorbar;
colormap(jet);
ylim([10 100]);

nexttile
plot(fast_time_s,fast,'b-')
xlim([Time(1) Time(end)])
xlabel('Time (s)');
title(sprintf('fast oscillation — burst %d', b));

nexttile
plot(Fr,PSD,'r-')
xlabel('Frequency (Hz)');
ylabel('Power')
title(sprintf('Power Spectrum Density — burst %d', b));

%% Do it yourself: put the code above in a loop to obtain all spectrograms

[~,Fr,Time] = spectrogram(fastOscillations(1,:),windowLength,overlap,nfft,Fs,'yaxis');

nsignals = height(bursts);

spectrograms = zeros(length(Fr),length(Time),nsignals);
PSDs = zeros(nsignals,length(Fr));
dominantFrequencies = zeros(nsignals,1);

for k = 1:nsignals
    
    [~,~,~,power] = spectrogram(fastOscillations(k,:),windowLength,overlap,nfft,fs,'yaxis');
    spectrograms(:,:,k) = power;

    PSDs(k,:) = mean(power,2);
    [~,dominantFrequency_idx] = max(PSDs(k,:));
    dominantFrequencies(k) = Fr(dominantFrequency_idx);

    oscillationTime_s = (0:binsize:length(fastOscillations(k,:))*binsize-1)/1000;

    tiledlayout(3,1,'TileSpacing','tight')

    nexttile
    imagesc(Time, Fr, power);
    axis xy;
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title(sprintf('Spectrogram (power) — burst %d', k));
    colorbar;
    colormap(jet);
    ylim([10 100]);
    
    nexttile
    plot(oscillationTime_s,fastOscillations(k,:),'b-')
    xlim([Time(1) Time(end)])
    xlabel('Time (s)');
    title(sprintf('fast oscillation — burst %d', k));

    nexttile
    plot(Fr,PSDs(k,:),'r-')
    xlabel('Frequency (Hz)');
    ylabel('Power')
    title(sprintf('Power Spectrum Density — burst %d', k));

    pause
    
end


%% (optional) Do it yourself find dominant frequency in the two datasets and compare
% usual two strategies: hard coding or loop

%% =±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±
% Session 4
%  =±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±

%% Intro to spatial analysis, FR per electrode in space

figure
FRMAP = NaN(12,12);

for i=1:numel(activeElectrodes)
    [elec_x, elec_y] = get_electrode_xy(activeElectrodes(i),metadata);
    FRMAP(elec_y, elec_x) = FRperElectrode(i); % DIY
end

%
imagesc(FRMAP, 'AlphaData',~isnan(FRMAP))
colormap(parula)
colorbar
set(gca,'PlotBoxAspectRatio', [1 1 1])

%% Example: CAT analysis for the onset of one population burst
whichBurst = 15;

singleCAT = computeCAT(BURST.burst_start_ms(whichBurst), BURST.burst_peak_ms(whichBurst), allspks, metadata, 5);

%
tiledlayout(3,1)
nexttile
idx1 = find(STH.time == BURST.burst_start_ms(whichBurst));
idx2 = find(STH.time == BURST.burst_peak_ms(whichBurst));
idx3 = find(STH.time == BURST.burst_end_ms(whichBurst));

plot(STH.time(idx1:idx3) - STH.time(idx2), STH.network(idx1:idx3),'k');
hold on
plot(STH.time(idx1:idx2) - STH.time(idx2), STH.network(idx1:idx2), 'r', 'LineWidth',2);

nexttile([2 1])
x_lim = [0.05 1.25];
y_lim = x_lim;

plot_x = [singleCAT.x nan];
plot_y = [singleCAT.y nan];
d = singleCAT.t - singleCAT.t(end);
plot_d = [d nan];

patch(plot_x,plot_y,plot_d, 'EdgeColor', 'interp','LineWidth',1)
hold on
scatter(singleCAT.x, singleCAT.y, 3,'w', 'filled')
colorbar
colormap(jet)
set(gca,...
    'XLim',x_lim,'YLim',y_lim,...
    'Color',[.1 .1 .1],'PlotBoxAspectRatio',[1 1 1],...
    'XColor','none','YColor','none')


%% Do it yourself: calculate the CAT for all the detected bursts
for i=1:height(BURST)
    CAT(i) = computeCAT(BURST.burst_start_ms(i), BURST.burst_peak_ms(i), allspks, metadata, 5);
end

%% Look at initiation and peak locations
nsamples = 3; % Number of frames to average to find the initiation location.
% Setting to 1 will consider only the first non NaN frame

burstInit_x = NaN(1,numel(CAT));
burstInit_y = NaN(1,numel(CAT));
burstPeak_x = NaN(1,numel(CAT));
burstPeak_y = NaN(1,numel(CAT));

for i = 1:numel(CAT)
    if numel(CAT(i).x) > nsamples
        xxs_idx = ~isnan(CAT(i).x);
        xxs = CAT(i).x(find(xxs_idx == 1,nsamples,'first'));
        if numel(xxs) > 1
            xxs = mean(xxs);
        end
        yys_idx = ~isnan(CAT(i).y);
        yys = CAT(i).y(find(xxs_idx == 1,nsamples,'first'));
        if numel(yys) > 1
            yys = mean(yys);
        end

        burstInit_x(i) = xxs;
        burstInit_y(i) = yys;

        % look at the peak CA
        endxxs = CAT(i).x(end);
        endyys = CAT(i).y(end);
        burstPeak_x(i) = endxxs;
        burstPeak_y(i) = endyys; 
    end
end


%%
figure
scatter(burstInit_x, burstInit_y,[],'g','filled')
hold on
scatter(burstPeak_x, burstPeak_y,[],'r','filled')

for i=1:numel(CAT)
    plot([burstInit_x(i) burstPeak_x(i)],[burstInit_y(i) burstPeak_y(i)],'k')
end
xlim([0.1 1.2])
ylim([0.1 1.2])
set(gca,'DataAspectRatio',[1 1 1])


%% 2D density distribution
figure
% Define grid
xgrid=0.1:0.1:1.2;
ygrid=0.1:0.1:1.2;
[x1,y1] = meshgrid(xgrid, ygrid);
% Perform kernel density estimate
% [x y] is actual data, xi is the desired grid points to evaluate
% f is an estimate of the density, ep(:,1) is the X location, ep(:,2) is the y location
xi = [x1(:) y1(:)];
[f,ep]=ksdensity([burstInit_x' burstInit_y'],xi,'Support','positive','BoundaryCorrection','reflection'); % remove the outputs to see a 3D plot of the distribution
% format data in matrix for contourf and plot
X = reshape(ep(:,1),length(xgrid),length(ygrid));
Y = reshape(ep(:,2),length(xgrid),length(ygrid));
Z = reshape(f,length(xgrid),length(ygrid));

imagesc(Z)
hold on
scatter(burstInit_x*10, burstInit_y*10,[],'g','filled','XJitter','density','XJitterWidth',0.3)
colorbar
set(gca,'YDir', 'normal',...
    'DataAspectRatio', [1 1 1],...
    'XColor','none',...
    'YColor','none')


%% Find the peaks
figure
thresh = multithresh(Z(:));
[pks, locs_y, locs_x] = peaks2(Z, 'MinPeakHeight', thresh);

imagesc(Z)
hold on
scatter(locs_x, locs_y,[], 'w','filled')
colorbar
set(gca,'YDir', 'normal',...
    'DataAspectRatio', [1 1 1],...
    'XColor','none',...
    'YColor','none')


%% Isolate the bursts from the top initiation regions
% Create a list of the burst idx from each initiation region, and one list
% for the bursts from other spatial origins

burstSubset = [];
for pk = 1:numel(pks)
    x_range = [locs_x(pk)-1 locs_x(pk)+1];
    y_range = [locs_y(pk)-1 locs_y(pk)+1];

    xs = burstInit_x*10 >= x_range(1) & burstInit_x*10 <= x_range(2);
    ys = burstInit_y*10 >= y_range(1) & burstInit_y*10 <= y_range(2);

    xy_combined = and(xs,ys);
    burstSubset(pk).x = locs_x(pk);
    burstSubset(pk).y = locs_y(pk);
    burstSubset(pk).bursts = find(xy_combined);
end

burstSubset(pk+1).bursts = setdiff([1:height(BURST)], horzcat(burstSubset.bursts));



%% Burst characteristics depending on origin
tiledlayout('flow')

nexttile
for i = 1:numel(burstSubset)
    violinplot(i,BURST.N_spikes(burstSubset(i).bursts))
    hold on
end
ylabel('Nr. spikes')

nexttile
for i = 1:numel(burstSubset)
    violinplot(i,BURST.burst_durations_ms(burstSubset(i).bursts))
    hold on
end
ylabel('Burst duration (ms)')

nexttile
for i = 1:numel(burstSubset)
    violinplot(i,BURST.TimeToPeak_ms(burstSubset(i).bursts))
    hold on
end
ylabel('Time to peak (ms)')

nexttile
for i = 1:numel(burstSubset)
    violinplot(i,BURST.halfDecay_ms(burstSubset(i).bursts))
    hold on
end
ylabel('Half decay (ms)')

nexttile
for i = 1:numel(burstSubset)
    violinplot(i,BURST.FWHM_ms(burstSubset(i).bursts))
    hold on
end

ylabel('FWHM (ms)')

ngroups = numel(burstSubset);
lbls = {};
for i = 1:ngroups-1
    lbls{i} = ['cluster ' num2str(i)];
end
lbls{end+1} = 'outside';
set(findobj(gcf,'type','axes'),'XTick',[1 2 3],'XTickLabel',lbls)
%% Do it yourself
%       Apply the CAT analysis to the burst decay

for i=1:height(BURST)
    CAT(i) = computeCAT(BURST.burst_peak_ms(i), BURST.burst_end_ms(i), allspks, metadata, 5);
end

%% (optional) Do it yourself:
%       Plot the CAT onset and decay together for one burst
whichBurst = 10;
onCAT = computeCAT(BURST.burst_start_ms(whichBurst), BURST.burst_peak_ms(whichBurst), allspks, metadata, 5);
offCAT = computeCAT(BURST.burst_peak_ms(whichBurst), BURST.burst_end_ms(whichBurst), allspks, metadata, 10);

%
tiledlayout(3,2)
nexttile([1 2])
idx1 = find(STH.time == BURST.burst_start_ms(whichBurst));
idx2 = find(STH.time == BURST.burst_peak_ms(whichBurst));
idx3 = find(STH.time == BURST.burst_end_ms(whichBurst));

plot(STH.time(idx1:idx3) - STH.time(idx2), STH.network(idx1:idx3),'k');
hold on
plot(STH.time(idx1:idx2) - STH.time(idx2), STH.network(idx1:idx2), 'r', 'LineWidth',2);
plot(STH.time(idx2:idx3) - STH.time(idx2), STH.network(idx2:idx3), 'b', 'LineWidth',2);

nexttile([2 1])
xlim = [0.05 1.25];
ylim = xlim;
plot_x = [onCAT.x nan];
plot_y = [onCAT.y nan];
d = onCAT.t - onCAT.t(end);
plot_d = [d nan];

patch(plot_x,plot_y,plot_d, 'EdgeColor', 'interp','LineWidth',1)
colorbar
hold on
scatter(onCAT.x, onCAT.y, 3,'w', 'filled')
colormap(jet)
set(gca,...
    'XLim',xlim,'YLim',ylim,...
    'Color',[.1 .1 .1],'PlotBoxAspectRatio',[1 1 1],...
    'XColor','none','YColor','none')

nexttile([2 1])
xlim = [0.05 1.25];
ylim = xlim;
plot_x = [offCAT.x nan];
plot_y = [offCAT.y nan];
d = offCAT.t - offCAT.t(1);
plot_d = [d nan];

patch(plot_x,plot_y,plot_d, 'EdgeColor', 'interp','LineWidth',1)
colorbar
hold on
scatter(offCAT.x, offCAT.y, 3,'w', 'filled')
colormap(jet)
set(gca,...
    'XLim',xlim,'YLim',ylim,...
    'Color',[.1 .1 .1],'PlotBoxAspectRatio',[1 1 1],...
    'XColor','none','YColor','none')

