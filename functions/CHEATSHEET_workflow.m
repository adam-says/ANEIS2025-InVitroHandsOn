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

%% Burst detection
detectThr = 0.01;
peakThr = 0.15;
minDuration = 100;
wait_time = 50;

BURST = detectBursts(STH,peakThr,detectThr,wait_time,minDuration,NactiveElectrodes);

%% Burst profilling

%% =±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±
% Session 3
%  =±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±=±

%% Burst slope

%% Burst oscillations

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



%% Extra: See video of activity

whichBurst = 17;
binningWindow = 5;
fullCAT = computeCAT(BURST.burst_start_ms(whichBurst), BURST.burst_end_ms(whichBurst), allspks, metadata, binningWindow);

% Prepare the frames
STHspace = NaN(12, 12, BURST.idx_StartPeakEnd(whichBurst,3)-BURST.idx_StartPeakEnd(whichBurst,1));
timebin = BURST.burst_start_ms(whichBurst):binningWindow:BURST.burst_end_ms(whichBurst);

for i=1:NactiveElectrodes
    spktimes = allspks(abs(allspks(:,2)) == activeElectrodes(i),1);
    [row col] = find(fullCAT.electrodeSpace == activeElectrodes(i));
    STHspace(row, col, :) = histcounts(spktimes, timebin);
end
%

wait = 0.5; % time between frames, in sec
step = 1;
h = figure;
tiledlayout(h, 2,1,'TileSpacing','tight')
sth_ax = nexttile;
idx1 = find(STH.time == BURST.burst_start_ms(whichBurst));
idx2 = find(STH.time == BURST.burst_end_ms(whichBurst));

plot(STH.time(idx1:idx2) - BURST.burst_peak_ms(whichBurst), STH.network(idx1:idx2),'k');
hold on
pre_xl = xline(sth_ax,STH.time(idx1) - BURST.burst_peak_ms(whichBurst));

set(sth_ax,'PlotBoxAspectRatio',[2 1 1])

cat_ax = nexttile;
hold(cat_ax,'on')
set(cat_ax,'Color',[0 0 0],'XLim',[0.5 12.5],'YLim',[0.5 12.5],...
    'PlotBoxAspectRatio',[1 1 1])

frame_lim = max(STHspace(:));

for i=1:step:size(STHspace,3)
    if ishghandle(h)
        imagesc(cat_ax,STHspace(:,:,i))
        hold on
        scatter(cat_ax,fullCAT.x(i)*10,fullCAT.y(i)*10,[],'r','filled')
        plot(cat_ax,fullCAT.x(1:i)*10,fullCAT.y(1:i)*10,'r')
        pre_xl.Value = STH.time(idx1+i-1) - BURST.burst_peak_ms(whichBurst);

        title(cat_ax,['Time: ' num2str(STH.time(idx1+i-1) - BURST.burst_peak_ms(whichBurst)) ' ms'])
        colormap(cat_ax,gray)
        caxis([0 frame_lim])
        set(cat_ax,'PlotBoxAspectRatio', [1 1 1],...
            'XColor','none','YColor','none')

        pause(wait)
    else
        break

    end
end
