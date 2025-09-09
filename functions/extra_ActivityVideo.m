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
