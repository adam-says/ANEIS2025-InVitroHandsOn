% Generate an array with all the bursts profiles (spike time histograms) from a single recording.
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



function varargout = burstProfiling(allspks, bursts, binsize, alignment, show)

    spiketrain = allspks(:,1);
    plot_extBefore = 100;     %--number of ms include before
    plot_extAfter = 150;      %--and after
    burst_max_ext = (round(max(bursts.burst_durations_ms/binsize))+ plot_extBefore/binsize + plot_extAfter/binsize)+1;   %--define maximal length of burst arrays (burst_ext included)

    N_bursts = height(bursts);
    burst_start = bursts.burst_start_ms;
    burst_end = bursts.burst_end_ms;
    burstProfiles = zeros(N_bursts,burst_max_ext);   

    burstsplot = figure("Visible","off");

    if show == 1
        set(burstsplot, 'visible', 'on')
    end

    
    for i = 1:height(bursts)

        %disp(i)

        if burst_start(i) - plot_extBefore <= 0
            ext_start = burst_start(i);
            ext_end = burst_end(i)+plot_extAfter;

        elseif burst_end(i) + plot_extAfter >= max(spiketrain)
            ext_start = burst_start(i)-plot_extBefore;
            ext_end = burst_end(i);

        else
            ext_start = burst_start(i)-plot_extBefore;
            ext_end = burst_end(i)+plot_extAfter;

        end

        burst_spikes = spiketrain(spiketrain>ext_start & spiketrain<ext_end);

        edges = ext_start:binsize:ext_end;
        burst_profile = histcounts(burst_spikes,edges);     
        
        burst_profile(end+1:burst_max_ext) = 0;

        burstProfiles(i,:) = burst_profile;

    end

    if strcmp(alignment,"start") == 1
            varargout{1} = burstProfiles;
        
    elseif strcmp(alignment,"peak") == 1
            burstProfiles = alignToPeak(burstProfiles);
            varargout{1} = burstProfiles;
    else
        disp('chose either "start" or "peak" for input variable "alignment"')
    end

    for j = 1:height(bursts)

        time_ms = (0:binsize:(size(burstProfiles,2)*binsize)-1);

        plot(time_ms, burstProfiles(j,:),'-', 'color',[0.6 0.6 0.6])
        hold on

        if show == 1
           pause
           disp('press any key to continue')
        end                
     end

     if nargout >= 2      
        burstProfiling_parameters = table(plot_extBefore,plot_extAfter,binsize);
        varargout{2} = burstProfiling_parameters;        
     end
        
     if nargout >= 3    
        mean_burst = mean(burstProfiles);
        plot(time_ms,mean_burst,'r-','linewidth',2)
        varargout{3} = mean_burst;    
     end
        
     set(burstsplot, 'visible', 'on')
     ax = gca;
     ax.TickDir = "out";
     ax.FontWeight = "bold";
     ax.FontSize = 12;
     set(gca,'TickLabelInterpreter', 'none');
     ylabel('# Spikes','FontSize',14)
     xlabel('Time (ms)','FontSize',14)

end

