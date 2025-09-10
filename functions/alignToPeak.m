% Function to align all bursts from a single recording to their peak.
% Function with 2 arguments:
% 1) burst_profiles -> the spike times histograms of all the burst detected
% in that recording and stored row-wise in an array
% 2) show -> if show == 0 don't show figures, if show == 1 show figure and
% pause
% Can be called with one output argument:
% alignedBursts = burstalign(burst_profiles, show)
% In this case the output will be the array with all the bursts profiles
% aligned to their peak and the relative plots
% Or two:
% [alignedBursts, mean_burst] = burstalign(burst_profiles, show)
% In this case the output will be the array with all the bursts profiles
% aligned to their peak and an array with the average burst calculated as point-bypoint average
% of all the bursts aligned at their peak. Their relative plot will show
% the average one in red


function alignedBursts = alignToPeak(burst_profiles)

    [burst_peaks, peaks_idx] = max(burst_profiles,[],2);
    
    [~, biggerBurst] = max(burst_peaks);
    biggerPeak_idx = peaks_idx(biggerBurst);
    
    allshifts = biggerPeak_idx-peaks_idx;
    max_shift = max(abs(allshifts));
    allburststoalign = padarray(burst_profiles,[0 max_shift],'both');
    
    alignedBursts = zeros(size(allburststoalign));    
        
    for h = 1:size(burst_profiles,1)
            
        shift = allshifts(h);
        burst = allburststoalign(h,:);
        alignedburst = zeros(1,length(burst));
        
        for z = 1:length(alignedburst)
            
           if z<=shift
               alignedburst(1,z) = 0;
           elseif z-shift>length(alignedburst)
               alignedburst(1,z) = 0;
           else
              alignedburst(1,z) = burst(z-shift);
           end
            
        end
                     
        alignedBursts(h,:) = alignedburst;
                      
    end
     
end
