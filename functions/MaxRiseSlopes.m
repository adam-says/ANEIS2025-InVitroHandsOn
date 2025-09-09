% Generate a vector with the max rising slopes of the input bursts
% INPUT:
% 1) smoothedBursts -> a cell array of burst profiles smoothed to avoid
% noisy oscillations
% 2) binsize -> size of the temporal bins used to generate the burst
% profiles
% 3) normalize -> 'yes' or 'no' depending if you want to calculate the
% rising slopes on burst normalized by their amplitude or not
%
% OUTPUT:
% a vector with the max rising slopes of each burst (relative or absolute)

function MaxRiseSlopes = MaxRiseSlopes(smoothedBursts,binsize,normalize)

    MaxRiseSlopes = zeros(size(smoothedBursts,1),1);
    
    for b = 1:size(smoothedBursts,1)
    
        burst_profile = smoothedBursts(b,:);
        [peak_ampl, peak_idx] = max(burst_profile);
        
        if strcmp(normalize,'yes')
            burst_profile = burst_profile/peak_ampl;
        end
        Burst_rise = burst_profile(1:peak_idx);
    
        RiseSlope_dt = diff(Burst_rise)/binsize^2;
        maxSlope = max(RiseSlope_dt);

        if isempty(RiseSlope_dt) == 1
           MaxRiseSlopes(b) = nan;
           continue
        end

        MaxRiseSlopes(b) = maxSlope; 

    end


end