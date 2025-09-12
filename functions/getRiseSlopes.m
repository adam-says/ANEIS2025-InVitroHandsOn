% Generate a matrix with the rising slopes of the input bursts
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

function RiseSlopes = getRiseSlopes(smoothedbursts,binsize,normalize)

    RiseSlopes = cell(size(smoothedbursts,1),1);
    
    for b = 1:size(smoothedbursts,1)
    
        burst_profile = smoothedbursts(b,:);
        [peak_ampl, peak_idx] = max(burst_profile);
        
        if strcmp(normalize,'yes')
            burst_profile = burst_profile/peak_ampl;
        end

        Burst_rise = burst_profile(1:peak_idx);
    
        RiseSlope_dt = diff(Burst_rise)/binsize;

        if isempty(RiseSlope_dt) == 1
           RiseSlopes(b) = nan;
           continue
        end

        RiseSlopes{b} = RiseSlope_dt; 

    end


end