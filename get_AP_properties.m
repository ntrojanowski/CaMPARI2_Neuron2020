function [amplitude, AHP, rise_slope, decay_slope, half_width, ap_peak_up, ap_peak_down, up_down_ratio] = get_AP_properties(data, dV_thresh)
%This function takes a small section of a current clamp recording and dV threshold, then uses
%that threshold to identify an action potential. It then returns various
%properties of the action potential (spike and afterhyperpolarization amplitude, 
%half-width, mean and peak upslope and downslope values, and ratio of up-to-down peak slopes).

%find all the points where dV is greater than a threshold
all_dV = find(diff(data) > dV_thresh);

%initialize, so that if a parameter measurement fails, the function still
%runs
amplitude = 0;
AHP = 0;
rise_slope = 0;
decay_slope = 0;
half_width = 0;
ap_peak_up = 0;
ap_peak_down = 0;
up_down_ratio = 0;

%if the dV never reaches threshold, plot and return
if isempty(all_dV)
    figure; plot(data);
    return
end

i_idx = 2;

spiketimes_t(1) = all_dV(1); %if there's only one spike, just take the first point where it crosses threshold
if max(diff(all_dV)) > 50 %look for spikes that are at least 5 ms apart
    for i = 2:numel(all_dV) %
        if all_dV(i) - all_dV(i-1) < 50
            continue
        else %if threshold crossings are more than 5 ms apart, record the time of the crossing
            spiketimes_t(i_idx) = all_dV(i); 
            i_idx = i_idx + 1;
        end
    end
end

%look through the spiketimes to find the spike to use for analysis
for k = 1:numel(spiketimes_t)
    if k == numel(spiketimes_t) %if it's the first one, use the whole trace for analysis
        test_spike = data(spiketimes_t(k):10000); %if last spike, end at the end of the pulse
    else %if there are multiple, only use the time between the first pulse and the second one
        test_spike = data(spiketimes_t(k):spiketimes_t(k+1));
        
        if (10000 - spiketimes_t(k+1)) < 100 %if it's too close to the end of the trace, don't use it. This happens (very rarely) if the spike occurs very close to the end of the current injection
            return
        end
    end
    
    %caluculate the amplitude of the spike and AHP
    amplitude(k) = max(test_spike) - test_spike(1);
    AHP(k) = test_spike(1) - min(test_spike);
    
    %find the time of the max and min points and calculate the rise and
    %decay slopes
    maxpoint = find(test_spike == max(test_spike));
    rise_slope(k) = amplitude(k)/maxpoint(1);
    
    minpoint = find(test_spike == min(test_spike));
    [~, startpoint_decay] = min(abs(test_spike(maxpoint:minpoint) - test_spike(1)));
    decay_slope(k) = amplitude(k)/startpoint_decay(1); 
    
    
    %find the time at which the amplitude is half of the maximum in the
    %rising and falling phases. The difference between these times is the half-width 
    %These times are calculated by finding the point at which the voltage
    %is closest to the calculated voltage halfway between the min and the
    %max
    [half_max_up, half_max_up_index] = min(abs(test_spike(1:maxpoint) - (test_spike(1) + amplitude(k)/2)));
    [~, half_max_down_index] = min(abs(test_spike(maxpoint:minpoint) - half_max_up));
    
    half_width(k) = half_max_down_index(1) - half_max_up_index(1);
    
    %find the maximum slopes in the rising and falling phases, and the
    %ratio of these values
    ap_peak_up(k) = max(diff(test_spike(1:maxpoint)));
    ap_peak_down(k) = - min(diff(test_spike(maxpoint:minpoint)));
    up_down_ratio(k) = ap_peak_up(k) / ap_peak_down(k);    
end

