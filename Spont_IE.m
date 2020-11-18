%this script runs on a folder of .h5 files recorded in voltage clamps. It
%first takes a histogram of the traces, then uses that histogram to measure
%the total charge transfer during the trace. Briefly, it assumes that the 
%noise is symmetrical around the peak of the histogram and that the 
%recordings were taken at the reversal potential for either excitation or 
%inhibition, and therefore subtracts the values from the side of the 
%histogram with the shorter tail from the side with the longer tail. The
%sum of the remaining values represents the total charge transfer. 

close all;
clear;

%initialize variables
total_time = zeros(numel(all_files) - 2, 1);
weighted_current = zeros(numel(all_files) - 2, 1);

%loop through all the files. See FI_curve_and_PP for filehandling details
all_files = dir;
for f = 3:numel(all_files)
    total = numel(all_files(f).name);
    fnum = str2double(all_files(f).name(total-6:total-3));
    if fnum < 10
        sname = strcat('000',num2str(fnum));
    elseif fnum >= 10 && fnum < 100
        sname = strcat('00',num2str(fnum));
    elseif fnum >= 100 && fnum < 1000
        sname = strcat('0',num2str(fnum));
    else
        sname = num2str(fnum);
    end
    
    hcall_location = strcat('untitled_',sname,'.h5');
    s = ws.loadDataFile(hcall_location);
    hcall_prefix = strcat('s.sweep_',sname,'.analogScans');
    full_data = eval(hcall_prefix);
    
    %discard regions at the beginning that tend to have substantial drift
    %after shifting holding potentials
    data = full_data(201001:end-10000); 
    time = (1:numel(data))/10000;
    
    %pick ROIs that don't have dramatic drift, save time points
    figure('units','normalized','outerposition',[0 0 1 1])
    if median(data) < mean(data)
        ylim([min(data), min(data) + 1000]);
    else
        ylim([max(data) - 1000, max(data)]);
    end
    title('pick regions to analyze');
    ylabel('current (pA)');
    xlabel('time (s)');
        
    hold on;
    plot(time(1:100000), data(320001:420000), '-k'); % using non-normalized to better remove drifting baselines
    plot([1,time(end)], [median(data),  median(data)], '-r'); %plot a line for visualization of drift
    
    [X,Y] = ginput;
    %X = [40000,  1000000];
    excised_points = X;
    
    
    %use excised timestamps to find the data to analyze
    yy = [];
    if mod(numel(excised_points), 2) ~= 0
        disp(strcat('Does not contain an appropriate number of excision points, keeping all'))
    else
        current = zeros(numel(excised_points)/2, 1);
        interval_regions = zeros(numel(excised_points)/2, 1);
        weighted_charge = zeros(numel(excised_points)/2, 1);
        smooth_ct = [];
        
        for m = 1:2:numel(excised_points-1)
            
            %work through the data one ROI at a time
            if excised_points(m+1) > numel(data)
                excised_points(m+1) = numel(data);
            end
            
            %detrend, zero  each chunk separately, then weighted-add them
            %based on their lengths
            warning('off','all')  %%%%% for some reason, some selected excision points are not integers...
            chunk = (data(excised_points(m):excised_points(m+1)));

            yy(numel(yy)+1:numel(yy)+numel(data(excised_points(m):excised_points(m+1)))) = chunk;
            warning('on','all')
            
            %since we are looking at currents only in one direction, all
            %values towards the opposite side of the center of the
            %histogram are noise. Since the noise should be symmetrical,
            %we can subtract the equivalent value from the side of the
            %currents we are measuring to remove the noise, then the
            %remainder of the area under that curve will be the signal we
            %want
            
            %find the max of the smoothed histogram, then fit a polynomial
            %to it so the peak detection is more accurate
            [ct, bin] = hist(chunk,numel(chunk)/500); %500 - 2000 gives the cleanest curve
            smooth_ct = smooth(ct, 15, 'lowess');
            [pt, idx] = max(smooth_ct);
            max_zero = bin(idx);
              
            weighted_charge((m+1)/2) = sum(ct .* (bin-max_zero));
            interval_regions((m+1)/2) = numel(chunk);        
        end        
    end
    
    %weight each contributing region to each cell's measurement by dividing by the length of that
    %region
    total_time(f-2) = sum(interval_regions)/10000;
    weighted_current(f-2) =  sum(weighted_charge);
end


output = horzcat(weighted_current, total_time);