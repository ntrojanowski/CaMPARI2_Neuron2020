%this script loops through all files in a folder to find the number of
%threshold crossings (action potentials) that occur in a voltage trace. 
%Using a user-defined threshold set on a trace-by-trace basis, it
%identifies each threshold crossing, and records the timestamp, then
%calculates the number of action potentials, as well as the CV of
%inter-event intervals and the Fano factor. 


close all;
clear;

all_files = dir;

f_idx = 0; %this is used for tracking the loop and storing the values

for f = 3:numel(all_files)
    total = numel(all_files(f).name);
    
    f_idx = f_idx + 1; %increment each time through the loop
    
    filename = str2double(all_files(f).name(total - 6:total - 3));
    
    %formatting file names to read them in
    if filename < 10
        start = strcat('000',num2str(filename));
    elseif filename >= 10 && range(1) < 100
        start = strcat('00',num2str(filename));
    elseif filename >= 100 && range(1) < 1000
        start = strcat('0',num2str(filename));
    else
        start = num2str(filename);
    end
    
    hcall_location = strcat('untitled_',start,'.h5');
    s = ws.loadDataFile(hcall_location);
    hcall_prefix = strcat('s.sweep_',start,'.analogScans');
    data = eval(hcall_prefix);
    
    samp_rate = 10000; %in Hz
    
    %initialize variables
    total_time = numel(data(:,1))/samp_rate;
    cell_1 = data(:,1);
    cell_2 = data(:,2);
    
    clear data_final
    clear times
    clear cv_events
    
    %loop through both channels (both electrodes)
    for j = 1:2
        ct = 0;
        spiketimes = [0,0];
        
        data_loop_cat = [];
        
        %loop 10000 samples at a tame (1 second)
        interval = 10000;
        for n = 1:floor(numel(data(:,1))/interval)%drop leftovers at the end (removes the last fraction of a second, out of a trace many minutes long)
            data_loop = detrend(data(((interval*(n-1))+1):(interval*n),j)); %detrend within each section, to make it easier to define an effective AP detector
            data_loop_cat = vertcat(data_loop_cat, data_loop); %concatenate all the segments after detrending
        end
        
        %plot all the segments, then select a spike threshold manually
        figure(1);
        plot(data_loop_cat)
        [X,Y] = ginput;
        
        thresh(j) = Y;%initially defined as [std(data_loop_cat) * 3.5], this proved unreliable;
        
        %loop through all the data, point by point, to find spike times
        for m = 2:numel(data_loop_cat)
            
            %first find the spike start time, add to the spike counter
            if data_loop_cat(m-1) <= thresh(j) && data_loop_cat(m) > thresh(j)
                ct = ct + 1;
                spiketimes(ct, 1) = m;
                
            %then find the spike end time
            elseif data_loop_cat(m-1) >= thresh(j) && data_loop_cat(m) < thresh(j)
                %need the index to be at least 1; traces with only one
                %spike are not reliable anyway, and are not used
                if ct < 1
                    ct = 1;
                end
                spiketimes(ct, 2) = m;
        
                %used to be some incorrect half-width calculations here
                %if the width is greater than 10 ms, discard the spike                
                if (spiketimes(ct,2) - spiketimes(ct,1)) > 15
                    ct = ct - 1;
                end
            end
        end
        
        %rename values for ease of transfer to excel
        data_final = data_loop_cat;
        numspikes(f_idx,j) = ct-1; %starting at 1
        times(1,j) = mean(spiketimes(2:end,2) - spiketimes(2:end,1)); %this calculates approximate spike widths
        dst = diff(spiketimes(:,1));
        cv_events = nanstd(dst)/nanmean(dst); %to look at the variation in the spike times
        fanofactor = (nanstd(dst))^2/nanmean(dst);
        
        numspikes(f_idx,3) = n * interval / samp_rate;
        numspikes(f_idx,j+3) = cv_events;
        numspikes(f_idx,j+5) = fanofactor;
    end
end