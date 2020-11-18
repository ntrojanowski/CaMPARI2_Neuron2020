%this script is used for measuring the strengths of connections between simulatneously
%recorded neurons. It runs on a folder of .h5 files. It is best run on
%voltage clamp data - it produces nonsensical results when run on current
%clamp traces. It finds the 3 point average current around the peak, and compares
%this to the baseline current value used to hold cells near -80. It also
%calculates the latency from the current injection in one cell to the
%reponse in another. 

close all
clear
% this script is written to be run from a folder of .h5 files

pulsetimes = [5000, 5500, 6000, 6500, 7000]; %these are the times at which pulses are delivered

%loop through all the files in the folder
all_files = dir;
f_idx = 0;
for f = 3:numel(all_files)
    total = numel(all_files(f).name);
    
    %for each file, first find the names of the traces, then convert the
    %name strings to manipulable numbers. See FI_curve_and_PP for file
    %input details. 
    if strfind(all_files(f).name, '-')
        f_idx = f_idx + 1;
        
        range(1) = str2double(all_files(f).name(total - 11:total - 8));
        range(2) = str2double(all_files(f).name(total - 6:total - 3));
        
        %formatting file names
        if range(1) < 10
            start = strcat('000',num2str(range(1)));
        elseif range(1) >= 10 && range(1) < 100
            start = strcat('00',num2str(range(1)));
        elseif range(1) >= 100 && range(1) < 1000
            start = strcat('0',num2str(range(1)));
        else
            start = num2str(range(1));
        end
        
        if range(2) < 10
            finish = strcat('000',num2str(range(2)));
        elseif range(2) >= 10 && range(2) < 100
            finish = strcat('00',num2str(range(2)));
        elseif range(2) >= 100 && range(2) < 1000
            finish = strcat('0',num2str(range(2)));
        else
            finish = num2str(range(2));
        end
        
        k_idx = 0;
        for k = range(1):range(2)
            k_idx = k_idx + 1;
            if k < 10
                sname = strcat('000',num2str(k));
            elseif k >= 10 && k < 100
                sname = strcat('00',num2str(k));
            elseif k >= 100 && k < 1000
                sname = strcat('0',num2str(k));
            else
                sname = num2str(k);
            end
            
            %use the file and trace names to read in the files
            hcall_location = strcat('untitled_',start,'-',finish,'.h5');
            s = ws.loadDataFile(hcall_location);
            hcall_prefix = strcat('s.sweep_',sname,'.analogScans');
            data = eval(hcall_prefix);
            
          
            %collect the raw data from each .h5 file into a one matrix
            for d=1:4
                raw_data(:,k_idx,d) = data(:,d);
            end
        end
        
        %average each sweep across all trials
        avg_trace(f,:,:) = squeeze(mean(raw_data,2));
        
        %calculate the mean amplitudes of each pulse for each cell
        for d = 1:4
            
            %loop through each pulse
            for p = 1:5
                
                %smooth the traces, find the baseline and trough of the
                %region surrounding each pulse
                pulsechunk = smooth(avg_trace(f,(pulsetimes(p)-100):(pulsetimes(p)+400),d)); %the chunk starts 100 points before the stimulus
                baseline(p,d) = mean(pulsechunk(1:75)); %baseline preceding pulse
                
                %avoid stimulus artifact and channel crosstalk by starting 0.8 ms after stimulus
                offset = 7;                
                [minnum, minval] = min(pulsechunk(1+100+offset):(1+200+offset)); 

                minadj = minval + 100 + offset;
                minlat(p,d) = minadj;
                
                %find the 3 point average surrounding the peak (and deal
                %with cases if the peak is near the beginning or end of the trace)
                if minadj == 1
                    minmean(p,d) = mean(pulsechunk(minval:minval+2));
                elseif minadj == 139
                    minmean(p,d) = mean(pulsechunk(minval-2:minval));
                else
                    minmean(p,d) = mean(pulsechunk(minadj-1:minadj+1));
                end
            end
        end
        
        %calculate the amplitude and latency
        amplitude(:,:,f-2) = baseline - minmean;
        latency(:,:,f-2) = minlat; %the latency is the time of the trough

        %plot the average traces for each cell
        figure(f_idx);
        for plot_ID = 1:4            
            subplot(4,1,plot_ID);  
            plot(smooth(avg_trace(f,4000:8000,plot_ID)));
        end
    end
end