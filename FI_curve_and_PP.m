%this script analyzes trains of action potentials resulting from current
%injections (to produce FI curves), obtaining AP measurements from the
%get_AP_properties function. It also measures passive properties. this 
%script runs on a folder of .h5 files. This script relies on three custom
%functions: get_spike_threshold, get_Ra_doubleexp_h5, and get_AP_properties

close all
clear

%if file name contains a dash, it's for FI curves because it has multiple
%sweeps. if not, it should have a seal test pulse in it.
all_files = dir;
num_FI = 0;

%all file indexes start at 3 due to hidden files
%this loop counts the number of files containing FI data
for f = 3:numel(all_files)
    if strfind(all_files(f).name, '-')
        num_FI = num_FI + 1;
    end
end

num_notFI = numel(all_files) - num_FI - 2; %this is the number of non-FI traces
%initialize a bunch of variables

Rs_calc = zeros(num_notFI, 1);
Rt = zeros(num_notFI, 1);
Cm = zeros(num_notFI, 1);
Vm = zeros(num_notFI, 1);
tauF = zeros(num_notFI, 1);
tauS = zeros(num_notFI, 1);

%initialize indices to be used in for loops
seal_test_index = 0;
fi_index = -1;%unorthodox, but it works
fr_index = -1;

for f = 3:numel(all_files)
    %first, read in the file names
    %for each file, the last 4 characters of the file name (before the
    %suffix) contain the number of the file. For example,
    %untitled_0001-0026.h5
    
    
    total = numel(all_files(f).name);
    
    if strfind(all_files(f).name, '-') %this is where the FI curves are calculated
        
        %find the numbers of the first and last traces
        range(1) = str2double(all_files(f).name(total - 11:total - 8));
        range(2) = str2double(all_files(f).name(total - 6:total - 3));
        
        %format the file names to strings file names
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
        
        
        fi_index = fi_index + 2; %to return cells on consecutive rows
        
        for d=1%:2 %change to 2 if recording from two cells simultaneously
            %return in sequential rows to avoid 3rd dimension. If only
            %recording from one cell, this makes alternating blank ros
            %that can later be deleted
            if d==1
                return_idx = fi_index;
            else
                return_idx = fi_index + 1;
            end
            
            
            %initialize a bunch of return variables
            lat_rheo(return_idx) = 0;
            lat_5(return_idx) = 0;
            lat_10(return_idx) = 0;
            lat_15(return_idx) = 0;
            adaptation_5(return_idx) = 0;
            adaptation_10(return_idx) = 0;
            adaptation_15(return_idx) = 0;
            adaptation_IFR_15_index_F2L(return_idx) = 0;
            adaptation_IFR_30_index_F2L(return_idx) = 0;
            adaptation_IFR_45_index_F2L(return_idx) = 0;
            adapt_ratio_cont = 0;
            vT(return_idx) = 0;
            
            %this opens and loops through the individual sweeps within
            %the .h5 files in the same way as above, by converting the
            %file names to numerical values
            for k = range(1):range(2)
                if k < 10
                    sname = strcat('000',num2str(k));
                elseif k >= 10 && k < 100
                    sname = strcat('00',num2str(k));
                elseif k >= 100 && k < 1000
                    sname = strcat('0',num2str(k));
                else
                    sname = num2str(k);
                end
                
                hcall_location = strcat('untitled_',start,'-',finish,'.h5');
                s = ws.loadDataFile(hcall_location);
                hcall_prefix = strcat('s.sweep_',sname,'.analogScans');
                data = eval(hcall_prefix);
                
                trial_idx = k - range(1) + 1;
                threshold = 10; %crossing this Vm counts as a spike for counting purposes
                dV_thresh = 20; %10 works better for finding picturesque AP waveforms
                %figure; hold on; plot(data(:,1));
                
                %identify spike times and count spikes by counting threshold crossings
                counter = 0;
                clear spiketimes;
                clear individual_adaptation
                %end at 15000 to calculate for half a second
                for j=10000:15000
                    if data(j,d) <= threshold && data(j + 1,d) > threshold
                        counter = counter + 1;
                        spiketimes(counter) = j;
                    end
                end
                
                MFR(return_idx,trial_idx) = counter * 2; %mean FR in Hz
                
                
                %then, if spikes are detected, calculate a whole bunch of parameters for the AP
                %trains
                
                if counter > 0
                    lat(return_idx,trial_idx) = (spiketimes(1) - 10000) / 10; %latency to first spike in ms
                    
                    %first, measure parameters at rheobase. These two
                    %loops find the first time in a series of traces a
                    %spike has occured
                    if counter >= 1
                        if lat_rheo(return_idx) == 0 %calculations at rheobase
                            lat_rheo(return_idx) = (spiketimes(1) - 10000) / 10; %latency at rheobase in ms
                            r_vT_ramp = get_spike_threshold(data(:,d), dV_thresh); %use get_spike_threshold to find the AP threshold
                            vT(return_idx) = r_vT_ramp(1);
                            %get AP parameters
                            [amplitude, AHP, rise_slope, decay_slope, half_width, ap_peak_up, ap_peak_down, up_down_ratio] = get_AP_properties(data(10000:20000,d), dV_thresh);
                            
                            %transfer parameters to return arrays
                            amp_first_rheo(return_idx) = amplitude(1);
                            AHP_first_rheo(return_idx) = AHP(1);
                            rise_slope_first_rheo(return_idx) = rise_slope(1);
                            decay_slope_first_rheo(return_idx) = decay_slope(1);
                            half_width_first_rheo(return_idx) = half_width(1);
                            ap_peak_up_first_rheo(return_idx) = ap_peak_up(1);
                            ap_peak_down_first_rheo(return_idx) = ap_peak_down(1);
                            up_down_ratio_first_rheo(return_idx) = up_down_ratio(1);
                        end
                        
                        d_ST = diff(spiketimes); %interspike intervals
                        
                        %if there are multiple spikes, calculate the
                        %initial and average instantaneous firing rates
                        if counter > 1
                            IFR(return_idx,trial_idx) = 1./(spiketimes(2) - spiketimes(1))*10000; %instantaneous FR in Hz
                            IFR(isnan(IFR)) = 0; %convert NaNs to zeros
                            IFR_avg(return_idx,trial_idx) = mean((1./diff(spiketimes)))*10000;
                            IFR_avg(isnan(IFR_avg)) = 0;
                            IFR_CV(return_idx,trial_idx) = var(diff(spiketimes),1);
                        end
                        
                        %if there are more than 2 spikes, look at
                        %adaptation
                        if counter > 2
                            for q = 2:numel(d_ST)
                                adaptation_index(q-1) = (d_ST(q) - d_ST(q-1))/(d_ST(q) + d_ST(q-1)); %ISI diff over sum
                            end
                            
                            summed_adaptation_index(return_idx, trial_idx) = sum(adaptation_index) / (numel(adaptation_index)-1);
                            summed_adaptation_index(~isfinite(summed_adaptation_index)) = 0;
                            summed_adaptation_index(summed_adaptation_index == 0) = NaN;
                            
                            %calculate the adaptation between the first
                            %and last spike pairs. This was never used
                            %for analysis
                            adaptation_F2L(return_idx, trial_idx) = (d_ST(end) - d_ST(1)) / (d_ST(end) + d_ST(1));
                            
                            %calculate adaptation at specific ISIs or
                            %spike numbers
                            
                            %first calculate the adaptation index for
                            %the whole trace sequentially, then the
                            %adaptation between the first and last
                            %spike pairs
                            IFR_idx_15 = [];
                            if IFR(return_idx,trial_idx) >= 15 && adaptation_IFR_15_index_F2L(return_idx) == 0
                                for  q = 1:numel(d_ST)-1
                                    IFR_idx_15(q) = ((d_ST(q+1) - d_ST(q))/(d_ST(q+1)+d_ST(q)));
                                end
                                IFR_idx_15_sum(return_idx) = sum(IFR_idx_15)/(numel(IFR_idx_15)-1);
                                adaptation_IFR_15(return_idx) = IFR_idx_15(1)/IFR_idx_15(end);
                                adaptation_IFR_15_index_F2L(return_idx) = (d_ST(end) - d_ST(1))/(d_ST(end) + d_ST(1));
                            end
                            
                            IFR_idx_30 = [];
                            if IFR(return_idx,trial_idx) >= 30 && adaptation_IFR_30_index_F2L(return_idx) == 0
                                for  q = 1:numel(d_ST)-1
                                    IFR_idx_30(q) = ((d_ST(q+1) - d_ST(q))/(d_ST(q+1)+d_ST(q)));
                                end
                                IFR_idx_30_sum(return_idx) = sum(IFR_idx_30)/(numel(IFR_idx_30)-1);
                                adaptation_IFR_30(return_idx) = IFR_idx_30(1)/IFR_idx_30(end);
                                adaptation_IFR_30_index_F2L(return_idx) = (d_ST(end) - d_ST(1))/(d_ST(end) + d_ST(1));
                            end
                            
                            IFR_idx_45 = [];
                            if IFR(return_idx,trial_idx) >= 45 && adaptation_IFR_45_index_F2L(return_idx) == 0
                                for  q = 1:numel(d_ST)-1
                                    IFR_idx_45(q) = ((d_ST(q+1) - d_ST(q))/(d_ST(q+1)+d_ST(q)));
                                end
                                IFR_idx_45_sum(return_idx) = sum(IFR_idx_45)/(numel(IFR_idx_45)-1);
                                adaptation_IFR_45(return_idx) = IFR_idx_45(1)/IFR_idx_45(end);
                                adaptation_IFR_45_index_F2L(return_idx) = (d_ST(end) - d_ST(1))/(d_ST(end) + d_ST(1));
                                
                                %sometimes depolarization block made ISI measurements inaccurate, so we
                                %ignored traces in which a high IFR (short initial ISI) did not produce more than 5 spikes
                                d_ST_45 = d_ST;
                                
                                if numel(d_ST_45) > 5
                                    IFR_45(return_idx, 1:6) = d_ST_45(1:6);
                                else
                                    IFR_45(return_idx, 1:6) = 0;
                                end
                            end
                            
                        end
                        
                        %we also wanted to measure adaptation at points where
                        %neurons were firing the same number of APs.
                        
                        if counter >= 5 && adaptation_5(return_idx) == 0
                            ISI_5(return_idx, 1:4) = d_ST(1:4);
                            adaptation_5(return_idx) = d_ST(1)/d_ST(4);
                            %adapt_index_5(return_idx, 1:4) = adaptation_index(1:4);
                            %adapt_index_5_F2L(return_idx) = adaptation_index(1)/adaptation_index(4);
                            lat_5(return_idx) = (spiketimes(1) - 10000) / 10; %in ms
                            
                            %We also wanted to measure AP properties at a
                            %point where the neurons were firing spike
                            %trains, since our rheobase measurements where
                            %somewhat coarse. Sometimes spikes
                            %cross threshold, but not rapidly enough so
                            %they are not detected with the more
                            %sophisticated get_AP_properties function
                            [amplitude, AHP, rise_slope, decay_slope, half_width, ap_peak_up, ap_peak_down, up_down_ratio] = get_AP_properties(data(10000:20000,d), dV_thresh);
                            if numel(amplitude) < 5
                                amp_fifth(return_idx) = NaN;
                                AHP_fifth(return_idx) = NaN;
                                rise_slope_fifth(return_idx) = NaN;
                                decay_slope_fifth(return_idx) = NaN;
                                half_width_fifth(return_idx) = NaN;
                                ap_peak_up_fifth(return_idx) = NaN;
                                ap_peak_down_fifth(return_idx) = NaN;
                                up_down_ratio_fifth(return_idx) = NaN;
                            else
                                amp_fifth(return_idx) = amplitude(5);
                                AHP_fifth(return_idx) = AHP(5);
                                rise_slope_fifth(return_idx) = rise_slope(5);
                                decay_slope_fifth(return_idx) = decay_slope(5);
                                half_width_fifth(return_idx) = half_width(5);
                                ap_peak_up_fifth(return_idx) = ap_peak_up(5);
                                ap_peak_down_fifth(return_idx) = ap_peak_down(5);
                                up_down_ratio_fifth(return_idx) = up_down_ratio(5);
                            end
                        end
                        
                        %measure adaptation at 10 and 15 spikes per
                        %pulse
                        if counter >= 10 && adaptation_10(return_idx) == 0
                            ISI_10(return_idx, 1:9) = d_ST(1:9);
                            adaptation_10(return_idx) = d_ST(1)/d_ST(9);
                            %adapt_index_10(return_idx, 1:9) = adaptation_index(1:9);
                            %adapt_index_10_F2L(return_idx) = adaptation_index(1)/adaptation_index(9);
                            lat_10(return_idx) = (spiketimes(1) - 10000) / 10; %in ms
                        end
                        
                        if counter >= 15 && adaptation_15(return_idx) == 0
                            ISI_15(return_idx, 1:14) = d_ST(1:14);
                            adaptation_15(return_idx) = d_ST(1)/d_ST(14);
                            %adapt_index_15(return_idx, 1:14) = adaptation_index(1:14);
                            %adapt_index_15_F2L(return_idx) = adaptation_index(1)/adaptation_index(14);
                            lat_15(return_idx) = (spiketimes(1) - 10000) / 10; %in ms
                        end
                    end
                end
                
                %calculate sag magnitude and fraction
                if k == range(1)
                    minpoint = min(data(:,d)) - mean(data(8000:10000,d));
                    steadystate = mean(data(18000:20000,d)) - mean(data(8000:10000,d));
                    sag_fraction(return_idx) = steadystate / minpoint;
                end
            end    
        end
        
    else %this is where traces that have just a seal test are used to calculate passive properties
        %first, find the files that only contain one sweep - these will
        %have seal tests at the beginning
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
        data = eval(hcall_prefix);
        
        %initialize values for measuring passive properties
        step_start = 0.5; %in s
        vstep = -0.005; %in V
        samprate = 10000;
        figures_on = 0; %1 for on, 0 for off
        seal_test_index = seal_test_index + 1;
        
        for g=1%:2 change to 2 if using a second electrode
            if std(data(:,g)) > 5 && data(1,g) > -5000%prevents running analysis on channels with good recording
                
                %measure passive properties
                [Rs_calc_scaled, Rt_scaled, Cm_scaled, Vm_scaled, tau_fast, tau_slow] = get_Ra_doubleexp_h5(data(:,g), step_start, vstep, samprate, figures_on);
                
                %save properties
                
                Rs_calc(seal_test_index,g) = Rs_calc_scaled;
                Rt(seal_test_index,g) = Rt_scaled;
                Cm(seal_test_index,g) = Cm_scaled;
                Vm(seal_test_index,g) = Vm_scaled;
                tauF(seal_test_index,g) = tau_fast*10000;
                tauS(seal_test_index,g) = tau_slow*10000; %this is never used
            end
        end
    end
end