% this script is used to select folders to search for data, then loop
% though the data files and allow the user to visually identify stable
% recording regions for analysis. Typically used for voltage clamp
% recordings. 

% this script saves the x coordinates and file names for each trace
%% to be changed for each new run
excised_filename = {'excised_EI_200320_all.mat'};
%%%%%% where to save excised x coordinates

experiment = {'200219', '200221', '200222', '200224', '200225', '200311', '200312', '200313'};
%%%%%% folders containing the data to analyze

fp_minianalysis = 'C:/Users/ntrojanowski/Documents/MATLAB/minianalysis';
%%%%%% Filepath where you keep this script and accompanying SAVED files/scripts

fp_data = 'C:/Users/ntrojanowski/Documents/DATA/';
%%%%%% Filepath where you keep data folders

%initialize
excised_points = cell(1,numel(experiment));
excised_wavenames = cell(1,numel(experiment));

%loop through files
for jj = 1:numel(experiment)
    current_experiment = experiment{jj};
    cd(strcat(fp_data,experiment{jj}))
    
    all_files = dir;
    
    excised_points{1,jj} = cell(1,numel(all_files)-2);
    excised_wavenames{1,jj} = cell(1,numel(all_files)-2);
    
    for kk = 3:numel(all_files)
        %%%% avoid FI test traces
        if strfind(all_files(kk).name, '-')
            continue
        else
            clear resampled; 
            
            %this method of reading in files was never corrected to the wavesurfer-provided 
            %function for reading in .h5 files, but since it only uses the x (time) values, this is
            %not an issue,
            total = numel(all_files(kk).name);
            hcall_location = strcat(fp_data,current_experiment,'/',all_files(kk).name);
            hcall_prefix=strcat('/sweep_',all_files(kk).name(total-6:total-3),'/analogScans');
            vals = double(h5read(hcall_location,hcall_prefix));
            
            %It is important to keep track of whether
            %or not resampled values are being used, as this will affect
            %the x coordinates that are saved.
            rs = resample(vals,1,2);
            resampled = rs(:,1);

            
            
            %plot the resampled values to allow user to select the regions
            %of good recordings. 
            figure('units','normalized','outerposition',[0 0 1 1])
            if median(resampled) < mean(resampled)
                ylim([min(resampled), min(resampled) + 1000]);
            else
                ylim([max(resampled) - 1000, max(resampled)]);
            end
            title(hcall_prefix);
            ylabel('current (pA)');
            xlabel('time (s)');            
            hold on;
            plot(resampled, '-k'); % using non-normalized to better remove drifting baselines

            
            [X,Y] = ginput;
            excised_points{1,jj}{kk-2} = X;
            excised_wavenames{1,jj}{kk-2} = all_files(kk).name;
            close all
            fclose all;
        end
    end
end

%%

%%%%%%%%%%%%% SAVE  (careful with overwriting!)
cd(fp_minianalysis)
result = input('Save results? This will overwrite previous file unless renamed! (1 = y, 2 = n):','s');

if str2double(result) == 1
    save(excised_filename,'excised_points','excised_wavenames')
end