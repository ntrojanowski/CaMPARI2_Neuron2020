%this script reads in the previously generated event amplitudes and timings, and a list of
%traces with stable regions for analysis, then uses the find_iso_amps
%function to randomly select the same amount of measurements from each
%cell and compare the distribution of the measurements between high and
%low activity neurons


clear;
load('MINIANALYSIS_EI_200320_initial_2.mat') %provide the amplitude data for all events
load('green_red_EI_indices.mat') %provide the cells to analyze


%%
close all;

%run 1000 iterations of selecting the same number of events from each
%trace, and combine the amplitude and inter-event intervals of these events
%into single arrays
for z = 1:1000
    %this selects the points and calculates the relevant parameters
    [green_randsamp_e, green_first_events_e, green_rand_ISI_e, green_first_ISI_e] = find_iso_amps(AMP_ALL, TIME_INDICES, greentraces, 1);    
    g_ra_e(:,z) = green_randsamp_e(:);
    g_rI_e(:,z) = green_rand_ISI_e(:);
    
    [red_randsamp_e, red_first_events_e, red_rand_ISI_e, red_first_ISI_e] = find_iso_amps(AMP_ALL, TIME_INDICES, redtraces, 1);    
    r_ra_e(:,z) = red_randsamp_e(:);
    r_rI_e(:,z) = red_rand_ISI_e(:);        
    
    [green_randsamp_i, green_first_events_i, green_rand_ISI_i, green_first_ISI_i] = find_iso_amps(AMP_ALL, TIME_INDICES, greentraces, 2);    
    g_ra_i(:,z) = green_randsamp_i(:);
    g_rI_i(:,z) = green_rand_ISI_i(:);
    
    [red_randsamp_i, red_first_events_i, red_rand_ISI_i, red_first_ISI_i] = find_iso_amps(AMP_ALL, TIME_INDICES, redtraces, 2);   
    r_ra_i(:,z) = red_randsamp_i(:);
    r_rI_i(:,z) = red_rand_ISI_i(:);
end

%this performs a KS test on the selected values to determine if they differ
%between hign and low activity neurons
[a1, b1] = kstest2(mean(sort(g_ra_e), 2), mean(sort(r_ra_e), 2));
[c1, d1] = kstest2(mean(sort(g_rI_e), 2), mean(sort(r_rI_e), 2));
[e1, f1] = kstest2(mean(sort(g_ra_i), 2), mean(sort(r_ra_i), 2));
[g1, h1] = kstest2(mean(sort(g_rI_i), 2), mean(sort(r_rI_i), 2));

%plot the cumulative distributions for each measurement
figure(1); hold on; cdfplot(mean(sort(g_ra_e), 2)), cdfplot(mean(sort(r_ra_e), 2));
figure(3); hold on; cdfplot(mean(sort(g_rI_e), 2)/10), cdfplot(mean(sort(r_rI_e), 2)/10);
figure(5); hold on; cdfplot(mean(sort(g_ra_i), 2)), cdfplot(mean(sort(r_ra_i), 2));
figure(7); hold on; cdfplot(mean(sort(g_rI_i), 2)/10), cdfplot(mean(sort(r_rI_i), 2)/10);
