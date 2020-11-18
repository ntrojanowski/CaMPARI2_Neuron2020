function [rand_events, first_events, rand_ISIs, first_ISIs] = find_iso_amps(AMP_ALL, TIME_INDICES, traces, flag)
%this function, which is called in the Cum_amps scrpit, takes arrays containing event amplitudes and times and an
%array indiciating which traces to analyze, and returns the amplitude and
%inter-event intervals of a user-provided number of randomly selected events

%if flag is 2, add one to the identified excitatory trace index to look at inhibition
if flag == 2
    traces(:,2) = traces(:,2)+1;
end

%first loop through the traces to identify the minimum number of events per trace
for d = 1:length(traces)
    row_trace = AMP_ALL{1,traces(d,1)};
    amps_trace = row_trace{1,traces(d,2)};
    num_events(d) = length(amps_trace(~isnan(amps_trace)));
end
%%

%this should find the total number of good events in each cell by combining the 
%numbers from different traces, but requires yet-unwritten code for
%concatenating traces in order for its output to be useful.  

events_by_cell(1) = num_events(1);
for m = 2:length(num_events)
    if traces(m,3) == traces(m-1, 3)
        events_by_cell(traces(m,3)) = events_by_cell(traces(m,3)) + num_events(m);
    else
        events_by_cell(traces(m,3)) = num_events(m);
    end
end

%select the number of events to select per cell, either by finding the
%fewest points in one trace manually or automatically, 
%or by providing a different user-definied number
downsample_val = 120; %min(events_by_cell);
%%

%initialize variables for speed
row_test = AMP_ALL{1,traces(1,1)};
amps_test = row_test{1,traces(1,2)};
per_cell = amps_test(~isnan(amps_test));

row_times = TIME_INDICES{1,traces(1,1)};
times_test = row_times{1,traces(1,2)};
ISIs = diff(times_test(~isnan(amps_test)));

%loop through all the traces. Essentially, this loop identifies which
%arrays come from the same cell, then concatenates the arrays and randomly
%samples amplitude and inter-event intervals from the concatenated traces
for k = 1:length(traces)   
    row_test = AMP_ALL{1,traces(k,1)};
    amps_test = row_test{1,traces(k,2)};
    
    row_times = TIME_INDICES{1,traces(k,1)};
    times_test = row_times{1,traces(k,2)};

   %the way this is written, the first trace has to be dealt with slightly differently
   %from the other traces. For this trace, first check if there is a second trace for that cell by
   %comparing the values in the 3rd column of the traces array (manually provided in cum_amps
   %script). If no, randomly sample the amplitudes and inter-event
   %intervals from this trace. If yes, remove the NaNs from the array of
   %amplitudes and inter-event intervals from this trace and proceed to the
   %next trace
    
    if k == 1        
        if traces(k,3) < traces(k+1,3)
            rand_events(:,traces(1,3)) = randsample(per_cell, downsample_val);
            first_events(:,traces(1,3)) = per_cell(1:downsample_val);
            
            rand_ISIs(:,traces(1,3)) = randsample(ISIs, downsample_val);
            first_ISIs(:,traces(1,3)) = ISIs(1:downsample_val);            
        else
            per_cell = amps_test(~isnan(amps_test));
            ISIs = diff(times_test(~isnan(amps_test)));
        end
        
        %From here, we use a more general approach of testing if the
        %current trace is from the same cell as the previous trace. If yes,
        %concatenate this trace with the previous. Then test if the
        %following trace is also from the same cell. If no, take the random
        %sample measurements. 
        
    elseif traces(k,3) == traces(k-1, 3)
        per_cell = horzcat(per_cell, amps_test(~isnan(amps_test)));
        ISIs = horzcat(ISIs, diff(times_test(~isnan(amps_test))));
        
        
        if k == length(traces)
            rand_events(:,traces(k-1,3)) = randsample(per_cell, downsample_val);
            first_events(:,traces(k-1,3)) = per_cell(1:downsample_val);

            rand_ISIs(:,traces(k-1,3)) = randsample(ISIs, downsample_val);
            first_ISIs(:,traces(k-1,3)) = ISIs(1:downsample_val);

        elseif traces(k,3) < traces(k+1,3)
            rand_events(:,traces(k-1,3)) = randsample(per_cell, downsample_val);
            first_events(:,traces(k-1,3)) = per_cell(1:downsample_val);
            
            rand_ISIs(:,traces(k-1,3)) = randsample(ISIs, downsample_val);
            first_ISIs(:,traces(k-1,3)) = ISIs(1:downsample_val);
        end
        
        %If the current trace is not from the same cell as the previous
        %trace, then take the random sample measurements here      
    else       
        per_cell = amps_test(~isnan(amps_test));
        ISIs = diff(times_test(~isnan(amps_test)));
    end
end