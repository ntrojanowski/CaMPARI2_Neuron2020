function [V_thresh, time] = get_spike_threshold(data, dV_thresh, plot_on)

%this function takes a current clamp trace and dV threshold and finds the
%time and voltage at which the trace crosses the dV threshold

%inputs:
%data - a current clamp trace with action potentials
%dV_thresh - the dV/dt threshold for the initiation of an action potential
%plot_on - optional, plots results


%determine whether or not to plot
if nargin < 3
    plot_on = 0;
end

%initialize data
V = data;
dV = diff(V);
time = NaN;

% get dV/dt (V/s)
dV_sec = dV.*10;
% find when dV/dt crosses V/s threshold
third_sample = find(dV_sec > dV_thresh,3,'first'); %used to weed out intial dV spike caused by current injection
if numel(third_sample) < 3 %if there is no action potential, only a subthreshold depolarization, return NaN
    V_thresh = NaN; 
    V_xi = NaN;
else
    next_sample = third_sample(3);
    prev_sample = find(dV_sec(1:next_sample) < dV_thresh,1,'last');
    
    % linear interpolation to find exact sample point at which the cross
    % happened, 
    dV_xi = dV_thresh; 
    dV_x = [dV_sec(prev_sample) dV_sec(next_sample)]';
    dV_y = [prev_sample next_sample]';
    dV_yi = interp1q(dV_x,dV_y,dV_xi);
    
    % linear interpolation on V to find Vm at which dV/dt crosses the
    % V/s threshold 
    V_x = dV_y + 1; %shifted to account for the difference taken to calculate dV
    V_y = V(V_x);
    V_xi = dV_yi + 1;
    V_yi = interp1q(V_x,V_y,V_xi);
    
    % this is our AP threshold and latency (time)
    V_thresh = V_yi;
    time = V_x(1);
    
end
%%
if plot_on
    figure(); hold on;
    plot(V,'linewidth',2);
    plot(dV_sec./10,'linewidth',1.5);
    plot(V_xi,V_thresh,'o','markersize',10,'markerfacecolor','k');
    set(gca,'xlim',[8000 15000]);
end



