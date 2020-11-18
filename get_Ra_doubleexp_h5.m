function [Ra_calc_scaled, Rt_scaled, Cm_scaled, Vm_scaled, tau_fast, tau_slow] = get_Ra_doubleexp_h5(data, step_start, vstep, samprate, figures_on)

%this function takes a voltage clamp trace containing a hyperpolarizing
%"seal test" and measures and returns various passive properties (series resistance, 
% input resistance, membrane capacitance, resting membrane potential, and
% time constant).

%Function inputs:
% data - contains a vector with the seal test (must be current clamp)
% step start - this is where the sealtest starts, should be in seconds
% vstep - the voltage of the step
% the sampling rate
% figures_on = 1 for plotting


data = data*10^-12; %convert to pA
%figure; plot(data);

V_test = -vstep;
V_hold = -7e-2; %holding potential, can be adjusted

%the most challenging passive property to calculate is the series
%resistance, because it is difficult to accurately measure the peak of the
%transient. To get around this, we fit an exponential to the decay of the
%transient to approximate the peak

t=step_start; %starting point for fit
tracevals = data(t*samprate:(t*samprate+7000));

peak = max(tracevals); %this is the measured peak
pulse_end = 0.5*samprate; %this is the pulse duration
startfit = find(tracevals==peak,1,'last');% offset to avoid the effects of pipette capacitance

ss_Start = 400; %timepoints after the peak used to define the steady state region
ss_End = 900; %timepoints after the peak used to define the steady state region

I_test = nanmedian(tracevals(pulse_end - 600:pulse_end - 100)); %calculate the average baseline current
I_baseline = nanmedian(tracevals(pulse_end + ss_Start:pulse_end + ss_End)); %calculate the average steady state current during the test
dI = I_test - I_baseline; %this is the measured current change during the depolarization
Rt = abs(V_test/dI); %total resistance

Vm = -(I_baseline * Rt) + V_hold; %calculate the resting membrane potential based on holding current and total resistance

%select the regions of the trace to use for fitting
whole_seal_fit_x = (startfit:startfit + ss_End)';
transient_seal_x = (pulse_end+1:pulse_end + ss_Start)';
ss_seal_fit_x = (startfit + ss_Start:startfit + ss_End)';  %steady state for BL subtract

%normalize the values of the transient to the baseline
baseline = nanmedian(tracevals(ss_seal_fit_x));
transient_seal_vals = tracevals(transient_seal_x) - baseline;
transient_vals = tracevals(whole_seal_fit_x) - baseline;
transient_vals = transient_vals*10^12; %converting units, though this isn't really necessary


%estimate the exponential fit, to provide a starting point
tau_est = (find(transient_vals<(transient_vals(1)*0.37), 1) - 1)/(samprate);

%select the region of the transient values to fit, and the values to use
transient_st = 4; 
transient_end = 90;
transient_vals = transient_vals(transient_st:transient_end); 
transient_vals = transient_vals - mean(transient_vals(end-5:end)); %normalize to the end of the trace (baseline)
time = (0:(1/samprate):(length(transient_vals)/samprate)-(1/samprate)).';


%filter out garbage traces, or those that have transients in the fitting
%region
if std(data)>1e-10
    Ra_calc_scaled = 0;
    Rt_scaled = 0;
    Cm_scaled = 0;
    Vm_scaled = 0;
    tau_fast = 0;
    tau_slow = 0;
    
    %if the data is good, fit it to a double exponential. The fast decay is
    %the one we want to use
else
    s = fitoptions('Method', 'NonlinearLeastSquares', ...
    'StartPoint', [transient_vals(1), tau_est*0.5, transient_vals(10), tau_est*1.5],...
    'Lower', [mean(transient_vals(1:2)*.9), tau_est*0.001, mean(transient_vals(1:5))*.1, tau_est*0.01],...
    'Upper', [transient_vals(1)*1.1, tau_est*5, transient_vals(1)*1.25, tau_est*15]);
    f = fittype('a*exp(-x/b) + c*exp(-x/d)','options',s);  
    [exp_fit,~] = fit(time,transient_vals,f);
    cval = coeffvalues(exp_fit);
    
    %Combine the area under the curve to calculate the capacitance, then
    %use this to find the series resistance
    fast_Tau = cval(2); %the first (fastest) exponent is the one we want
    Q1 = sum(abs(transient_seal_vals(1:200).*(1/samprate))); %this is the steady state charge    
    Q2 = -dI * fast_Tau; %this is the transient charge
    Qt = Q1 + Q2; %total charge
    Cm = Qt/V_test; 
    Ra_calc = fast_Tau/Cm;
    
    %Originally used, but replaced with above code that proved more
    %accurate when using a model cell
%     num_pts_back_for_extrap = ((11+transient_st)-find(tracevals(startfit-10:startfit)>I_test+5*I_test_std,1));
%     %%%%need to compensate above line for the shift in the first transient
%     %%%%value above
%     
%     %disp(num_pts_back_for_extrap)
%     [vals] = find(tracevals==peak);
%     if length(vals) < transient_st + 1
%         if num_pts_back_for_extrap > transient_st + 1
%             num_pts_back_for_extrap = transient_st + 1;
%         end
%     end
%     
%     
%     Rs_est = V_test/((exp_fit((-(num_pts_back_for_extrap/samprate)))*10^-12) + abs(I_test)); %usually closest off by around 2-3 Mohms for model
    
    if figures_on == 1
        figure;
        plot(tracevals)
        figure;
        plot(exp_fit,time,transient_vals)
        %text(10,50,num2str(Rs*1e-6))
        %axis([0 50 min(transient_seal_vals) 100])
    end
    
    Ra_calc_scaled = Ra_calc*1e-6; %MOhm
    Rt_scaled = Rt*1e-6; %MOhm
    Cm_scaled = Cm*1e12; %pF
    Vm_scaled = Vm*1e3; %mV
    tau_fast = cval(2);
    tau_slow = cval(4);
end