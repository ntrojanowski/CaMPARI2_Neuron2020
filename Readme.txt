Many scripts require ws.loadDataFile from WaveSurfer

Scripts:

Raw_FR: %this script loops through all files in a folder to find the number of
%threshold crossings (action potentials) that occur in a voltage trace. 
%Using a user-defined threshold set on a trace-by-trace basis, it
%identifies each threshold crossing, and records the timestamp, then
%calculates the number of action potentials, as well as the CV of
%inter-event intervals and the Fano factor.
 
Spont_IE: %this script runs on a folder of .h5 files recorded in voltage clamps. It
%first takes a histogram of the traces, then uses that histogram to measure
%the total charge transfer during the trace. Briefly, it assumes that the 
%noise is symmetrical around the peak of the histogram and that the 
%recordings were taken at the reversal potential for either excitation or 
%inhibition, and therefore subtracts the values from the side of the 
%histogram with the shorter tail from the side with the longer tail. The
%sum of the remaining values represents the total charge transfer.

TestConnections:%this script is used for measuring the strengths of connections between simulatneously
%recorded neurons. It runs on a folder of .h5 files. It is best run on
%voltage clamp data - it produces nonsensical results when run on current
%clamp traces. It finds the 3 point average current around the peak, and compares
%this to the baseline current value used to hold cells near -80. It also
%calculates the latency from the current injection in one cell to the
%reponse in another. 

Excise_X_Points: % this script is used to select folders to search for data, then loop
% though the data files and allow the user to visually identify stable
% recording regions for analysis. Typically used for voltage clamp
% recordings. this script saves the x coordinates and file names for each trace

Cum_amps: %this script reads in the previously generated event amplitudes and timings, and a list of
%traces with stable regions for analysis, then uses the find_iso_amps
%function to randomly select the same amount of measurements from each
%cell and compare the distribution of the measurements between high and
%low activity neurons

FI_curve_and_PP: %this script analyzes trains of action potentials resulting from current
%injections (to produce FI curves), obtaining AP measurements from the
%get_AP_properties function. It also measures passive properties. this 
%script runs on a folder of .h5 files. This script relies on three custom
%functions: get_spike_threshold, get_Ra_doubleexp_h5, and get_AP_properties


Functions:

find_iso_amps: %this function, which is called in the Cum_amps scrpit, takes arrays containing event amplitudes and times and an
%array indiciating which traces to analyze, and returns the amplitude and
%inter-event intervals of a user-provided number of randomly selected events

get_spike_threshold: %this function takes a current clamp trace and dV threshold and finds the
%time and voltage at which the trace crosses the dV threshold

get_Ra_doubleexp_h5: %this function takes a voltage clamp trace containing a hyperpolarizing
%"seal test" and measures and returns various passive properties (series resistance, 
% input resistance, membrane capacitance, resting membrane potential, and
% time constant).

get_AP_properties: %This function takes a small section of a current clamp recording and dV threshold, then uses
%that threshold to identify an action potential. It then returns various
%properties of the action potential (spike and afterhyperpolarization amplitude, 
%half-width, mean and peak upslope and downslope values, and ratio of up-to-down peak slopes).



