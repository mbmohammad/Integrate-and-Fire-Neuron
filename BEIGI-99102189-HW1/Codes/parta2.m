%% parta
clc;clear all; close all;
% Simulate LIF model
[v, sp] = run_LIF(20, 0)

% Visualize
figure
plot(0.001:0.001:0.1, v)
xlabel("time(msec)")
ylabel("voltage(mV)")
title("time-course of the membrane potential")
%% part b
clc;clear all; close all;
tau_m = 30;
[ spikeMat , tVec ] = poissonSpikeGen ( 3000 , 1/10000, 0.1 , 1 );
for i = 1 : length(spikeMat)
    realisticCurrent(i) = 4.5*spikeMat(i) * i/10000 *exp(-i/100);
end
x = 1 : 1000;
y = exp(-x/10000/tau_m);
    voltages = 4000/tau_m*conv(realisticCurrent, y);
figure
plot(0.0001 : 0.0001 : 0.1999,voltages(1:1999))
title(["Stimulated neuron membrane potential by",...
    "a time-varying input current without thr"])
ylabel("Membrane potential(mV)")
xlabel("time(sec)")
thr = 25;
for i = 1 : length(voltages)
   if(voltages(i)>=thr) 
       voltages(i) = 40;
       voltages(i + 1) = -5;
       voltages(i + 2 : length(voltages)) = voltages(i + 2 : length(voltages))...
           - 20*ones(1, length(i + 2 : length(voltages)));
      i = i + 1;
   end
end
figure
plot(0.0001 : 0.0001 : 0.1,voltages(1:1000))
title("Stimulated neuron membrane potential by a time-varying input current")
ylabel("Membrane potential(mV)")
xlabel("time(sec)")
% figure
% plot(realisticCurrent)
% title("rc")
%% part c,d
close all

for p = 50 : 50 : 300 
    spikeinter = zeros(100000, 1000);
for l = 1 : 100000
A = ones(1, 1000);
r = randperm(1000, p);
A(r) = -1;
newRealisticCurrent = A.*realisticCurrent;
x = 1 : 1000;
y = exp(-x/10000/tau_m);
    newVoltages = 7000/tau_m*conv(newRealisticCurrent, y);
% figure
% plot(0.0001 : 0.0001 : 0.1999,newVoltages(1:1999))
% title(["Stimulated neuron membrane potential by",...
%     "a time-varying input current without thr",...
%     "considering IPSPs as same as EPSPs"])
% ylabel("Membrane potential(mV)")
% xlabel("time(sec)")
for i = 1 : length(newVoltages)
   if(newVoltages(i)>=15) 
       newVoltages(i) = 40;
       newVoltages(i + 1) = -5;
       newVoltages(i + 2 : length(newVoltages)) = newVoltages(i + 2 : length(newVoltages))...
           - 20*ones(1, length(i + 2 : length(newVoltages)));
      i = i + 1;
   end
end
newVoltages(newVoltages~= 40) = 0;
newVoltages(newVoltages==40) = 1;
spikeinter(l, 1:1000) = newVoltages(1:1000);
end
p
spikeIntervals = ISIHistogramPlot(spikeinter);
meanOfISIs = mean(spikeIntervals)
CV = sqrt(var(spikeIntervals))/mean(spikeIntervals)

end
% figure
% plot(0.0001 : 0.0001 : 0.1,newVoltages(1:1000))
% title(["Stimulated neuron membrane potential by a time-varying input current",...
% sprintf("considering %d percent IPSPs and others EPSPs", p/10)])
% ylabel("Membrane potential(mV)")
% xlabel("time(sec)")

%% part d
close all

for p = 50 : 50 : 300 
    spikeinter = zeros(100000, 1000);
for l = 1 : 100000
A = ones(1, 1000);
r = randperm(1000, p);
A(r) = -1;
newRealisticCurrent = A.*realisticCurrent;
x = 1 : 1000;
y = exp(-x/10000/tau_m);
    newVoltages = 7000/tau_m*conv(newRealisticCurrent, y);
% figure
% plot(0.0001 : 0.0001 : 0.1999,newVoltages(1:1999))
% title(["Stimulated neuron membrane potential by",...
%     "a time-varying input current without thr",...
%     "considering IPSPs as same as EPSPs"])
% ylabel("Membrane potential(mV)")
% xlabel("time(sec)")
for i = 1 : length(newVoltages)
   if(newVoltages(i)>=thr) 
       newVoltages(i) = 40;
       newVoltages(i + 1) = -5;
       newVoltages(i + 2 : length(newVoltages)) = newVoltages(i + 2 : length(newVoltages))...
           - 20*ones(1, length(i + 2 : length(newVoltages)));
      i = i + 1;
   end
end
newVoltages(newVoltages~= 40) = 0;
newVoltages(newVoltages==40) = 1;
spikeinter(l, 1:1000) = newVoltages(1:1000);
end
p
spikeIntervals = ISIHistogramPlot(spikeinter);
meanOfISIs = mean(spikeIntervals)
CV = sqrt(var(spikeIntervals))/mean(spikeIntervals)

end


%%
function [rec_v, rec_sp] = run_LIF( Iinj, stop)
  % Simulate the LIF dynamics with external input current
  % Set parameters
  V_th = 15
  V_reset = 0
  tau_m = 10
  g_L = 1000
  V_init = 0
  E_L = 0
  dt = 0.1
  range_t = 100
  Lt = range_t
  tref = 0.001
  % Initialize voltage
  v = zeros(Lt, 1);
  v(1) = V_init
  % Set current time course
  Iinj = Iinj * ones(Lt, 1)
  % If current pulse, set beginning and end to 0
  if stop == 1
    Iinj(1:floor(length(Iinj) / 2) - 1000) = 0;
    Iinj(floor(length(Iinj) / 2) + 1000:end) = 0;
  end
  % Loop over time
  rec_spikes = []; % record spike times
  tr = 0; % the count for refractory duration
  for it = 1:Lt - 1
      it
    if tr > 0 % check if in refractory period
      v(it) = V_reset; % set voltage to reset
      tr = tr - 1; % reduce running counter of refractory period
    elseif v(it) >= V_th % if voltage over threshold
      rec_spikes(end+1) = it; % record spike event
      v(it) = V_reset; % reset voltage
      tr = tref / dt; % set refractory time
    end
    % Calculate the increment of the membrane potential
    dv = (-g_L * (v(it) - E_L) + Iinj(it)) / tau_m * dt;
    % Update the membrane potential
    v(it + 1) = v(it) + dv;
  end
  % Get spike times in ms
  rec_spikes = rec_spikes * dt;
  rec_v = v;
  rec_sp = rec_spikes;
end


function [ spikeMat , tVec ] = poissonSpikeGen ( fr , dt, tSim , nTrials )
nBins = floor ( tSim / dt ) ;
spikeMat = rand ( nTrials , nBins ) < fr * dt ;
tVec = 0: dt : tSim - dt ;
end

function spikeIntervals = ISIHistogramPlot(spikeMat)
[r, c] = size(spikeMat);
    idx = 1;
for i = 1 : r
    temp = 0;
    for j = 1 : c
        if(spikeMat(i, j) == 1 )
            spikeIntervals(idx) = temp + 1;
            idx = idx + 1;
            temp = 0;
        elseif(spikeMat(i, j) == 0 )
            temp = temp + 1;
        end
    end
end
end

