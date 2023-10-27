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
plot(voltages)
title("volt")
thr = 15;
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
