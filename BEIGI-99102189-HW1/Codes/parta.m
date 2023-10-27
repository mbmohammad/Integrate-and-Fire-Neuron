%% part abc poisson
clc; 
clear all;
close all;
[ spikeMat , tVec ] = poissonSpikeGen ( 100 , 1/1000, 1 , 30000 );
numOfOnes = histogramPlot(spikeMat);
spikeIntervals = ISIHistogramPlot(spikeMat);
sqrt(var(spikeIntervals))/mean(spikeIntervals)
sqrt(var(numOfOnes))/mean(numOfOnes)
%% functions
function [ spikeMat , tVec ] = poissonSpikeGen ( fr , dt, tSim , nTrials )
nBins = floor ( tSim / dt ) ;
spikeMat = rand ( nTrials , nBins ) < fr * dt ;
tVec = 0: dt : tSim - dt ;
end

function numOfOnes = histogramPlot(spikeMat)
numOfOnes = sum(transpose(spikeMat(:, :)) == 1);
figure; 
histogram(numOfOnes,'Normalization','probability', 'FaceColor','g');
hold on
x = 0:150;
y = poisspdf(x,100);
plot(60:150, y(60:150),"r")
xlabel("number of spikes in a trial")
ylabel("probability")
title(["normalized histogram of number of spikes in each trial",...
    "superimposed with poisson mean 100"])
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
figure; 
histogram(spikeIntervals/1000, [0.000001:0.001:0.070],'Normalization',...
    'probability', 'FaceColor','g');
hold on
x = 0.001:0.001:0.07;
y = exppdf(x*1000,10);
plot(x, y,"r")
xlabel("time to observe next spike")
ylabel("probability")
title(["normalized histogram of ISIs",...
    "superimposed with exponential mean 0.1"])
end
