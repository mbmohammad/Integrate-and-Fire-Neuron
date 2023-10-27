%% part abc poisson
clc; 
clear all;
close all;
[ spikeMat , tVec ] = erlangSpikeGen ( 100 , 1/1000, 1 , 100000, 5 );
numOfOnes = GammaHistogramPlot(spikeMat, 5);
spikeIntervals = ISIGammaHistogramPlot(spikeMat, 5);
gamma5 = sqrt(var(spikeIntervals))/mean(spikeIntervals);
secondSpikeGenerating5 = sqrt(var(numOfOnes))/mean(numOfOnes)
[ spikeMat , tVec ] = erlangSpikeGen ( 100 , 1/1000, 1 , 100000, 9 );
numOfOnes = GammaHistogramPlot(spikeMat, 9);
spikeIntervals = ISIGammaHistogramPlot(spikeMat, 9);
gamma9 = sqrt(var(spikeIntervals))/mean(spikeIntervals);
secondSpikeGenerating9 = sqrt(var(numOfOnes))/mean(numOfOnes)
[ spikeMat , tVec ] = erlangSpikeGen ( 100 , 1/1000, 1 , 100000, 20 );
numOfOnes = GammaHistogramPlot(spikeMat, 20);
spikeIntervals = ISIGammaHistogramPlot(spikeMat, 20);
gamma20 = sqrt(var(spikeIntervals))/mean(spikeIntervals);
secondSpikeGenerating20 = sqrt(var(numOfOnes))/mean(numOfOnes)
%% functions
function [ spikeMat , tVec ] = erlangSpikeGen ( fr , dt, tSim , nTrials, k )
nBins = floor ( tSim / dt ) ;
spikeMat = rand ( nTrials , nBins ) < fr * dt ;
tVec = 0: dt : tSim - dt ;
[r, c] = size(spikeMat);
for i = 1 : r
    temp = 0;
    for j = 1 : c
        if(spikeMat(i, j) == 1)
            temp = temp + 1;
            if(temp ~= k)
                spikeMat(i, j) = 0;
            elseif(temp == k)
                temp = 0;
            end
        end
    end
end
end

function numOfOnes = GammaHistogramPlot(spikeMat, k)
numOfOnes = sum(transpose(spikeMat(:, :)) == 1);
figure; 
histogram(numOfOnes,15,'Normalization','probability', 'FaceColor','g');
hold on
% figure
x = 1:1:30;
for i = 1 : length(x)
   y(i) = k*exp(-100)*(100)^(k*x(i))/factorial(k*x(i));
end
plot(x, y,"r")
xlabel("number of spikes in a trial")
ylabel("pdf")
title(["histogram of number of spikes in each trial"...
    "superimposed with poisson(5x) with fr = 100"])
end

function spikeIntervals = ISIGammaHistogramPlot(spikeMat, k)
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
histogram(spikeIntervals/1000, 45, 'Normalization',...
    'probability', 'FaceColor','g');
 hold on
x = 0:0.001:0.2;
y = 0.0055*gampdf(x,k, 1/100);
plot(x, y,"r")
xlabel(sprintf("time to observe next %d th spikes", k))
ylabel("pdf")
title(sprintf("histogram of ISIs superimposed with gamma(%d , 100)", k))
end







