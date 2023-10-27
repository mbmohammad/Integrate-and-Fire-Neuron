%%
% clc; 
% clear all;
% close all;
% CVs = zeros(3, 20);
% means = zeros(3, 20);
%%
k = 1;
deltaT = 1/1000;
for i = 1 : 18
[ spikeMat , tVec ] = erlangSpikeGenRef ( 1000 - 50*(i - 1) , 1/1000, 5 , 5000, k, deltaT );
spikeIntervals = ISIGammaHistogramPlotRef(spikeMat, k);
i
CVs(1, i) = sqrt(var(spikeIntervals))/mean(spikeIntervals)
means(1, i) = mean(spikeIntervals)
end
%%
for i = 19 : 35
[ spikeMat , tVec ] = erlangSpikeGenRef ( 100 - 5*(i - 18) , 1/1000, 5 , 5000, k, deltaT );
spikeIntervals = ISIGammaHistogramPlotRef(spikeMat, k);
i
CVs(1, i) = sqrt(var(spikeIntervals))/mean(spikeIntervals)
means(1, i) = mean(spikeIntervals)
end
for i = 36 : 40
[ spikeMat , tVec ] = erlangSpikeGenRef ( 15 - 2*(i - 35) , 1/1000, 5 , 5000, k, deltaT );
spikeIntervals = ISIGammaHistogramPlotRef(spikeMat, k);
i
CVs(1, i) = sqrt(var(spikeIntervals))/mean(spikeIntervals)
means(1, i) = mean(spikeIntervals)
end
%%
k = 4;
for i = 1 : 30
[ spikeMat , tVec ] = erlangSpikeGenRef ( 1000 - 30*(i - 1) , 1/10000, 5 , 1000, k, deltaT );
spikeIntervals = ISIGammaHistogramPlotRef(spikeMat, k);
i
CVs(2, i) = sqrt(var(spikeIntervals))/mean(spikeIntervals)
means(2, i) = mean(spikeIntervals)
end
for i = 31 : 48
[ spikeMat , tVec ] = erlangSpikeGenRef ( 100 - 5*(i - 30) , 1/10000, 5 , 1000, k, deltaT );
spikeIntervals = ISIGammaHistogramPlotRef(spikeMat, k);
i
CVs(2, i) = sqrt(var(spikeIntervals))/mean(spikeIntervals)
means(2, i) = mean(spikeIntervals)
end
%%
k = 16;
for i = 1 : 40
[ spikeMat , tVec ] = erlangSpikeGenRef ( 1000 - 23*(i - 1) , 1/10000, 5 , 2000, k, deltaT );
spikeIntervals = ISIGammaHistogramPlotRef(spikeMat, k);
CVs(3, i) = sqrt(var(spikeIntervals))/mean(spikeIntervals)
means(3, i) = mean(spikeIntervals)
end
for i = 41 : 45
[ spikeMat , tVec ] = erlangSpikeGenRef ( 80 - 5*(i - 40) , 1/10000, 5 , 800, k, deltaT );
spikeIntervals = ISIGammaHistogramPlotRef(spikeMat, k);
CVs(3, i) = sqrt(var(spikeIntervals))/mean(spikeIntervals)
means(3, i) = mean(spikeIntervals)
end
%%

figure;
plot(means(1, 1:39), CVs(1, 1:39))
hold on
plot(means(2, :)/10, CVs(2, :))
hold on
plot(means(3, 1:45)/10, CVs(3, 1:45))
hold on
yline(1)
hold on
yline(0.5)
hold on
yline(0.25)
xlabel("mean of ISIs(msec)")
ylabel("CV")
title("CVs for different Nth(num of spikes to fire) and t0(refractory period)")
legend('Nth = 1, t0 = 1ms','Nth = 4, t0 = 1ms', 'Nth = 16, t0 = 1ms', 'Nth = 1, t0 = 0ms','Nth = 4, t0 = 0ms', 'Nth = 16, t0 = 0ms')
%% functions
function [ spikeMat , tVec ] = erlangSpikeGenRef ( fr , dt, tSim , nTrials, k , t0)
nBins = floor ( tSim / dt ) ;
refractoryPeriod = floor ( t0/dt ) ;
for i = 1 : nTrials
    tmp = -inf;
    for j = 1 : nBins
        if(j < tmp + refractoryPeriod)
            spikeMat(i, j) = 0;
        elseif(j >= tmp + refractoryPeriod)
            spikeMat(i, j) = rand (1) < fr * dt ;
        end
        if(spikeMat(i, j) == 1)
            tmp = j;
        end
    end 
end
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

function spikeIntervals = ISIGammaHistogramPlotRef(spikeMat, k)
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