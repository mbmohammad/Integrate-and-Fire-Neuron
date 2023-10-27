clc;clear all; close all;
dt  = 1/10000;
[ spikeMat , tVec ] = erlangSpikeGenRef ( 2500 , dt, 0.1 , 2000, 1 , 5/10000 );
p = 100;
A = ones(1, 1000);
r = randperm(1000, p);
A(r) = -1;
spikeMat = A.*spikeMat;
for D = 1 : 10
    for N = 2 : 8
        i = D
        N
        resultSpikes = coincidenceDetector(spikeMat, dt, i, N);
        spikeIntervals = ISIGammaHistogramPlotRef(resultSpikes, 1);
        CV(D , N - 1)  = sqrt(var(spikeIntervals))/mean(spikeIntervals);
    end
end
CV
%% plots
close all
figure
plot(2:20,CV(1, :))
hold on
plot(2:20,CV(2, :))
hold on
plot(2:20,CV(3, :))
hold on
plot(2:20,CV(4, :))
hold on
plot(2:20,CV(5, :))
hold on
plot(2:20,CV(6, :))
legend('D = 15','D = 16','D = 17','D = 18','D = 19', 'D = 20', 'Location','northwest')
xlabel("N")
ylabel("CV")
title("Values of CV per N for 7 different values of D")

%%
function [ spikeMat , tVec ] = poissonSpikeGen ( fr , dt, tSim , nTrials )
nBins = floor ( tSim / dt ) ;
spikeMat = rand ( nTrials , nBins ) < fr * dt ;
tVec = 0: dt : tSim - dt ;
end

function resultSpikes = coincidenceDetector(spikeMat, dt, D, N)
[r, c] = size(spikeMat);
resultSpikes = zeros(r, c);
D = D/1000/dt;
for i = 1 : r
    for j = 1 : c - D
        temp = 0;
        for k = 0 : D - 1
            if(spikeMat(i, j + k) == 1)
                temp = temp + 1;
            elseif(spikeMat(i, j + k) == -1)
                temp = temp - 1;
            end
        end
        if(temp >= N)
            resultSpikes(i, floor(j + D)) = 1;
        end
        
    end
end
for i = 1 : r
    for j = 1 : D
        temp = 0;
        for k = 1 : j - 1
            if(spikeMat(i, k) == 1)
                temp = temp + 1;
            end
        end
        if(temp >= N)
            resultSpikes(i, j) = 1;
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