function [ spikeMat , tVec ] = poissonSpikeGen ( fr , dt, tSim , nTrials )
nBins = floor ( tSim / dt ) ;
spikeMat = rand ( nTrials , nBins ) < fr * dt ;
tVec = 0: dt : tSim - dt ;
end