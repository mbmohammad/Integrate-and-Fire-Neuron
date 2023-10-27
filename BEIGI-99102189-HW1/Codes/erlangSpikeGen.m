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