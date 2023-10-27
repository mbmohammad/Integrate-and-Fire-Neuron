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