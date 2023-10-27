function numOfOnes = GammaHistogramPlot(spikeMat, k)
numOfOnes = sum(transpose(spikeMat(:, :)) == 1);
figure; 
histogram(numOfOnes,15,'Normalization','probability', 'FaceColor','g');
hold on
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