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