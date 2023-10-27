%A = [1 1 1 1 1  0 0 1 0 1 0 1]
resultSpikes = coincidenceDetector(A, 1/1000, 5, 2)
%%
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
            end
        end
        if(temp >= N)
            resultSpikes(i, j + D) = 1;
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