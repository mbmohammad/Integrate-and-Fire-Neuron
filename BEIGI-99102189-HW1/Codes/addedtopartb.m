A = ones(1, 1000);
r = randperm(1000, 300);
A(r) = -1;
newRealisticCurrent = A.*realisticCurrent;