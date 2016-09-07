function SI = computeSIarray(n, numElements);
SI = zeros(1,numElements);
SI(1) = computeSI(n);
for i = 2:numElements
  SI(i) = SI(1);
end