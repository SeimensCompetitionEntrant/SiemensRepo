function MI = computeMIarray(len,separation,beta,numElements)
MI = zeros(numElements,numElements);

for i = 1:numElements
  for k = 1:numElements
    if i != k
      MI(i,k) = computeMIparallel(len,abs(k-i)*separation,beta);
    end
  end
end
    