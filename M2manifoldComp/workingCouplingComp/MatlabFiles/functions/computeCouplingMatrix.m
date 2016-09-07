function cM = computeCouplingMatrix(SI, MI, ZL, numElements)

cM = zeros(numElements, numElements);
for outCount = 1:length(cM)
  for inCount = 1:length(cM)
    if outCount == inCount
      cM(outCount,inCount) = 1;
    else
      cM(outCount,inCount) = MI(outCount,inCount)/(SI(inCount) + ZL);
    end
  end
end
    