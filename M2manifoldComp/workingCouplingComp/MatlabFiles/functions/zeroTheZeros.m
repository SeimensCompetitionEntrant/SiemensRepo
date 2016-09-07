function [editedManifold] = zeroTheZeros(manifold,threshold)
%This function goes through the values given in "manifold" and makes those
%that are ~zero... actually zero.  This is done for purposes of avoiding
%misnomer phase errors in plots.  A const to decide what is zero is done
%based on average values and a user stipulated threshold.
editedManifold = manifold;
const = threshold;

%Assumes manifold has different angle pairs in the columns 
    %[thetaGain phiGain]

numPts = length(manifold(:,1));

for i = 1:numPts
    curTheta = abs( manifold(i,1) );
    curPhi = abs( manifold(i,2) );

    if(curTheta < const)
        editedManifold(i,1) = 0;
    end
    if(curPhi < const)
        editedManifold(i,2) = 0;
    end    
end
