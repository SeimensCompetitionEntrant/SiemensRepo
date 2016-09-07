function Z = computeSI(n)
%Computes the self impedance of antenna of length n*lambda/2 where n is odd
    %NOTE: This only works for odd multiples of lambda/2 lengths of the
    %antenna


Z = 30*( 0.577 + log(2*pi*n) - cosint(2*pi*n) + 1j*sinint(2*pi*n) );