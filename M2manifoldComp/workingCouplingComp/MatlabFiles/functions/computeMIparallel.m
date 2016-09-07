function Z = computeMIparallel(length,dist,beta)
%Computes the mutual impedance between two parallel antennas side by side
%given an antenna length (length) and distance apart (dist).
    %NOTE: This only works for odd multiples of lambda/2 lengths of the
    %antenna
L = length;
d = dist;

%Real Part
R = 30*( 2*cosint(beta*d) - cosint(beta*(sqrt(d^2+L^2)+L)) - ...
    cosint(beta*(sqrt(d^2+L^2)-L)));

%Imaginary Part
X = -30*( 2*sinint(beta*d) - sinint(beta*(sqrt(d^2+L^2)+L)) - ...
    sinint(beta*(sqrt(d^2+L^2)-L)));

Z = R + 1j*X;