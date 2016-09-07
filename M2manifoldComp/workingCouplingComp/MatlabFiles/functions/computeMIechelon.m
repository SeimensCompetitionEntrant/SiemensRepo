function Z = computeMIechelon(length,dist,height,beta)
%Computes the mutual impedance between two parallel antennas in echelon
%given an antenna length (length) and distance apart (dist).
    %NOTE: This only works for odd multiples of lambda/2 lengths of the
    %antenna
L = length;
d = dist;
h = height;


%Scalers
A = beta*( sqrt(d^2 + h^2) + h );
Ap = beta*( sqrt(d^2 + h^2) - h );
B = beta*( sqrt(d^2 + (h-L)^2) + (h-L) );
Bp = beta*( sqrt(d^2 + (h-L)^2) - (h-L) );
C = beta*( sqrt(d^2 + (h+L)^2) + (h+L) );
Cp = beta*( sqrt(d^2 + (h+L)^2) - (h+L) );

%Real Part
R = -15*cos(beta*h)*( -2*cosint(A)-2*cosint(Ap)+cosint(B)+cosint(Bp) + ...
    cosint(C)+cosint(Cp) ) + 15*sin(beta*h)*( 2*sinint(A)-2*sinint(Ap) - ...
    sinint(B)+sinint(Bp)-sinint(C)+sinint(Cp) );

%Imaginary Part
X = -15*cos(beta*h)*( 2*sinint(A)+2*sinint(Ap)-sinint(B)-sinint(Bp) - ...
    sinint(C)-sinint(Cp) ) + 15*sin(beta*h)*( 2*cosint(A)-2*cosint(Ap) - ...
    cosint(B)+cosint(Bp)-cosint(C)+cosint(Cp) );

Z = R + 1j*X;