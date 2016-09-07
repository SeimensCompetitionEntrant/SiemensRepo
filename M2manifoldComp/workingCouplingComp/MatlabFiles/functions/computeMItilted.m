function Z = computeMItilted(antennaDim)
%Computes the mutual impedance between two antennas - one vertical and the
%other tilted and in echelon form, displaced on the y and z axis
    %antennaDim is 6x2 and contains all of the dimensions (x1 y1 z1 x2 y2
    %z2) in a column of one of the two antennas
    
    %This is based on the paper by H.C. Baker and H. LaGrone 1962 "Digital
    %Computation of the Mutual Impedance Between Thin Dipoles"
    
    %x1,y1,z1 should be lower than x2,y2,z2 in value
    
    %antennaDim should be in wavelength units

%Degrees to Radians
    d2r = pi/180;    
%unit cosine in direction of propogation
    u = @(theta_k, phi_k) [cos(pi/180*phi_k)*sin(pi/180*theta_k); sin(pi/180*phi_k)*sin(pi/180*theta_k); cos(pi/180*theta_k)];     
%Rotation Matrices
    rotX = @(beta) [1 0 0; 0 cos(d2r*beta) sin(d2r*beta); 0 -sin(d2r*beta) cos(d2r*beta)];
%Rect to Sph Transformations
    azi = @(x,y) atan2(y,x);
    elev = @(x,y,z) atan2(sqrt(x^2 + y^2),z);

%% Constants (vectorized for numerical integration)        
ant1 = antennaDim(:,1);
ant2 = antennaDim(:,2);

    %Extract Antenna Coordinates
    x1 = ant1(1:3:end);
    y1 = ant1(2:3:end);
    z1 = ant1(3:3:end);  
    
    x2 = ant2(1:3:end);
    y2 = ant2(2:3:end);
    z2 = ant2(3:3:end);
    
    %Find Antenna Lengths
    L1 = sqrt( (x1(2)-x1(1))^2 + (y1(2)-y1(1))^2 + (z1(2)-z1(1))^2  );
    L2 = sqrt( (x2(2)-x2(1))^2 + (y2(2)-y2(1))^2 + (z2(2)-z2(1))^2  );

    %defining the integrating variable (s) & Continuity of Derivative
    N = 1000;
    s = linspace(-L2/2,L2/2,N);    
    
    %Find Antenna centers
    xAve1 = sum(x1)/2;
    yAve1 = sum(y1)/2;
    zAve1 = sum(z1)/2;
    
    xAve2 = sum(x2)/2;
    yAve2 = sum(y2)/2;
    zAve2 = sum(z2)/2;   
    
    %Find Tilt    
    x2new = x2 - xAve2;
    y2new = y2 - yAve2;
    z2new = z2 - zAve2;
    
        theta = elev(x2new(2),y2new(2),z2new(2));
        phi = azi(x2new(2),y2new(2));
        
    %Find cartesian components of the vector s made by tilted dipole (to be
    %later multiplied by s)
        sx = s*sin(theta)*cos(phi);
        sy = s*sin(theta)*sin(phi);
        sz = s*cos(theta);
        
    %Get Antenna Displacement
    displacement = [xAve2 yAve2 zAve2] - [xAve1 yAve1 zAve1];
        xo = displacement(1);
        yo = displacement(2);
        zo = displacement(3);
    
    %calculate radials
        p = sqrt( (xo+sx).^2 + (yo+sy).^2  );
        r = sqrt(p.^2 + (zo+sz).^2);
        r1 = sqrt( p.^2 + (zo+sz+L1/2).^2);
        r2 = sqrt( p.^2 + (zo+sz-L1/2).^2);

%% Numerical Integration
R21ds = -30*( (1./(p.^2).*( sin(2*pi*r1).*((sz+zo+L1/2)./r1 ) + sin(2*pi*r2).*((sz+zo-L1/2)./r2 ) - 2*cos(pi*L1).*sin(2*pi.*r).*((sz+zo)./r) ).*(sx.^2+yo*sy+sy.^2)) + ...
    ((2*(sin(2*pi*r).*cos(pi*L1))./r - sin(2*pi*r1)./r1 - sin(2*pi*r2)./r2).*sz) ) .* ( sin(2*pi*(L2/2 - abs(s)))./s  );
R21 = simps(s,R21ds);

X21ds = -30*( (1./(p.^2).*( cos(2*pi*r1).*((sz+zo+L1/2)./r1 ) + cos(2*pi*r2).*((sz+zo-L1/2)./r2 ) - 2*cos(pi*L1).*cos(2*pi.*r).*((sz+zo)./r) ).*(sx.^2+yo*sy+sy.^2)) + ...
    ((2*(cos(2*pi*r).*cos(pi*L1))./r - cos(2*pi*r1)./r1 - cos(2*pi*r2)./r2).*sz) ) .* ( sin(2*pi*(L2/2 - abs(s)))./s  );
X21 = simps(s,X21ds);

Z = R21 + X21*1j;













