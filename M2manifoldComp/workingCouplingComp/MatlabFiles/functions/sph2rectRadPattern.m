function [ThetaAll,PhiAll,thetaGain,phiGain,Truth,radTotalGain,p_m] = sph2rectRadPattern(modelSpecs,AntNum)
%Conversions
    d2r = pi/180;   %degrees to radians
    c = 299792458;        %speed of light [meters per second]
%For graphing
    colors = [1 0 0;
              0 1 0;
              0 0 1];
%Model Functions     
    %unit cosine in direction of propogation
        u = @(theta_k, phi_k) [cos(pi/180*phi_k)*sin(pi/180*theta_k); sin(pi/180*phi_k)*sin(pi/180*theta_k); cos(pi/180*theta_k)];
    %v is a vertical unit cosine perpendicular to u and h
        v = @(theta_k, phi_k) [cos(pi/180*phi_k)*cos(pi/180*theta_k); sin(pi/180*phi_k)*cps(pi/180*theta_k); -sin(pi/180*theta_k)];
    %h is a horizontal unit cosine perpendicular to u and v
        h = @(theta_k) [-sin(pi/180*phi_k); cos(pi/180*phi_k); 0];
    %Directinoal Gain Vector of a Vertical Dipole Antenna
        vDipole = @(phi_k) [0 -sin(pi/180*theta_k)];  
        
    res = abs( modelSpecs.radPattern(1,1) - modelSpecs.radPattern(2,1) );



%Extract Antenna Characteristics
x1 = modelSpecs.X1(AntNum);
x2 = modelSpecs.X2(AntNum);
y1 = modelSpecs.Y1(AntNum);
y2 = modelSpecs.Y2(AntNum);
z1 = modelSpecs.Z1(AntNum);
z2 = modelSpecs.Z2(AntNum);

xAve = (x1 + x2)/2;
yAve = (y1 + y2)/2;
zAve = (z1 + z2)/2;
    p_m = [xAve yAve zAve]';
    
%Extract Radiation Power Gain
ThetaAll = modelSpecs.radPattern(:,1)*pi/180;  
PhiAll = modelSpecs.radPattern(:,2)*pi/180;     
PwrRadGainAll = 10.^(modelSpecs.radPattern(:,3)/20);        %not in dB

%Extract Radiation Pattern (V/m)
thetaGain = modelSpecs.radPattern(:,4);
phiGain = modelSpecs.radPattern(:,5);
radTotalGain = abs( thetaGain + phiGain );

%Normalize
maxGain = max(radTotalGain);
radTotalGain = radTotalGain./maxGain;

%Produce Radiation 3d Plot    
thetaCounter = 180 / res + 1;
phiCounter = 360 / res + 1;
Truth = [];
for k = 1:phiCounter
    for m = 1:thetaCounter
        counter = (k-1)*thetaCounter + m;
        r = radTotalGain(counter);

        Truth = [Truth; r*u(ThetaAll(counter)*180/pi,PhiAll(counter)*180/pi)'];
    end
end    
    