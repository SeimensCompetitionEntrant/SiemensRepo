function [gainDb,coordinates] = plotGainPattern(modelSpecs,analyticModel,choice)
%Plots the Gain (in dB w.r.t. iso souce)
    %automatically turns "hold" on
    %"choice" is whether the function is being fed
        % 1)  A Pattern generated and extracted from NEC, a struct modelSpecs
        % 2)  An Analytic Pattern calculated, an array analyticModel with the
        % first column as thetas, 2nd as phis, 3rd as radiation values
    
%Constants
    d2r = pi/180;   %degrees to radians
    c = 299792458;        %speed of light [meters per second]
    eta = 376.8194; %intrinsic impedance of medium = sqrt(mu/eps)
%Model Functions     
    %unit cosine in direction of propogation
        u = @(theta_k, phi_k) [cos(pi/180*phi_k)*sin(pi/180*theta_k); sin(pi/180*phi_k)*sin(pi/180*theta_k); cos(pi/180*theta_k)];
    %v is a vertical unit cosine perpendicular to u and h
        v = @(theta_k, phi_k) [cos(pi/180*phi_k)*cos(pi/180*theta_k); sin(pi/180*phi_k)*cps(pi/180*theta_k); -sin(pi/180*theta_k)];
    %h is a horizontal unit cosine perpendicular to u and v
        h = @(theta_k) [-sin(pi/180*phi_k); cos(pi/180*phi_k); 0];
    %Find Resolution
        res = abs( modelSpecs.radPattern(1,1) - modelSpecs.radPattern(2,1) );
        
if(choice == 1)
    %Acquire Pattern
        [ThetaAll,PhiAll,thetaGain1,phiGain1,truth1,radTotalGain1,p_1] = sph2rectRadPattern(modelSpecs,1);
        manifoldTruth = sqrt( thetaGain1.^2 + phiGain1.^2 );  
    %Store manifold in generic name
        values = manifoldTruth;
end

if(choice == 2)
    %Store manifold in generic name
        values = analyticModel(:,3);
end


%Find ideal isolated point source power
Siso = modelSpecs.RadPower / (4*pi*modelSpecs.FarFieldDist^2);


%Stipulate Axes
    thetaAxis = -180:res:0;
    phiAxis = 0:res:360;
    
%Graph it
thetaCounter = 180 / res + 1;
phiCounter = 360 / res + 1;

coordinates = [];
gainDb = [];
    for k = 1:length(phiAxis)
        for m = 1:length(thetaAxis)
            counter = (k-1)*thetaCounter + m;
            
            %Note the values here come as:  sqrt( thetaGain1.^2 + phiGain1.^2 );
            E = values(counter);
            Sr = 0.5 * E^2 / eta;
            Srelative = abs(Sr)/Siso;
            Gdb = 10*log10(Srelative);
            gainDb = [gainDb; Gdb];
            
            
            newCoord = abs(Sr)*u(thetaAxis(m),phiAxis(k));
            coordinates = [coordinates newCoord];

        end
    end     
    
hold on



