%% Comparison of Manifolds of 4 element array, parallel in echelon
    %Patrick Ellis
    %April 8th, 2015       
    %Assumes excitation in each antenna of 1V
    
addpath(genpath('/Users/PatrickEllis/Desktop/Research/current/workingCouplingComp')) 


clear all
close all
clc

% tic

%Constants
    d2r = pi/180;   %degrees to radians
    c = 299792458;        %speed of light [meters per second]
    eta = 376.8194; %intrinsic impedance of medium = sqrt(mu/eps)
%Functions     
    %unit cosine in direction of propogation
        u = @(theta_k, phi_k) [cos(pi/180*phi_k)*sin(pi/180*theta_k); sin(pi/180*phi_k)*sin(pi/180*theta_k); cos(pi/180*theta_k)];       
    %Etheta
        eThetaDist = @(Io,dist,beta,len,theta) 1j*60*(Io)*exp(-1j*beta*dist)/dist * ...
                    (cos(beta*len*cos(theta*pi/180)/2) - cos(beta*len/2)) / sin(theta*pi/180+0.000001);
        eThetaNoDist = @(Io,dist,beta,len,theta) 1j*60*(Io) * ...
                    (cos(beta*len*cos(theta*pi/180)/2) - cos(beta*len/2)) / sin(theta*pi/180+0.000001); 
%% Load NEC Manifold with both turned on
name = 'vDipoleArrayM4staggered';
fileName = strcat(name,'.txt');
modelSpecs = readNecOutfile(fileName);   
    %Acquire final Manifold
        [ThetaAll,PhiAll,thetaGainAllOn,phiGainAllOn,truthAllOn,radTotalGainAllOn,p_AllOn] = sph2rectRadPattern(modelSpecs,1);
        manifoldTruth = [thetaGainAllOn phiGainAllOn];
        %Takes away misnomer phase errors on zero values
        threshold = 10e-10;
        [manifoldTruth] = zeroTheZeros(manifoldTruth,threshold); 
        %Find Terminal Current
        Io = modelSpecs.current; 
        %Find Resolution of Recording 
        res = abs( modelSpecs.radPattern(1,1) - modelSpecs.radPattern(2,1) );  
            %Define Axis
            thetaAxis = -180:res:0;
            phiAxis = 0:res:360;
%See Orientation
p_all = plotArrayConfig(modelSpecs);

 
%% Analytic Manifold
    %Load Wave Parameters
        fc = modelSpecs.FREQUENCY;
        lambda = modelSpecs.WAVELENGTH;                  %wavelength [meters]
        beta = 2*pi/lambda;
    %Load Antenna Position and Length
        len = abs(modelSpecs.Z2(1) - modelSpecs.Z1(1));
    %Find Resolution of Recording 
        res = abs( modelSpecs.radPattern(1,1) - modelSpecs.radPattern(2,1) );
    %Load FarField Distance    
        rDist = modelSpecs.FarFieldDist;       
    %Produce Radiation 3d Plot    
        thetaAxis = -180:res:0;
        phiAxis = 0:res:360;
    %Get Isolated Current
        name = 'arrayVant1iso';
        fileName = strcat(name,'.txt');
        isoModelSpecs = readNecOutfile(fileName);  
        Iiso = isoModelSpecs.current;
        
%Compute 4x4 Analytic Coupling Matrix
    %Self and Mutual Impedance
    n = 1; %since length of antenna is n* lambda/2
    Zself = computeSI(n); 
        Z11 = Zself; Z22 = Zself; Z33 = Zself; Z44 = Zself;
        
    Z12 = computeMIechelon(len,max( abs(p_all(1,1:2) - p_all(2,1:2)) ), ...
                            abs(p_all(1,3) - p_all(2,3)),beta);
        Z21 = Z12;
    Z13 = computeMIechelon(len,max( abs(p_all(1,1:2) - p_all(3,1:2)) ), ...
                            abs(p_all(1,3) - p_all(3,3)),beta);
        Z31 = Z13;
    Z14 = computeMIechelon(len,max( abs(p_all(1,1:2) - p_all(4,1:2)) ), ...
                            abs(p_all(1,3) - p_all(4,3)),beta);
        Z41 = Z14;
    Z23 = computeMIechelon(len,max( abs(p_all(2,1:2) - p_all(3,1:2)) ), ...
                            abs(p_all(2,3) - p_all(3,3)),beta);
        Z32 = Z23;
    Z24 = computeMIechelon(len,max( abs(p_all(2,1:2) - p_all(4,1:2)) ), ...
                            abs(p_all(2,3) - p_all(4,3)),beta);
        Z42 = Z24;
    Z34 = computeMIechelon(len,max( abs(p_all(3,1:2) - p_all(4,1:2)) ), ...
                            abs(p_all(3,3) - p_all(4,3)),beta);
        Z43 = Z34;   

    %The Matrix
    C1a = [ 1 Z12/Z22 Z13/Z33 Z14/Z44; ...
            Z21/Z11 1 Z23/Z33 Z24/Z44; ...
            Z31/Z11 Z32/Z22 1 Z34/Z44; ...
            Z41/Z11 Z42/Z22 Z43/Z33 1];
    Ca = C1a^-1;
        
%Compute Isolated Efields
EfieldVertical = [];
Efield1iso = [];
Efield2iso = [];
Efield3iso = [];
Efield4iso = [];
    for k = 1:length(phiAxis)
        for m = 1:length(thetaAxis)
            
            ShapeGain1 = exp(1j*2*pi/lambda*p_all(1,:)*u(thetaAxis(m),phiAxis(k)));
                IsoAnt1 = ShapeGain1*eThetaNoDist(Iiso,rDist,beta,len,thetaAxis(m));
                Efield1iso = [Efield1iso; IsoAnt1];
            ShapeGain2 = exp(1j*2*pi/lambda*p_all(2,:)*u(thetaAxis(m),phiAxis(k)));
                IsoAnt2 = ShapeGain2*eThetaNoDist(Iiso,rDist,beta,len,thetaAxis(m));  
                Efield2iso = [Efield2iso; IsoAnt2];
            ShapeGain3 = exp(1j*2*pi/lambda*p_all(3,:)*u(thetaAxis(m),phiAxis(k)));
                IsoAnt3 = ShapeGain3*eThetaNoDist(Iiso,rDist,beta,len,thetaAxis(m));
                Efield3iso = [Efield3iso; IsoAnt3];
            ShapeGain4 = exp(1j*2*pi/lambda*p_all(4,:)*u(thetaAxis(m),phiAxis(k)));
                IsoAnt4 = ShapeGain4*eThetaNoDist(Iiso,rDist,beta,len,thetaAxis(m));  
                Efield4iso = [Efield4iso; IsoAnt4];                
                
                
            holder = Ca*[IsoAnt1 ; IsoAnt2 ; IsoAnt3 ; IsoAnt4];
            curMani = sum(holder);

            EfieldVertical = [EfieldVertical; curMani];   
        end
    end 
EfieldHorizontal = zeros(size(EfieldVertical)); 

analyticArrayManifold = [EfieldVertical EfieldHorizontal];


%Do the same to the analytic manifold as we did to the truth    
[analyticArrayManifold] = zeroTheZeros(analyticArrayManifold,threshold); 

%% Analyze & Plot Results
const = 10e-8;


%First do the truth and the nec computations
%Change shape of manifolds
lenTheta = length(thetaAxis);
lenPhi = length(phiAxis);

imageAnalyticManifold = reshape(analyticArrayManifold(:,1),lenTheta,lenPhi);
imageTruthManifold = reshape(manifoldTruth(:,1),lenTheta,lenPhi);

figure;
subplot(321)
    imagesc(phiAxis,thetaAxis,abs(imageTruthManifold))
    xlabel('\phi Azimuth'),ylabel('\theta Elevation')
    title('| Truth |'),grid,colorbar
subplot(322)
    imagesc(phiAxis,thetaAxis,abs(imageAnalyticManifold))
    xlabel('\phi Azimuth'),ylabel('\theta Elevation')
    title('| Analytic |'),grid,colorbar
subplot(323)
    imagesc(phiAxis,thetaAxis,angle(imageTruthManifold))
    xlabel('\phi Azimuth'),ylabel('\theta Elevation')
    title('< Truth '),grid,colorbar
subplot(324)
    imagesc(phiAxis,thetaAxis,angle(imageAnalyticManifold))
    xlabel('\phi Azimuth'),ylabel('\theta Elevation')
    title('< Analytic '),grid,colorbar    
subplot(325)
    truthAbs = abs(imageTruthManifold) + const;
    analyticAbs = abs(imageAnalyticManifold) + const;
    toPlot = truthAbs./analyticAbs;
    imagesc(phiAxis,thetaAxis(1:end),toPlot(1:end,:))
    xlabel('\phi Azimuth'),ylabel('\theta Elevation')
    title('| Truth | / | Analytic |'),grid,colorbar
subplot(326)
    truthAngle = angle(imageTruthManifold) + const;
    analyticAngle = angle(imageAnalyticManifold) + const;
    toPlot = abs( truthAngle ./ analyticAngle );
    imagesc(phiAxis,thetaAxis(1:end),toPlot(1:end,:))
    xlabel('\phi Azimuth'),ylabel('\theta Elevation')
    title('<Truth  /  <Analytic '),grid,colorbar

