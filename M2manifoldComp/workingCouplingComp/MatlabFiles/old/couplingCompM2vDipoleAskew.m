%% Comparison of Manifolds of 2 element array, askew and coplanar
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
    %Rotation Matrices
        rotX = @(beta) [1 0 0; 0 cos(d2r*beta) sin(d2r*beta); 0 -sin(d2r*beta) cos(d2r*beta)];
    %Rect to Sph Transformations
        azi = @(x,y) atan2(y,x);
        elev = @(x,y,z) atan2(sqrt(x^2 + y^2),z);                
%% Load NEC Manifold with both turned on
name = 'hwDipoleVaskewM2';
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
axis([-4 4 0 15 0 20])
 
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
        Z11a = Zself; Z22a = Zself;
        
        one = [modelSpecs.X1(1) modelSpecs.Y1(1) modelSpecs.Z1(1) modelSpecs.X2(1) modelSpecs.Y2(1) modelSpecs.Z2(1)]';
        two = [modelSpecs.X1(2) modelSpecs.Y1(2) modelSpecs.Z1(2) modelSpecs.X2(2) modelSpecs.Y2(2) modelSpecs.Z2(2)]';
    antennaDim = [one two];    
    Z12a = computeMItilted(antennaDim);
        Z21a = Z12a; 

    %The Matrix
    C1a = [1 Z12a/Z22a; Z21a/Z11a 1];
    Ca = C1a^-1;
        
%Compute Isolated Efields
EfieldVertical = [];
Efield1iso = [];
Efield2iso = [];
    for k = 1:length(phiAxis)
        for m = 1:length(thetaAxis)
            
            ShapeGain1 = exp(1j*2*pi/lambda*p_all(1,:)*u(thetaAxis(m),phiAxis(k)));
                IsoAnt1 = ShapeGain1*eThetaNoDist(Iiso,rDist,beta,len,thetaAxis(m));
                Efield1iso = [Efield1iso; IsoAnt1];
            ShapeGain2 = exp(1j*2*pi/lambda*p_all(2,:)*u(thetaAxis(m),phiAxis(k)));
                IsoAnt2 = ShapeGain2*eThetaNoDist(Iiso,rDist,beta,len,thetaAxis(m));  
                Efield2iso = [Efield2iso; IsoAnt2];
                
            holder = Ca*[IsoAnt1 ; IsoAnt2 ];
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

