%% Scratch work
    
    
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
    %Vertical Dipole Manifold
        vA = @(beta,len,theta)  (cos(beta*len*cos(theta*pi/180)/2) - cos(beta*len/2)) / sin(theta*pi/180+0.000001);    
  

%% Load the NEC isolated responses and multiply with NEC computed coupling matrix
%Takes away misnomer phase errors on zero values
threshold = 10e-10;  
        
%Load Active Element Manifolds
name = 'ant1';
fileName = strcat(name,'.txt');
modelSpecs = readNecOutfile(fileName);
    %Acquire final Isolated Manifold
        [ThetaAll,PhiAll,thetaGainAnt1,phiGainAnt1,truthAnt1,radTotalGainAnt1,p1] = sph2rectRadPattern(modelSpecs,1);
        EthetaAnt1 = [thetaGainAnt1 phiGainAnt1];
        %Takes away misnomer phase errors on zero values
        [EthetaAnt1] = zeroTheZeros(EthetaAnt1,threshold);        
name = 'ant2';
fileName = strcat(name,'.txt');
modelSpecs = readNecOutfile(fileName);
    %Acquire final Isolated Manifold
        [ThetaAll,PhiAll,thetaGainAnt2,phiGainAnt2,truthAnt2,radTotalGainAnt2,p2] = sph2rectRadPattern(modelSpecs,2);
        EthetaAnt2 = [thetaGainAnt2 phiGainAnt2];
        %Takes away misnomer phase errors on zero values
        [EthetaAnt2] = zeroTheZeros(EthetaAnt2,threshold);        
             
%Load all Elements On Manifold
name = 'ant1ant2';
fileName = strcat(name,'.txt');
modelSpecs = readNecOutfile(fileName);
    %Acquire final Isolated Manifold
        [ThetaAll,PhiAll,thetaGainAnt1Ant2,phiGainAnt1Ant2,truthAnt1Ant2,radTotalGainAnt1Ant2,p] = sph2rectRadPattern(modelSpecs,1);
        EthetaAnt1Ant2 = [thetaGainAnt1Ant2 phiGainAnt1Ant2];
        %Takes away misnomer phase errors on zero values
        [EthetaAnt1Ant2] = zeroTheZeros(EthetaAnt1Ant2,threshold);     
%% Analytic Manifold

%Load some needed Constants
    %Antenna Length
    len = abs(modelSpecs.Z2(1) - modelSpecs.Z1(1));
    %Get Wave Parameters
    fc = modelSpecs.FREQUENCY;
    lambda = modelSpecs.WAVELENGTH/100;     %MUST FIX FOR NEGATIVE EXPONENTS!!!!!!!!
    beta = 2*pi/lambda;
    %Get Axis and Resolution   
    res = abs( modelSpecs.radPattern(1,1) - modelSpecs.radPattern(2,1) );  
    thetaAxis = -180:res:0;
    phiAxis = 0:res:360;   
    
%Compute Analytic Coupling Matrix
    %The Load
    ZL = 50;
    
    %Self Impedance
    n = 1; %since length of antenna is n* lambda/2
    Z11a = computeSI(n);
    Z22a = Z11a;
    
    %Mutual Impedance
    Z12a = computeMIparallel(len,lambda*0.2,beta);
    Z21a = computeMIparallel(len,lambda*0.2,beta);

    %The Matrix
    C1a = [1 Z12a/(Z22a+ZL); Z21a/(Z11a+ZL) 1];
    Ca = C1a^-1;
    
%Compute the Manifold    
manifoldV1iso = []; manifoldV2iso = [];    
manifoldV1 = []; manifoldV2 = [];
    for k = 1:length(phiAxis)
        for m = 1:length(thetaAxis)
            
            %Calculate Isolated Responses
            ShapeGain1 = exp(1j*2*pi/lambda*p1'*u(thetaAxis(m),phiAxis(k)));
                curManiV1 = ShapeGain1*vA(beta,len,thetaAxis(m));
                manifoldV1iso = [manifoldV1iso; curManiV1];
                curManiH1 = 0;
                curMani1 = [curManiV1 curManiH1];
                
            ShapeGain2 = exp(1j*2*pi/lambda*p2'*u(thetaAxis(m),phiAxis(k)));
            curManiV2 = ShapeGain2*vA(beta,len,thetaAxis(m));
                manifoldV2iso = [manifoldV2iso; curManiV2];
                curManiH2 = 0;
                curMani2 = [curManiV2 curManiH2];                   
                
            %Incorporate Coupling
            curMani = Ca*[curMani1 ; curMani2 ];
            curManifold = sum(curMani);
                manifoldV1 = [manifoldV1; curMani(1,1)];
                manifoldV2 = [manifoldV2; curMani(2,1)]; 
            
        end
    end 
%Incorporate the Horizontal Manifold    
manifoldH1 = zeros(size(manifoldV1));   manifoldH2 = manifoldH1;  
aManifold1 = [manifoldV1 manifoldH1];
aManifold2 = [manifoldV2 manifoldH2];  
  
%Find the Complex Constant for the coupled Manifold
[aManifold1] = zeroTheZeros(aManifold1,threshold); 
    constC1f = constManifoldMin(EthetaAnt1(:,1),aManifold1(:,1));
    aManifold1Scaled = constC1f*aManifold1;
[aManifold2] = zeroTheZeros(aManifold2,threshold); 
    constC2f = constManifoldMin(EthetaAnt2(:,1),aManifold2(:,1));
    aManifold2Scaled = constC2f*aManifold2;
    
    
%% Compare          

%All On vs. NEC generated Sum of Active Elements   
ActiveSum = EthetaAnt1 + EthetaAnt2;
plotManifolds(EthetaAnt1Ant2(:,1),ActiveSum(:,1),thetaAxis,phiAxis,' Total vs. Sum of Active Elements (Truths)')

%All On vs. Analytic Sum of Active Elements   
ActiveSum = aManifold1Scaled + aManifold2Scaled;
plotManifolds(EthetaAnt1Ant2(:,1),ActiveSum(:,1),thetaAxis,phiAxis,' Sum of Active Elements')

%Antenna 1:  NEC active Element vs. Analytic active element
plotManifolds(EthetaAnt1(:,1),aManifold1Scaled(:,1),thetaAxis,phiAxis,' Active Element Ant 1')

%Antenna 2:  NEC active Element vs. Analytic active element
plotManifolds(EthetaAnt2(:,1),aManifold2Scaled(:,1),thetaAxis,phiAxis,' Active Element Ant 2')

