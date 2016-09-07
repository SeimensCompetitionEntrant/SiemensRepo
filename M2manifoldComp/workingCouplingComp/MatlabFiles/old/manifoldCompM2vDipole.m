%% Compares manifolds in a 2 element parallel side by side array
    %Patrick Ellis
    %April 14th, 2015
    %Assumes 1 volt excitation of antennas
    
    
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

%Load Isolated Manifolds and get Self Impedances
name = 'arrayVant1iso';
fileName = strcat(name,'.txt');
modelSpecs = readNecOutfile(fileName);   
    %Acquire final Isolated Manifold
        [ThetaAll,PhiAll,thetaGainAnt1iso,phiGainAnt1iso,truthAnt1iso,radTotalGainAnt1iso,p1] = sph2rectRadPattern(modelSpecs,1);
        EthetaAnt1iso = [thetaGainAnt1iso phiGainAnt1iso];
        %Takes away misnomer phase errors on zero values
        [EthetaAnt1iso] = zeroTheZeros(EthetaAnt1iso,threshold);
    %Find Self Impedance
        Iiso1 = modelSpecs.current; Viso1 = 1;
        Z11 = Viso1 / Iiso1;
    %Antenna Length
        len1 = abs(modelSpecs.Z2(1) - modelSpecs.Z1(1));
name = 'arrayVant2iso';
fileName = strcat(name,'.txt');
modelSpecs = readNecOutfile(fileName);   
    %Acquire final Isolated Manifold
        [ThetaAll,PhiAll,thetaGainAnt2iso,phiGainAnt2iso,truthAnt2iso,radTotalGainAnt2iso,p2] = sph2rectRadPattern(modelSpecs,1);
        EthetaAnt2iso = [thetaGainAnt2iso phiGainAnt2iso];
        %Takes away misnomer phase errors on zero values
        [EthetaAnt2iso] = zeroTheZeros(EthetaAnt2iso,threshold);
    %Find Self Impedance
        Iiso2 = modelSpecs.current; Viso2 = 1;
        Z22 = Viso2 / Iiso2;  
    %Antenna Length
        len2 = abs(modelSpecs.Z2(1) - modelSpecs.Z1(1));        
%Get Wave Parameters
fc = modelSpecs.FREQUENCY;
lambda = modelSpecs.WAVELENGTH;                  %wavelength [meters]
beta = 2*pi/lambda;

%Get Axis and Resolution   
res = abs( modelSpecs.radPattern(1,1) - modelSpecs.radPattern(2,1) );  
    thetaAxis = -180:res:0;
    phiAxis = 0:res:360;
    
%Acquire terminal current with both antennas turned on
name = 'arrayVant1ant2';
fileName = strcat(name,'.txt');
modelSpecs = readNecOutfile(fileName); 
    Io = modelSpecs.current; 
        
%Find Mutual Impedances
V12 = Viso1 - Z11*Io;
V21 = Viso2 - Z22*Io;
Z12 = V21/Io;
Z21 = V12/Io;        
        
%Compute Coupling Matrix
C1 = [1 Z12/Z22; Z21/Z11 1];
C = C1^-1;
        
        
%Compute Final Coupled Manifold
%Now Computations with Coupling Matrix and Isolated E Fields
thetaCounter = 180 / res + 1;
phiCounter = 360 / res + 1;
manifoldComputed = [];

EfieldTotalV = [];
EfieldAnt1isoV = []; EfieldAnt2isoV = [];
EfieldAnt1V = [];    EfieldAnt2V = [];
    for k = 1:length(phiAxis)
        for m = 1:length(thetaAxis)
            counter = (k-1)*thetaCounter + m;
            
            cur_Iso1Efield = EthetaAnt1iso(counter,:);
            cur_Iso2Efield = EthetaAnt2iso(counter,:);
            
            holder = C*[cur_Iso1Efield ; cur_Iso2Efield ];
            curEfieldTotalV = [sum(holder(:,1)) sum(holder(:,2))];
            EfieldAnt1V = [EfieldAnt1V; holder(1,1)];
            EfieldAnt2V = [EfieldAnt2V; holder(2,1)];

            EfieldTotalV = [EfieldTotalV; curEfieldTotalV];
        end
    end    
EfieldHorizontal = zeros(size(EfieldTotalV));
EfieldTotal = [EfieldTotalV EfieldHorizontal];
    EfieldAnt1 = [EfieldAnt1V EfieldHorizontal];
    EfieldAnt2 = [EfieldAnt2V EfieldHorizontal];


%Do the same to the analytic manifold as we did to the truth    
[EfieldTotal] = zeroTheZeros(EfieldTotal,threshold); 
[EfieldAnt1] = zeroTheZeros(EfieldAnt1,threshold);
[EfieldAnt2] = zeroTheZeros(EfieldAnt2,threshold);


%% Analytic Manifold

%Compute Analytic Coupling Matrix
    %Self Impedance
    n = 1; %since length of antenna is n* lambda/2
    Z11a = computeSI(n);
    Z22a = Z11a;
    
    %Mutual Impedance
    Z12a = computeMIparallel(len1,lambda/2,beta);
    Z21a = computeMIparallel(len2,lambda/2,beta);

    %The Matrix
    C1a = [1 Z12a/Z22a; Z21a/Z11a 1];
    Ca = C1a^-1;
    
%Compute the Manifold    
manifoldV1iso = []; manifoldV2iso = [];    
manifoldV1 = []; manifoldV2 = [];
manifoldTotal = [];
    for k = 1:length(phiAxis)
        for m = 1:length(thetaAxis)
            
            %Calculate Isolated Responses
            ShapeGain1 = exp(1j*2*pi/lambda*p1'*u(thetaAxis(m),phiAxis(k)));
                curManiV1 = ShapeGain1*vA(beta,len1,thetaAxis(m));
                manifoldV1iso = [manifoldV1iso; curManiV1];
                curManiH1 = 0;
                curMani1 = [curManiV1 curManiH1];
                
            ShapeGain2 = exp(1j*2*pi/lambda*p2'*u(thetaAxis(m),phiAxis(k)));
            curManiV2 = ShapeGain2*vA(beta,len2,thetaAxis(m));
                manifoldV2iso = [manifoldV2iso; curManiV2];
                curManiH2 = 0;
                curMani2 = [curManiV2 curManiH2];                   
                
            %Incorporate Coupling
            curMani = Ca*[curMani1 ; curMani2 ];
            curManifoldV = sum(holder);  wrong
                manifoldV1 = [manifoldV1; curMani(1)];
                manifoldV2 = [manifoldV2; curMani(2)]; 
                manifoldTotal = [manifoldTotal; curManifoldV];
            
        end
    end 
%Incorporate the Horizontal Manifold    
manifoldH1 = zeros(size(manifoldV1));   manifoldH2 = manifoldH1;  
aManifold1 = [manifoldV1 manifoldH1];
aManifold2 = [manifoldV2 manifoldH2];  
  
%Find the Complex Constant for the coupled Manifold
[aManifold1] = zeroTheZeros(aManifold1,threshold); 
    constC1f = constManifoldMin(EfieldAnt1(:,1),aManifold1(:,1));
    aManifold1Scaled = constC1f*aManifold1;
[aManifold2] = zeroTheZeros(aManifold2,threshold); 
    constC2f = constManifoldMin(EfieldAnt2(:,1),aManifold2(:,1));
    aManifold2Scaled = constC2f*aManifold2;

%% Analyze & Plot Results

%Antenna 1:  CA Truth vs. Analytical CA   
plotManifolds(EfieldAnt1(:,1),aManifold1Scaled(:,1),thetaAxis,phiAxis,' Antenna 1 Vertical')

%Antenna 2:  CA Truth vs. Analytical CA   
plotManifolds(EfieldAnt2(:,1),aManifold2Scaled(:,1),thetaAxis,phiAxis,' Antenna 2 Vertical')




