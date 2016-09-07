%% Compares manifolds in a 4 element ULA
    %Patrick Ellis
    %April 10th, 2015
    
    
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
    
%% Load NEC Manifold with both turned on
name = 'ant1';
fileName = strcat(name,'.txt');
modelSpecs = readNecOutfile(fileName);   
    %Acquire final Manifold
        [ThetaAll,PhiAll,thetaGainAllOn,phiGainAllOn,truthAllOn,radTotalGainAllOn,p_AllOn] = sph2rectRadPattern(modelSpecs,1);
        manifoldTruthAnt1 = [thetaGainAllOn phiGainAllOn];
        %Takes away misnomer phase errors on zero values
        threshold = 10e-10;
        [manifoldTruthAnt1] = zeroTheZeros(manifoldTruthAnt1,threshold); 
name = 'ant2';
fileName = strcat(name,'.txt');
modelSpecs = readNecOutfile(fileName);   
    %Acquire final Manifold
        [ThetaAll,PhiAll,thetaGainAllOn,phiGainAllOn,truthAllOn,radTotalGainAllOn,p_AllOn] = sph2rectRadPattern(modelSpecs,1);
        manifoldTruthAnt2 = [thetaGainAllOn phiGainAllOn];
        %Takes away misnomer phase errors on zero values
        [manifoldTruthAnt2] = zeroTheZeros(manifoldTruthAnt2,threshold); 
name = 'ant3';
fileName = strcat(name,'.txt');
modelSpecs = readNecOutfile(fileName);   
    %Acquire final Manifold
        [ThetaAll,PhiAll,thetaGainAllOn,phiGainAllOn,truthAllOn,radTotalGainAllOn,p_AllOn] = sph2rectRadPattern(modelSpecs,1);
        manifoldTruthAnt3 = [thetaGainAllOn phiGainAllOn];
        %Takes away misnomer phase errors on zero values
        [manifoldTruthAnt3] = zeroTheZeros(manifoldTruthAnt3,threshold);         
name = 'ant4';
fileName = strcat(name,'.txt');
modelSpecs = readNecOutfile(fileName);   
    %Acquire final Manifold
        [ThetaAll,PhiAll,thetaGainAllOn,phiGainAllOn,truthAllOn,radTotalGainAllOn,p_AllOn] = sph2rectRadPattern(modelSpecs,1);
        manifoldTruthAnt4 = [thetaGainAllOn phiGainAllOn];
        %Takes away misnomer phase errors on zero values
        [manifoldTruthAnt4] = zeroTheZeros(manifoldTruthAnt4,threshold);         
        
        
        
        
%Find Wave Parameters
	lambda = modelSpecs.WAVELENGTH;                  %wavelength [meters]
	beta = 2*pi/lambda;
%Find Resolution of Recording 
res = abs( modelSpecs.radPattern(1,1) - modelSpecs.radPattern(2,1) );  
%Define Axis
    thetaAxis = -180:res:0;
    phiAxis = 0:res:360;
    
%See Orientation and get Antenna Length
p_all = plotArrayConfig(modelSpecs);    
    axis([-2 42 -4 4 -0.55*modelSpecs.WAVELENGTH 0.55*modelSpecs.WAVELENGTH])
len = abs(modelSpecs.Z2(1) - modelSpecs.Z1(1));
%% Analytic Manifold

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

manifoldV1 = []; manifoldV2 = [];
    for k = 1:length(phiAxis)
        for m = 1:length(thetaAxis)
            
            %Calculate Isolated Responses
            %Ant1
            ShapeGain1 = exp(1j*2*pi/lambda*p_all(1,:)*u(thetaAxis(m),phiAxis(k)));
                curManiV1 = ShapeGain1*vA(beta,len,thetaAxis(m));
                curManiH1 = 0;
                curMani1 = [curManiV1 curManiH1];
            %Ant2    
            ShapeGain2 = exp(1j*2*pi/lambda*p_all(2,:)*u(thetaAxis(m),phiAxis(k)));
            curManiV2 = ShapeGain2*vA(beta,len,thetaAxis(m));
                curManiH2 = 0;
                curMani2 = [curManiV2 curManiH2];  
            %Ant3    
            ShapeGain3 = exp(1j*2*pi/lambda*p_all(3,:)*u(thetaAxis(m),phiAxis(k)));
            curManiV3 = ShapeGain3*vA(beta,len,thetaAxis(m));
                curManiH3 = 0;
                curMani3 = [curManiV3 curManiH3];   
            %Ant4
            ShapeGain4 = exp(1j*2*pi/lambda*p_all(4,:)*u(thetaAxis(m),phiAxis(k)));
            curManiV4 = ShapeGain4*vA(beta,len,thetaAxis(m));
                curManiH2 = 0;
                curMani4 = [curManiV4 curManiH4];    
                
%             %Incorporate Coupling
%             curMani = [curMani1 ; curMani2 ];
%                 manifoldV1 = [manifoldV1; curMani(1,1)];
%                 manifoldV2 = [manifoldV2; curMani(2,1)];
        end
    end 


% manifoldH1 = zeros(size(manifoldV1));   manifoldH2 = manifoldH1;  
% analyticArrayManifold1 = [manifoldV1 manifoldH1];
%     [analyticArrayManifold1] = zeroTheZeros(analyticArrayManifold1,threshold); 
% analyticArrayManifold2 = [manifoldV2 manifoldH2];
%     [analyticArrayManifold2] = zeroTheZeros(analyticArrayManifold2,threshold);
%     
% A = [analyticArrayManifold1(:,1); analyticArrayManifold2(:,1)];
% B = [manifoldTruthAnt1(:,1); manifoldTruthAnt2(:,1)];
% constC = constManifoldMin(B,A);
%     manifoldAnalyticAnt1Scaled = constC*analyticArrayManifold1;
%     [manifoldAnalyticAnt1Scaled] = zeroTheZeros(manifoldAnalyticAnt1Scaled,threshold);
%     manifoldAnalyticAnt2Scaled = constC*analyticArrayManifold2;
%     [manifoldAnalyticAnt2Scaled] = zeroTheZeros(manifoldAnalyticAnt2Scaled,threshold); 


%% Analyze & Plot Results
% const = 10e-8; %to avoid dividing by zero
% 
% %Change shape of manifolds
% lenTheta = length(thetaAxis);
% lenPhi = length(phiAxis);
% 
% %Check the response of Antenna 1 vs. Truth
% imageAnalyticManifold = reshape(manifoldAnalyticAnt1Scaled(:,1),lenTheta,lenPhi);
% imageTruthManifold = reshape(manifoldTruthAnt1(:,1),lenTheta,lenPhi);
% 
% figure;
% subplot(321)
%     imagesc(phiAxis,thetaAxis,abs(imageTruthManifold))
%     xlabel('\phi Azimuth'),ylabel('\theta Elevation')
%     title('| Truth |'),grid,colorbar
% subplot(322)
%     imagesc(phiAxis,thetaAxis,abs(imageAnalyticManifold))
%     xlabel('\phi Azimuth'),ylabel('\theta Elevation')
%     title('| Analytic |'),grid,colorbar
% subplot(323)
%     imagesc(phiAxis,thetaAxis,angle(imageTruthManifold))
%     xlabel('\phi Azimuth'),ylabel('\theta Elevation')
%     title('< Truth '),grid,colorbar
% subplot(324)
%     imagesc(phiAxis,thetaAxis,angle(imageAnalyticManifold))
%     xlabel('\phi Azimuth'),ylabel('\theta Elevation')
%     title('< Analytic '),grid,colorbar    
% subplot(325)
%     truthAbs = abs(imageTruthManifold) + const;
%     analyticAbs = abs(imageAnalyticManifold) + const;
%     toPlot = truthAbs./analyticAbs;
%     imagesc(phiAxis,thetaAxis(1:end),toPlot(1:end,:))
%     xlabel('\phi Azimuth'),ylabel('\theta Elevation')
%     title('| Truth | / | Analytic |'),grid,colorbar
% subplot(326)
%     truthAngle = abs(angle(imageTruthManifold));
%     analyticAngle = abs(angle(imageAnalyticManifold));
%     toPlot =  truthAngle - analyticAngle ;
%     imagesc(phiAxis,thetaAxis(1:end),toPlot(1:end,:))
%     xlabel('\phi Azimuth'),ylabel('\theta Elevation')
%     title('<Truth  /  <Analytic '),grid,colorbar
% 
% %Check the response of Antenna 2 vs. Truth
% imageAnalyticManifold = reshape(manifoldAnalyticAnt2Scaled(:,1),lenTheta,lenPhi);
% imageTruthManifold = reshape(manifoldTruthAnt2(:,1),lenTheta,lenPhi);
% 
% figure;
% subplot(321)
%     imagesc(phiAxis,thetaAxis,abs(imageTruthManifold))
%     xlabel('\phi Azimuth'),ylabel('\theta Elevation')
%     title('| Truth |'),grid,colorbar
% subplot(322)
%     imagesc(phiAxis,thetaAxis,abs(imageAnalyticManifold))
%     xlabel('\phi Azimuth'),ylabel('\theta Elevation')
%     title('| Analytic |'),grid,colorbar
% subplot(323)
%     imagesc(phiAxis,thetaAxis,angle(imageTruthManifold))
%     xlabel('\phi Azimuth'),ylabel('\theta Elevation')
%     title('< Truth '),grid,colorbar
% subplot(324)
%     imagesc(phiAxis,thetaAxis,angle(imageAnalyticManifold))
%     xlabel('\phi Azimuth'),ylabel('\theta Elevation')
%     title('< Analytic '),grid,colorbar    
% subplot(325)
%     truthAbs = abs(imageTruthManifold) + const;
%     analyticAbs = abs(imageAnalyticManifold) + const;
%     toPlot = truthAbs./analyticAbs;
%     imagesc(phiAxis,thetaAxis(1:end),toPlot(1:end,:))
%     xlabel('\phi Azimuth'),ylabel('\theta Elevation')
%     title('| Truth | / | Analytic |'),grid,colorbar
% subplot(326)
%     truthAngle = abs(angle(imageTruthManifold));
%     analyticAngle = abs(angle(imageAnalyticManifold));
%     toPlot =  truthAngle - analyticAngle ;
%     imagesc(phiAxis,thetaAxis(1:end),toPlot(1:end,:))
%     xlabel('\phi Azimuth'),ylabel('\theta Elevation')
%     title('<Truth  /  <Analytic '),grid,colorbar
%     
    
    
    