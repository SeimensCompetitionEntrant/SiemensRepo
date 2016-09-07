%% A comparison of antenna responses for a 2 dipole array
function returnNum = M2vArrayFinal(numElements, separation, folderPath, thetaRes, phiRes)

%numElements = 16;
%separation = 0.5*(299792458/430000000);
%folderPath = 'D:\10th Grade\SIP\M2manifoldComp';
%thetaRes = 5;
%phiRes = 5;

pkg load specfun;   
addpath(genpath('D:\10th Grade\SIP\M2manifoldComp\workingCouplingComp'));

  
%Constants
    d2r = pi/180;   %degrees to radians
    c = 299792458;  %speed of light [meters per second]
    eta = 3;        %intrinsic impedance of medium = sqrt(mu/eps) (Medium is free space)
%Functions     
    %unit cosine in direction of propogation
        u = @(theta_k, phi_k) [cos(pi/180*phi_k)*sin(pi/180*theta_k); sin(pi/180*phi_k)*sin(pi/180*theta_k); cos(pi/180*theta_k)];       
    %Vertical Dipole Manifold
        vA = @(beta,len,theta)  (cos(beta*len*cos(theta*pi/180)/2) - cos(beta*len/2)) / sin(theta*pi/180+0.000001);

tic;

%% Load the NEC generated responses
%Takes away misnomer phase errors on zero values
threshold = 10e-10;  

antNames = zeros(numElements,4+floor(log10(numElements)));
simulAntName = '';
numGainValues = ((360/thetaRes)+1)*((180/phiRes)+1);
thetaGainAnt = zeros(numGainValues,1);
phiGainAnt = zeros(numGainValues,1);
EthetaAnt = zeros(numGainValues, 2*numElements);
radTotalGainAnt = zeros(numGainValues,1);
p = zeros(3,numElements);
truthAnt = zeros(numGainValues, 3*numElements);
for i = 1:numElements
  antNames(i,1:4+floor(log10(i))) = strcat('ant', num2str(i));
  simulAntName = strcat(simulAntName, antNames(i,1:4+floor(log10(i))));
end
%Load Active Element Manifolds
%Antenna i Transmitting
for i = 1:numElements
  name = antNames(i,1:4+floor(log10(i)));
  fileName = strcat(name,'.txt');
  modelSpecs = readNecOutfile(fileName);
      %Acquire final Isolated Manifold
          [ThetaAll,PhiAll,thetaGainAnt(:,i),phiGainAnt(:,i),truthAnt(:,(3*(i-1))+1:3*i),radTotalGainAnt(:,i),p(:,i)] = sph2rectRadPattern(modelSpecs,i);
          EthetaAnt(:,(2*(i-1))+1:2*i) = [thetaGainAnt(:,i) phiGainAnt(:,i)];
          %Takes away misnomer phase errors on zero values
          [EthetaAnt(:,(2*(i-1))+1:2*i)] = zeroTheZeros(EthetaAnt(:,(2*(i-1))+1:2*i),threshold); 
end
  %Antenna 2 Transmitting      
%name = 'ant2';
%fileName = strcat(name,'.txt');
%modelSpecs = readNecOutfile(fileName);
    %Acquire final Isolated Manifold
%        [ThetaAll,PhiAll,thetaGainAnt2,phiGainAnt2,truthAnt2,radTotalGainAnt2,p2] = sph2rectRadPattern(modelSpecs,2);
%        EthetaAnt2 = [thetaGainAnt2 phiGainAnt2];
        %Takes away misnomer phase errors on zero values
%        [EthetaAnt2] = zeroTheZeros(EthetaAnt2,threshold);

%Load all Elements Transmitting Manifold
name = simulAntName;
fileName = strcat(name,'.txt');
modelSpecs = readNecOutfile(fileName);
    %Acquire final Isolated Manifold
        [ThetaAll,PhiAll,thetaGainAnt1Ant2,phiGainAnt1Ant2,truthAnt1Ant2,radTotalGainAnt1Ant2,psimul] = sph2rectRadPattern(modelSpecs,1);
        EthetaAnt1Ant2 = [thetaGainAnt1Ant2 phiGainAnt1Ant2];
        %Takes away misnomer phase errors on zero values
        [EthetaAnt1Ant2] = zeroTheZeros(EthetaAnt1Ant2,threshold);     
%% Analytic Manifold

timerVal = toc;
fprintf('%s%f\n','Parsing Time: ', timerVal);
tic

%Load some needed Constants
    %Antenna Length
    len = abs(modelSpecs.Z2(1) - modelSpecs.Z1(1));
    %Get Wave Parameters
    fc = modelSpecs.FREQUENCY;
    lambda = modelSpecs.WAVELENGTH/100;     %HARCODED - MUST FIX EXTRACTION TO CORRECTLY INTERPRET NEGATIVE EXPONENTS!!!!!!!!
    beta = 2*pi/lambda;
    %Get Axis and Resolution   
    res = abs(modelSpecs.radPattern(1,1) - modelSpecs.radPattern(2,1));  
    thetaAxis = -180:res:0;
    phiAxis = 0:res:360;   
    
%Compute Analytic Coupling Matrix
    %The Load
    ZL = 50;
    
    %Self Impedance
    n = 1; %since length of antenna is n* lambda/2
    SI = computeSIarray(n,numElements);
    
    %Mutual Impedance
    MI = computeMIarray(len, separation, beta, numElements);

    %The Matrix
    C1a = computeCouplingMatrix(SI, MI, ZL, numElements);
    Ca = C1a^-1;
    
%Compute the Manifold    
manifoldViso = zeros(numGainValues,numElements);
manifoldV = zeros(numGainValues,numElements);
ShapeGain = zeros(1,numElements);
curManiV = zeros(1,numElements);
curManiH = zeros(1,numElements);
curManiIndiv = zeros(numElements,2);
    for k = 1:length(phiAxis)
        for m = 1:length(thetaAxis)
            for i = 1:numElements
              %Calculate Isolated Responses
              ShapeGain(i) = exp(1j*2*pi/lambda*p(:,i)'*u(thetaAxis(m),phiAxis(k)));
                 curManiV(i) = ShapeGain(i)*vA(beta,len,thetaAxis(m));
                 if k == 1 && m == 1
                  manifoldViso(1,i) = [curManiV(i)];
                 else
                  manifoldViso(1:((k-1)*(length(thetaAxis)))+m,i) = [manifoldViso(1:((k-1)*(length(thetaAxis)))+m-1,i); curManiV(i)];
                 end
                 curManiH(i) = 0;
                 curManiIndiv(i,:) = [curManiV(i) curManiH(i)];
                
          %  ShapeGain2 = exp(1j*2*pi/lambda*p2'*u(thetaAxis(m),phiAxis(k)));
          %  curManiV2 = ShapeGain2*vA(beta,len,thetaAxis(m));
          %      manifoldV2iso = [manifoldV2iso; curManiV2];
          %      curManiH2 = 0;
          %      curMani2 = [curManiV2 curManiH2];
            end    
            %Incorporate Coupling
            curMani = Ca*curManiIndiv;
            curManifold = sum(curMani);
                for i = 1:numElements
                  if k == 1 && m == 1
                    manifoldV(1,i) = [curMani(i,1)];
                  else
                    manifoldV(1:((k-1)*(length(thetaAxis)))+m,i) = [manifoldV(1:((k-1)*(length(thetaAxis)))+m-1,i); curMani(i,1)];
                  end
                end
            
        end
    end 
    
%Incorporate the Horizontal Manifold    
manifoldH = zeros(numGainValues,numElements);
aManifold = zeros(numGainValues, 2*numElements);
sumMatrix = zeros(numGainValues,2);
constCf = zeros(1, numElements);
for i = 1:numElements
  aManifold(:,(2*(i-1))+1:2*i) = [manifoldV(:,i) manifoldH(:,i)];
  sumMatrix = sumMatrix + aManifold(:,(2*(i-1))+1:2*i);
  
%Find the Complex Constant Betwenn Manifold and NEC Output (Constants
%should be nearly identical)
[aManifold(:,(2*(i-1))+1:2*i)] = zeroTheZeros(aManifold(:,(2*(i-1))+1:2*i),threshold); 
    constCf(i) = constManifoldMin(EthetaAnt(:,(2*(i-1))+1),aManifold(:,(2*(i-1))+1));
end
sumConstCf = constManifoldMin(EthetaAnt1Ant2(:,1), sumMatrix(:,1));


timerVal = toc;
fprintf('%s%f\n','Calculation Time: ', timerVal);

%% Compare          

%Theory Check
% %NEC All On vs. NEC generated Sum of Active Elements   
% ActiveSum = EthetaAnt1 + EthetaAnt2;
% plotManifolds(EthetaAnt1Ant2(:,1),ActiveSum(:,1),thetaAxis,phiAxis,' Total vs. Sum of Active Elements (Truths)')
% 
% %NEC All On vs. Analytic Sum of Active Elements   
% ActiveSum = aManifold1*constC1f + aManifold2*constC2f;
% plotManifolds(EthetaAnt1Ant2(:,1),ActiveSum(:,1),thetaAxis,phiAxis,' Sum of Active Elements')
%display(EthetaAnt);
%display(abs(EthetaAnt));
%printf("%s%f",'Truth for 5000 percent error look direction: ',EthetaAnt(:))
%printf("%s%f",'Analytic for 5000 percent error look direction: ',EthetaAnt(:))
for i = 1:numElements
  plotManifolds(EthetaAnt(:,(2*(i-1)+1))/constCf(i),aManifold(:,(2*(i-1)+1)),thetaAxis,phiAxis,strcat(' Active Element Ant ',num2str(i)),numElements,modelSpecs.FREQUENCY,separation/lambda,i)
end
%plotManifolds(EthetaAnt1Ant2(:,1)/sumConstCf,sumMatrix(:,1)./numElements,thetaAxis,phiAxis,strcat(' Active Element Ant ',num2str(0)),numElements,modelSpecs.FREQUENCY,separation/lambda,0)


%counter=1;
%Manifold Check
%Antenna 1:  NEC active Element vs. Analytic active element
%plotManifolds(EthetaAnt1(:,1)/constC1f,aManifold1(:,1),thetaAxis,phiAxis,' Active Element Ant 1',numElements,modelSpecs.FREQUENCY,separation/lambda,counter)
%counter=2;
%Antenna 2:  NEC active Element vs. Analytic active element
%plotManifolds(EthetaAnt2(:,1)/constC2f,aManifold2(:,1),thetaAxis,phiAxis,' Active Element Ant 2',numElements,modelSpecs.FREQUENCY,separation/lambda,counter)
returnNum = 0;