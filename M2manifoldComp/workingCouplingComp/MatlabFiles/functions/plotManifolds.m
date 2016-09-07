function plotManifolds(manifoldTruth,manifoldAnalytic,thetaAxis,phiAxis,titleStr,numberOfElements,frequency,spacing,antNum)
%plots a variety of comparisons of gain and phase of the two manifolds,
%manifoldTruth and manifoldAnalytic with the string "title" within the
%titles of the figures.  The manifold vector is transformed to an (mxn)
%image as perscribed by the passed in axes
const = 10e-8;

%average = sum(sum(abs(abs(manifoldAnalytic)-abs(manifoldTruth))))/numel(manifoldAnalytic)
%display(abs(manifoldTruth));

specificLookTheta = 19;
specificLookPhi = 19;

lenTheta = length(thetaAxis);
lenPhi = length(phiAxis);

imageAnalyticManifold = reshape(manifoldAnalytic,lenTheta,lenPhi);
imageTruthManifold = reshape(manifoldTruth,lenTheta,lenPhi);

a=figure;
%subplot(321)
%    imagesc(phiAxis,thetaAxis,abs(imageTruthManifold))
%    xlabel('\phi Azimuth'),ylabel('\theta Elevation')
%    title(strcat('| Truth | for ',titleStr)),grid,colorbar
%subplot(322)
%   imagesc(phiAxis,thetaAxis,abs(imageAnalyticManifold))
%    xlabel('\phi Azimuth'),ylabel('\theta Elevation')
%    title(strcat('| Analytic | for ',titleStr)),grid,colorbar
%subplot(323)
%    imagesc(phiAxis,thetaAxis,angle(imageTruthManifold))
%    xlabel('\phi Azimuth'),ylabel('\theta Elevation')
%    title(strcat('< Truth for ',titleStr)),grid,colorbar
%subplot(324)
%    imagesc(phiAxis,thetaAxis,angle(imageAnalyticManifold))
%    xlabel('\phi Azimuth'),ylabel('\theta Elevation')
%    title(strcat('< Analytic for ',titleStr)),grid,colorbar    
subplot(121)
    truthAbs = abs(imageTruthManifold) + const;
    analyticAbs = abs(imageAnalyticManifold) + const;
    toPlot = truthAbs./analyticAbs;
    magError = toPlot;
    magSliceYValues = toPlot(19,:);
    specificMag = toPlot(specificLookTheta,specificLookPhi);
    imagesc(phiAxis,thetaAxis,(truthAbs-analyticAbs)./truthAbs);
    xlabel('\phi Azimuth'),ylabel('\theta Elevation')
    title(strcat('| Truth | / | Analytic | for ',titleStr)),grid,colorbar
subplot(122)
    truthAngle = abs(angle(imageTruthManifold))+ const;
    analyticAngle = abs(angle(imageAnalyticManifold))+ const;
    toPlot =  truthAngle - analyticAngle;
    phaseError = toPlot;
    phaseSliceYValues = toPlot(19,:);
    specificPhase = toPlot(specificLookTheta,specificLookPhi);
    imagesc(phiAxis,thetaAxis,(truthAngle-analyticAngle)./truthAngle);
    xlabel('\phi Azimuth'),ylabel('\theta Elevation')
    title(strcat('<Truth  -  <Analytic for ',titleStr)),grid,colorbar
    
%bslice=figure;
%  plot(phiAxis, magSliceYValues,'r');
%  xlabel('\phi Azimuth'), ylabel('Value'), grid on
%  hold on
%  plot(phiAxis, phaseSliceYValues,'b');
%  hold on
%  magSliceMean = mean(magSliceYValues);
%  %magSliceMeanArray = repmat(magSliceMean,73);
%  line([0,360],[magSliceMean, magSliceMean]);
%  hold on
  %phaseSliceMean = mean(phaseSliceYValues);
%  %phaseSliceMeanArray = repmat(phaseSliceMean,73);
%  line([0,360],[phaseSliceMean, phaseSliceMean]);
%  legend(strcat('|Truth|  -  |Analytic| for ',titleStr),strcat('<Truth  -  <Analytic for ',titleStr),strcat('Mean |Truth|  -  |Analytic| for ',titleStr),strcat('Mean <Truth  -  <Analytic for ',titleStr))

%fprintf('%s%d%s%f\n','Mean Magnitude percent error for ant ',antNum,': ',(sum(sum(abs(truthAbs-analyticAbs)./truthAbs)))/numel(truthAbs));
%fprintf('%s%d%s%f\n','Mean phase error for ant ',antNum,': ',(sum(sum(abs((truthAngle-analyticAngle)./truthAngle)))/numel(truthAbs)));
%fprintf('%s%d%s%f\n','Variance of magnitude in ant ',antNum,': ',var(reshape(magError,[numel(magError),1])));
%fprintf('%s%d%s%f\n','Variance of phase in ant ',antNum,': ',var(reshape(magError,[numel(phaseError),1])));
%fprintf('%s%d%s%f\n','Mean Magnitude percent error in slice of ant ',antNum,': ',(sum(sum(abs(truthAbs(19,:)-analyticAbs(19,:)))./truthAbs(19,:)))/numel(truthAbs));
%fprintf('%s%d%s%f\n','Mean Phase error in slice of ant ',antNum,': ',abs(phaseSliceMean));
%fprintf('%s%d%s%f\n','Variance of magnitude in slice of ant ',antNum,': ',var(magSliceYValues));
%fprintf('%s%d%s%f\n','Variance of phase in slice of ant ',antNum,': ',var(phaseSliceYValues));
%fprintf('%s%d%s%d%s%d%s%f\n','Error in Magnitude of Ant ',antNum,' Theta at ',thetaAxis(specificLookTheta),' and Phi at ',phiAxis(specificLookPhi),': ',specificMag);
%fprintf('%s%d%s%d%s%d%s%f\n','Error in Phase of ',antNum,' Theta at ',thetaAxis(specificLookTheta),' and Phi at ',phiAxis(specificLookPhi),': ',specificPhase);
 
 
 %uncomment and change path to save the graphs
 % pathh =strcat( 'C:\Users\atvtyt\Desktop\2016_SIP_Project_ELE-03\octave-4.0.2\M2manifoldComp\workingCouplingComp\NECuserMadeFiles\vDipoleArray\test\', num2str(numberOfElements),"ant_",num2str(frequency),"hz_",num2str(spacing),"wvl_",num2str(antNum),".png")
 %saveas(a, pathh,'png');
 %pathhh =strcat( 'C:\Users\atvtyt\Desktop\2016_SIP_Project_ELE-03\octave-4.0.2\M2manifoldComp\workingCouplingComp\NECuserMadeFiles\vDipoleArray\test\', num2str(numberOfElements),"ant_",num2str(frequency),"hz_",num2str(spacing),"wvl_",num2str(antNum),'slicegraph',".png")
 %saveas(bslice, pathhh,'png');
 
