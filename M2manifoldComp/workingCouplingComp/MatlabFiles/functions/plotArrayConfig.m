function p = plotArrayConfig(modelSpecs)
%Plots the current antenna configuration (must be wire antenna, dipole or
%monopole)
    %only works with vertical dipoles at the moment

figure;
flagX = 0; flagY = 0; flagZ = 0;
radius = 0.025; %of antenna

%Get number of antennas
M = modelSpecs.Num(end);

%Get Wavelength
lambda = modelSpecs.WAVELENGTH;

%Get positions of antennas
p = [];
for m = 1:M
    %Extract Antenna Coordinates
    x1 = modelSpecs.X1(m)/lambda;
    x2 = modelSpecs.X2(m)/lambda;
    y1 = modelSpecs.Y1(m)/lambda;
    y2 = modelSpecs.Y2(m)/lambda;
    z1 = modelSpecs.Z1(m)/lambda;
    z2 = modelSpecs.Z2(m)/lambda;

    xAve = (x1 + x2)/2;
    yAve = (y1 + y2)/2;
    zAve = (z1 + z2)/2;
    
    p(m,:) = [xAve yAve zAve]';
    
    %Find orientation of Antenna and Length
    if(x1 ~= x2)        %Dipole laying on x axis
        len = abs(modelSpecs.X2(1) - modelSpecs.X1(1))/lambda;
        flagX = 1;
    elseif(y1 ~= y2)    %Dipole laying on y axis
        len = abs(modelSpecs.Y2(1) - modelSpecs.Y1(1))/lambda;
        flagY = 1;
    elseif(z1 ~= y2)    %Dipole laying on z axis
        len = abs(modelSpecs.Z2(1) - modelSpecs.Z1(1))/lambda;
        flagZ = 1;
    else
        disp('error with positioning')
    end
    
    %Plot it
    hold on
    [Xcur,Ycur,Zcur] = cylinder(radius);
    pts = length(Zcur);
        Zcur(1,:) = ones(1,pts)*-len/2;
        Zcur(2,:) = ones(1,pts)*len/2;
    surf(Xcur+p(m,1),Ycur+p(m,2),Zcur+p(m,3))
    
    %Add Annotation
    text(xAve,yAve,zAve,['      Antenna #' ...
        num2str(m)],'HorizontalAlignment','left','FontSize',15);
    
end


grid
title('Array Configuration')
xlabel('X Axis (Wavelengths)'),ylabel('Y Axis (Wavelengths)'),zlabel('Z Axis (Wavelengths)')
view(3)


p = p*lambda;



