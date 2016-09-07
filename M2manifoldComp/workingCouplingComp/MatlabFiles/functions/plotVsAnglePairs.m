function plotVsAnglePairs(modelSpecs,values,format)
%Will Produce a stem plot with the values to be plotted as given, plus a
%custom x-axis showing all the angle pairs
    %format is a string such as format.line = 'b-';  format.linewidth = 3;

d2r = pi/180;
line = format.line;
lineWidth = format.linewidth;

%Extraction
[ThetaAll,PhiAll,thetaGain1,phiGain1,truth1,radTotalGain1,p_1] = sph2rectRadPattern(modelSpecs,1);
    %Resolution
        res = abs( modelSpecs.radPattern(1,1) - modelSpecs.radPattern(2,1) );
    %Create Viewer Friendly Error Plot - For Graph Annotation (xlabel)
    degTheta = round( ThetaAll/d2r );
    degPhi = round( PhiAll/d2r );
    degAll = [degTheta degPhi];
    strDegAll = [];
    for i = 1:length(degPhi)

        aa =  strcat ( num2str(degTheta(i)),',',num2str(degPhi(i)) );
        strDegAll = [strDegAll ' ' aa];
    end
   
    thetaCounter = 180 / res + 1;
    phiCounter = 360 / res + 1;

    strLen = 8;
    strDegAllFinal = [];
    others = strDegAll;
    for k = 1:phiCounter
        for m = 1:thetaCounter

            [val,others] = strtok(others); 
            curLen = length(val);
            if curLen < strLen
                filler = blanks(strLen - curLen);
               val = [val filler];
            end
            strDegAllFinal = [strDegAllFinal; val];

        end
    end  
    %Add space to fill in the zero spot where there should be nothing
    strDegAllFinal = [blanks(8) ; strDegAllFinal ];
 
%Now do the actual plotting    
stem(values,line,'linewidth',lineWidth)
axis( [0 length(degPhi)  min(values)+0.25*min(values) max(values)+0.25*max(values) ]);
set(gca,'xtick',[0:length(degPhi)-1],'xticklabel',...
    {strDegAllFinal});

        
        
        
        