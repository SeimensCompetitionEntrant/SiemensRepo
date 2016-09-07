function modelSpecs = readNecOutfile(fileName)

warning('off');
fileID = fopen(fileName);                   % open function
%Find Number of Lines in the File
numLines = 0;
tline = fgetl(fileID);
while ischar(tline)
    tline = fgetl(fileID);
    numLines = numLines+1;
end
fclose(fileID);

%Initializations
antennaSpecs_lineStart = [];    check1 = 1;
antennaSpecs_lineEnd = [];      check2 = 1;
freqANDwave_lineStart = [];     check3 = 1;
CurANDVol_lineStart = [];        check35 = 1;
CurAndVol_lineEnd = [];         check355 = 1;
radPattern_lineStart = [];      check4 = 1;
radPattern_lineEnd = [];        check5 = 1;
    radPatternArray = [];       counter1 = 0;
    radPatternAdder = 0;        check6 = 1;     holder6 = 0;

fileID = fopen(fileName);
%Loop Through all Lines and find Relevant starting & Ending point Lines
for curLine = 1:numLines
    tline = fgetl(fileID);                     % read a line
    if ischar(tline)                        % if the line is string

        %Find Start Point of Antenna Physical Characteristics
        if check1 == 1
            aa1 = strfind(tline, 'STRUCTURE SPECIFICATION'); % where the string start (if at all)
            if isfinite(aa1) == 1;                % if it is a number actually
               antennaSpecs_lineStart = curLine;     % found
               check1 = 0;
            end
        end

        %Find Endpoint of Antenna Physical Characteristics
        if check2 == 1
            aa2 = strfind(tline, 'MULTIPLE WIRE JUNCTIONS'); % where the string start (if at all)
            if isfinite(aa2) == 1;                % if it is a number actually
               antennaSpecs_lineEnd = curLine;     % found
               check2 = 0;
            end
        end

        %Find Start Point of Antenna Physical Characteristics
        if check3 == 1
            aa3 = strfind(tline, 'FREQUENCY'); % where the string start (if at all)
            if isfinite(aa3) == 1;                % if it is a number actually
               freqANDwave_lineStart = curLine;     % found
               check3 = 0;
            end
        end

        %Find Start point of ANTENNA INPUT PARAMETERS
        if check35 == 1
            aa35 = strfind(tline, 'ANTENNA INPUT PARAMETERS'); % where the string start (if at all)
            if isfinite(aa35) == 1;                % if it is a number actually
               CurANDVol_lineStart = curLine;     % found
               check35 = 0;
            end
        end

        %Find End point of ANTENNA INPUT PARAMETERS
        if check355 == 1
            aa355 = strfind(tline, 'CURRENTS AND LOCATION');
            if isfinite(aa355) == 1;
                CurAndVol_lineEnd = curLine;
                check355 = 0;
            end
        end

        %Find Start Point of Antenna Radiation Patterns
        if check4 == 1
            aa4 = strfind(tline, 'RADIATION PATTERNS'); % where the string start (if at all)
            if isfinite(aa4) == 1;                % if it is a number actually
                radPattern_lineStart = curLine;     % found
                holder6 = curLine + 1;
               check4 = 0;
            end
        end

        %Check to see if Distance to record Efield at is stipulated
        if check6 == 1
            %Check to see if NEC user has stipulated a distance to
            %extract EField at
            if curLine == holder6 && check4 == 0
                [a,rem1] = strtok(tline);
                tf = strcmp(a,'');
                if(tf == true(1))
                    radPatternAdder = 5;    %no stipulation

                else
                    radPatternAdder = 8;    %Stipulated Distance
                end
                check6 = 0;
            end
        end


        %Find End Point of Antenna Radiation Patterns
        if check5 == 1
            aa5 = strfind(tline, 'AVERAGE POWER GAIN'); % where the string start (if at all)
            if isfinite(aa5) == 1;                % if it is a number actually
               radPattern_lineEnd = curLine;     % found
               check5 = 0;
            end
        end


    end
end
fclose(fileID);

%Now Extract the Information Using the Found Hotpoints
fileID = fopen(fileName);

    %inits
    modelSpecs = [];
    allFieldNames = [];

%Loop Through all Lines and find Relevant starting & Ending point Lines
for curLine = 1:numLines
    tline = fgetl(fileID);                     % read a line

    %Antenna Physical Characteristics
        %Field Names
        if curLine == antennaSpecs_lineStart + 8
            remain = tline;
            for vals = 1:9
                [token,remain] = strtok(remain);
                %Take care of some issues with extracting names with periods
                %and naming fields
                if(vals == 1)
                    token = 'Num';
                end
                if(vals == 9)
                    token = 'SegNum';
                end
                modelSpecs.(token) = [];
                allFieldNames = [allFieldNames ' ' token];
            end
        end

        %Number of Antennas
        NumAntennas = (antennaSpecs_lineEnd - 5) - (antennaSpecs_lineStart + 8);
            flagger = NumAntennas;
        %Field Values
        remain = tline; remain2 = allFieldNames;
        if curLine >= antennaSpecs_lineStart + 9 && curLine <= (antennaSpecs_lineEnd - 5)
            if flagger > 0
                for vals = 1:9
                    %Get Value
                    [token,remain] = strtok(remain);
                    %Get FieldName to Store the Value
                    [token2,remain2] = strtok(remain2);

                    modelSpecs.(token2) = [modelSpecs.(token2) str2double(token)];
                end
                flagger = flagger - 1;
            end
        end

    %Wave Characteristics
        %Frequency
        if curLine == freqANDwave_lineStart + 2
            [token3,remain3] = strtok(tline);
            pos1 = strfind(remain3, 'E');
                decimal = remain3(1:pos1-1);
            [DC,DNC] = strtok(remain3);
                tensPwr = DC(pos1+1:end);

            %field name
            modelSpecs.(token3(1:end-1)) = str2double(decimal)*10^(str2double(tensPwr)) * 10^6;
        end
        %Frequency
        if curLine == freqANDwave_lineStart + 3
            [token4,remain4] = strtok(tline);
            pos1 = strfind(remain4, 'E');
                decimal = remain4(1:pos1-1);
            [DC,DNC] = strtok(remain4);
                tensPwr = DC(pos1+1:end);

            %field name
            modelSpecs.(token4(1:end-1)) = str2double(decimal)*10^(str2double(tensPwr));
        end

    %Antenna Current and Voltages
        if curLine >= CurANDVol_lineStart + 4 && curLine <= CurAndVol_lineEnd - 4
            [a,rem1] = strtok(tline);
            [b,rem2] = strtok(rem1);
            [c,rem3] = strtok(rem2);
            [d,rem4] = strtok(rem3);
            [e,rem5] = strtok(rem4);

            tester = strfind(e, 'E');
            if length(tester) > 1       %sometimes the output file concatentates the real and imaginary parts
                ptOne = tester(1)+3;
                realPart = e(1:ptOne);
                    %Do Extraction and Conversion
                    pos1 = strfind(realPart, 'E');
                    decimal = realPart(1:pos1-1);
                    tensPwr = realPart(pos1+1:end);
                    curReal = str2double(decimal)*10^(str2double(tensPwr));
                imagPart = e(ptOne+1:end);
                    %Do Extraction and Conversion
                    pos1 = strfind(imagPart, 'E');
                    decimal = imagPart(1:pos1-1);
                    tensPwr = imagPart(pos1+1:end);
                    curImag = str2double(decimal)*10^(str2double(tensPwr));
            else
                [f,rem6] = strtok(rem5);
                    %Do Extraction and Conversion
                    pos1 = strfind(e, 'E');
                    decimal = e(1:pos1-1);
                    tensPwr = e(pos1+1:end);
                    curReal = str2double(decimal)*10^(str2double(tensPwr));
                    %Do Extraction and Conversion
                    pos1 = strfind(f, 'E');
                    decimal = f(1:pos1-1);
                    tensPwr = f(pos1+1:end);
                    curImag = str2double(decimal)*10^(str2double(tensPwr));
            end
            modelSpecs.current = curReal + 1j*curImag;
        end




        %Radiation Patterns
        if radPatternAdder == 8     %So file has an Efield distance that must be extracted
            if curLine == radPattern_lineStart + 1
                [a,rem1] = strtok(tline);
                [b,rem2] = strtok(rem1);
                pos1 = strfind(b, 'E');
                    decimal = b(1:pos1-1);
                    tensPwr = b(pos1+1:end);
                    FFdist = str2double(decimal)*10^(str2double(tensPwr));
                modelSpecs.FarFieldDist = FFdist;
            end
        else
            modelSpecs.FarFieldDist = 0;
        end


        if curLine >= radPattern_lineStart + radPatternAdder && curLine <= radPattern_lineEnd - 3  %plus 5 if radius not stipulated
            counter1 = counter1 + 1;
            [thetaCur,rem1] = strtok(tline);
            [psiCur,rem2] = strtok(rem1);
            [dnc,rem3] = strtok(rem2);
            [dnc,rem4] = strtok(rem3);
            [gainTot,rem5] = strtok(rem4);

            [nothing,rem6] = strtok(rem5);
            [nothing,rem7] = strtok(rem6);

            [thetaMagRem,rem8] = strtok(rem7);
            if(strcmp(thetaMagRem,'LINEAR') == 1)
                [thetaMagRem2,rem9] = strtok(rem8);
                    pos1 = strfind(thetaMagRem2, 'E');
                        decimal = thetaMagRem2(1:pos1-1);
                    [DC,DNC] = strtok(thetaMagRem2);
                        tensPwr = DC(pos1+1:end);
                    thetaMagCur = str2double(decimal)*10^(str2double(tensPwr));

                    [thetaPhaseCur,rem10] = strtok(rem9);
                    [phiMagRem2,rem10] = strtok(rem10);
                        pos1 = strfind(phiMagRem2, 'E');
                            decimal = phiMagRem2(1:pos1-1);
                        [DC,DNC] = strtok(phiMagRem2);
                            tensPwr = DC(pos1+1:end);
                        phiMagCur = str2double(decimal)*10^(str2double(tensPwr));
                    [phiPhaseCur,rem11] = strtok(rem10);
            elseif(strcmp(thetaMagRem,'LINEAR') == 0)
                pos1 = strfind(thetaMagRem, 'E');
                    decimal = thetaMagRem(1:pos1-1);
                [DC,DNC] = strtok(thetaMagRem);
                    tensPwr = DC(pos1+1:end);
                thetaMagCur = str2double(decimal)*10^(str2double(tensPwr));

                [thetaPhaseCur,rem9] = strtok(rem8);
                [phiMagRem,rem10] = strtok(rem9);
                    pos1 = strfind(phiMagRem, 'E');
                        decimal = phiMagRem(1:pos1-1);
                    [DC,DNC] = strtok(phiMagRem);
                        tensPwr = DC(pos1+1:end);
                    phiMagCur = str2double(decimal)*10^(str2double(tensPwr));
                [phiPhaseCur,rem11] = strtok(rem10);
            else
                ;
            end
            radPatternArray(counter1,1) = str2double(thetaCur);
            radPatternArray(counter1,2) = str2double(psiCur);
            radPatternArray(counter1,3) = str2double(gainTot);
            radPatternArray(counter1,4) = thetaMagCur*exp(1j*pi/180*str2double(thetaPhaseCur));
            radPatternArray(counter1,5) = phiMagCur*exp(1j*pi/180*str2double(phiPhaseCur));

        end

        if curLine == radPattern_lineEnd
            modelSpecs.radPattern = radPatternArray;
        end


end
fclose(fileID);
