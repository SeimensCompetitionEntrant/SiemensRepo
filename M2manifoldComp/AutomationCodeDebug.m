%AutomationCode.m
%Gautham Gujjula, Yash Tandon
%This m file automates the process of generating 4nec2 files,
%parsing through it to get the data we care about, and creating,
%a manifold with just the steering vector, coupling matrix, etc.

%Close windows, clear variables, clear command line
clear
close all
clc

[good, folderPath] = system('cd');
folderPath = strrep(folderPath,'/','\');
folderPath = strrep(folderPath,sprintf('\n'),'');

necPath = 'C:\4nec2';
simulFileName = '';
antFileName = '';
numAnts = 20;
pattern = 1;
freq = 430.0;
wvl = 0.0;
spacing = 0.5;
dist = 0.0;
c = 299792458.0;

%necPath = input(sprintf('Please enter the full path to the 4nec2 installation.\nPlease make sure that the installation of 4nec2 and\nthis folder are on the same drive: '), 's');
%necPath = strrep(necPath,sprintf('\n'),'')
%numAnts = input('Number of antennas: ');
%pattern = input('Pattern? 1 for ULA, 2 for UCA: ');
%freq = input('Carrier frequency, in MHz: ');
%spacing = input('Spacing between antennas in wavelengths: ');
wvl = c/(freq * 1000000.0);
dist = spacing*wvl;

tic;

  antFileName = strcat('ant',num2str(11));
  simulFileName = strcat(simulFileName,antFileName);
  system(strcat({'COPY NUL '},antFileName,'.nec'));
  fileID = fopen(strcat(antFileName,'.nec'),'wt');
  fprintf(fileID,'%s\r%s\r','CM','CE');
  for k = 1:numAnts
    fprintf(fileID, '%4s%d%6s%f%1s%f%3s%f%1s%f%6s\r', 'GW  ', k, '  9 0 ', (k-1)*dist, ' ', -wvl/4, ' 0 ', (k-1)*dist, ' ', wvl/4, ' .0001');
  end
  fprintf(fileID, '%5s\r', 'GE  0');
  for k = 1:numAnts
    fprintf(fileID,'%6s%d%7s\r', 'LD  4 ', k, ' 5 5 50');
  end
  fprintf(fileID, '%5s\r', 'GN -1');
  fprintf(fileID, '%2s\r', 'EK');
  fprintf(fileID, '%6s%d%10s\r', 'EX  0 ', 11 , ' 5 0 1 0 0');
  fprintf(fileID, '%12s%3d%2s\r', 'FR  0 0 0 0 ', freq, ' 0');
  fprintf(fileID, '%51s\r', 'RP  0 37  73    1003 -180     0         5         5');
  fprintf(fileID, '%2s', 'EN');
  fclose(fileID);


  antFileName = strcat('ant', num2str(11));
  movefile(strcat(antFileName,'.nec'), strcat('.\workingCouplingComp\NECuserMadeFiles\vDipoleArray\test\', antFileName, '.nec'));
  system(strcat({'4nec2.exe '}, folderPath, '\workingCouplingComp\NECuserMadeFiles\vDipoleArray\test\', antFileName, '.nec -I'));
  if exist(strcat('.\workingCouplingComp\NECuserMadeFiles\vDipoleArray\test\', antFileName, '.txt'),'file') != 0
    delete(strcat('.\workingCouplingComp\NECuserMadeFiles\vDipoleArray\test\', antFileName, '.txt'));
  end
  movefile(strcat(necPath, '\out\', antFileName, '.out'),strcat('.\workingCouplingComp\NECuserMadeFiles\vDipoleArray\test\', antFileName, '.txt'));

timeElapsed = toc;
fprintf('%s%f','4nec2 Time: ',timeElapsed);
cd('.\workingCouplingComp\MatlabFiles');
M2vArrayFinalDebug(numAnts, dist, folderPath, 5, 5);
cd('..\..');