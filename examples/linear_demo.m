% Copyright (C) 2021 Josiah W. Smith
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

%% Include Necessary Directories
addpath(genpath("../"))

%% Create the Objects
fmcw = fmcwChirpParameters();
ant = sarAntennaArray(fmcw);
sar = sarScenario(ant);
target = sarTarget(fmcw,ant,sar);
im = sarImage(fmcw,ant,sar,target);

%% Set FMCW Parameters
% When the parameters of an fmcwChirpParameters object are changed by the
% user, the object automatically updates itself, namely the property 'k'
% and other dependencies of the changed parameters.
fmcw.ADCSamples = 50;
fmcw.f0 = 300*1e9;
fmcw.RampEndTime_s = 50.1*1e-6;
fmcw.fS = 1000*1e3;
fmcw.K = 200*1e12;
fmcw.IdleTime_s = 0*1e-6;
fmcw.TXStartTime_s = 0*1e-6;
fmcw.ADCStartTime_s = 0*1e-6;
fmcw.fC = 305*1e9;

%% Set Antenna Array Properties
% When the parameters of an antennaArray object are changed by the user, 
% the object automatically updates itself
ant.isEPC = false;
ant.z0_m = 0;
% Small MIMO Array
ant.tableTx = [
    0   0   2   0   1
    0   0   4   0   1];
ant.tableRx = [
    0   0   0   0   1
    0   0   0.5 0   1
    0   0   1   0   1
    0   0   1.5 0   1];

% Display the Antenna Array
ant.displayAntennaArray();

%% Set SAR Scenario Parameters
% When the parameters of a sarScenario object are changed by the user, the
% object automatically updates itself
sar.scanMethod = 'Linear';
sar.yStep_m = fmcw.lambda_m*2;
sar.numY = 25;

% Display the SAR Scenario
sar.displaySarScenario();

%% Set Target Parameters
% When the parameters of a sarTarget object are changed by the user, the
% object automatically updates itself
target.isAmplitudeFactor = false;

target.tableTarget = [
    0   0       0.1    1
    0   0.01    0.1    1];

target.rp.numTargets = 16;
target.rp.xMin_m = 0;
target.rp.xMax_m = 0;
target.rp.yMin_m = -0.05;
target.rp.yMax_m = 0.05;
target.rp.zMin_m = 0.05;
target.rp.zMax_m = 0.2;
target.rp.ampMin = 0.5;
target.rp.ampMax = 1;

% Which to use
target.isTable = true;
target.isRandomPoints = false;

% Display the target
target.displayTarget();

%% Compute Beat Signal
target.isGPU = true;
target.computeTarget();

%% Set Image Reconstruction Parameters and Create sarImage Object
% When the parameters of a sarImage object are changed by the user, the
% object automatically updates itself
im.nFFTy = 512;
im.nFFTz = 512;

im.yMin_m = -0.05;
im.yMax_m = 0.05;

im.zMin_m = 0;
im.zMax_m = 0.2;

im.numY = 512;
im.numZ = 128;

im.isGPU = false;
im.zSlice_m = 0.1; % Use if reconstructing a 1-D image
im.method = "Uniform 1-D SAR 1-D FFT";

im.isMult2Mono = true;
im.zRef_m = 0.1;

%% Reconstruct the Image
im.computeImage();
im.displayImage();

%% Display the Image with Different Parameters
im.dBMin = -30;
im.fontSize = 12;
im.displayImage();