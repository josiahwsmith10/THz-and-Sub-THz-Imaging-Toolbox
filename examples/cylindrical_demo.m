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
ant.z0_m = 0.1;
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
sar.scanMethod = 'Cylindrical';
sar.numY = 25;
sar.yStep_m = fmcw.lambda_m*2;
sar.numTheta = 1024;
sar.thetaMax_deg = 360;

% Display the SAR Scenario
sar.displaySarScenario();

%% Set Target Parameters
% When the parameters of a sarTarget object are changed by the user, the
% object automatically updates itself
target.isAmplitudeFactor = false;

target.tableTarget = [
    0   0       0.01    1
    0   0.05    -0.01   1];

target.png.fileName = 'circle.png';
target.png.xStep_m = 2e-4;
target.png.yStep_m = 2e-4;
target.png.xOffset_m = -0.0025;
target.png.yOffset_m = 0.005;
target.png.zOffset_m = 0;
target.png.reflectivity = 1;
target.png.downsampleFactor = 4;

target.stl.fileName = 'ar15.stl';
target.stl.zCrop_m = 0.11;
target.stl.xOffset_m = 0;
target.stl.yOffset_m = 0;
target.stl.zOffset_m = 0;
target.stl.reflectivity = 1;
target.stl.downsampleFactor = 40;

target.rp.numTargets = 16;
target.rp.xMin_m = -0.05;
target.rp.xMax_m = 0.05;
target.rp.yMin_m = -0.05;
target.rp.yMax_m = 0.05;
target.rp.zMin_m = 0.05;
target.rp.zMax_m = 0.2;
target.rp.ampMin = 0.5;
target.rp.ampMax = 1;

% Which to use
target.isTable = true;
target.isPNG = false;
target.isSTL = false;
target.isRandomPoints = false;

% Display the target
target.displayTarget();

%% Compute Beat Signal
target.isGPU = true;
target.computeTarget();

%% Set Image Reconstruction Parameters and Create sarImage Object
% When the parameters of a sarImage object are changed by the user, the
% object automatically updates itself
im.nFFTx = 4096;
im.nFFTy = 512;
im.nFFTz = 512;

im.xMin_m = -0.05;
im.xMax_m = 0.05;

im.yMin_m = -0.05;
im.yMax_m = 0.05;

im.zMin_m = -0.05;
im.zMax_m = 0.05;

im.numX = 128;
im.numY = 128;
im.numZ = 128;

im.isGPU = false;
im.thetaUpsampleFactor = 4;
im.method = "Uniform 2-D CSAR 3-D PFA";

im.isMult2Mono = true;
im.zRef_m = 0.1;

%% Reconstruct the Image
im.computeImage();
im.displayImage();

%% Display the Image with Different Parameters
im.dBMin = -25;
im.fontSize = 12;
im.displayImage();