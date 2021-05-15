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
fmcw.f0 = 77*1e9;
fmcw.K = 100.036*1e12;
fmcw.IdleTime_s = 0*1e-6;
fmcw.TXStartTime_s = 0*1e-6;
fmcw.ADCStartTime_s = 0*1e-6;
fmcw.ADCSamples = 79;
fmcw.fS = 2000*1e3;
fmcw.RampEndTime_s = 39.98*1e-6;
fmcw.fC = 79*1e9;

%% Set Antenna Array Properties
% When the parameters of an antennaArray object are changed by the user,
% the object automatically updates itself
ant.isEPC = false;
ant.z0_m = 0;
% Antenna array from xWR1x43
ant.tableTx = [
    0   0   1.5 5   1
    0.5 0   2.5 5   0
    0   0   3.5 5   1];
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
    0   0   0.25    1
    0   0.1 0.25    1];

target.rp.numTargets = 16;
target.rp.xMin_m = 0;
target.rp.xMax_m = 0;
target.rp.yMin_m = -0.1;
target.rp.yMax_m = 0.1;
target.rp.zMin_m = 0.1;
target.rp.zMax_m = 0.5;
target.rp.ampMin = 0.5;
target.rp.ampMax = 1;

% Which to use
target.isTable = false;
target.isRandomPoints = true;

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

im.yMin_m = -0.2;
im.yMax_m = 0.2;

im.zMin_m = 0;
im.zMax_m = 0.5;

im.numY = 256;
im.numZ = 256;

im.isGPU = true;
im.method = "Uniform 1-D SAR 2-D RMA";

im.isMult2Mono = true;
im.zRef_m = 0.25;

%% Reconstruct the Image
im.computeImage();
im.displayImage();

%% Compute the Ideal Image from the YZ Imaging Scene and Display
imIdeal2D = computeIdeal2D(target,im,"YZ",3e-3,3e-3);
figure
plotXYdB(imIdeal2D,im.y_m,im.z_m,-25,"y (m)","z (m)","Ideal Image",12)

%% Setup Dataset Scenario
numSamples = 32;
numRandomPointsMax = 64;
numRandomPointsMin = 4;

radarImages = single(zeros(im.numY,im.numZ,numSamples));
idealImages = single(zeros(im.numY,im.numZ,numSamples));

% Set objects to silent
target.isSilent = true;
im.isSilent = true;

%% Loop through Image Generation
for indSample = 1:numSamples
    target.rp.numTargets = randi([numRandomPointsMin,numRandomPointsMax]);
    target.computeTarget();
    
    im.computeImage();
    
    radarImages(:,:,indSample) = im.imXYZ;
    idealImages(:,:,indSample) = computeIdeal2D(target,im,"YZ",3e-3,3e-3);
    disp("Iteration " + indSample + "/" + numSamples + " Done");
end

%% Compare Samples Side by Side
indRand = randi(numSamples);

figure
subplot(121)
plotXYdB(radarImages(:,:,indRand),im.y_m,im.z_m,-25,"y (m)","z (m)","Radar Image #" + indRand,12)
subplot(122)
plotXYdB(idealImages(:,:,indRand),im.y_m,im.z_m,-25,"y (m)","z (m)","Ideal Image #" + indRand,12)