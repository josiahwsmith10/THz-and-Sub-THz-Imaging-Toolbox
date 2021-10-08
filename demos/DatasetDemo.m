% AntennaArray_demo - A short demonstration of dataset generation
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
wav = THzWaveformParameters();

ant = THzAntennaArray(wav);

scanner = THzScanner(ant);

target = THzTarget(wav,ant,scanner);

im = THzImageReconstruction(wav,ant,scanner,target);

%% Set Waveform Parameters
wav.Nk = 50;
wav.f0 = 300*1e9;
wav.fS = 1000*1e3;
wav.fC = 305*1e9;
wav.B = 10e9;

%% Compute the Waveform Parameters
% Computes the wavenumber vector etc.
wav.Compute();

%% Set Antenna Array Properties
ant.isEPC = false;
ant.z0_m = 0;
% Large MIMO Array
ant.tableTx = [
    0   0   4   0   1
    0   0   8   0   1];
ant.tableRx = [
    0   0   0   0   1
    0   0   0.5 0   1
    0   0   1   0   1
    0   0   1.5 0   1
    0   0   2   0   1
    0   0   2.5 0   1
    0   0   3   0   1
    0   0   3.5 0   1];
ant.Compute();

% Display the Antenna Array
ant.Display();

%% Set Scanner Parameters
scanner.method = "Linear";
scanner.yStep_m = wav.lambda_m*4;
scanner.numY = 16;

scanner.Compute();

% Display the Synthetic Array
scanner.Display();

%% Set Target Parameters
target.isAmplitudeFactor = false;

target.rp.numTargets = 16;
target.rp.xMin_m = 0;
target.rp.xMax_m = 0;
target.rp.yMin_m = -0.03;
target.rp.yMax_m = 0.03;
target.rp.zMin_m = 0.02;
target.rp.zMax_m = 0.18;
target.rp.ampMin = 0.5;
target.rp.ampMax = 1;

% Which to use
target.isTable = false;
target.isRandomPoints = true;

target.Get();

% Display the target
target.Display();

%% Compute the Beat Signal
target.isGPU = true;
target.Compute();

%% Set Image Reconstruction Parameters and Create THzImageReconstruction Object
im.nFFTy = 512;
im.nFFTz = 512;

im.yMin_m = -0.05;
im.yMax_m = 0.05;

im.zMin_m = 0;
im.zMax_m = 0.2;

im.numY = 256;
im.numZ = 256;

im.isGPU = false;
% im.zSlice_m = 0.1; % Use if reconstructing a 1-D image
im.method = "Uniform 1-D SAR 2-D RMA";

im.isMult2Mono = true;
im.zRef_m = 0.1;

%% Reconstruct the Image
im.Compute();
im.Display();

%% Compute the Ideal Image from the YZ Imaging Scene and Display
imIdeal2D = computeIdeal2D(target,im,"YZ",0.8e-3,0.8e-3);
figure
plotXYdB(imIdeal2D,im.y_m,im.z_m,-25,"y (m)","z (m)","Ideal Image",12)

%% Setup Dataset Scenario
numSamples = 32;
numRandomPointsMax = 16;
numRandomPointsMin = 4;

radarImages = single(zeros(im.numY,im.numZ,numSamples));
idealImages = single(zeros(im.numY,im.numZ,numSamples));

% Set objects to silent
target.isSilent = true;
im.isSilent = true;

%% Loop through Image Generation
for indSample = 1:numSamples
    target.rp.numTargets = randi([numRandomPointsMin,numRandomPointsMax]);
    target.Compute();
    
    im.Compute();
    
    radarImages(:,:,indSample) = im.imXYZ;
    idealImages(:,:,indSample) = computeIdeal2D(target,im,"YZ",0.8e-3,0.8e-3);
    disp("Iteration " + indSample + "/" + numSamples + " Done");
end

%% Compare Samples Side by Side
indRand = randi(numSamples);

figure
subplot(121)
plotXYdB(radarImages(:,:,indRand),im.y_m,im.z_m,-25,"y (m)","z (m)","Radar Image #" + indRand,12)
subplot(122)
plotXYdB(idealImages(:,:,indRand),im.y_m,im.z_m,-25,"y (m)","z (m)","Ideal Image #" + indRand,12)