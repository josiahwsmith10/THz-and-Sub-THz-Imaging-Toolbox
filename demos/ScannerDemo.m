% AntennaArray_demo - A short demonstration of the THzScanner class
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

doc THzScanner

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
    0   0   3.5 5   1
    0   0   7.5 5   1];
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
scanner.method = "Rectilinear";
scanner.xStep_m = wav.lambda_m/4;
scanner.numX = 256;
scanner.yStep_m = wav.lambda_m*4;
scanner.numY = 16;

scanner.Compute()

% Display the Synthetic Array
scanner.Display()