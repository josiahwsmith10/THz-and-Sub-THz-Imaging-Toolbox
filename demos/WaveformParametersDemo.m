% WaveformParameters_demo - A short demonstration of the 
% THzWaveformParameters class
%
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

doc THzWaveformParameters

%% Set Waveform Parameters
wav.Nk = 50;
wav.f0 = 300*1e9;
wav.fS = 1000*1e3;
wav.fC = 305*1e9;
wav.B = 10e9;

%% Compute the Waveform Parameters
% Computes the wavenumber vector etc.
wav.Compute();

%% Load in Saved Waveform Parameters
wav.Load();

%% Save New Waveform Parameters
wav.Nk = 79;
wav.f0 = 77*1e9;
wav.fS = 2000*1e3;
wav.fC = 79*1e9;
wav.B = 4e9;
wav.Compute();

wav.Save();