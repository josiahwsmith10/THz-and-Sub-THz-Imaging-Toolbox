% fmcwChirpParameters_demo - A Short Demonstration of the fmcwChirpParameters class
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

doc fmcwChirpParameters

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

%% Load in Saved FMCW Parameters
fmcw.loadChirpParameters("fmcw_THz");

%% Save New FMCW Parameters
fmcw.ADCSamples = 79;
fmcw.f0 = 77*1e9;
fmcw.RampEndTime_s = 39.98*1e-6;
fmcw.fS = 2000*1e3;
fmcw.K = 100.036*1e12;
fmcw.IdleTime_s = 0*1e-6;
fmcw.TXStartTime_s = 0*1e-6;
fmcw.ADCStartTime_s = 0*1e-6;
fmcw.fC = 79*1e9;
fmcw.saveChirpParameters("fmcw_77GHz");