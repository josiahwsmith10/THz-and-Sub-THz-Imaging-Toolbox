% fmcwChirpParameters_app The chirp parameters of the FMCW radar under
% test, adhering loosely to the TI mmWave Studio definitions
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

classdef fmcwChirpParameters_app < handle
    properties
        f0 = 77e9                   % Starting frequency in Hz
        K = 100.036e12              % Chirp slope in Hz/s
        IdleTime_s = 0              % Duration of time before starting the chirp/ramp after the previous chirp
        TXStartTime_s = 0           % Duration of time before chirp/ramp starts wherein transmitter is ON
        ADCStartTime_s = 0          % Duration of time after transmitter is ON and TXStartTime_s has passed wherein chirp/ramp is active, but samples are not being collected (acts as a guard against initial chirp non-linearities)
        RampEndTime_s = 39.98e-6    % Duration of time after transmitter is ON and TXStartTime_s has passed wherein chirp/ramp is active
        ADCSamples = 79             % Number of ADC samples
        fS = 2000e3                 % Sampling frequency in Hz
        fC = 79e9                   % Center frequency in Hz
        ADCSampleTime_time = 0      % Duration of time during which ADC samples are taken, dictated by time parameters
        ADCSampleTime_sample = 0    % Duration of time during which ADC samples are taken, dictated by sampling parameters (ADCSamples & fS)
        B_total = 0                 % Total bandwidth used
        B_useful = 0                % Useful/actual bandwidth covered by ADC samples
        c = 3e8                     % Speed of light
        rangeMax_m = 0              % Maximum range of FMCW chirp parameters in meters
        rangeResolution_m = 0       % Range resolution in meters
        k                           % Instantaneous wavenumber vector
        lambda_m = 0                % Wavelength of center frequency
        B_max = 0                   % Maximum allowable bandwidth, as specified by the user in the app
    end
    
    methods
        function obj = fmcwChirpParameters_app(app)
            % Attaches the listener to the fmcwChirpParameters object and
            % updates the properties from the app
            obj = update(obj,app);
        end
        
        function obj = update(obj,app)
            % Update the FMCW chirp parameters from the app
            
            obj = getChirpParameters(obj,app);
            obj = computeChirpParameters(obj,app);
            displayChirpParameters(obj,app);
        end
        
        function obj = getChirpParameters(obj,app)
            % Get the parameters from the app
            
            obj.f0 = app.StartingFrequencyGHzEditField.Value*1e9;
            obj.K = app.FreqSlopeMHzusEditField.Value*1e12;
            obj.IdleTime_s = app.IdleTimeusEditField.Value*1e-6;
            obj.TXStartTime_s = app.TXStartTimeusEditField.Value*1e-6;
            obj.ADCStartTime_s = app.ADCStartTimeusEditField.Value*1e-6;
            obj.ADCSamples = app.ADCSamplesEditField.Value;
            obj.fS = app.SamplingFrequencykspsEditField.Value*1e3;
            obj.RampEndTime_s = app.RampEndTimeusEditField.Value*1e-6;
            obj.fC = app.CenterFrequencyGHzEditField.Value*1e9;
            obj.B_max = app.MaximumBandwidthGHzEditField.Value*1e9;
        end
        
        function obj = computeChirpParameters(obj,app)
            % Compute the parameters given the current values in the app
            
            obj.ADCSampleTime_time = obj.RampEndTime_s - obj.ADCStartTime_s;
            obj.ADCSampleTime_sample = obj.ADCSamples/obj.fS;
            
            obj.c = physconst('lightspeed');
            
            if obj.ADCSampleTime_sample > obj.ADCSampleTime_time
                uiconfirm(app.UIFigure,"Not enough time to collect " + obj.ADCSamples + " samples at " + obj.fS*1e-3 + " ksps",'Not Enough Time!',...
                    "Options",{'OK'},'Icon','warning');
            end
            obj.B_total = obj.RampEndTime_s*obj.K;
            obj.B_useful = obj.ADCSampleTime_sample*obj.K;
            
            if obj.B_total > obj.B_max || obj.B_useful > obj.B_max
                warning("Bandwidth exceeds maximum!")
                uiconfirm(app.UIFigure,"Bandwidth exceeds maximum! Reduce ADC Samples and/or Ramp End Time",'Not Enough Bandwidth!',...
                    "Options",{'OK'},'Icon','warning');
            end
            
            obj.rangeMax_m = obj.fS*obj.c/(2*obj.K);
            obj.rangeResolution_m = obj.c/(2*obj.B_useful);
            
            f0_temp = obj.f0 + obj.ADCStartTime_s*obj.K; % This is for ADC sampling offset
            f = f0_temp + (0:obj.ADCSamples-1)*obj.K/obj.fS; % wideband frequency
            
            obj.k = 2*pi*f/obj.c;
            obj.lambda_m = obj.c/(obj.fC);
        end
        
        function displayChirpParameters(obj,app)
            % Display the new calculated parameters in the app
            
            app.TotalBandwidthGHzEditField.Value = obj.B_total*1e-9;
            app.UsefulBandwidthGHzEditField.Value = obj.B_useful*1e-9;
            app.WavelengthmmEditField.Value = obj.lambda_m*1e3;
            app.MaximumRangemEditField.Value = obj.rangeMax_m;
            app.RangeResolutionmmEditField.Value = obj.rangeResolution_m*1e3;
        end
        
        function obj = loadChirpParameters(obj,app)
            % Load the chirp parameters from a file
            
            loadPathFull = "./saved/chirpParameters/" + app.FMCWLoadNameEditField.Value + ".mat";
            if ~exist(loadPathFull,'file')
                uiconfirm(app.UIFigure,"No file called " + app.FMCWLoadNameEditField.Value + ".mat to load",'Cannot Load',...
                    "Options",{'OK'},'Icon','warning');
                warning("Parameters not loaded!");
                return;
            end
            
            load(loadPathFull,"savedfmcw");
            
            app.StartingFrequencyGHzEditField.Value = savedfmcw.f0*1e-9;
            app.FreqSlopeMHzusEditField.Value = savedfmcw.K*1e-12;
            app.IdleTimeusEditField.Value = savedfmcw.IdleTime_s*1e6;
            app.TXStartTimeusEditField.Value = savedfmcw.TXStartTime_s*1e6;
            app.ADCStartTimeusEditField.Value = savedfmcw.ADCStartTime_s*1e6;
            app.ADCSamplesEditField.Value = savedfmcw.ADCSamples;
            app.SamplingFrequencykspsEditField.Value = savedfmcw.fS*1e-3;
            app.RampEndTimeusEditField.Value = savedfmcw.RampEndTime_s*1e6;
            app.CenterFrequencyGHzEditField.Value = savedfmcw.fC*1e-9;
            app.MaximumBandwidthGHzEditField.Value = savedfmcw.B_max*1e-9;
            
            obj = update(obj,app);
        end
        
        function saveChirpParameters(obj,app)
            % Save the chirp parameters to a file
            
            savePathFull = "./saved/chirpParameters/" + app.FMCWSaveNameEditField.Value + ".mat";
            if exist(savePathFull,'file')
                selection = uiconfirm(app.UIFigure,'Are you sure you want to overwrite?','Confirm Overwrite',...
                    'Icon','warning');
                if string(selection) == "Cancel"
                    warning("Parameters not saved!");
                    return;
                end
            end
            
            savedfmcw = obj;
            save(savePathFull,"savedfmcw");
        end
    end
end