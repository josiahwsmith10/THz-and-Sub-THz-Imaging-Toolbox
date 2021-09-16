% THzWaveformParameters The chirp parameters of the THz sensors under test
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

classdef THzWaveformParameters < handle
    properties
        f0 = 300e9                  % Starting frequency in Hz
        K = 200e12                  % Chirp slope in Hz/s
        Nk = 50                     % Number of frequency/wavenumber samples
        fS = 1000e3                 % Sampling frequency in Hz
        fC = 305e9                  % Center frequency in Hz
        B = 10e9                    % Bandwidth in Hz
        c = 299792458               % Speed of light in m/s
        rangeMax_m = 0              % Unambiguous range in m
        rangeResolution_m = 0       % Range resolution in meters
        k                           % Instantaneous wavenumber vector
        f                           % Instantaneous frequency vector
        lambda_m = 0                % Wavelength of center frequency
        
        app = struct("UIFigure",[]) % App handle
        isApp                       % Boolean whether or not to use the app functionality
    end
    
    methods
        function obj = THzWaveformParameters(app)
            % Intializes the WaveformParameters object taking the option
            % input argument of the app
            
            if nargin == 0
                obj.isApp = false;
            elseif nargin == 1
                obj.isApp = true;
                obj.app = app;
            end
            
            obj = Update(obj);
        end
        
        function obj = Update(obj)
            % Update the FMCW chirp parameters from the app
            
            if obj.isApp
                obj = Get(obj);
            end
            obj = Compute(obj);
            obj = Display(obj);
        end
        
        function obj = Get(obj)
            % Get the parameters from the app
            
            if ~obj.isApp
                warning("Must be using app to call Get method");
                return;
            end
            
            obj.f0 = obj.app.StartingFrequencyGHzEditField.Value*1e9;
            obj.Nk = obj.app.ADCSamplesEditField.Value;
            obj.fS = obj.app.SamplingFrequencykspsEditField.Value*1e3;
            obj.fC = obj.app.CenterFrequencyGHzEditField.Value*1e9;
            obj.B = obj.app.BandwidthGHzEditField.Value*1e9;
        end
        
        function obj = Compute(obj)
            % Compute the parameters given the current values in the app
            
            if obj.Nk > 1
                obj.K = obj.fS*obj.B/(obj.Nk-1);
            else
                obj.K = 0;
            end
            
            obj.rangeMax_m = obj.fS*obj.c/(2*obj.K);
            obj.rangeResolution_m = obj.c/(2*obj.B);
            obj.f = obj.f0 + (0:obj.Nk-1)*obj.K/obj.fS;
            
            obj.k = 2*pi*obj.f/obj.c;
            obj.lambda_m = obj.c/(obj.fC);
        end
        
        function obj = Display(obj)
            % Display the new calculated parameters in the app
            
            obj.app.WavelengthmmEditField.Value = obj.lambda_m*1e3;
            obj.app.MaximumRangemEditField.Value = obj.rangeMax_m;
            obj.app.RangeResolutionmmEditField.Value = obj.rangeResolution_m*1e3;
        end
        
        function obj = Load(obj)
            % Load the waveform parameters from a file
            
            if obj.isApp
                loadPathFull = "./saved/waveformParameters/" + obj.app.WaveformLoadNameEditField.Value + ".mat";
            else
                loadPathFull = [];
            end
            if ~exist(loadPathFull,'file')
                [filename,pathname] = uigetfile("./saved/waveformParamters/*.mat","Select Desired File to Load");
                
                if filename == 0
                    warning("Waveform parameters file not loaded!");
                    return;
                else
                    loadPathFull = string(pathname) + string(filename);
                end
            end
            
            load(loadPathFull,"savedwav");
            
            obj = getFields(obj,savedwav,["app","isApp"]);
            obj = Update(obj);
        end
        
        function Save(obj)
            % Save the chirp parameters to a file
            
            if obj.isApp
                savePathFull = "./saved/waveformParameters/" + obj.app.WaveformSaveNameEditField.Value + ".mat";
            else
                savePathFull = [];
            end
            if ~exist(savePathFull,'file')
                [filename,pathname] = uiputfile("./saved/waveformParamters/*.mat","Select Desired File Location + Name for Save");
                
                if filename == 0
                    warning("Waveform parameters file not saved!");
                    return;
                else
                    savePathFull = string(pathname) + string(filename);
                end
            end
            
            savedwav = struct();
            savedwav = getFields(savedwav,obj,["app","isApp"]);
            save(savePathFull,"savedwav");
        end
    end
end