% THzAntennaArray object holds the properties and methods of the MIMO
% antenna array as specified by the user
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

classdef THzAntennaArray < handle
    properties
        tx = struct('xy_m',[],'xyz_m',[])   % Structure containing the parameters (location, etc.) of the transmitter antennas
        rx = struct('xy_m',[],'xyz_m',[])   % Structure containing the parameters (location, etc.) of the receiver antennas
        vx = struct('xy_m','xyz_m')         % Structure containing the parameters (location, etc.) of the virtual element antennas
        
        a = struct                          % Struct holding the antenna element properties
        p = struct                          % Struct holding the antenna pattern properties
        
        wav                                 % A THzWaveformParameters handle
        fig = struct("f",[],"h",[])         % Structure containing the figure and handle used for showing the antenna array
        
        isEPC = false                       % Boolean whether or not to use the equivalent phase center virtual element locations instead of the MIMO physical element locations. Set by the user, not the app
        
        app = struct("UIFigure",[])         % App handle
        isApp                               % Boolean whether or not to use the app functionality
        z0_m = 0                            % Location of the antenna array in the z-plane
        
        % tableTx - Array of transmitter antenna locations and states in the
        % following form:
        %   x (m)   |   x (lambda)  |   y (m)   |   y (lambda)  |   state
        %   A       |   B           |   C       |   D           |   1
        % Where A, B, C, and D are doubles and the x-location of the
        % transmit element is computed by x = A + B*lambda, using lambda as
        % the lambda_m property of the THzWaveformParameters object wav
        tableTx = [
            0   0   3.5 1   1
            0   0   7.5 1   1]
        
        % tableRx - Array of receiver antenna locations and states in the
        % following form:
        %   x (m)   |   x (lambda)  |   y (m)   |   y (lambda)  |   state
        %   A       |   B           |   C       |   D           |   1
        % Where A, B, C, and D are doubles and the x-location of the
        % receive element is computed by x = A + B*lambda, using lambda as
        % the lambda_m property of the THzWaveformParameters object wav
        tableRx = [
            0   0   0   0   1
            0   0   0.5 0   1
            0   0   1   0   1
            0   0   1.5 0   1
            0   0   2   0   1
            0   0   2.5 0   1
            0   0   3   0   1
            0   0   3.5 0   1]
    end
    
    methods
        function obj = THzAntennaArray(wav,app)
            % Sets up the wav and app properties
            
            obj.wav = wav;
            if nargin == 1
                obj.isApp = false;
            elseif nargin == 2
                obj.isApp = true;
                obj.app = app;
            end
            
            obj.a.typeOld = "";
            obj.a.type = "Isotropic";
            obj.p.isConfigured = false;
            obj.p.isWidebandConfigured = false;
            obj.p.isWideband = false;
        end
        
        function obj = Update(obj)
            % Updates the antenna array and displays in the figures
            
            if obj.isApp
                obj = Get(obj);
            end
            obj = Compute(obj);
            
            Display(obj);
        end
        
        function obj = Get(obj)
            % Gets the antenna array values from the app
            
            if ~obj.isApp
                warning("Must be using app to call Get method");
                return;
            end
            
            obj.wav = obj.app.wav;
            obj.tx.z0_m = obj.app.TransmitterZmEditField.Value;
            obj.rx.z0_m = obj.app.ReceiverZmEditField.Value;
            obj.z0_m = obj.app.ReceiverZmEditField.Value;
            obj.tx.xy_m = table2array(obj.app.TxTable.Data);
            obj.rx.xy_m = table2array(obj.app.RxTable.Data);
        end
        
        function obj = Compute(obj)
            % Computes the physical and virtual antenna locations
            
            if ~obj.isApp && (isempty(obj.tableTx) || isempty(obj.tableRx))
                warning("Must specify tableTx and tableRx if not using the app. See ""doc AntennaArray"" for more help on the formatting");
            elseif ~obj.isApp && ~isempty(obj.tableTx) && ~isempty(obj.tableRx)
                obj.tx.xy_m = single(obj.tableTx);
                obj.rx.xy_m = single(obj.tableRx);
                obj.tx.z0_m = obj.z0_m;
                obj.rx.z0_m = obj.z0_m;
            end
            
            % Get number or Tx, Rx, and Vx
            obj.tx.numTx = sum(obj.tx.xy_m(:,5));
            obj.rx.numRx = sum(obj.rx.xy_m(:,5));
            obj.vx.numVx = obj.tx.numTx * obj.rx.numRx;
            
            if obj.vx.numVx == 0
                return;
            end
            
            % Get only the enabled elements
            if obj.tx.numTx > 0
                obj.tx.xy_m = obj.tx.xy_m(logical(obj.tx.xy_m(:,5)),[1,3])*obj.wav.lambda_m + obj.tx.xy_m(logical(obj.tx.xy_m(:,5)),[2,4])*1e-3;
            else
                obj.tx.xy_m = [];
            end
            
            if obj.rx.numRx > 0
                obj.rx.xy_m = obj.rx.xy_m(logical(obj.rx.xy_m(:,5)),[1,3])*obj.wav.lambda_m + obj.rx.xy_m(logical(obj.rx.xy_m(:,5)),[2,4])*1e-3;
            else
                obj.rx.xy_m = [];
            end
            
            obj.vx.xy_m = [];
            obj.vx.dxy = [];
            
            for indTx = 1:obj.tx.numTx
                obj.vx.xy_m = cat(1,obj.vx.xy_m,(obj.tx.xy_m(indTx,:) + obj.rx.xy_m)/2);
                obj.vx.dxy = cat(1,obj.vx.dxy,obj.tx.xy_m(indTx,:) - obj.rx.xy_m);
            end
            
            % Setup Tx/Rx xyz spacing
            obj.tx.xyz_m = [obj.tx.xy_m,obj.tx.z0_m*ones(obj.tx.numTx,1)];
            obj.rx.xyz_m = [obj.rx.xy_m,obj.rx.z0_m*ones(obj.rx.numRx,1)];
            obj.tx.xyz_m = repmat(obj.tx.xyz_m,obj.rx.numRx,1);
            obj.tx.xyz_m = reshape(obj.tx.xyz_m,obj.tx.numTx,obj.rx.numRx,3);
            obj.tx.xyz_m = permute(obj.tx.xyz_m,[2,1,3]);
            obj.tx.xyz_m = reshape(obj.tx.xyz_m,obj.vx.numVx,3);
            obj.rx.xyz_m = repmat(obj.rx.xyz_m,obj.tx.numTx,1);
            
            obj.tx.xyz_m = reshape(obj.tx.xyz_m,obj.vx.numVx,1,3);
            obj.rx.xyz_m = reshape(obj.rx.xyz_m,obj.vx.numVx,1,3);
            obj.vx.xyz_m = reshape([obj.vx.xy_m,obj.tx.z0_m*ones(obj.vx.numVx,1)],obj.vx.numVx,1,3);
            
            temp = mean(obj.vx.xyz_m,1);
            temp(3) = 0;
            
            obj.tx.xyz_m = obj.tx.xyz_m - temp;
            obj.rx.xyz_m = obj.rx.xyz_m - temp;
            obj.vx.xyz_m = obj.vx.xyz_m - temp;
        end
        
        function obj = GetAntenna(obj)
            % Gets the antenna element properties from the app
            
            if ~obj.isApp
                warning("Must be using app to call GetAntenna method");
                return;
            end
            
            obj.p.precision = obj.app.AngularPrecisionEditField.Value;
            obj.p.isWideband = obj.app.UsewidebandantennapatternCheckBox.Value;
            obj.a.type = string(obj.app.TypeDropDown.Value);
            obj.a.length_m = obj.app.LengthmmEditField.Value*1e-3;
            obj.a.width_m = obj.app.WidthmmEditField.Value*1e-3;
            
            switch obj.a.type
                case "Isotropic"
                    obj.app.LengthmmEditField.Enable = false;
                    obj.app.WidthmmEditField.Enable = false;
                    obj.app.ConfigureAntennaPatternButton.Enable = false;
                    obj.app.ShowAntennaPatternButton.Enable = false;
                    obj.app.UsewidebandantennapatternCheckBox.Value = true;
                    obj.app.UsewidebandantennapatternCheckBox.Enable = false;
                case "HFSS"
                    obj.app.LengthmmEditField.Enable = false;
                    obj.app.WidthmmEditField.Enable = false;
                    obj.app.ShowAntennaPatternButton.Enable = false;
                    obj.app.ConfigureAntennaPatternButton.Enable = true;
                    obj.app.UsewidebandantennapatternCheckBox.Value = false;
                    obj.app.UsewidebandantennapatternCheckBox.Enable = false;
                case "Patch"
                    obj.app.LengthmmEditField.Enable = true;
                    obj.app.WidthmmEditField.Enable = true;
                    obj.app.ConfigureAntennaPatternButton.Enable = true;
                    obj.app.ShowAntennaPatternButton.Enable = true;
                    obj.app.UsewidebandantennapatternCheckBox.Enable = true;
                case "Dipole"
                    obj.app.LengthmmEditField.Enable = true;
                    obj.app.WidthmmEditField.Enable = true;
                    obj.app.ConfigureAntennaPatternButton.Enable = true;
                    obj.app.ShowAntennaPatternButton.Enable = true;
                    obj.app.UsewidebandantennapatternCheckBox.Enable = true;
                otherwise
            end
        end
        
        function obj = UpdateAntenna(obj)
            % Updates the antenna property
            
            if obj.isApp
                obj.app.AntennaArrayUptoDateLamp.Color = 'yellow';
                drawnow
                
                obj = GetAntenna(obj);
            end
            
            try
                obj = ConfigureAntenna(obj);
            catch ME
                if obj.isApp
                    uialert(obj.app.UIFigure,ME.message,"Invalid Frequency");
                    obj.app.AntennaArrayUptoDateLamp.Color = 'red';
                end
                obj.p.isConfigured = false;
                return;
            end
            
            if obj.isApp
                obj.app.AntennaArrayUptoDateLamp.Color = 'green';
            end
        end
        
        function obj = ConfigureAntenna(obj)
            % Computes the values of the antenna element
            
            switch obj.a.type
                case "Isotropic"
                    obj.p.isConfigured = true;
                case "HFSS"
                    % Do nothing
                case "Patch"
                    if obj.wav.fC + obj.wav.B > 200e9
                        error("Maximum frequency MATLAB Antenna Toolbox can tolerate is 200 GHz!");
                    end
                    obj.a.ant = design(patchMicrostrip,obj.wav.fC);
                    obj.a.ant.Length = obj.a.length_m;
                    obj.a.ant.Width = obj.a.width_m;
                    obj.p.isConfigured = false;
                    obj.p.isWidebandConfigured = false;
                case "Dipole"
                    if obj.wav.fC + obj.wav.B > 200e9
                        error("Maximum frequency MATLAB Antenna Toolbox can tolerate is 200 GHz!");
                    end
                    con = metal('Name','Copper','Conductivity', 5.96e7,'Thickness',17e-6);
                    obj.a.ant = dipole('Length',obj.a.length_m,'Width',obj.a.width_m,...
                        'Conductor',con,'Tilt',-90,'TiltAxis',[1 0 0]);
                    obj.p.isConfigured = false;
                    obj.p.isWidebandConfigured = false;
                otherwise
                    
            end
        end
        
        function obj = DisplayAntennaPattern(obj)
            % Displays the antenna pattern
            
            figure;
            pattern(obj.a.ant,obj.wav.fC);
        end
        
        function obj = ConfigurePattern(obj)
            % Creates the patternFunc property of the THzAntennaArray
            % object which will return the normalized gain of the antenna
            % pattern with the specified azimuth and elevation angle and
            % frequency
            
            if obj.isApp
                obj = GetAntenna(obj);
            end
                
            % Only execute this command if the type is new
            if obj.p.isConfigured || obj.a.type == "Isotropic"
                obj.p.isConfigured = true;
                return;
            end
            
            % Bypass if HFSS mode
            if obj.a.type == "HFSS"
                obj = ConfigureHFSSAntenna(obj);
                return;
            end
            
            if obj.wav.fC + obj.wav.B > 200e9
                if obj.isApp
                    uialert(obj.app.UIFigure,"Maximum frequency MATLAB Antenna Toolbox can tolerate is 200 GHz!","Invalid Frequency");
                    obj.app.AntennaArrayUptoDateLamp.Color = 'red';
                end
                obj.p.isConfigured = false;
                return;
            end
            
            obj = UpdateAntenna(obj);
            
            % Create the progress dialog
            if obj.isApp
                d = uiprogressdlg(obj.app.UIFigure,'Title','Computing Antenna Pattern',...
                    'Message',"Estimated Time Remaining: 0:0:0","Cancelable","on");
            else
                d = waitbar(0,'1','Name',' Computing Antenna Pattern...',...
                    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
            end
            
            az = 0:obj.p.precision:180;
            el = 0:obj.p.precision:90;
            
            if ~obj.p.isWideband
                temp = pattern(obj.a.ant,obj.wav.fC,az,el);
                temp = 10.^(temp/10);
                obj.p.pattern = single(temp);
                
                obj.p.pattern = obj.p.pattern/max(obj.p.pattern(:));
                obj.p.patternVec = obj.p.pattern(:);
                obj.p.isWidebandConfigured = false;
            elseif obj.p.isWideband
                obj.p.pattern = zeros(length(el),length(az),obj.wav.Nk,'single');
                tocs = zeros(1,obj.wav.Nk,'single');
                
                for indK = 1:obj.wav.Nk
                    if obj.isApp
                        % Check for cancelation
                        if d.CancelRequested
                            warning("Antenna Pattern not Computed!")
                            obj.p.pattern = ones(length(el),length(az),obj.wav.Nk);
                            obj.p.isConfigured = false;
                            return;
                        end
                    else
                        if getappdata(d,'canceling')
                            warning("Antenna Pattern not Computed!")
                            obj.p.pattern = ones(length(el),length(az),obj.wav.Nk);
                            obj.p.isConfigured = false;
                            delete(d);
                            return;
                        end
                    end
                    
                    tic
                    temp = pattern(obj.a.ant,obj.wav.f(indK),az,el);
                    temp = 10.^(temp/10);
                    obj.p.pattern(:,:,indK) = single(temp);
                    tocs(indK) = toc;
                    
                    if obj.isApp
                        d.Value = indK/obj.wav.Nk;
                        d.Message = "Estimated Time Remaining: " + getEstTime(obj,tocs,indK,obj.wav.Nk);
                    else
                        waitbar(indK/obj.wav.Nk,d,"Computing Antenna Pattern. Estimated Time Remaining: " + getEstTime(obj,tocs,indK,obj.wav.Nk));
                    end
                end
                
                obj.p.pattern = obj.p.pattern/max(obj.p.pattern(:));
                obj.p.patternVec = reshape(obj.p.pattern,[],obj.wav.Nk);
                obj.p.isWidebandConfigured = true;
            end
            
            delete(d);
            
            obj.p.isConfigured = true;
            
            if obj.isApp
                obj.app.ConfigureAntennaPatternButton.FontWeight = 'normal';
            end
        end
        
        function gain = computePattern(obj,az,el)
            % Returns the gain for a given look direction across frequency
            % as a column vector
            
            s = size(az);
            az = az(:);
            el = el(:);
            
            if obj.a.type == "Isotropic"
                gain = 1;
                return;
            end
            
            if ~obj.p.isConfigured
                gain = 1;
                return;
            end
            
            if obj.p.isWideband
                if ~obj.p.isWidebandConfigured
                    obj.p.isConfigured = false;
                    obj = ConfigurePattern(obj);
                end
            end
            
            indAz = round(az/obj.p.precision + 1);
            az = [];
            indEl = round(el/obj.p.precision + 1);
            el = [];
            ind = sub2ind(size(obj.p.pattern(:,:,1)),indEl,indAz);
            indAz = [];
            indEl = [];
            gain = obj.p.patternVec(ind,:);
            gain = reshape(gain,s(1),s(2),[]);
        end
        
        function obj = ConfigureHFSSAntenna(obj)
            % Configures the antenna when importing an HFSS antenna pattern
            % from a csv file
            
            % Get the user's file selection
            [filename,pathname] = uigetfile("./saved/antennas/*.csv","Select Desired HFSS File to Load");
            if filename == 0
                warning("Antenna file not loaded!");
                obj.p.isConfigured = false;
                return;
            else
                loadPathFull = string(pathname) + string(filename);
            end
            
            az = 0:obj.p.precision:180;
            el = 0:obj.p.precision:90;
            
            hfss.t = readCSV_HFSS(loadPathFull);
            hfss.az = unique(hfss.t.Phideg);
            hfss.el = 90 - unique(hfss.t.Thetadeg);
            hfss.data = reshape(hfss.t.dBDirTotal,length(hfss.az),length(hfss.el)).';
            [EL,AZ] = ndgrid(el,az);
            
            temp = interpn(hfss.el(:),hfss.az(:),hfss.data,EL,AZ);
            temp = 10.^(temp/10);
            obj.p.pattern = single(temp);
            
            obj.p.pattern = obj.p.pattern/max(obj.p.pattern(:));
            obj.p.patternVec = obj.p.pattern(:);
            obj.p.isWidebandConfigured = false;
            obj.p.isConfigured = true;
            
            function hfssTable = readCSV_HFSS(filename)
                % Set up the Import Options and import the data
                opts = delimitedTextImportOptions("NumVariables", 3);
                
                % Specify range and delimiter
                opts.DataLines = [2, Inf];
                opts.Delimiter = ",";
                
                % Specify column names and types
                opts.VariableNames = ["Phideg", "Thetadeg", "dBDirTotal"];
                opts.VariableTypes = ["double", "double", "double"];
                
                % Specify file level properties
                opts.ExtraColumnsRule = "ignore";
                opts.EmptyLineRule = "read";
                
                % Import the data
                hfssTable = readtable(filename, opts);
            end
        end
        
        function obj = InitializeFigures(obj)
            % Initializes the figures
            
            CloseFigures(obj)
            set(0,'DefaultFigureWindowStyle','docked')
            
            if obj.isApp
                try
                    obj.fig.f = obj.app.sar.fig.f;
                    obj.fig.anth = obj.app.sar.fig.anth;
                    obj.fig.sarh = obj.app.sar.fig.sarh;
                    if ~(ishandle(obj.fig.anth) && ishandle(obj.fig.sarh))
                        error("error");
                    end
                catch
                    obj.fig.f = figure;
                    obj.fig.anth = handle(subplot(121));
                    obj.fig.sarh = handle(subplot(122));
                end
            else
                obj.fig.f = figure;
                obj.fig.h = handle(axes);
            end
        end
        
        function CloseFigures(obj)
            % Attempts to close the figures
            
            try
                close(obj.fig.f)
            catch
            end
        end
        
        function Display(obj)
            % Plots either the MIMO or EPC array
            
            if obj.isApp
                if isempty(obj.fig.f) || ~isvalid(obj.fig.anth) || ~isvalid(obj.fig.sarh)
                    obj = InitializeFigures(obj);
                end
            else
                if isempty(obj.fig.f) || ~isvalid(obj.fig.h)
                    obj = InitializeFigures(obj);
                end
            end
            
            if isAppEPC(obj.app) || obj.isEPC
                DisplayEPC(obj);
            else
                DisplayMIMO(obj);
            end
        end
        
        function DisplayMIMO(obj)
            % Plots the MIMO antenna array
            
            if ~obj.tx.numTx || ~obj.rx.numRx
                return;
            end
            
            if obj.isApp
                h = obj.fig.anth;
            else
                h = obj.fig.h;
            end
            hold(h,'off')
            scatter(h,obj.tx.xy_m(:,1)/obj.wav.lambda_m,obj.tx.xy_m(:,2)/obj.wav.lambda_m,'xr');
            hold(h,'on')
            scatter(h,obj.rx.xy_m(:,1)/obj.wav.lambda_m,obj.rx.xy_m(:,2)/obj.wav.lambda_m,'ob');
            xlim(h,[min(min(obj.tx.xy_m(:,1)/obj.wav.lambda_m),min(obj.rx.xy_m(:,1)/obj.wav.lambda_m))-1,max(max(obj.tx.xy_m(:,1)/obj.wav.lambda_m),max(obj.rx.xy_m(:,1)/obj.wav.lambda_m))+1])
            ylim(h,[min(min(obj.tx.xy_m(:,2)/obj.wav.lambda_m),min(obj.rx.xy_m(:,2)/obj.wav.lambda_m))-1,max(max(obj.tx.xy_m(:,2)/obj.wav.lambda_m),max(obj.rx.xy_m(:,2)/obj.wav.lambda_m))+1])
            xlabel(h,"x ($\lambda$ m)",'Interpreter','latex')
            ylabel(h,"y ($\lambda$ m)",'Interpreter','latex')
            legend(h,"Tx","Rx",'Interpreter','latex')
            title(h,"Physical Array (x-y)",'Interpreter','latex')
        end
        
        function DisplayEPC(obj)
            % Plots the virtual array
            
            if ~obj.tx.numTx || ~obj.rx.numRx
                return;
            end
            
            if obj.isApp
                h = obj.fig.anth;
            else
                h = obj.fig.h;
            end
            hold(h,'off')
            scatter(h,obj.vx.xy_m(:,1)/obj.wav.lambda_m,obj.vx.xy_m(:,2)/obj.wav.lambda_m,'.k');
            hold(h,'on')
            xlim(h,[min(min(obj.tx.xy_m(:,1)/obj.wav.lambda_m),min(obj.rx.xy_m(:,1)/obj.wav.lambda_m))-1,max(max(obj.tx.xy_m(:,1)/obj.wav.lambda_m),max(obj.rx.xy_m(:,1)/obj.wav.lambda_m))+1])
            ylim(h,[min(min(obj.tx.xy_m(:,2)/obj.wav.lambda_m),min(obj.rx.xy_m(:,2)/obj.wav.lambda_m))-1,max(max(obj.tx.xy_m(:,2)/obj.wav.lambda_m),max(obj.rx.xy_m(:,2)/obj.wav.lambda_m))+1])
            xlabel(h,"x ($\lambda$ m)",'Interpreter','latex')
            ylabel(h,"y ($\lambda$ m)",'Interpreter','latex')
            legend(h,"Vx",'Interpreter','latex')
            title(h,"Virtual Array (x-y)",'Interpreter','latex')
        end
        
        function obj = Load(obj,wav)
            % Loads the antenna array from a file. If second argument is
            % present, sets obj.wav to wav so the THzWaveformParameters
            % object handle is synced
            
            if obj.isApp
                loadPathFull = "./saved/antennaArrays/" + obj.app.ArrayLoadNameEditField.Value + ".mat";
            else
                loadPathFull = [];
            end
            if ~exist(loadPathFull,'file')
                [filename,pathname] = uigetfile("./saved/antennaArrays/*.mat","Select Desired File to Load");
                
                if filename == 0
                    warning("Antenna array file not loaded!");
                    return;
                else
                    loadPathFull = string(pathname) + string(filename);
                end
            end
            
            load(loadPathFull,"savedant","savedwav");
            
            obj = getFields(obj,savedant,["app","fig","wav","isApp"]);
            
            if nargin == 1
                showErrorMessage(obj,"Using savedwav from file, might not be properly linked","Warning")
                % Load the saved wav
                if obj.isApp
                    obj.wav = THzWaveformParameters(obj.app);
                else
                    obj.wav = THzWaveformParameters();
                end
                obj.wav = getFields(obj.wav,savedwav,["app","isApp"]);
                obj.wav = obj.wav.Compute();
                obj.wav = obj.wav.Display();
            else
                obj.wav = wav;
            end
            
            % Compute the new values
            obj.tx.xy_m = single(obj.tableTx);
            obj.rx.xy_m = single(obj.tableRx);
            obj.app.TxTable.Data = array2table(obj.tx.xy_m);
            obj.app.TxTable.Data.Var5 = logical(obj.tx.xy_m(:,end));
            obj.app.RxTable.Data = array2table(obj.rx.xy_m);
            obj.app.RxTable.Data.Var5 = logical(obj.rx.xy_m(:,end));
            obj = Compute(obj);
        end
        
        function Save(obj)
            % Saves the THzAntennaArray object handle as a struct for
            % reloading
            
            if obj.isApp
                savePathFull = "./saved/antennaArrays/" + obj.app.ArraySaveNameEditField.Value + ".mat";
            else
                savePathFull = [];
            end
            if exist(savePathFull,'file')
                [filename,pathname] = uiputfile("./saved/antennaArrays/*.mat","Select Desired File Location + Name for Save");
                
                if filename == 0
                    warning("Antenna array file not saved!");
                    return;
                else
                    savePathFull = string(pathname) + string(filename);
                end
            end
            
            savedant = struct();
            savedant = getFields(savedant,obj,["app","fig","wav","isApp"]);
            
            savedwav = struct();
            savedwav = getFields(savedwav,obj.wav,["app","isApp"]);
            save(savePathFull,"savedant","savedwav");
        end
        
        function outstr = getEstTime(obj,tocs,currentInd,totalInd)
            % Estimates the time until completion
            
            avgtoc = mean(tocs(1:currentInd))*(totalInd - currentInd);
            hrrem = floor(avgtoc/3600);
            avgtoc = avgtoc - floor(avgtoc/3600)*3600;
            minrem = floor(avgtoc/60);
            avgtoc = avgtoc - floor(avgtoc/60)*60;
            secrem = round(avgtoc);
            
            if hrrem < 10
                hrrem = "0" + hrrem;
            end
            if minrem < 10
                minrem = "0" + minrem;
            end
            if secrem < 10
                secrem = "0" + secrem;
            end
            outstr = hrrem + ":" + minrem + ":" + secrem;
        end
    end
end