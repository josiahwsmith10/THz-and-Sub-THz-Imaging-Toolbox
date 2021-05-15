% sarAntennaArray object holds the properties and methods
% of the MIMO antenna array as specified by the user
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

classdef sarAntennaArray < handle
    properties
        
        tx = struct('xy_m',[],'xyz_m',[])   % Structure containing the parameters (location, etc.) of the transmitter antennas
        rx = struct('xy_m',[],'xyz_m',[])   % Structure containing the parameters (location, etc.) of the receiver antennas
        vx = struct('xy_m','xyz_m')         % Structure containing the parameters (location, etc.) of the virtual element antennas
        
        fig = struct('f',[],'h',[])         % Structure containing the figure and handle used for showing the antenna array
    end
    
    properties(SetObservable)
        z0_m = 0                            % Location of the antenna array in the z-plane
        fmcw                                % An fmcwChirpParameters object
        
        % tableTx - Array of transmitter antenna locations and states in the
        % following form:
        %   x (m)   |   x (lambda)  |   y (m)   |   y (lambda)  |   state
        %   A       |   B           |   C       |   D           |   1
        % Where A, B, C, and D are doubles and the x-location of the
        % transmit element is computed by x = A + B*lambda, using lambda as
        % the lambda_m property of the fmcwChirpParameters object fmcw
        tableTx = [
            0   0   1.5 5   1
            0.5 0   2.5 5   0
            0   0   3.5 5   1]
        
        % tableRx - Array of receiver antenna locations and states in the
        % following form:
        %   x (m)   |   x (lambda)  |   y (m)   |   y (lambda)  |   state
        %   A       |   B           |   C       |   D           |   1
        % Where A, B, C, and D are doubles and the x-location of the
        % receive element is computed by x = A + B*lambda, using lambda as
        % the lambda_m property of the fmcwChirpParameters object fmcw
        tableRx = [
            0   0   0   0   1
            0   0   0.5 0   1
            0   0   1   0   1
            0   0   1.5 0   1]
        
        isEPC = false                       % Boolean whether or not to use the equivalent phase center virtual element locations instead of the MIMO physical element locations
    end
    
    methods
        function obj = sarAntennaArray(fmcw)
            % Attaches the listener to the sarAntennaArray object so
            % observable properties can be watched for changes and sets the
            % fmcw parameter to the fmcwChirpParameters object input
            
            attachListener(obj);
            obj.fmcw = fmcw;
        end
        
        function obj = computeAntennaArray(obj)
            % Computes the physical and virtual antenna locations
            
            obj.tx.xy_m = single(obj.tableTx);
            obj.rx.xy_m = single(obj.tableRx);
            
            % Get number or Tx, Rx, and Vx
            obj.tx.numTx = sum(obj.tx.xy_m(:,5));
            obj.rx.numRx = sum(obj.rx.xy_m(:,5));
            obj.vx.numVx = obj.tx.numTx * obj.rx.numRx;
            
            % Get only the enabled elements
            if obj.tx.numTx > 0
                obj.tx.xy_m = obj.tx.xy_m(logical(obj.tx.xy_m(:,5)),[1,3])*obj.fmcw.lambda_m + obj.tx.xy_m(logical(obj.tx.xy_m(:,5)),[2,4])*1e-3;
            else
                obj.tx.xy_m = [];
            end
            
            if obj.rx.numRx > 0
                obj.rx.xy_m = obj.rx.xy_m(logical(obj.rx.xy_m(:,5)),[1,3])*obj.fmcw.lambda_m + obj.rx.xy_m(logical(obj.rx.xy_m(:,5)),[2,4])*1e-3;
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
            obj.tx.xyz_m = [obj.tx.xy_m,obj.z0_m*ones(obj.tx.numTx,1)];
            obj.rx.xyz_m = [obj.rx.xy_m,obj.z0_m*ones(obj.rx.numRx,1)];
            obj.tx.xyz_m = repmat(obj.tx.xyz_m,obj.rx.numRx,1);
            obj.tx.xyz_m = reshape(obj.tx.xyz_m,obj.tx.numTx,obj.rx.numRx,3);
            obj.tx.xyz_m = permute(obj.tx.xyz_m,[2,1,3]);
            obj.tx.xyz_m = reshape(obj.tx.xyz_m,obj.vx.numVx,3);
            obj.rx.xyz_m = repmat(obj.rx.xyz_m,obj.tx.numTx,1);
            
            obj.tx.xyz_m = reshape(obj.tx.xyz_m,obj.vx.numVx,1,3);
            obj.rx.xyz_m = reshape(obj.rx.xyz_m,obj.vx.numVx,1,3);
            obj.vx.xyz_m = reshape([obj.vx.xy_m,obj.z0_m*ones(obj.vx.numVx,1)],obj.vx.numVx,1,3);
            
            temp = mean(obj.vx.xyz_m,1);
            temp(3) = 0;
            
            obj.tx.xyz_m = obj.tx.xyz_m - temp;
            obj.rx.xyz_m = obj.rx.xyz_m - temp;
            obj.vx.xyz_m = obj.vx.xyz_m - temp;
        end
        
        function initializeFigures(obj)
            % Initializes the figures
            
            closeFigures(obj)
            set(0,'DefaultFigureWindowStyle','docked')
            
            obj.fig.f = figure;
            obj.fig.h = handle(axes);
        end
        
        function closeFigures(obj)
            % Attempts to close the figures
            
            try
                close(obj.fig.f)
            catch
            end
        end
        
        function displayAntennaArray(obj)
            % Displays the antenna arrays in the figure
            
            if isempty(obj.fig.f) || ~isvalid(obj.fig.h)
                initializeFigures(obj);
            end
            
            if obj.isEPC
                displayAntennaArrayEPC(obj);
            else
                displayAntennaArrayMIMO(obj);
            end
        end
        
        function displayAntennaArrayMIMO(obj)
            % Plots the MIMO array
            
            if ~obj.tx.numTx || ~obj.rx.numRx
                return;
            end
            
            h = obj.fig.h;
            hold(h,'off')
            scatter(h,obj.tx.xy_m(:,1)/obj.fmcw.lambda_m,obj.tx.xy_m(:,2)/obj.fmcw.lambda_m,'xr');
            hold(h,'on')
            scatter(h,obj.rx.xy_m(:,1)/obj.fmcw.lambda_m,obj.rx.xy_m(:,2)/obj.fmcw.lambda_m,'ob');
            legend(h,"Tx","Rx")
            xlabel(h,"x (\lambda m)")
            xlim(h,[min(min(obj.tx.xy_m(:,1)/obj.fmcw.lambda_m),min(obj.rx.xy_m(:,1)/obj.fmcw.lambda_m))-1,max(max(obj.tx.xy_m(:,1)/obj.fmcw.lambda_m),max(obj.rx.xy_m(:,1)/obj.fmcw.lambda_m))+1])
            ylim(h,[min(min(obj.tx.xy_m(:,2)/obj.fmcw.lambda_m),min(obj.rx.xy_m(:,2)/obj.fmcw.lambda_m))-1,max(max(obj.tx.xy_m(:,2)/obj.fmcw.lambda_m),max(obj.rx.xy_m(:,2)/obj.fmcw.lambda_m))+1])
            ylabel(h,"y (\lambda m)")
            title(h,"Physical Array (x-y)")
        end
        
        function displayAntennaArrayEPC(obj)
            % Plots the equivalent-phase-center (EPC) virtual array
            
            if ~obj.tx.numTx || ~obj.rx.numRx
                return;
            end
            
            h = obj.fig.h;
            hold(h,'off')
            scatter(h,obj.vx.xy_m(:,1)/obj.fmcw.lambda_m,obj.vx.xy_m(:,2)/obj.fmcw.lambda_m,'.k');
            hold(h,'on')
            legend(h,"Vx")
            xlabel(h,"x (\lambda m)")
            xlim(h,[min(min(obj.tx.xy_m(:,1)/obj.fmcw.lambda_m),min(obj.rx.xy_m(:,1)/obj.fmcw.lambda_m))-1,max(max(obj.tx.xy_m(:,1)/obj.fmcw.lambda_m),max(obj.rx.xy_m(:,1)/obj.fmcw.lambda_m))+1])
            ylim(h,[min(min(obj.tx.xy_m(:,2)/obj.fmcw.lambda_m),min(obj.rx.xy_m(:,2)/obj.fmcw.lambda_m))-1,max(max(obj.tx.xy_m(:,2)/obj.fmcw.lambda_m),max(obj.rx.xy_m(:,2)/obj.fmcw.lambda_m))+1])
            ylabel(h,"y (\lambda m)")
            title(h,"Virtual Array (x-y)")
        end
        
        function saveAntennaArray(obj,saveName)
            % Saves the sarAntennaArray object
            
            if exist(saveName + ".mat",'file')
                str = input('Are you sure you want to overwrite? Y/N: ','s');
                if str ~= 'Y'
                    warning("Antenna array not saved!");
                    return;
                end
            end
            
            savedant = obj;
            savedant.fig = [];
            save(saveName,"savedant");
            disp("Anteanna array saved to: " + saveName);
        end
        
        function loadAntennaArray(obj,loadName)
            % Loads the sarAntennaArray object
            
            if ~exist(loadName + ".mat",'file')
                warning("No file called " + loadName + ".mat to load. Antenna array not loaded!");
                return;
            end
            
            load(loadName,"savedant");
            
            obj.isEPC = savedant.isEPC;
            obj.z0_m = savedant.z0_m;
            obj.tableTx = savedant.tableTx;
            obj.tableRx = savedant.tableRx;
        end
        
        function attachListener(obj)
            % Attaches the listener to the object handle
            
            addlistener(obj,{'isEPC','z0_m','tableTx','tableRx','fmcw'},'PostSet',@sarAntennaArray.propChange);
        end
    end
    
    methods(Static)
        function propChange(metaProp,eventData)
            % Recomputes the array positions if the watched properties change
            
            computeAntennaArray(eventData.AffectedObject);
        end
    end
end