% sarScenario object holds the properties and methods
% used for simulating a MIMO-SAR scanning scenario in four modes:
% linear, rectilinear (planar), circular, or cylindrical
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

classdef sarScenario < handle
    properties(SetObservable)
        scanMethod                  % Scanning mode: "Linear", "Rectilinear", "Circular", or "Cylindrical"
        xStep_m = 1e-3              % Step size along the x-dimension to move the antenna array in meters
        yStep_m = 8e-3              % Step size along the y-dimension to move the antenna array in meters
        thetaMax_deg = 360          % Maximum angle for circular or cylindrical scans in degrees
        numX = 1                    % Number of steps along the x-dimension
        numY = 1                    % Number of steps along the y-dimension
        numTheta = 1                % Number of angular steps
        
        ant                         % sarAntennaArray object
    end
    
    properties
        X_m                         % Temporary storage of the entire aperture's x-coordinates
        Y_m                         % Temporary storage of the entire aperture's y-coordinates
        Z_m                         % Temporary storage of the entire aperture's z-coordinates
        xyz_m                       % Location in x-y-z of the synthetic aperture, an Nx3 array, where N is the total number of SAR elements and each row contains the [x,y,z] coordinates of each element
        tx                          % Structure containing the parameters (location, etc.) of the transmitter antennas
        rx                          % Structure containing the parameters (location, etc.) of the receiver antennas
        vx                          % Structure containing the parameters (location, etc.) of the virtual element antennas
        x_m                         % Scanning points along which the antenna array is scanned in the x-dimension
        y_m                         % Scanning points along which the antenna array is scanned in the y-dimension
        z_m                         % Scanning points along which the antenna array is scanned in the z-dimension
        theta_rad                   % Scanning points along which the antenna array is scanned in the angular dimension
        xSize_m                     % Size of the synthetic aperture along the x-dimension
        ySize_m                     % Size of the synthetic aperture along the y-dimension
        
        sarSize                     % Size of the SAR data after simulation (dimensions of the scan)
        
        fig = struct('f',[],'h',[]) % Structure containing the figure and handle used for showing the target
    end
    
    methods
        function obj = sarScenario(ant)
            % Attaches the listener to the sarAntennaArray object so
            % observable properties can be watched for changes and sets the
            % ant parameter to the sarAntennaArray object input
            
            obj.ant = ant;
            attachListener(obj);
        end
        
        function computeSarScenario(obj)
            % Computes the SAR scenario depending on scanMethod
            
            switch obj.scanMethod
                case "Linear"
                    % Verify linearity of MIMO array
                    if verifyLinearity(obj) || verifyLinear(obj)
                        return
                    end
                    computeLinearScan(obj);
                    
                case "Rectilinear"
                    % Verify linearity of MIMO array
                    if verifyLinearity(obj) || verifyRectilinear(obj)
                        return
                    end
                    computeRectilinearScan(obj);
                    
                case "Circular"
                    if verifyCircular(obj)
                        return
                    end
                    computeCircularScan(obj);
                    
                case "Cylindrical"
                    % Verify linearity of MIMO array
                    if verifyLinearity(obj) || verifyCylindrical(obj)
                        return
                    end
                    computeCylindricalScan(obj);
            end
            
            obj.tx.xyz_m = single(obj.tx.xyz_m);
            obj.rx.xyz_m = single(obj.rx.xyz_m);
            obj.vx.xyz_m = single(obj.vx.xyz_m);
        end
        
        function computeLinearScan(obj)
            % Computes the SAR scenario for a "Linear" scan
            
            getsarAxes(obj);
            
            [obj.X_m,obj.Y_m,obj.Z_m] = ndgrid(obj.x_m,obj.y_m,obj.z_m);
            obj.xyz_m = reshape(cat(3,obj.X_m,obj.Y_m,obj.Z_m),1,[],3);
            obj.xyz_m = repmat(obj.xyz_m,obj.ant.vx.numVx,1,1);
            
            obj.tx.xyz_m = single(obj.xyz_m + obj.ant.tx.xyz_m);
            obj.rx.xyz_m = single(obj.xyz_m + obj.ant.rx.xyz_m);
            obj.vx.xyz_m = single(obj.xyz_m + obj.ant.vx.xyz_m);
            
            obj.tx.xyz_m = reshape(obj.tx.xyz_m,[],3);
            obj.rx.xyz_m = reshape(obj.rx.xyz_m,[],3);
            obj.vx.xyz_m = reshape(obj.vx.xyz_m,[],3);
            
            if ~obj.ant.isEPC
                % Unwrap obj.tx.xyz_m & obj.rx.xyz_m as [numRx,numTx,numY,3]
                obj.sarSize = [obj.ant.rx.numRx,obj.ant.tx.numTx,obj.numY];
            else
                obj.sarSize = obj.ant.vx.numVx*obj.numY;
            end
        end
        
        function computeRectilinearScan(obj)
            % Computes the SAR scenario for a "Rectilinear" scan
            
            getsarAxes(obj);
            
            [obj.Y_m,obj.X_m,obj.Z_m] = ndgrid(obj.y_m,obj.x_m,obj.z_m);
            obj.xyz_m = reshape(cat(3,obj.X_m,obj.Y_m,obj.Z_m),1,[],3);
            obj.xyz_m = repmat(obj.xyz_m,obj.ant.vx.numVx,1,1);
            
            obj.tx.xyz_m = single(obj.xyz_m + obj.ant.tx.xyz_m);
            obj.rx.xyz_m = single(obj.xyz_m + obj.ant.rx.xyz_m);
            obj.vx.xyz_m = single(obj.xyz_m + obj.ant.vx.xyz_m);
            
            obj.tx.xyz_m = reshape(obj.tx.xyz_m,[],3);
            obj.rx.xyz_m = reshape(obj.rx.xyz_m,[],3);
            obj.vx.xyz_m = reshape(obj.vx.xyz_m,[],3);
            
            if ~obj.ant.isEPC
                % Unwrap obj.tx.xyz_m & obj.rx.xyz_m as [numRx,numTx,numY,numX,3]
                obj.sarSize = [obj.ant.rx.numRx,obj.ant.tx.numTx,obj.numY,obj.numX];
            else
                obj.sarSize = [obj.ant.vx.numVx*obj.numY,obj.numX];
            end
        end
        
        function computeCircularScan(obj)
            % Computes SAR scenario for "circular" scan
            
            % Verify single element array
            if obj.ant.tx.numTx ~= 1 || obj.ant.rx.numRx ~= 1
                error("Array must have only 1 Tx and 1 Rx. Please disable necessary elements or change the antenna array.")
            end
            
            getsarAxes(obj);
            
            obj.x_m = obj.ant.z0_m*cos(obj.theta_rad);
            obj.y_m = zeros(size(obj.theta_rad));
            obj.z_m = obj.ant.z0_m*sin(obj.theta_rad) - obj.ant.z0_m;
            
            obj.xyz_m = reshape([obj.x_m(:),obj.y_m(:),obj.z_m(:)],1,[],3);
            
            obj.tx.xyz_m = single(obj.xyz_m + obj.ant.tx.xyz_m);
            obj.rx.xyz_m = single(obj.xyz_m + obj.ant.rx.xyz_m);
            obj.vx.xyz_m = single(obj.xyz_m + obj.ant.vx.xyz_m);
            
            obj.tx.xyz_m = reshape(obj.tx.xyz_m,[],3);
            obj.rx.xyz_m = reshape(obj.rx.xyz_m,[],3);
            obj.vx.xyz_m = reshape(obj.vx.xyz_m,[],3);
            
            % Unwrap obj.tx.xyz_m & obj.rx.xyz_m as [numTheta,3]
            obj.sarSize = obj.numTheta;
        end
        
        function computeCylindricalScan(obj)
            % Computes SAR scenario for "Cylindrical" scan
            
            getsarAxes(obj);
            
            obj.x_m = obj.ant.z0_m*cos(obj.theta_rad);
            obj.z_m = obj.ant.z0_m*sin(obj.theta_rad) - obj.ant.z0_m;
            
            % Use y' as first dimension -> so we obtain s(y',theta)
            obj.X_m = repmat(obj.x_m(:).',obj.numY,1);
            obj.Z_m = repmat(obj.z_m(:).',obj.numY,1);
            obj.Y_m = repmat(obj.y_m(:),1,obj.numTheta);
            
            obj.xyz_m = reshape(cat(3,obj.X_m,obj.Y_m,obj.Z_m),1,[],3);
            
            obj.tx.xyz_m = single(obj.xyz_m + obj.ant.tx.xyz_m);
            obj.rx.xyz_m = single(obj.xyz_m + obj.ant.rx.xyz_m);
            obj.vx.xyz_m = single(obj.xyz_m + obj.ant.vx.xyz_m);
            
            obj.tx.xyz_m = reshape(obj.tx.xyz_m,[],3);
            obj.rx.xyz_m = reshape(obj.rx.xyz_m,[],3);
            obj.vx.xyz_m = reshape(obj.vx.xyz_m,[],3);
            
            if ~obj.ant.isEPC
                % Unwrap obj.tx.xyz_m & obj.rx.xyz_m as [numRx,numTx,numY,numTheta,3]
                obj.sarSize = [obj.ant.rx.numRx,obj.ant.tx.numTx,obj.numY,obj.numTheta];
            else
                obj.sarSize = [obj.ant.vx.numVx*obj.numY,obj.numTheta];
            end
        end
        
        function getsarAxes(obj)
            % Computes the scanning position axes
            
            % Update number of steps
            obj.xSize_m = obj.numX * obj.xStep_m;
            obj.ySize_m = obj.numY * obj.yStep_m;
            
            % Create synthetic aperture step axes
            obj.x_m = (-(obj.numX - 1)/2 : (obj.numX - 1)/2) * obj.xStep_m;
            obj.y_m = (-(obj.numY - 1)/2 : (obj.numY - 1)/2) * obj.yStep_m;
            obj.theta_rad = linspace(0,obj.thetaMax_deg - obj.thetaMax_deg/obj.numTheta,obj.numTheta)*2*pi/360;
            obj.z_m = 0;
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
        
        function displaySarScenario(obj)
            % Display the synthetic aperture
            
            if isempty(obj.fig.f) || ~isvalid(obj.fig.h)
                initializeFigures(obj);
            end
            
            if obj.ant.isEPC
                displaySarScenarioEPC(obj);
            else
                displaySarScenarioMIMO(obj);
            end
        end
        
        function displaySarScenarioMIMO(obj)
            % Plots the MIMO-SAR array
            
            if obj.scanMethod == "-"
                return;
            end
            
            h = obj.fig.h;
            hold(h,'off')
            temp = obj.tx.xyz_m;
            scatter3(h,temp(:,1),temp(:,3),temp(:,2),'.r')
            hold(h,'on')
            temp = obj.rx.xyz_m;
            scatter3(h,temp(:,1),temp(:,3),temp(:,2),'.b')
            xlabel(h,"x (m)")
            temp1 = obj.tx.xyz_m(:,1);
            temp2 = obj.rx.xyz_m(:,1);
            xlim(h,[min(min(temp1),min(temp2))-0.01,max(max(temp1),max(temp2))+0.01])
            ylabel(h,"z (m)")
            temp1 = obj.tx.xyz_m(:,3);
            temp2 = obj.rx.xyz_m(:,3);
            ylim(h,[min(min(temp1),min(temp2))-0.01,max(max(temp1),max(temp2))+0.01])
            zlabel(h,"y (m)")
            temp1 = obj.tx.xyz_m(:,2);
            temp2 = obj.rx.xyz_m(:,2);
            zlim(h,[min(min(temp1),min(temp2))-0.01,max(max(temp1),max(temp2))+0.01])
            title(h,"MIMO Synthetic Aperture")
            legend(h,"Tx","Rx")
            view(h,3)
            daspect(h,[1 1 1])
        end
        
        function displaySarScenarioEPC(obj)
            % Plots the equivalent-phase-center (EPC) SAR array
            
            if obj.scanMethod == "-"
                return;
            end
            
            h = obj.fig.h;
            hold(h,'off')
            temp = obj.vx.xyz_m;
            scatter3(h,temp(:,1),temp(:,3),temp(:,2),'.k')
            xlabel(h,"x (m)")
            temp1 = obj.tx.xyz_m(:,1);
            temp2 = obj.rx.xyz_m(:,1);
            xlim(h,[min(min(temp1),min(temp2))-0.01,max(max(temp1),max(temp2))+0.01])
            ylabel(h,"z (m)")
            temp1 = obj.tx.xyz_m(:,3);
            temp2 = obj.rx.xyz_m(:,3);
            ylim(h,[min(min(temp1),min(temp2))-0.01,max(max(temp1),max(temp2))+0.01])
            zlabel(h,"y (m)")
            temp1 = obj.tx.xyz_m(:,2);
            temp2 = obj.rx.xyz_m(:,2);
            zlim(h,[min(min(temp1),min(temp2))-0.01,max(max(temp1),max(temp2))+0.01])
            title(h,"Virtual Synthetic Aperture")
            legend(h,"Vx")
            view(h,3)
            daspect(h,[1 1 1])
        end
        
        function isFail = verifyLinearity(obj)
            % Verifies that the MIMO antenna array is colinear
            
            if max(diff([obj.ant.tx.xy_m(:,1);obj.ant.rx.xy_m(:,1)])) > 8*eps
                warning("MIMO array must be colinear. Please disable necessary elements.")
                obj.scanMethod = "-";
                isFail = true;
            else
                isFail = false;
            end
        end
        
        function isFail = verifyLinear(obj)
            % Checks parameters to verify if linear scan can be computed
            
            isFail = false;
            if obj.numX ~= 1
                warning("numX must be 1 for linear scan!")
                isFail = true;
            end
            if obj.numY < 1
                warning("numY must be an integer greater than 1 for linear scan!")
                isFail = true;
            end
            if obj.numTheta ~= 1
                warning("numTheta must be 1 for linear scan!")
                isFail = true;
            end
        end
        
        function isFail = verifyRectilinear(obj)
            % Checks parameters to verify if rectilinear scan can be
            % computed
            
            isFail = false;
            if obj.numX < 1
                warning("numX must be an integer greater than 1 for rectilinear scan!")
                isFail = true;
            end
            if obj.numY < 1
                warning("numY must be an integer greater than 1 for rectilinear scan!")
                isFail = true;
            end
            if obj.numTheta ~= 1
                warning("numTheta must be 1 for rectilinear scan!")
                isFail = true;
            end
        end
        
        function isFail = verifyCircular(obj)
            % Checks parameters to verify if circular scan can be computed
            
            isFail = false;
            if obj.numX ~= 1
                warning("numX must be 1 for circular scan!")
                isFail = true;
            end
            if obj.numY ~= 1
                warning("numY must be 1 for circular scan!")
                isFail = true;
            end
            if obj.numTheta < 1
                warning("numTheta must be an integer greater than 1 for circular scan!")
                isFail = true;
            end
        end
        
        function isFail = verifyCylindrical(obj)
            % Verify proper antenna and SAR scenario parameters for a
            % cylindrical scan
            
            isFail = false;
            if obj.numX ~= 1
                warning("numX must be 1 for cylindrical scan!")
                isFail = true;
            end
            if obj.numY < 1
                warning("numY must be an integer greater than 1 for cylindrical scan!")
                isFail = true;
            end
            if obj.numTheta < 1
                warning("numTheta must be an integer greater than 1 for cylindrical scan!")
                isFail = true;
            end
        end
        
        function attachListener(obj)
            % Attaches the listener to the object handle
            
            addlistener(obj,{'scanMethod','xStep_m','numX','yStep_m','numY','thetaMax_deg','numTheta','ant'},'PostSet',@sarScenario.propChange);
        end
    end
    
    methods(Static)
        function propChange(metaProp,eventData)
            % Recomputes the SAR antenna positions if watched properties are changed
            
            computeSarScenario(eventData.AffectedObject);
        end
    end
end