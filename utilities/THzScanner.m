% THzScanner object holds the properties and methods used for simulating a
% THz imaging scanning scenario in four modes: linear, rectilinear
% (planar), circular, or cylindrical
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

classdef THzScanner < handle
    properties
        method                      % Scanning mode: "Linear", "Rectilinear", "Circular", or "Cylindrical"
        xStep_m = 1e-3              % Step size along the x-dimension to move the antenna array in meters
        yStep_m = 16e-3             % Step size along the y-dimension to move the antenna array in meters
        thetaMax_deg = 360          % Maximum angle for circular or cylindrical scans in degrees
        numX = 1                    % Number of steps along the x-dimension
        numY = 1                    % Number of steps along the y-dimension
        numTheta = 1                % Number of angular steps
        
        X_m                         % Temporary storage of the entire aperture's x-coordinates
        Y_m                         % Temporary storage of the entire aperture's y-coordinates
        Z_m                         % Temporary storage of the entire aperture's z-coordinates
        xyz_m                       % Location in x-y-z of the synthetic aperture, an Nx3 array, where N is the total number of synthetic elements and each row contains the [x,y,z] coordinates of each element
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
        ant                         % A THzAntennaArray handle
        
        app = struct("UIFigure",[]) % App handle
        isApp                       % Boolean whether or not to use the app functionality
    end
    
    methods
        function obj = THzScanner(ant,app)
            % Sets up the ant and app properties
            
            obj.ant = ant;
            if nargin == 1
                obj.isApp = false;
            elseif nargin == 2
                obj.isApp = true;
                obj.app = app;
            end
        end
        
        function obj = Update(obj)
            % Update the scanner and displays in the figures
            
            if obj.isApp
                obj = Get(obj);
            end
            obj = Compute(obj);
            
            Display(obj);
        end
        
        function obj = Get(obj)
            % Get parameters of THz scanning pattern from the app
            
            if ~obj.isApp
                warning("Must be using app to call Get method");
                return;
            end
            
            obj.method = obj.app.SARMethodDropDown.Value;
            
            % Update step sizes and theta aperture size
            if obj.app.XStepSwitch.Value == "λ"
                obj.xStep_m = obj.app.XStepSizeEditField.Value * obj.app.wav.lambda_m;
            elseif obj.app.XStepSwitch.Value == "mm"
                obj.xStep_m = obj.app.XStepSizeEditField.Value * 1e-3;
            end
            
            if obj.app.YStepSwitch.Value == "λ"
                obj.yStep_m = obj.app.YStepSizeEditField.Value * obj.app.wav.lambda_m;
            elseif obj.app.YStepSwitch.Value == "mm"
                obj.yStep_m = obj.app.YStepSizeEditField.Value * 1e-3;
            end
            
            obj.thetaMax_deg = obj.app.ThetaSizedegEditField.Value;
            
            obj.numX = obj.app.NumXStepsEditField.Value;
            obj.numY = obj.app.NumYStepsEditField.Value;
            obj.numTheta = obj.app.NumThetaStepsEditField.Value;
            
            obj.app.XSizemEditField.Value = (obj.numX-1) * obj.xStep_m;
            obj.app.YSizemEditField.Value = (obj.numY-1) * obj.yStep_m;
        end
        
        function obj = Compute(obj)
            % Compute positions of the THz scanning synthetic array
            
            switch obj.method
                case "Linear"
                    if verifyLinearity(obj) && verifyLinear(obj)
                        obj = ComputeLinear(obj);
                    else
                        return;
                    end
                case "Rectilinear"
                    if verifyLinearity(obj) && verifyRectilinear(obj)
                        obj = ComputeRectilinear(obj);
                    else
                        return;
                    end
                case "Circular"
                    if verifyCircular(obj)
                        obj = ComputeCircular(obj);
                    else
                        return;
                    end
                case "Cylindrical"
                    if verifyLinearity(obj) && verifyCylindrical(obj)
                        obj = ComputeCylindrical(obj);
                    else
                        return;
                    end
            end
            
            obj.tx.xyz_m = single(obj.tx.xyz_m);
            obj.rx.xyz_m = single(obj.rx.xyz_m);
            obj.vx.xyz_m = single(obj.vx.xyz_m);
            
            % Set SAR method
            if obj.isApp
                obj.method = obj.app.SARMethodDropDown.Value;
            end
        end
                
        function obj = ComputeLinear(obj)
            % Computes the linear scan type
            
            if obj.isApp
                % Enable necessary fields
                obj.app.XStepSizeEditField.Enable = false;
                obj.app.YStepSizeEditField.Enable = true;
                obj.app.ThetaSizedegEditField.Enable = false;
                obj.app.NumXStepsEditField.Enable = false;
                obj.app.NumYStepsEditField.Enable = true;
                obj.app.NumThetaStepsEditField.Enable = false;
                obj = Get(obj);
            end
            
            obj = getsarAxes(obj);
            
            [obj.X_m,obj.Y_m,obj.Z_m] = ndgrid(obj.x_m,obj.y_m,obj.z_m);
            obj.xyz_m = reshape(cat(3,obj.X_m,obj.Y_m,obj.Z_m),1,[],3);
            obj.xyz_m = repmat(obj.xyz_m,obj.ant.vx.numVx,1,1);
            
            obj.tx.xyz_m = single(obj.xyz_m + obj.ant.tx.xyz_m);
            obj.rx.xyz_m = single(obj.xyz_m + obj.ant.rx.xyz_m);
            obj.vx.xyz_m = single(obj.xyz_m + obj.ant.vx.xyz_m);
            
            obj.tx.xyz_m = reshape(obj.tx.xyz_m,[],3);
            obj.rx.xyz_m = reshape(obj.rx.xyz_m,[],3);
            obj.vx.xyz_m = reshape(obj.vx.xyz_m,[],3);
            
            if isAppEPC(obj.app) || obj.ant.isEPC
                obj.sarSize = obj.ant.vx.numVx*obj.numY;
            else
                % Unwrap obj.tx.xyz_m & obj.rx.xyz_m as [numRx,numTx,numY,3]
                obj.sarSize = [obj.ant.rx.numRx,obj.ant.tx.numTx,obj.numY];
            end
        end
        
        function obj = ComputeRectilinear(obj)
            % Computes the rectilinear scan type
            
            % Enable necessary fields
            if obj.isApp
                obj.app.XStepSizeEditField.Enable = true;
                obj.app.YStepSizeEditField.Enable = true;
                obj.app.ThetaSizedegEditField.Enable = false;
                obj.app.NumXStepsEditField.Enable = true;
                obj.app.NumYStepsEditField.Enable = true;
                obj.app.NumThetaStepsEditField.Enable = false;
                obj = Get(obj);
            end
            
            obj = getsarAxes(obj);
            
            [obj.Y_m,obj.X_m,obj.Z_m] = ndgrid(obj.y_m,obj.x_m,obj.z_m);
            obj.xyz_m = reshape(cat(3,obj.X_m,obj.Y_m,obj.Z_m),1,[],3);
            obj.xyz_m = repmat(obj.xyz_m,obj.ant.vx.numVx,1,1);
            
            obj.tx.xyz_m = single(obj.xyz_m + obj.ant.tx.xyz_m);
            obj.rx.xyz_m = single(obj.xyz_m + obj.ant.rx.xyz_m);
            obj.vx.xyz_m = single(obj.xyz_m + obj.ant.vx.xyz_m);
            
            obj.tx.xyz_m = reshape(obj.tx.xyz_m,[],3);
            obj.rx.xyz_m = reshape(obj.rx.xyz_m,[],3);
            obj.vx.xyz_m = reshape(obj.vx.xyz_m,[],3);
            
            if isAppEPC(obj.app) || obj.ant.isEPC
                obj.sarSize = [obj.ant.vx.numVx*obj.numY,obj.numX];
            else
                % Unwrap obj.tx.xyz_m & obj.rx.xyz_m as [numRx,numTx,numY,numX,3]
                obj.sarSize = [obj.ant.rx.numRx,obj.ant.tx.numTx,obj.numY,obj.numX];
            end
        end
        
        function obj = ComputeCircular(obj)
            % Computes the circular scan type
            
            % Enable necessary fields
            if obj.isApp
                obj.app.XStepSizeEditField.Enable = false;
                obj.app.YStepSizeEditField.Enable = false;
                obj.app.ThetaSizedegEditField.Enable = true;
                obj.app.NumXStepsEditField.Enable = false;
                obj.app.NumYStepsEditField.Enable = false;
                obj.app.NumThetaStepsEditField.Enable = true;
                obj = Get(obj);
            end
            
            obj = getsarAxes(obj);
            
            obj.x_m = obj.ant.tx.z0_m*cos(obj.theta_rad);
            obj.y_m = zeros(size(obj.theta_rad));
            obj.z_m = obj.ant.tx.z0_m*sin(obj.theta_rad) - obj.ant.tx.z0_m;
            
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
        
        function obj = ComputeCylindrical(obj)
            % Computes the cylindrical scan type
            
            % Enable necessary fields
            if obj.isApp
                obj.app.XStepSizeEditField.Enable = false;
                obj.app.YStepSizeEditField.Enable = true;
                obj.app.ThetaSizedegEditField.Enable = true;
                obj.app.NumXStepsEditField.Enable = false;
                obj.app.NumYStepsEditField.Enable = true;
                obj.app.NumThetaStepsEditField.Enable = true;
                obj = Get(obj);
            end
            
            obj = getsarAxes(obj);
            
            obj.x_m = obj.ant.tx.z0_m*cos(obj.theta_rad);
            obj.z_m = obj.ant.tx.z0_m*sin(obj.theta_rad) - obj.ant.tx.z0_m;
            
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
            
            if isAppEPC(obj.app) || obj.ant.isEPC
                obj.sarSize = [obj.ant.vx.numVx*obj.numY,obj.numTheta];
            else
                % Unwrap obj.tx.xyz_m & obj.rx.xyz_m as [numRx,numTx,numY,numTheta,3]
                obj.sarSize = [obj.ant.rx.numRx,obj.ant.tx.numTx,obj.numY,obj.numTheta];
            end
        end
        
        function obj = getsarAxes(obj)
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
        
        function obj = InitializeFigures(obj)
            % Initilizes the figures
            
            CloseFigures(obj)
            set(0,'DefaultFigureWindowStyle','docked')
            
            if obj.isApp
                try
                    obj.fig.f = obj.app.ant.fig.f;
                    obj.fig.anth = obj.app.ant.fig.anth;
                    obj.fig.sarh = obj.app.ant.fig.sarh;
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
        
        function tf = verifyLinear(obj)
            % Verifies the linear scan type
            
            tf = true;
            if obj.numX ~= 1
                showErrorMessage(obj,"numX must be 1 for linear scan!","Linear scan error")
                obj.method = "-";
                tf = false;
            end
            if obj.numY < 1
                showErrorMessage(obj,"numY must be an integer greater than 1 for linear scan!","Linear scan error")
                obj.method = "-";
                tf = false;
            end
            if obj.numTheta ~= 1
                showErrorMessage(obj,"numTheta must be 1 for linear scan!","Linear scan error")
                obj.method = "-";
                tf = false;
            end
        end
        
        function tf = verifyRectilinear(obj)
            % Verifies the rectilinear scan type
            
            tf = true;
            if obj.numX < 1
                showErrorMessage(obj,"numX must be an integer greater than 1 for rectilinear scan!","Rectilinear scan error")
                obj.method = "-";
                tf = false;
            end
            if obj.numY < 1
                showErrorMessage(obj,"numY must be an integer greater than 1 for rectilinear scan!","Rectilinear scan error")
                obj.method = "-";
                tf = false;
            end
            if obj.numTheta ~= 1
                showErrorMessage(obj,"numTheta must be 1 for rectilinear scan!","Rectilinear scan error")
                obj.method = "-";
                tf = false;
            end
        end
        
        function tf = verifyCircular(obj)
            % Verifies the circular scan type
            
            tf = true;
            % Verify single element array
            if obj.ant.tx.numTx ~= 1 || obj.ant.rx.numRx ~= 1
                showErrorMessage(obj,"Array must have only 1 Tx and 1 Rx. Please disable necessary elements. Setting method to ""-""","Array topology error");
                obj.app.SARMethodDropDown.Value = "-";
                tf = false;
            end
            
            if obj.numX ~= 1
                showErrorMessage(obj,"numX must be 1 for circular scan!","Circular scan error")
                obj.method = "-";
                tf = false;
            end
            if obj.numY ~= 1
                showErrorMessage(obj,"numY must be 1 for circular scan!","Circular scan error")
                obj.method = "-";
                tf = false;
            end
            if obj.numTheta < 1
                showErrorMessage(obj,"numTheta must be an integer greater than 1 for circular scan!","Circular scan error")
                obj.method = "-";
                tf = false;
            end
        end
        
        function tf = verifyCylindrical(obj)
            % Verifies the cylindrical scan type
            
            tf = true;
            if obj.numX ~= 1
                showErrorMessage(obj,"numX must be 1 for cylindrical scan!","Cylindrical scan error")
                obj.method = "-";
                tf = false;
            end
            if obj.numY < 1
                showErrorMessage(obj,"numY must be an integer greater than 1 for cylindrical scan!","Cylindrical scan error")
                obj.method = "-";
                tf = false;
            end
            if obj.numTheta < 1
                showErrorMessage(obj,"numTheta must be an integer greater than 1 for cylindrical scan!","Cylindrical scan error")
                obj.method = "-";
                tf = false;
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
            
            if isAppEPC(obj.app) || obj.ant.isEPC
                DisplayEPC(obj);
            else
                DisplayMIMO(obj);
            end
        end
        
        function DisplayMIMO(obj)
            % Plots the MIMO synthetic array
            
            if obj.method == "-"
                return;
            end
            
            if obj.isApp
                h = obj.fig.sarh;
            else
                h = obj.fig.h;
            end
            hold(h,'off')
            temp = obj.tx.xyz_m;
            scatter3(h,temp(:,1),temp(:,3),temp(:,2),'.r')
            hold(h,'on')
            temp = obj.rx.xyz_m;
            scatter3(h,temp(:,1),temp(:,3),temp(:,2),'.b')
            temp1 = obj.tx.xyz_m(:,1);
            temp2 = obj.rx.xyz_m(:,1);
            xlim(h,[min(min(temp1),min(temp2))-0.01,max(max(temp1),max(temp2))+0.01])
            temp1 = obj.tx.xyz_m(:,3);
            temp2 = obj.rx.xyz_m(:,3);
            ylim(h,[min(min(temp1),min(temp2))-0.01,max(max(temp1),max(temp2))+0.01])
            temp1 = obj.tx.xyz_m(:,2);
            temp2 = obj.rx.xyz_m(:,2);
            zlim(h,[min(min(temp1),min(temp2))-0.01,max(max(temp1),max(temp2))+0.01])
            xlabel(h,"x (m)",'interpreter','latex')
            ylabel(h,"z (m)",'interpreter','latex')
            zlabel(h,"y (m)",'interpreter','latex')
            legend(h,"Tx","Rx",'interpreter','latex')
            title(h,"MIMO Synthetic Aperture",'interpreter','latex')
            view(h,3)
            daspect(h,[1 1 1])
        end
        
        function DisplayEPC(obj)
            % Plots the equivalent-phase-center (EPC) synthetic array
            
            if obj.method == "-"
                return;
            end
            
            if obj.isApp
                h = obj.fig.sarh;
            else
                h = obj.fig.h;
            end
            
            hold(h,'off')
            temp = obj.vx.xyz_m;
            scatter3(h,temp(:,1),temp(:,3),temp(:,2),'.k')
            temp1 = obj.tx.xyz_m(:,1);
            temp2 = obj.rx.xyz_m(:,1);
            xlim(h,[min(min(temp1),min(temp2))-0.01,max(max(temp1),max(temp2))+0.01])
            temp1 = obj.tx.xyz_m(:,3);
            temp2 = obj.rx.xyz_m(:,3);
            ylim(h,[min(min(temp1),min(temp2))-0.01,max(max(temp1),max(temp2))+0.01])
            temp1 = obj.tx.xyz_m(:,2);
            temp2 = obj.rx.xyz_m(:,2);
            zlim(h,[min(min(temp1),min(temp2))-0.01,max(max(temp1),max(temp2))+0.01])
            xlabel(h,"x (m)",'interpreter','latex')
            ylabel(h,"z (m)",'interpreter','latex')
            zlabel(h,"y (m)",'interpreter','latex')
            legend(h,"Vx",'interpreter','latex')
            title(h,"Virtual Synthetic Aperture",'interpreter','latex')
            view(h,3)
            daspect(h,[1 1 1])
        end
        
        function tf = verifyLinearity(obj)
            % Verifies that the MIMO antenna array is colinear
            %   tf = true - antenna array is colinear
            %   tf = false - antenna array is not colinear
            
            if max(diff([obj.ant.tx.xy_m(:,1);obj.ant.rx.xy_m(:,1)])) > 8*eps
                showErrorMessage(obj,"MIMO array must be colinear. Please disable necessary elements. Setting method to ""-""","Array topology error");
                obj.app.SARMethodDropDown.Value = "-";
                obj.method = "-";
                tf = false;
            else
                tf = true;
            end
        end
    end
end