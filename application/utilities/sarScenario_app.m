% sarScenario_app object holds the properties and methods
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

classdef sarScenario_app < handle
    properties
        method                      % Scanning mode: "Linear", "Rectilinear", "Circular", or "Cylindrical"
        xStep_m                     % Step size along the x-dimension to move the antenna array in meters
        yStep_m                     % Step size along the y-dimension to move the antenna array in meters
        thetaMax_deg                % Maximum angle for circular or cylindrical scans in degrees
        numX                        % Number of steps along the x-dimension
        numY                        % Number of steps along the y-dimension
        numTheta                    % Number of angular steps
        
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
        
        isMIMO                      % Boolean whether or not to use the MIMO physical element locations instead of the equivalent phase center virtual element locations
    end
    
    methods
        function obj = sarScenario_app(app)
            % Initialize the figures and update the SAR scenario
            
            obj = initializeFigures(obj,app);
            obj = update(obj,app);
        end
        
        function obj = update(obj,app)
            % Update the SAR scenario by getting values from app and
            % recomputing SAR scenario
            
            obj = getSarScenario(obj,app);
            obj = computeSarScenario(obj,app);
            if verifyMIMO(obj,app)
                displaySarScenario(obj);
            else
                displayVirtualSarScenario(obj);
            end
        end
        
        function obj = getSarScenario(obj,app)
            % Get parameters of SAR scenario from app
            
            obj.isMIMO = verifyMIMO(obj,app);
            
            obj.method = app.SARMethodDropDown.Value;
            
            % Update step sizes and theta aperture size
            if app.XStepSwitch.Value == "λ"
                obj.xStep_m = app.XStepSizeEditField.Value * app.fmcw.lambda_m;
            elseif app.XStepSwitch.Value == "mm"
                obj.xStep_m = app.XStepSizeEditField.Value * 1e-3;
            end
            
            if app.YStepSwitch.Value == "λ"
                obj.yStep_m = app.YStepSizeEditField.Value * app.fmcw.lambda_m;
            elseif app.YStepSwitch.Value == "mm"
                obj.yStep_m = app.YStepSizeEditField.Value * 1e-3;
            end
            
            obj.thetaMax_deg = app.ThetaSizedegEditField.Value;
            
            obj.numX = app.NumXStepsEditField.Value;
            obj.numY = app.NumYStepsEditField.Value;
            obj.numTheta = app.NumThetaStepsEditField.Value;
            
            app.XSizemEditField.Value = obj.numX * obj.xStep_m;
            app.YSizemEditField.Value = obj.numY * obj.yStep_m;
        end
        
        function obj = computeSarScenario(obj,app)
            % Compute the SAR scenario, namely the positions of the SAR
            % aperture
            
            switch obj.method
                case "Linear"
                    % Enable necessary fields
                    app.XStepSizeEditField.Enable = false;
                    app.YStepSizeEditField.Enable = true;
                    app.ThetaSizedegEditField.Enable = false;
                    app.NumXStepsEditField.Enable = false;
                    app.NumYStepsEditField.Enable = true;
                    app.NumThetaStepsEditField.Enable = false;
                    
                    % Verify linearity of MIMO array
                    if verifyLinearity(obj,app)
                        return
                    end
                    
                    obj = getSarScenario(obj,app);
                    obj = getsarAxes(obj);
                    
                    [obj.X_m,obj.Y_m,obj.Z_m] = ndgrid(obj.x_m,obj.y_m,obj.z_m);
                    obj.xyz_m = reshape(cat(3,obj.X_m,obj.Y_m,obj.Z_m),1,[],3);
                    obj.xyz_m = repmat(obj.xyz_m,app.ant.vx.numVx,1,1);
                    
                    obj.tx.xyz_m = single(obj.xyz_m + app.ant.tx.xyz_m);
                    obj.rx.xyz_m = single(obj.xyz_m + app.ant.rx.xyz_m);
                    obj.vx.xyz_m = single(obj.xyz_m + app.ant.vx.xyz_m);
                    
                    obj.tx.xyz_m = reshape(obj.tx.xyz_m,[],3);
                    obj.rx.xyz_m = reshape(obj.rx.xyz_m,[],3);
                    obj.vx.xyz_m = reshape(obj.vx.xyz_m,[],3);
                    
                    if obj.isMIMO
                        % Unwrap obj.tx.xyz_m & obj.rx.xyz_m as [numRx,numTx,numY,3]
                        obj.sarSize = [app.ant.rx.numRx,app.ant.tx.numTx,obj.numY];
                    else
                        obj.sarSize = app.ant.vx.numVx*obj.numY;
                    end
                    
                case "Rectilinear"
                    % Enable necessary fields
                    app.XStepSizeEditField.Enable = true;
                    app.YStepSizeEditField.Enable = true;
                    app.ThetaSizedegEditField.Enable = false;
                    app.NumXStepsEditField.Enable = true;
                    app.NumYStepsEditField.Enable = true;
                    app.NumThetaStepsEditField.Enable = false;
                    
                    obj = getSarScenario(obj,app);
                    obj = getsarAxes(obj);
                    
                    [obj.Y_m,obj.X_m,obj.Z_m] = ndgrid(obj.y_m,obj.x_m,obj.z_m);
                    obj.xyz_m = reshape(cat(3,obj.X_m,obj.Y_m,obj.Z_m),1,[],3);
                    obj.xyz_m = repmat(obj.xyz_m,app.ant.vx.numVx,1,1);
                    
                    obj.tx.xyz_m = single(obj.xyz_m + app.ant.tx.xyz_m);
                    obj.rx.xyz_m = single(obj.xyz_m + app.ant.rx.xyz_m);
                    obj.vx.xyz_m = single(obj.xyz_m + app.ant.vx.xyz_m);
                    
                    obj.tx.xyz_m = reshape(obj.tx.xyz_m,[],3);
                    obj.rx.xyz_m = reshape(obj.rx.xyz_m,[],3);
                    obj.vx.xyz_m = reshape(obj.vx.xyz_m,[],3);
                    
                    if obj.isMIMO
                        % Unwrap obj.tx.xyz_m & obj.rx.xyz_m as [numRx,numTx,numY,numX,3]
                        obj.sarSize = [app.ant.rx.numRx,app.ant.tx.numTx,obj.numY,obj.numX];
                    else
                        % Verify linearity of MIMO array
                        if verifyLinearity(obj,app)
                            return
                        end
                        obj.sarSize = [app.ant.vx.numVx*obj.numY,obj.numX];
                    end
                    
                case "Circular"
                    % Enable necessary fields
                    app.XStepSizeEditField.Enable = false;
                    app.YStepSizeEditField.Enable = false;
                    app.ThetaSizedegEditField.Enable = true;
                    app.NumXStepsEditField.Enable = false;
                    app.NumYStepsEditField.Enable = false;
                    app.NumThetaStepsEditField.Enable = true;
                    
                    % Verify single element array
                    if app.ant.tx.numTx ~= 1 || app.ant.rx.numRx ~= 1
                        uiconfirm(app.UIFigure,"Array must have only 1 Tx and 1 Rx. Please disable necessary elements.",'Array Topology Error!',...
                            "Options",{'OK'},'Icon','warning');
                        app.SARMethodDropDown.Value = "-";
                        return
                    end
                    
                    obj = getSarScenario(obj,app);
                    obj = getsarAxes(obj);
                    
                    obj.x_m = app.ant.tx.z0_m*cos(obj.theta_rad);
                    obj.y_m = zeros(size(obj.theta_rad));
                    obj.z_m = app.ant.tx.z0_m*sin(obj.theta_rad) - app.ant.tx.z0_m;
                    
                    obj.xyz_m = reshape([obj.x_m(:),obj.y_m(:),obj.z_m(:)],1,[],3);
                    
                    obj.tx.xyz_m = single(obj.xyz_m + app.ant.tx.xyz_m);
                    obj.rx.xyz_m = single(obj.xyz_m + app.ant.rx.xyz_m);
                    obj.vx.xyz_m = single(obj.xyz_m + app.ant.vx.xyz_m);
                    
                    obj.tx.xyz_m = reshape(obj.tx.xyz_m,[],3);
                    obj.rx.xyz_m = reshape(obj.rx.xyz_m,[],3);
                    obj.vx.xyz_m = reshape(obj.vx.xyz_m,[],3);
                    
                    % Unwrap obj.tx.xyz_m & obj.rx.xyz_m as [numTheta,3]
                    obj.sarSize = obj.numTheta;
                    
                case "Cylindrical"
                    % Enable necessary fields
                    app.XStepSizeEditField.Enable = false;
                    app.YStepSizeEditField.Enable = true;
                    app.ThetaSizedegEditField.Enable = true;
                    app.NumXStepsEditField.Enable = false;
                    app.NumYStepsEditField.Enable = true;
                    app.NumThetaStepsEditField.Enable = true;
                    
                    % Verify linearity of MIMO array
                    if verifyLinearity(obj,app)
                        return
                    end
                    
                    obj = getSarScenario(obj,app);
                    obj = getsarAxes(obj);
                    
                    obj.x_m = app.ant.tx.z0_m*cos(obj.theta_rad);
                    obj.z_m = app.ant.tx.z0_m*sin(obj.theta_rad) - app.ant.tx.z0_m;
                    
                    % Use y' as first dimension -> so we obtain s(y',theta)
                    obj.X_m = repmat(obj.x_m(:).',obj.numY,1);
                    obj.Z_m = repmat(obj.z_m(:).',obj.numY,1);
                    obj.Y_m = repmat(obj.y_m(:),1,obj.numTheta);
                    
                    obj.xyz_m = reshape(cat(3,obj.X_m,obj.Y_m,obj.Z_m),1,[],3);
                    
                    obj.tx.xyz_m = single(obj.xyz_m + app.ant.tx.xyz_m);
                    obj.rx.xyz_m = single(obj.xyz_m + app.ant.rx.xyz_m);
                    obj.vx.xyz_m = single(obj.xyz_m + app.ant.vx.xyz_m);
                    
                    obj.tx.xyz_m = reshape(obj.tx.xyz_m,[],3);
                    obj.rx.xyz_m = reshape(obj.rx.xyz_m,[],3);
                    obj.vx.xyz_m = reshape(obj.vx.xyz_m,[],3);
                    
                    if obj.isMIMO
                        % Unwrap obj.tx.xyz_m & obj.rx.xyz_m as [numRx,numTx,numY,numTheta,3]
                        obj.sarSize = [app.ant.rx.numRx,app.ant.tx.numTx,obj.numY,obj.numTheta];
                    else
                        obj.sarSize = [app.ant.vx.numVx*obj.numY,obj.numTheta];
                    end
            end
            
            obj.tx.xyz_m = single(obj.tx.xyz_m);
            obj.rx.xyz_m = single(obj.rx.xyz_m);
            obj.vx.xyz_m = single(obj.vx.xyz_m);
            
            % Set SAR method
            obj.method = app.SARMethodDropDown.Value;
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
        
        function obj = initializeFigures(obj,app)
            % Initilizes the figures
            
            closeFigures(obj)
            set(0,'DefaultFigureWindowStyle','docked')
            
            try
                obj.fig.f = app.ant.fig.f;
                obj.fig.anth = app.ant.fig.anth;
                obj.fig.sarh = app.ant.fig.sarh;
                if ~(ishandle(obj.fig.anth) && ishandle(obj.fig.sarh))
                    error("error");
                end
            catch
                obj.fig.f = figure;
                obj.fig.anth = handle(subplot(121));
                obj.fig.sarh = handle(subplot(122));
            end
        end
        
        function closeFigures(obj)
            % Attempts to close the figures
            
            try
                close(obj.fig.f)
            catch
            end
        end
        
        function displaySarScenario(obj)
            % Plots the MIMO-SAR array
            
            if obj.method == "-"
                return;
            end
            
            h = obj.fig.sarh;
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
        
        function displayVirtualSarScenario(obj)
            % Plots the equivalent-phase-center (EPC) SAR array
            
            if obj.method == "-"
                return;
            end
            
            h = obj.fig.sarh;
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
        
        function failtf = verifyLinearity(obj,app)
            % Verifies that the MIMO antenna array is colinear
            
            if max(diff([app.ant.tx.xy_m(:,1);app.ant.rx.xy_m(:,1)])) > 8*eps
                uiconfirm(app.UIFigure,"MIMO array must be colinear. Please disable necessary elements.",'Array Topology Error!',...
                    "Options",{'OK'},'Icon','warning');
                app.SARMethodDropDown.Value = "-";
                failtf = true;
            else
                failtf = false;
            end
        end
        
        function tf = verifyMIMO(obj,app)
            % Verifies if the user has specified MIMO or EPC in the app
            
            if app.MIMOSwitch.Value == "Use MIMO Array"
                tf = true;
            elseif app.MIMOSwitch.Value == "Use EPC Virtual Elements"
                tf = false;
            end
        end
    end
end