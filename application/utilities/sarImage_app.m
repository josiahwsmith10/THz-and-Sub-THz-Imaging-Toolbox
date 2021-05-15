% sarImage_app object holds the properties and methods used for
% reconstructing the image from the simulated SAR scenario
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

classdef sarImage_app < handle
    properties
        nFFTx = 512                 % Number of FFT points along the x-dimension, when using FFT-based reconstruction algorithms
        nFFTy = 512                 % Number of FFT points along the y-dimension, when using FFT-based reconstruction algorithms
        nFFTz = 512                 % Number of FFT points along the z-dimension, when using FFT-based reconstruction algorithms
        
        numX = 128                  % Number of voxels in the reconstructed image along the x-dimension
        numY = 128                  % Number of voxels in the reconstructed image along the y-dimension
        numZ = 128                  % Number of voxels in the reconstructed image along the z-dimension
        
        x_m                         % Reconstructed image x axis
        y_m                         % Reconstructed image y axis
        z_m                         % Reconstructed image z axis
        
        xMin_m = -0.2               % Minimum value of reconstructed image along x-dimension
        xMax_m = 0.2                % Maximum value of reconstructed image along x-dimension
        
        yMin_m = -0.2               % Minimum value of reconstructed image along y-dimension
        yMax_m = 0.2                % Maximum value of reconstructed image along y-dimension
        
        zMin_m = 0                  % Minimum value of reconstructed image along z-dimension
        zMax_m = 0.5                % Maximum value of reconstructed image along z-dimension
        
        method = "-"                % Reconstruction algorithm to use
        
        imXYZ                       % Reconstructed image
        
        fig = struct('f',[],'h',[]) % Structure containing the figure and handle used for showing the target
        dBMin = -25                 % Minimum value displayed, in decibells
        fontSize = 12               % Font size of displayed image
        vSliceIndex                 % Slices of the 3-D image along the z-dimension to use (default: use all)
        
        reconstructor               % Reconstructor object: depends on which reconstruction algorithm is being used
        isGPU                       % Boolean whether or not to use the GPU for image reconstruction (results vary depending on imaging scenario and parameters)
    end
    methods
        function obj = sarImage_app(app)
            obj = initializeFigures(obj);
        end
        
        function obj = update(obj,app)
            % Update the imaging parameters from the app and verify the
            % parameters
            
            obj = getImagingParameters(obj,app);
            if app.ReconstructionAlgorithmDropDown.Value ~= "-"
                obj.reconstructor = obj.reconstructor.update(app,obj);
            else
                obj.reconstructor.isFail = true;
            end
            
            if obj.reconstructor.isFail
                app.ReconstructionAlgorithmDropDown.Value = "-";
            end
        end
        
        function obj = computeImage(obj,app)
            % Attempt to reconstruct the image
            
            if isempty(app.target.sarData)
                uiconfirm(app.UIFigure,"Must compute beat signal before image reconstruction!",'Beat Signal Error!',...
                    "Options",{'OK'},'Icon','warning');
            end
            
            obj = update(obj,app);
            
            if app.ReconstructionAlgorithmDropDown.Value ~= "-"
                [obj.reconstructor,obj.imXYZ] = obj.reconstructor.computeReconstruction(app,obj);
                
                displayImage(obj,app);
            end
            
            if isempty(obj.reconstructor.isFail) || obj.reconstructor.isFail
                app.ImageReconstructionCompleteLamp.Color = "red";
            else
                app.ImageReconstructionCompleteLamp.Color = "green";
            end
        end
        
        function obj = getImagingParameters(obj,app)
            % Get the imaging parameters from the app, get the imaging
            % axes, and construct the image reconstruction algorithm object
            % for use later
            
            obj.method = app.ReconstructionAlgorithmDropDown.Value;
            
            obj.dBMin = app.MindBEditField.Value;
            obj.fontSize = app.FontSizeEditField.Value;
            
            obj.nFFTx = app.NumberofXFFTPointsEditField.Value;
            obj.nFFTy = app.NumberofYFFTPointsEditField.Value;
            obj.nFFTz = app.NumberofZFFTPointsEditField.Value;
            
            obj.xMin_m = app.XMinmEditField_im.Value;
            obj.xMax_m = app.XMaxmEditField_im.Value;
            
            obj.yMin_m = app.YMinmEditField_im.Value;
            obj.yMax_m = app.YMaxmEditField_im.Value;
            
            obj.zMin_m = app.ZMinmEditField_im.Value;
            obj.zMax_m = app.ZMaxmEditField_im.Value;
            
            obj.numX = app.NumXVoxelsEditField.Value;
            obj.numY = app.NumYVoxelsEditField.Value;
            obj.numZ = app.NumZVoxelsEditField.Value;
            
            obj = generateAxes(obj);
            
            obj.isGPU = app.UseGPUCheckBox_2.Value;
            
            switch obj.method
                case "Uniform 2-D SAR 3-D RMA"
                    obj.reconstructor = uniform_XY_SAR_XYZ_RMA_app(app,obj);
                    
                case "Uniform 1-D SAR 2-D RMA"
                    obj.reconstructor = uniform_Y_SAR_YZ_RMA_app(app,obj);
                    
                case "Uniform 2-D SAR 2-D FFT"
                    obj.reconstructor = uniform_XY_SAR_XY_FFT_app(app,obj);
                    
                case "Uniform 1-D SAR 1-D FFT"
                    obj.reconstructor = uniform_Y_SAR_Y_FFT_app(app,obj);
                    
                case "2-D SAR 3-D BPA"
                    obj.reconstructor = nonuniform_XY_SAR_XYZ_BPA_app(app,obj);
                    
                case "2-D SAR 2-D BPA"
                    obj.reconstructor = nonuniform_XY_SAR_XY_BPA_app(app,obj);
                    
                case "Uniform 2-D CSAR 3-D PFA"
                    obj.reconstructor = uniform_thetaY_CSAR_XYZ_PFA_app(app,obj);
                    
                case "2-D CSAR 3-D BPA"
                    obj.reconstructor = nonuniform_thetaY_CSAR_XYZ_BPA_app(app,obj);
                    
                case "Uniform 1-D CSAR 2-D PFA"
                    obj.reconstructor = uniform_theta_CSAR_XZ_PFA_app(app,obj);
                    
                case "1-D CSAR 2-D BPA"
                    obj.reconstructor = nonuniform_theta_CSAR_XZ_BPA_app(app,obj);
                    
                case "Custom"
                    obj.reconstructor = reconstructionAlgorithmTemplate_app(app,obj);
            end
        end
        
        function obj = generateAxes(obj)
            % Generate the imaging axes
            
            obj.x_m = single(linspace(obj.xMin_m,obj.xMax_m-(obj.xMax_m-obj.xMin_m)/obj.numX,obj.numX));
            obj.y_m = single(linspace(obj.yMin_m,obj.yMax_m-(obj.yMax_m-obj.yMin_m)/obj.numY,obj.numY));
            obj.z_m = single(linspace(obj.zMin_m,obj.zMax_m-(obj.zMax_m-obj.zMin_m)/obj.numZ,obj.numZ));
        end
        
        function obj = initializeFigures(obj)
            % Initialize the figures
            
            closeFigures(obj);
            set(0,'DefaultFigureWindowStyle','docked')
            
            obj.fig.f = figure;
            obj.fig.h = handle(axes);
        end
        
        function closeFigures(obj)
            % Attempt to close the figures
            
            try
                close(obj.fig.f)
            catch
            end
        end
        
        function displayImage(obj,app)
            % Displays the reconstructed image
            
            obj.dBMin = app.MindBEditField.Value;
            obj.fontSize = app.FontSizeEditField.Value;
            
            obj.reconstructor.displayImage(obj);
        end
        
        function openInVolumeViewer(obj)
            % If the reconstructed image is 3-D, open in MATLAB
            % volumeViewer
            
            if ismatrix(obj.imXYZ)
                warning("Cannot open 2D image in volume viewer!");
            end
            
            try
                vvim = abs(permute(obj.imXYZ,[2,1,3]));
                volumeViewer(vvim);
            catch
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