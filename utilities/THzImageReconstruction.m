% THzImageReconstruction object holds the properties and methods used for
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

classdef THzImageReconstruction < handle
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
        
        reconstructor = struct("isFail",false) % Reconstructor object: depends on which reconstruction algorithm is being used
        isGPU = true                % Boolean whether or not to use the GPU for image reconstruction (results vary depending on imaging scenario and parameters)
        isGPUVerified = false       % Boolean whether or not the GPU is has been verified
        isMult2Mono = true          % Boolean whether or not to use the multistatic-to-monostatic approximation
        zRef_m = 0.25               % z location of reference plane for multistatic-to-monostatic approximation
        zSlice_m                    % z slice of interest when reconstructing a 2-D x-y image, in meters
        thetaUpsampleFactor = 1     % Upsample factor on the theta dimension of the SAR data operated in circular or cylindrical mode
        
        isSilent = false            % Boolean whether or not to display updates to the user
        
        wav                         % A THzWaveformParameters object handle
        ant                         % A THzAntennaArray object handle
        scanner                     % A THzScanner object handle
        target                      % A THzTarget object handle
        
        app = struct("UIFigure",[]) % App handle
        isApp                       % Boolean whether or not to use the app functionality
    end
    
    methods
        function obj = THzImageReconstruction(wav,ant,scanner,target,app)
            % Sets wav, ant, scanner, target, and app properties
            
            obj.wav = wav;
            obj.ant = ant;
            obj.scanner = scanner;
            obj.target = target;
            if nargin == 4
                obj.isApp = false;
            elseif nargin == 5
                obj.isApp = true;
                obj.app = app;
            else
                showErrorMessage("Must input wav, ant, scanner, and target at least to THzImageReconstruction");
            end
            
            VerifyGPU(obj);
        end
        
        function obj = Update(obj)
            % Update the imaging parameters from the app and verify the
            % parameters
            
            obj = Get(obj);
            if obj.method == "-" || (obj.isApp && obj.app.ReconstructionAlgorithmDropDown.Value == "-")
                obj.reconstructor.isFail = true;
            end
            
            if obj.isApp && obj.reconstructor.isFail
                obj.app.ReconstructionAlgorithmDropDown.Value = "-";
            elseif obj.isApp
                obj.reconstructor.update();
                if obj.reconstructor.isFail
                    obj.app.ReconstructionAlgorithmDropDown.Value = "-";
                    obj.method = "-";
                end
            end
            
            if obj.isApp && obj.reconstructor.isFail
                obj.app.ReconstructionAlgorithmDropDown.Value = "-";
            elseif obj.reconstructor.isFail
                obj.method = "-";
            end
        end
        
        function obj = Compute(obj)
            % Attempt to reconstruct the image
            
            Get(obj);
            if obj.reconstructor.isFail
                showErrorMessage(obj,"Resolve errors before image reconstruction","Parameter error!");
                return;
            end
            
            if isempty(obj.target.sarData)
                showErrorMessage(obj,"Must compute beat signal before image reconstruction!","Beat Signal Error!");
                return;
            end
            
            obj = Update(obj);
            
            if ~obj.reconstructor.isFail && (obj.method ~= "-" || (obj.isApp && obj.app.ReconstructionAlgorithmDropDown.Value ~= "-"))
                
                if ~obj.isSilent
                    disp("Attempting image reconstruction using " + obj.method + " method.")
                end
                
                obj.imXYZ = obj.reconstructor.computeReconstruction();
                
                if obj.reconstructor.isFail
                    if obj.isApp
                        obj.app.ImageReconstructionCompleteLamp.Color = "red";
                    end
                    showErrorMessage(obj,"Reconstruction Failed","Reconstruction error");
                else
                    obj.app.ImageReconstructionCompleteLamp.Color = "green";
                    if ~obj.isSilent
                        disp("Successfully image reconstruction using " + obj.method + " method.")
                    end
                end
                
                if obj.reconstructor.isFail
                    disp(obj.method + " failed!!")
                end
            end
        end
        
        function obj = Get(obj)
            % Get the imaging parameters from the app, get the imaging
            % axes, and construct the image reconstruction algorithm object
            % for use later
            
            if obj.isApp
                obj.method = obj.app.ReconstructionAlgorithmDropDown.Value;
                
                obj.dBMin = obj.app.MindBEditField.Value;
                obj.fontSize = obj.app.FontSizeEditField.Value;
                
                obj.nFFTx = obj.app.NumberofXFFTPointsEditField.Value;
                obj.nFFTy = obj.app.NumberofYFFTPointsEditField.Value;
                obj.nFFTz = obj.app.NumberofZFFTPointsEditField.Value;
                
                obj.xMin_m = obj.app.XMinmEditField_im.Value;
                obj.xMax_m = obj.app.XMaxmEditField_im.Value;
                
                obj.yMin_m = obj.app.YMinmEditField_im.Value;
                obj.yMax_m = obj.app.YMaxmEditField_im.Value;
                
                obj.zMin_m = obj.app.ZMinmEditField_im.Value;
                obj.zMax_m = obj.app.ZMaxmEditField_im.Value;
                
                obj.numX = obj.app.NumXVoxelsEditField.Value;
                obj.numY = obj.app.NumYVoxelsEditField.Value;
                obj.numZ = obj.app.NumZVoxelsEditField.Value;
                
                obj.isGPU = obj.app.UseGPUCheckBox_2.Value;
            end
            
            obj = generateAxes(obj);
            
            switch obj.method
                case "Uniform 2-D SAR 3-D RMA"
                    obj.reconstructor = uniform_XY_SAR_XYZ_RMA(obj);
                    
                case "Uniform 1-D SAR 2-D RMA"
                    obj.reconstructor = uniform_Y_SAR_YZ_RMA(obj);
                    
                case "Uniform 2-D SAR 2-D FFT"
                    obj.reconstructor = uniform_XY_SAR_XY_FFT(obj);
                    
                case "Uniform 1-D SAR 1-D FFT"
                    obj.reconstructor = uniform_Y_SAR_Y_FFT(obj);
                    
                case "2-D SAR 3-D BPA"
                    obj.reconstructor = nonuniform_XY_SAR_XYZ_BPA(obj);
                    
                case "2-D SAR 2-D BPA"
                    obj.reconstructor = nonuniform_XY_SAR_XY_BPA(obj);
                    
                case "1-D SAR 2-D BPA"
                    obj.reconstructor = nonuniform_Y_SAR_YZ_BPA(obj);
                    
                case "Uniform 2-D CSAR 3-D PFA"
                    obj.reconstructor = uniform_thetaY_CSAR_XYZ_PFA(obj);
                    
                case "2-D CSAR 3-D BPA"
                    obj.reconstructor = nonuniform_thetaY_CSAR_XYZ_BPA(obj);
                    
                case "Uniform 1-D CSAR 2-D PFA"
                    obj.reconstructor = uniform_theta_CSAR_XZ_PFA(obj);
                    
                case "1-D CSAR 2-D BPA"
                    obj.reconstructor = nonuniform_theta_CSAR_XZ_BPA(obj);
                case "custom"
                    obj.reconstructor = reconstructionAlgorithmTemplate(obj);
                case "-"
                    
                otherwise
                    showErrorMessage(obj,"Incorrect reconstruction algorithm method! See documentation for list of allowed reconstruction algorithms or call ""custom"" to use your custom algorithm","Reconstruction error");
            end
        end
        
        function obj = generateAxes(obj)
            % Generate the imaging axes
            
            obj.x_m = single(linspace(obj.xMin_m,obj.xMax_m-(obj.xMax_m-obj.xMin_m)/obj.numX,obj.numX));
            obj.y_m = single(linspace(obj.yMin_m,obj.yMax_m-(obj.yMax_m-obj.yMin_m)/obj.numY,obj.numY));
            obj.z_m = single(linspace(obj.zMin_m,obj.zMax_m-(obj.zMax_m-obj.zMin_m)/obj.numZ,obj.numZ));
        end
        
        function obj = InitializeFigures(obj)
            % Initialize the figures
            
            CloseFigures(obj);
            set(0,'DefaultFigureWindowStyle','docked')
            
            obj.fig.f = figure(3);
            obj.fig.h = handle(axes);
        end
        
        function CloseFigures(obj)
            % Attempt to close the figures
            
            try
                close(obj.fig.f)
            catch
            end
        end
        
        function Display(obj)
            % Displays the reconstructed image in a figure
            
            if isempty(obj.fig.f) || ~isvalid(obj.fig.h)
                obj = InitializeFigures(obj);
            end
            
            if obj.isApp
                obj.dBMin = obj.app.MindBEditField.Value;
                obj.fontSize = obj.app.FontSizeEditField.Value;
            end
            
            obj.reconstructor.displayImage();
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
        
        function obj = VerifyGPU(obj)
            % Verifies if the GPU can be used
            
            % If the checkbox in the app is unchecked
            if obj.isApp && ~obj.app.UseGPUCheckBox.Value
                obj.isGPU = false;
                return;
            end
            
            if obj.isGPU && ~obj.isGPUVerified
                try
                    reset(gpuDevice);
                    obj.isGPUVerified = true;
                catch
                    obj.isGPU = false;
                    obj.isGPUVerified = false;
                    showErrorMessage(obj,"Unable to locate Nvidia GPU","No GPU available");
                    if obj.isApp
                        obj.app.UseGPUCheckBox.Value = false;
                    end
                    return;
                end
            end
        end
    end
end