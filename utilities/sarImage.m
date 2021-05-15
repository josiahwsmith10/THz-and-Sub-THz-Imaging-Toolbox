% sarImage object holds the properties and methods used for
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

classdef sarImage < handle
    properties(SetObservable)
        fmcw                        % An fmcwChirpParameters object
        ant                         % A sarAntennaArray object
        sar                         % A sarScenario object
        target                      % A sarTarget object
        
        method = "-"                % Reconstruction algorithm to use
    end
    
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
        
        imXYZ                       % Reconstructed image
        
        fig = struct('f',[],'h',[]) % Structure containing the figure and handle used for showing the target
        dBMin = -25                 % Minimum value displayed, in decibells
        fontSize = 12               % Font size of displayed image
        vSliceIndex                 % Slices of the 3-D image along the z-dimension to use (default: use all)
        
        reconstructor = struct("isFail",false) % Reconstructor object: depends on which reconstruction algorithm is being used
        isGPU = true                % Boolean whether or not to use the GPU for image reconstruction (results vary depending on imaging scenario and parameters)
        isGPUVerified = false       % Boolean whether or not the GPU is has been verified
        isMult2Mono = false         % Boolean whether or not to use the multistatic-to-monostatic approximation
        zRef_m = 0.25               % z location of reference plane for multistatic-to-monostatic approximation
        zSlice_m                    % z slice of interest when reconstructing a 2-D x-y image, in meters
        thetaUpsampleFactor = 1     % Upsample factor on the theta dimension of the SAR data operated in circular or cylindrical mode
        
        isSilent = false            % Boolean whether or not to display updates to the user
    end
    
    methods
        function obj = sarImage(fmcw,ant,sar,target)
            % Attaches the listener to the sarTarget object so
            % observable properties can be watched for changes and sets the
            % necessary parameters from the other SAR and FMCW type
            % objects. Then verifies whether the GPU can be used, if
            % necessary
            
            attachListener(obj);
            obj.fmcw = fmcw;
            obj.ant = ant;
            obj.sar = sar;
            obj.target = target;
            
            verifyGPU(obj);
        end
        
        function update(obj)
            % Update the imaging parameters
            
            getImagingParameters(obj);
            if obj.method == "-"
                obj.reconstructor.isFail = true;
            end
        end
        
        function computeImage(obj)
            % Attempt to reconstruct the image
            
            if isempty(obj.target.sarData)
                warning("Must compute beat signal before image reconstruction!")
                return
            end
            
            %update(obj);
            
            if obj.method ~= "-" && ~obj.reconstructor.isFail
                
                if ~obj.isSilent
                    disp("Attempting image reconstruction using " + obj.method + " method.")
                end
                
                obj.imXYZ = obj.reconstructor.computeReconstruction();
                
                if ~obj.isSilent
                    disp("Done reconstructing image using " + obj.method + " method.")
                end
                
                if obj.reconstructor.isFail
                    disp(obj.method + " failed!!")
                end
            end
        end
        
        function getImagingParameters(obj)
            % Get the imaging axes and construct the image reconstruction
            % algorithm object for use later. During the construction
            % process, the imaging parameters are verified
            
            generateAxes(obj);
            
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
                    
                case "Uniform 2-D CSAR 3-D PFA"
                    obj.reconstructor = uniform_thetaY_CSAR_XYZ_PFA(obj);
                    
                case "2-D CSAR 3-D BPA"
                    obj.reconstructor = nonuniform_thetaY_CSAR_XYZ_BPA(obj);
                    
                case "Uniform 1-D CSAR 2-D PFA"
                    obj.reconstructor = uniform_theta_CSAR_XZ_PFA(obj);
                    
                case "1-D CSAR 2-D BPA"
                    obj.reconstructor = nonuniform_theta_CSAR_XZ_BPA(obj);                    
                case "custom"
%                     obj.reconstructor = reconstructionAlgorithmTemplate(obj);
                    obj.reconstructor = freehand_linear_RMA_jws(obj);
                case "-"
                    
                otherwise
                    error("Incorrect reconstruction algorithm method! See documentation for list of allowed reconstruction algorithms or call ""custom"" to use your custom algorithm");
            end
        end
        
        function generateAxes(obj)
            % Generate the imaging axes
            
            obj.x_m = single(linspace(obj.xMin_m,obj.xMax_m-(obj.xMax_m-obj.xMin_m)/obj.numX,obj.numX));
            obj.y_m = single(linspace(obj.yMin_m,obj.yMax_m-(obj.yMax_m-obj.yMin_m)/obj.numY,obj.numY));
            obj.z_m = single(linspace(obj.zMin_m,obj.zMax_m-(obj.zMax_m-obj.zMin_m)/obj.numZ,obj.numZ));
        end
        
        function initializeFigures(obj)
            % Initializes the figures
            
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
        
        function displayImage(obj)
            % Displays the reconstructed image
            
            if isempty(obj.fig.f) || ~isvalid(obj.fig.h)
                initializeFigures(obj);
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
        
        function attachListener(obj)
            % Attaches the listener to the object handle
            
            addlistener(obj,{'method','fmcw','ant','sar','target'},'PostSet',@sarImage.propChange);
        end
        
        function verifyGPU(obj)
            % Verifies if the GPU can be used
            
            if obj.isGPU && ~obj.isGPUVerified
                try
                    reset(gpuDevice);
                catch
                    obj.isGPU = false;
                    warning("Unable to locate Nvidia GPU");
                    return;
                end
                obj.isGPUVerified = true;
            end
        end
    end
    
    methods(Static)
        function propChange(metaProp,eventData)
            % Updates the image parameters if watched properties are changed
            
            update(eventData.AffectedObject);
        end
    end
end