% reconstructionAlgorithmTemplate_app see
% reconstructionAlgorithmTemplate documentation
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

classdef reconstructionAlgorithmTemplate_app < handle
    properties
        sarData             % Computed beat signal
        
        nFFTx = 512         % Number of FFT points along the x-dimension, when using FFT-based reconstruction algorithms
        nFFTy = 512         % Number of FFT points along the y-dimension, when using FFT-based reconstruction algorithms
        nFFTz = 512         % Number of FFT points along the z-dimension, when using FFT-based reconstruction algorithms
        
        x_m                 % Reconstructed image x axis
        y_m                 % Reconstructed image y axis
        z_m                 % Reconstructed image z axis
        
        imXYZ               % Reconstructed image
        
        isGPU               % Boolean whether or not to use the GPU for image reconstruction
        isAmplitudeFactor   % Boolean whether or not to include the amplitude factor in the image reconstruction process
        isFail = false      % Boolean whether or not the reconstruction has failed
        
        k_vec               % Instantaneous wavenumber vector
        z0_m                % Location of the antenna array in the z-plane
        xStep_m = 1e-3      % Step size along the x-dimension to move the antenna array in meters
        yStep_m = 8e-3      % Step size along the y-dimension to move the antenna array in meters
    end
    
    methods
        function obj = reconstructionAlgorithmTemplate_app(app,im)
            
        end
        function obj = update(obj,app,im)
            obj = getParameters(obj,app,im);
            obj = verifyParameters(obj,app);
            obj = verifyReconstruction(obj,app);
        end
        
        function obj = getParameters(obj,app,im)
            obj.nFFTx = im.nFFTx;
            obj.nFFTy = im.nFFTy;
            obj.nFFTz = im.nFFTz;
            
            obj.x_m = im.x_m;
            obj.y_m = im.y_m;
            obj.z_m = im.z_m;
            
            obj.sarData = app.target.sarData;
            
            obj.isGPU = im.isGPU;
            obj.isAmplitudeFactor = app.target.isAmplitudeFactor;
            
            obj.k_vec = app.fmcw.k;
            obj.z0_m = app.ant.tx.z0_m;
            obj.xStep_m = app.sar.xStep_m;
            obj.yStep_m = app.sar.yStep_m;
        end
        
        function [obj,imXYZ_out] = computeReconstruction(obj,app,im)
            obj = update(obj,app,im);
            if ~obj.isFail
                obj = reconstruct(obj,app);
            end
            imXYZ_out = obj.imXYZ;
        end
        
        function obj = reconstruct(obj,app)
            
        end
    end
end