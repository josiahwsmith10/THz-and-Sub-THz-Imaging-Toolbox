% nonuniform_XY_SAR_XYZ_BPA is a reconstructor class that performs a
% 3-D Backprojection Algorithm image reconstruction. The synthetic
% aperture must span an x-y plane and can be nonuniformly spaced or use
% a MIMO array, even overlapping. The target can be a 1-D, 2-D, or 3-D
% target in x-y-z space
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

classdef nonuniform_XY_SAR_XYZ_BPA < handle
    properties
        sarData             % Computed beat signal
        
        tx_xyz_m            % x-y-z coordinates of the transmitter antennas in the synthetic aperture
        rx_xyz_m            % x-y-z coordinates of the receiver antennas in the synthetic aperture
        vx_xyz_m            % x-y-z coordinates of the virtual elements in the synthetic aperture
        target_xyz_m        % x-y-z coordinates of the voxels in the target domain, as specified by the user
        sizeTarget          % Dimensions of the desired target, as specified by the user
        
        imXYZ               % Reconstructed image
        
        isGPU               % Boolean whether or not to use the GPU for image reconstruction
        isAmplitudeFactor   % Boolean whether or not to include the amplitude factor in the image reconstruction process
        isFail = false      % Boolean whether or not the reconstruction has failed
        
        k_vec               % Instantaneous wavenumber vector
        
        estTime             % Structure for holding the estimated time until completion parameters
        
        wav                 % A THzWaveformParameters object handle
        ant                 % A THzAntennaArray object handle
        scanner             % A THzScanner object handle
        target              % A THzTarget object handle
        im                  % A THzImageReconstruction object handle
    end
    
    methods
        function obj = nonuniform_XY_SAR_XYZ_BPA(im)
            % Set the properties corresponding to the object handles for
            % the imaging scenario and get the parameters from those object
            % handles
            
            obj.wav = im.wav;
            obj.ant = im.ant;
            obj.scanner = im.scanner;
            obj.target = im.target;
            obj.im = im;
            
            getParameters(obj);
            obj.estTime.old = inf;
            obj.estTime.count = 0;
        end
        
        function update(obj)
            % Update the reconstruction algorithm by getting the parameters
            % from the object handles and verifying the paramters
            
            getParameters(obj);
            verifyReconstruction(obj);
        end
        
        function getParameters(obj)
            % Get the parameters from the object handles
            
            obj.tx_xyz_m = obj.scanner.tx.xyz_m;
            obj.rx_xyz_m = obj.scanner.rx.xyz_m;
            obj.vx_xyz_m = obj.scanner.vx.xyz_m;
            
            [X,Y,Z] = ndgrid(obj.im.x_m(:),obj.im.y_m(:),obj.im.z_m(:));
            obj.sizeTarget = [obj.im.numX,obj.im.numY,obj.im.numZ];
            
            obj.target_xyz_m = [X(:),Y(:),Z(:)];
            
            obj.sarData = reshape(obj.target.sarData,[],obj.wav.Nk);
            
            obj.isGPU = obj.im.isGPU;
            obj.isAmplitudeFactor = obj.target.isAmplitudeFactor;
            
            obj.k_vec = obj.wav.k;
        end
        
        function verifyReconstruction(obj)
            % Verify that the reconstruction can continue. The scan must be
            % a rectilinear scan
            
            obj.isFail = false;
            
            if obj.scanner.method ~= "Rectilinear"
                showErrorMessage(obj.im,"Must use 2-D XY SAR scan to use 2-D SAR 3-D BPA image reconstruction method!","3D BPA Error");
                obj.isFail = true;
                return
            end
        end
        
        function imXYZ_out = computeReconstruction(obj)
            % Update the reconstruction algorithm and attempt the
            % reconstruction
            
            update(obj);
            
            if ~obj.isFail
                reconstruct(obj);
                imXYZ_out = obj.imXYZ;
            else
                imXYZ_out = single(zeros(obj.im.numX,obj.im.numY,obj.im.numZ));
            end
        end
        
        function reconstruct(obj)
            % Reconstruct the image using the 2-D Backprojection Algorithm
            
            if obj.isGPU
                reset(gpuDevice);
            end
            
            % Orient vx_xyz x target_xyz x k
            obj.sarData = reshape(obj.sarData,[],1,length(obj.k_vec));
            k = single(reshape(obj.k_vec,1,1,[]));
            
            try
                % Fast way
                fastBPA(obj,k);
            catch
                mediumBPA(obj,k);
            end
            
            try
                obj.imXYZ = reshape(obj.imXYZ,obj.sizeTarget);
            catch
                showErrorMessage(obj.im,"Error in BPA","BPA Error!");
                obj.isFail = true;
            end
        end
        
        function displayImage(obj)
            % Display the x-y-z reconstructed image
            
            displayImage3D(obj.im);
        end
    end
end