% nonuniform_theta_CSAR_XZ_BPA is a reconstructor class that performs
% a 2-D Backprojection Algorithm image reconstruction. The synthetic
% aperture must span a circular aperture and the target can be a
% 1-D or 2-D target in the x-z place at the y-coordinate y = 0
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

classdef nonuniform_theta_CSAR_XZ_BPA < handle
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
        
        fmcw                % fmcwChirpParameters object
        ant                 % sarAntennaArray object
        sar                 % sarScenario object
        target              % sarTarget object
        im                  % sarImage object
    end
    
    methods
        function obj = nonuniform_theta_CSAR_XZ_BPA(im)
            % Set the properties corresponding to the object handles for
            % the imaging scenario and get the parameters from those object
            % handles
            
            obj.fmcw = im.fmcw;
            obj.ant = im.ant;
            obj.sar = im.sar;
            obj.target = im.target;
            obj.im = im;
            
            getParameters(obj);
            obj.estTime.old = inf;
            obj.estTime.count = 0;
        end
        
        function update(obj)
            % Update the reconstruction algorithm by getting the parameters
            % from the object handles and verifying the parameters
            
            getParameters(obj);
            verifyReconstruction(obj);
        end
        
        function getParameters(obj)
            % Get the parameters from the object handles
            
            obj.tx_xyz_m = obj.sar.tx.xyz_m;
            obj.rx_xyz_m = obj.sar.rx.xyz_m;
            obj.vx_xyz_m = obj.sar.vx.xyz_m;
            
            [X,Y,Z] = ndgrid(obj.im.x_m(:),0,obj.im.z_m(:));
            obj.sizeTarget = [obj.im.numX,obj.im.numZ];
            
            obj.target_xyz_m = [X(:),Y(:),Z(:)];
            
            obj.sarData = reshape(obj.target.sarData,[],obj.fmcw.ADCSamples);
            
            obj.isGPU = obj.im.isGPU;
            obj.isAmplitudeFactor = obj.target.isAmplitudeFactor;
            
            obj.k_vec = obj.fmcw.k;
        end
        
        function verifyReconstruction(obj)
            % Verify that the reconstruction can continue. The scan must be
            % a circular scan
            
            obj.isFail = false;
            
            if obj.sar.scanMethod ~= "Circular"
                warning("Must use 1-D θ Circular CSAR scan to use 1-D CSAR 2-D BPA image reconstruction method!");
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
                imXYZ_out = single(zeros(obj.im.numX,obj.im.numZ));
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
                warning("WHAT!")
            end
        end
        
        function displayImage(obj)
            % Display the x-z reconstructed image
            
            displayImage2D(obj.im,obj.im.x_m,obj.im.z_m,"x (m)","z (m)");
        end
    end
end