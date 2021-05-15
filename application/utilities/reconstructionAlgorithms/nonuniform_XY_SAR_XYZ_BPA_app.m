% nonuniform_XY_SAR_XYZ_BPA_app see nonuniform_XY_SAR_XYZ_BPA
% documentation
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

classdef nonuniform_XY_SAR_XYZ_BPA_app < handle
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
        isMIMO              % Boolean whether or not to use the MIMO physical element locations instead of the equivalent phase center virtual element locations
        
        k_vec               % Instantaneous wavenumber vector
        
        estTime             % Structure for holding the estimated time until completion parameters
    end
    
    methods
        function obj = nonuniform_XY_SAR_XYZ_BPA_app(app,im)
            obj = update(obj,app,im);
            obj.estTime.old = inf;
            obj.estTime.count = 0;
        end
        function obj = update(obj,app,im)
            % Update the reconstruction algorithm by getting the parameters
            % from the object handles and verifying the paramters
            
            obj = getParameters(obj,app,im);
            obj = verifyReconstruction(obj,app);
        end
        
        function obj = getParameters(obj,app,im)
            obj.tx_xyz_m = app.sar.tx.xyz_m;
            obj.rx_xyz_m = app.sar.rx.xyz_m;
            obj.vx_xyz_m = app.sar.vx.xyz_m;
            
            [X,Y,Z] = ndgrid(im.x_m(:),im.y_m(:),im.z_m(:));
            obj.sizeTarget = [im.numX,im.numY,im.numZ];
            
            obj.target_xyz_m = [X(:),Y(:),Z(:)];
            
            obj.sarData = reshape(app.target.sarData,[],app.fmcw.ADCSamples);
            
            obj.isGPU = im.isGPU;
            obj.isAmplitudeFactor = app.target.isAmplitudeFactor;
            obj.isMIMO = app.sar.isMIMO;
            
            obj.k_vec = app.fmcw.k;
        end
        
        function obj = verifyReconstruction(obj,app)
            obj.isFail = false;
            
            if app.sar.method ~= "Rectilinear"
                uiconfirm(app.UIFigure,"Must use 2-D XY SAR scan to use 2-D SAR 3-D BPA image reconstruction method!",'SAR Scenario Error!',...
                    "Options",{'OK'},'Icon','warning');
                obj.isFail = true;
                return
            end
        end
        
        function [obj,imXYZ_out] = computeReconstruction(obj,app,im)
            obj = update(obj,app,im);
            obj = reconstruct(obj,app);
            imXYZ_out = obj.imXYZ;
        end
        
        function obj = reconstruct(obj,app)
            % Create the progress dialog
            d = uiprogressdlg(app.UIFigure,'Title','Performing XYZ BPA',...
                'Message',"Estimated Time Remaining: 0:0:0","Cancelable","on");
            
            if obj.isGPU
                reset(gpuDevice);
            end
            
            % Orient vx_xyz x target_xyz x k
            obj.sarData = reshape(obj.sarData,[],1,length(obj.k_vec));
            k = single(reshape(obj.k_vec,1,1,[]));
            
            try
                % Fast way
                obj = fastBPA_app(obj,k,d);
            catch
                d.Title = "Performing XYZ BPA";
                obj = mediumBPA_app(obj,k,d);
            end
            
            delete(d)
            try
                obj.imXYZ = reshape(obj.imXYZ,obj.sizeTarget);
            catch
                warning("WHAT!")
            end
        end
        
        function obj = displayImage(obj,im)
            displayImage3D_app(im);
        end
    end
end