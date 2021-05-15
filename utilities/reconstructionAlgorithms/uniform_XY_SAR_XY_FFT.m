% uniform_XY_SAR_XY_FFT is a reconstructor class that performs 2-D FFT
% image reconstruction. The synthetic aperture must span an x-y plane
% at the z-coordinate z0_m and the target must be a 2-D target at in
% x-y plane at the z-coordinate zSlice_m
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

classdef uniform_XY_SAR_XY_FFT < handle
    properties
        sarData             % Computed beat signal
        
        nFFTx = 512         % Number of FFT points along the x-dimension, when using FFT-based reconstruction algorithms
        nFFTy = 512         % Number of FFT points along the y-dimension, when using FFT-based reconstruction algorithms
        
        x_m                 % Reconstructed image x axis
        y_m                 % Reconstructed image y axis
        
        imXYZ               % Reconstructed image
        
        isGPU               % Boolean whether or not to use the GPU for image reconstruction
        isAmplitudeFactor   % Boolean whether or not to include the amplitude factor in the image reconstruction process
        isFail = false      % Boolean whether or not the reconstruction has failed
        isMult2Mono = false   % Boolean whether or not to use the multistatic-to-monostatic approximation
        
        zRef_m = 0.25       % z location of reference plane for multistatic-to-monostatic approximation
        zSlice_m            % z slice of interest when reconstructing a 2-D x-y image, in meters
        k_vec               % Instantaneous wavenumber vector
        z0_m                % Location of the antenna array in the z-plane
        xStep_m = 1e-3      % Step size along the x-dimension to move the antenna array in meters
        yStep_m = 8e-3      % Step size along the y-dimension to move the antenna array in meters
        
        fmcw                % fmcwChirpParameters object
        ant                 % sarAntennaArray object
        sar                 % sarScenario object
        target              % sarTarget object
        im                  % sarImage object
    end
    
    methods
        function obj = uniform_XY_SAR_XY_FFT(im)
            % Set the properties corresponding to the object handles for
            % the imaging scenario and get the parameters from those object
            % handles
            
            obj.fmcw = im.fmcw;
            obj.ant = im.ant;
            obj.sar = im.sar;
            obj.target = im.target;
            obj.im = im;
            
            getParameters(obj);
        end
        
        function obj = update(obj)
            % Update the reconstruction algorithm by getting the parameters
            % from the object handles and verifying the parameters
            
            getParameters(obj);
            verifyParameters(obj);
            verifyReconstruction(obj);
        end
        
        function getParameters(obj)
            % Get the parameters from the object handles
            
            obj.nFFTx = obj.im.nFFTx;
            obj.nFFTy = obj.im.nFFTy;
            
            obj.x_m = obj.im.x_m;
            obj.y_m = obj.im.y_m;
            
            obj.sarData = obj.target.sarData;
            
            obj.isGPU = obj.im.isGPU;
            obj.isAmplitudeFactor = obj.target.isAmplitudeFactor;
            obj.isMult2Mono = obj.im.isMult2Mono;
            
            obj.zRef_m = obj.im.zRef_m;
            obj.k_vec = obj.fmcw.k;
            obj.z0_m = obj.ant.z0_m;
            obj.zSlice_m = obj.im.zSlice_m;
            obj.xStep_m = obj.sar.xStep_m;
            obj.yStep_m = obj.sar.yStep_m;
        end
        
        function verifyParameters(obj)
            % Verify the parameters allow for imaging
            
            obj.isFail = false;
            
            x_m_temp = make_x(obj,obj.sar.xStep_m,obj.nFFTx);
            y_m_temp = make_x(obj,obj.sar.yStep_m,obj.nFFTy);
            
            if max(abs(obj.x_m)) > max(abs(x_m_temp))
                warning("xMax_m is too large for nFFTx. Decrease xMax_m or increase nFFTx")
                obj.isFail = true;
                return;
            end
            if max(abs(obj.y_m)) > max(abs(y_m_temp))
                warning("yMax_m is too large for nFFTy. Decrease yMax_m or increase nFFTy")
                obj.isFail = true;
                return;
            end
            if isempty(obj.zSlice_m)
                warning("Must specify zSlice_m property of sarImage!")
                obj.isFail = true;
                return;
            end
        end
        
        function verifyReconstruction(obj)
            % Verify the reconstruction can continue
            
            if obj.sar.scanMethod ~= "Rectilinear"
                warning("Must use 2-D XY SAR scan to use 2-D SAR 2-D FFT image reconstruction method!");
                obj.isFail = true;
                return
            end
            
            % Ensure array is colinear
            if max(diff([obj.ant.tx.xy_m(:,1);obj.ant.rx.xy_m(:,1)])) > 8*eps
                warning("MIMO array must be colinear. Please disable necessary elements.");
                obj.isFail = true;
                return
            end
            
            % Ensure virtual array is uniform
            if mean(diff(obj.ant.vx.xyz_m(:,2),2)) > eps
                warning("Virtual antenna array is nonuniform! Change antenna positions.");
                obj.isFail = true;
                return
            end
            
            % And sar step size is correct
            if obj.sar.yStep_m - mean(diff(obj.ant.vx.xyz_m(:,2)))*obj.ant.vx.numVx > 8*eps
                warning("SAR step size is incorrect!");
                obj.isFail = true;
                return
            end
            
            if ~obj.ant.isEPC && ~obj.isMult2Mono
                % If using MIMO Array
                % Ensure multistatic-to-monostatic approximation is employed
                warning("Must use multistatic-to-monostatic approximation to use uniform method!");
                obj.isFail = true;
                return
            end
            
            % Everything is okay to continue
            obj.yStep_m = obj.sar.yStep_m/obj.ant.vx.numVx;
        end
        
        function imXYZ_out = computeReconstruction(obj)
            % Update the reconstruction algorithm and attempt the
            % reconstruction
            
            update(obj);
            
            if ~obj.isFail
                if ~obj.ant.isEPC && obj.isMult2Mono
                    mult2mono(obj);
                end
                
                reconstruct(obj);
                imXYZ_out = obj.imXYZ;
            else
                imXYZ_out = single(zeros(obj.im.numX,obj.im.numY));
            end
        end
        
        function obj = reconstruct(obj)
            % Reconstruct the image using the 2-D FFT Method
            
            % sarData is of size (sar.numY, sar.numX, fmcw.ADCSamples)
            % Zero-Pad Data: s(y,x,k)
            sarDataPadded = obj.sarData;
            sarDataPadded = padarray(sarDataPadded,[floor((obj.nFFTy-size(obj.sarData,1))/2) 0],0,'pre');
            sarDataPadded = padarray(sarDataPadded,[0 floor((obj.nFFTx-size(obj.sarData,2))/2)],0,'pre');
            clear sarData
            
            % Compute Wavenumbers
            k = single(reshape(obj.k_vec,1,1,[]));
            L_x = obj.nFFTx * obj.xStep_m;
            dkX = 2*pi/L_x;
            kX = make_kX(obj,dkX,obj.nFFTx);
            
            L_y = obj.nFFTy * obj.yStep_m;
            dkY = 2*pi/L_y;
            kY = make_kX(obj,dkY,obj.nFFTy)';
            
            if obj.isGPU
                reset(gpuDevice)
                k = gpuArray(k);
                kX = gpuArray(kX);
                kY = gpuArray(kY);
                sarDataPadded = gpuArray(sarDataPadded);
            end
            
            % Compute FFT across Y & X Dimensions: S(kY,kX,k)
            sarDataFFT = fftshift(fftshift(fft(fft(sarDataPadded,obj.nFFTy,1),obj.nFFTx,2),1),2)/obj.nFFTx/obj.nFFTy;
            clear sarDataPadded
            
            kZ = single(sqrt((4 * k.^2 - kX.^2 - kY.^2) .* (4 * k.^2 > kX.^2 + kY.^2)));
            
            % Compute Focusing Filter
            focusingFilter = exp(-1j * kZ * (obj.zSlice_m - obj.z0_m));
            focusingFilter(4 * k.^2 < kX.^2 + kY.^2) = 0;
            
            if obj.isGPU
                focusingFilter = gpuArray(focusingFilter);
            end
            
            sarImageFFT = sum(sarDataFFT .* focusingFilter,3);
            clear sarDataFFT focusingFilter kY kX k kZ
            
            % Recover Image by IFT: p(y,x,z)
            sarImage = single(ifftn(sarImageFFT));
            clear sarImageFFT focusingFilter
            
            % Reorient Image: p(x,y,z)
            sarImage = permute(sarImage,[2,1,3]);
            
            % Declare Spatial Vectors
            x_m_temp = make_x(obj,obj.xStep_m,obj.nFFTx);
            y_m_temp = make_x(obj,obj.yStep_m,obj.nFFTy);
            
            [X,Y] = ndgrid(obj.x_m(:),obj.y_m(:));
            obj.imXYZ = single(gather(interpn(x_m_temp(:),y_m_temp(:),sarImage,X,Y,"nearest",0)));
        end
        
        function displayImage(obj)
            % Display the reconstructed x-y image
            
            displayImage2D(obj.im,obj.im.x_m,obj.im.y_m,"x (m)","y (m)");
        end
        
        function x = make_x(obj,xStep_m,nFFTx)
            % Make an imaging axis from the step size and number of FFT
            % points
            
            x = xStep_m * (-(nFFTx-1)/2 : (nFFTx-1)/2);
            x = single(x);
        end
        
        function kX = make_kX(obj,dkX,nFFTx)
            % Make a spatial wavenumber vector from the step size of the
            % spatial wavenumber and number of FFT points
            
            if mod(nFFTx,2)==0
                kX = dkX * ( -nFFTx/2 : nFFTx/2-1 );
            else
                kX = dkX * ( -(nFFTx-1)/2 : (nFFTx-1)/2 );
            end
            kX = single(kX);
        end
    end
end