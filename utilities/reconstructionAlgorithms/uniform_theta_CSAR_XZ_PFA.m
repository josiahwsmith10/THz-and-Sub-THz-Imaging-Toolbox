% uniform_theta_CSAR_XZ_PFA is a reconstructor class that performs 2-D
% Polar Formatting Algorithm image reconstruction. The synthetic
% aperture must span the theta dimension, forming a circular array, and
% the target can be a 1-D or 2-D target in the x-z plane at the
% y-coordinate y = 0
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

classdef uniform_theta_CSAR_XZ_PFA < handle
    properties
        sarData             % Computed beat signal
        
        nFFTx = 512         % Number of FFT points along the x-dimension, when using FFT-based reconstruction algorithms
        nFFTz = 512         % Number of FFT points along the z-dimension, when using FFT-based reconstruction algorithms
        thetaUpsampleFactor = 1 % Upsample factor on the theta dimension of the SAR data operated in circular or cylindrical mode
        
        x_m                 % Reconstructed image x axis
        z_m                 % Reconstructed image z axis
        
        imXYZ               % Reconstructed image
        
        isGPU               % Boolean whether or not to use the GPU for image reconstruction
        isAmplitudeFactor   % Boolean whether or not to include the amplitude factor in the image reconstruction process
        isFail = false      % Boolean whether or not the reconstruction has failed
        
        theta_rad_vec       % Angular scan locations vector
        k_vec               % Instantaneous wavenumber vector
        R0_m                % Radius of scan
        
        wav                 % A THzWaveformParameters object handle
        ant                 % A THzAntennaArray object handle
        scanner             % A THzScanner object handle
        target              % A THzTarget object handle
        im                  % A THzImageReconstruction object handle
    end
    
    methods
        function obj = uniform_theta_CSAR_XZ_PFA(im)
            % Set the properties corresponding to the object handles for
            % the imaging scenario and get the parameters from those object
            % handles
            
            obj.wav = im.wav;
            obj.ant = im.ant;
            obj.scanner = im.scanner;
            obj.target = im.target;
            obj.im = im;
            
            getParameters(obj);
        end
        
        function update(obj)
            % Update the reconstruction algorithm by getting the parameters
            % from the object handles and verifying the parameters
            
            getParameters(obj);
            verifyReconstruction(obj);
            verifyParameters(obj);
        end
        
        function getParameters(obj)
            % Get the parameters from the object handles
            
            obj.nFFTx = obj.im.nFFTx;
            obj.nFFTz = obj.im.nFFTz;
            
            obj.x_m = obj.im.x_m;
            obj.z_m = obj.im.z_m;
            
            obj.sarData = obj.target.sarData;
            
            obj.isGPU = obj.im.isGPU;
            obj.isAmplitudeFactor = obj.target.isAmplitudeFactor;
            
            obj.thetaUpsampleFactor = obj.im.thetaUpsampleFactor;
            
            obj.theta_rad_vec = obj.scanner.theta_rad;
            obj.k_vec = obj.wav.k;
            obj.R0_m = obj.ant.z0_m;
            
            if obj.im.isApp
                obj.thetaUpsampleFactor = obj.im.app.ThetaUpsampleFactorEditField_2.Value;
            end
        end
        
        function verifyParameters(obj)
            % Verify the parameters allow for imaging
            
            [x_m_temp,z_m_temp] = getTempXZ(obj);
            
            if obj.im.isApp
                app = obj.im.app;
                app.XMinmEditField_im_5.Value = min(x_m_temp);
                app.XMaxmEditField_im_5.Value = max(x_m_temp);
                
                app.ZMinmEditField_im_5.Value = min(z_m_temp);
                app.ZMaxmEditField_im_5.Value = max(z_m_temp);
            end
            
            if max(abs(obj.x_m)) > max(abs(x_m_temp))
                showErrorMessage(obj.im,"xMax_m is too large for nFFTx. Decrease xMax_m or increase nFFTx","2D PFA Error")
                obj.isFail = true;
                return;
            end
            
            if obj.thetaUpsampleFactor*obj.scanner.numTheta > obj.nFFTx
                showErrorMessage(obj.im,"thetaUpsampleFactor is too large for nFFTx. Decrease thetaUpsampleFactor or increase nFFTx","2D PFA Error")
                obj.isFail = true;
                return;
            end
            
            if max(obj.z_m) > max(z_m_temp)
                showErrorMessage(obj.im,"zMax_m is too large for nFFTz. Decrease zMax_m or increase nFFTz","2D PFA Error")
                obj.isFail = true;
                return;
            end
        end
        
        function verifyReconstruction(obj)
            % Verify the reconstruction can continue
            
            if obj.scanner.method ~= "Circular"
                showErrorMessage(obj.im,"Must use 1-D θ Circular CSAR scan to use 1-D CSAR 2-D PFA image reconstruction method!","2D PFA Error");
                obj.isFail = true;
                return
            end
            
            % Ensure single element array
            if obj.ant.tx.numTx ~= 1 || obj.ant.rx.numRx ~= 1
                showErrorMessage(obj.im,"Array must have only 1 Tx and 1 Rx. Please disable necessary elements.","2D PFA Error");
                obj.isFail = true;
                return
            end
            obj.isFail = false;
        end
        
        function [x_m_temp,z_m_temp] = getTempXZ(obj)
            % Compute the axes of the recovered image for parameter
            % verification
            
            theta_rad = single(reshape(obj.theta_rad_vec,1,[]));
            
            % Compute Wavenumbers
            k = single(reshape(obj.k_vec,1,1,[]));
            
            kR = 2*k;
            kX = kR.*cos(theta_rad);
            kZ = kR.*sin(theta_rad);
            
            kXmax = max(kX,[],'all');
            kZmax = max(kZ,[],'all');
            kXmin = min(kX,[],'all');
            kZmin = min(kZ,[],'all');
            clear kX kZ
            
            kXU = single(reshape(linspace(kXmin,kXmax,obj.nFFTx),1,[],1));
            kZU = single(reshape(linspace(kZmin,kZmax,obj.nFFTz),1,1,[]));
            dkXU = kXU(2) - kXU(1);
            dkZU = kZU(2) - kZU(1);
            
            x_m_temp = double(make_x(obj,2*pi/(dkXU*obj.nFFTx),obj.nFFTx));
            z_m_temp = double(make_x(obj,2*pi/(dkZU*obj.nFFTz),obj.nFFTz));
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
            % Reconstruct the image using the 2-D Polar Formatting
            % Algorithm
            
            % sarData is of size (sar.numTheta, fmcw.ADCSamples)
            % Zero-Pad Data: s(theta,k)
            sarDataPadded = obj.sarData;
            
            theta_rad = single(reshape(obj.theta_rad_vec,[],1));
            
            % Compute Wavenumbers
            k = single(reshape(obj.k_vec,1,[]));
            
            kR = 2*k;
            kX = kR.*cos(theta_rad);
            kZ = kR.*sin(theta_rad);
            
            kXmax = max(kX,[],'all');
            kZmax = max(kZ,[],'all');
            kXmin = min(kX,[],'all');
            kZmin = min(kZ,[],'all');
            clear kX kZ
            
            kXU = single(reshape(linspace(kXmin,kXmax,obj.nFFTx),[],1));
            kZU = single(reshape(linspace(kZmin,kZmax,obj.nFFTz),1,[]));
            dkXU = kXU(2) - kXU(1);
            dkZU = kZU(2) - kZU(1);
            
            % Declare Uniform Grid for Interpolation
            theta_radU = atan2(kZU,kXU);
            kU = 1/2*sqrt(kXU.^2 + kZU.^2);
            clear kXU kZU
            
            % Upsample data
            if obj.thetaUpsampleFactor > 1
                theta_radUp = single(reshape(linspace(theta_rad(1),theta_rad(end),length(theta_rad)*obj.thetaUpsampleFactor),[],1));
                [T,K] = ndgrid(theta_radUp(:),reshape(single(1:size(sarDataPadded,2)),1,[]));
                sarDataPadded = interpn(theta_rad(:),single(1:size(sarDataPadded,2))',sarDataPadded,T,K);
                clear T K
                theta_rad = theta_radUp;
            end
            
            % Compute Azimuth Filter H(Theta,k)
            azimuthFilterFFT = single(fft(exp(1j*kR*obj.R0_m.*cos(theta_rad)),[],1));
            clear kR
            
            % Compute Azimuth Filtered Data: p(theta,k) = IFT_Theta[ S(Theta,k) * H*(Theta,k)]
            azimuthFiltered = ifft(fft(sarDataPadded,[],1) .* conj(azimuthFilterFFT),[],1);
            clear azimuthFilterFFT sarDataPadded
            
            if obj.isGPU
                reset(gpuDevice)
                theta_rad = gpuArray(theta_rad);
                k = gpuArray(k);
                azimuthFiltered = gpuArray(azimuthFiltered);
                theta_radU = gpuArray(theta_radU);
                kU = gpuArray(kU);
            end
            
            % Stolt Interpolation
            sarImageFFT = interpn(theta_rad(:),k(:), azimuthFiltered ,theta_radU,kU,'linear',0);
            clear azimuthFiltered kY theta_rad k azimuthfiltered kYU theta_radU kU
            
            if obj.isGPU
                sarImageFFT = gather(sarImageFFT);
                reset(gpuDevice);
                sarImageFFT = gpuArray(sarImageFFT);
            end
            
            % Recover Image by IFT: p(x,z)
            sarImage = single(ifftshift(ifftshift(ifftn(sarImageFFT),1),2));
            clear sarImageFFT
            
            % Declare Spatial Vectors
            x_m_temp = make_x(obj,2*pi/(dkXU*obj.nFFTx),obj.nFFTx);
            z_m_temp = make_x(obj,2*pi/(dkZU*obj.nFFTz),obj.nFFTz);
            
            [X,Z] = ndgrid(obj.x_m(:),obj.z_m(:));
            obj.imXYZ = single(gather(interpn(x_m_temp(:),z_m_temp(:),sarImage,X,Z,"nearest",0)));
        end
        
        function displayImage(obj)
            % Display the reconstructed x-z image
            
            displayImage2D(obj.im,obj.im.x_m,obj.im.z_m,"x (m)","z (m)");
        end
        
        function x = make_x(obj,xStep_m,nFFTx)
            % Make an imaging axis from the step size and number of FFT
            % points
            
            x = xStep_m * (-(nFFTx-1)/2 : (nFFTx-1)/2);
            x = single(x);
        end
    end
end