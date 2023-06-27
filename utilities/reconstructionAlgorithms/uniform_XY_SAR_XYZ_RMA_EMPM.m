% uniform_XY_SAR_XYZ_RMA_EMPM is a reconstructor class that performs 3-D Range
% Migration Algorithm image reconstruction using EMPM algorithm introduced
% in "Efficient 3-D Near-Field MIMO-SAR Imaging for Irregular Scanning
% Geometries" by J. Smith and M. Torlak
%
% Copyright (C) 2021 Josiah W. Smith
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.

classdef uniform_XY_SAR_XYZ_RMA_EMPM < handle
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
        isMult2Mono = false   % Boolean whether or not to use the multistatic-to-monostatic approximation
        
        zRef_m = 0.25       % z location of reference plane for multistatic-to-monostatic approximation
        k_vec               % Instantaneous wavenumber vector
        z0_m                % Location of the antenna array in the z-plane
        xStep_m = 1e-3      % Step size along the x-dimension to move the antenna array in meters
        yStep_m = 8e-3      % Step size along the y-dimension to move the antenna array in meters
        
        im_method           % Interpolation method for image interpolation
        stolt_method        % Interpolation method for Stolt interpolation
        
        wav                 % A THzWaveformParameters object handle
        ant                 % A THzAntennaArray object handle
        scanner             % A THzScanner object handle
        target              % A THzTarget object handle
        im                  % A THzImageReconstruction object handle
    end
    
    methods
        function obj = uniform_XY_SAR_XYZ_RMA_EMPM(im)
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
            
            obj.im_method = obj.im.im_method;
            obj.stolt_method = obj.im.stolt_method;
            
            obj.nFFTx = obj.im.nFFTx;
            obj.nFFTy = obj.im.nFFTy;
            obj.nFFTz = obj.im.nFFTz;
            
            obj.x_m = obj.im.x_m;
            obj.y_m = obj.im.y_m;
            obj.z_m = obj.im.z_m;
            
            obj.sarData = obj.target.sarData;
            
            obj.isGPU = obj.im.isGPU;
            obj.isAmplitudeFactor = obj.target.isAmplitudeFactor;
            obj.isMult2Mono = obj.im.isMult2Mono;
            
            obj.zRef_m = obj.im.zRef_m;
            obj.k_vec = obj.wav.k;
            obj.z0_m = obj.ant.z0_m;
            obj.xStep_m = obj.scanner.xStep_m;
            obj.yStep_m = obj.scanner.yStep_m;
        end
        
        function verifyParameters(obj)
            % Verify the parameters allow for imaging
            
            kZU = single(reshape(linspace(0,2*max(obj.k_vec) - 2*max(obj.k_vec)/obj.nFFTz,obj.nFFTz),1,1,[]));
            dkZU = kZU(2) - kZU(1);
            x_m_temp = double(make_x(obj,obj.xStep_m,obj.nFFTx));
            y_m_temp = double(make_x(obj,obj.yStep_m,obj.nFFTy));
            z_m_temp = double(2*pi / (dkZU * obj.nFFTz) * (0:obj.nFFTz-1));
            
            if obj.im.isApp
                app = obj.im.app;
                app.XMinmEditField_im_2.Value = min(x_m_temp);
                app.XMaxmEditField_im_2.Value = max(x_m_temp);
                
                app.YMinmEditField_im_2.Value = min(y_m_temp);
                app.YMaxmEditField_im_2.Value = max(y_m_temp);
                
                app.ZMinmEditField_im_2.Value = min(z_m_temp);
                app.ZMaxmEditField_im_2.Value = max(z_m_temp);
            end
            
            if max(abs(obj.x_m)) > max(abs(x_m_temp))
                showErrorMessage(obj.im,"xMax_m is too large for nFFTx. Decrease xMax_m or increase nFFTx","3D RMA Error")
                obj.isFail = true;
                return;
            end
            if max(abs(obj.y_m)) > max(abs(y_m_temp))
                showErrorMessage(obj.im,"yMax_m is too large for nFFTy. Decrease yMax_m or increase nFFTy","3D RMA Error")
                obj.isFail = true;
                return;
            end
            if max(obj.z_m) > max(z_m_temp)
                showErrorMessage(obj.im,"zMax_m is too large for nFFTz. Decrease zMax_m or increase nFFTz","3D RMA Error")
                obj.isFail = true;
                return;
            end
        end
        
        function verifyReconstruction(obj)
            % Verify the reconstruction can continue
            
            if obj.scanner.method ~= "Rectilinear"
                showErrorMessage(obj.im,"Must use 2-D XY SAR scan to use 2-D SAR 3-D RMA image reconstruction method!","3D RMA Error");
                obj.isFail = true;
                return
            end
            
            % Ensure array is colinear
            if max(diff([obj.ant.tx.xy_m(:,1);obj.ant.rx.xy_m(:,1)])) > 8*sqrt(eps)
                showErrorMessage(obj.im,"MIMO array must be colinear. Please disable necessary elements.","3D RMA Error");
                obj.isFail = true;
                return
            end
            
            % Ensure virtual array is uniform
            if mean(diff(obj.ant.vx.xyz_m(:,2),2)) > sqrt(eps)
                showErrorMessage(obj.im,"Virtual antenna array is nonuniform! Change antenna positions.","3D RMA Error");
                obj.isFail = true;
                return
            end
            
            % And sar step size is correct
            if obj.scanner.yStep_m - mean(diff(obj.ant.vx.xyz_m(:,2)))*obj.ant.vx.numVx > 8*sqrt(eps)
                showErrorMessage(obj.im,"SAR step size is incorrect!","3D RMA Error");
                obj.isFail = true;
                return
            end
            
            if ~(isAppEPC(obj.im.app) || obj.ant.isEPC) && ~obj.isMult2Mono
                % If using MIMO Array Ensure multistatic-to-monostatic
                % approximation is employed
                showErrorMessage(obj.im,"Must use multistatic-to-monostatic approximation to use uniform method!","3D RMA Error");
                obj.isFail = true;
                return
            end
            
            % Everything is okay to continue
            obj.yStep_m = obj.scanner.yStep_m/obj.ant.vx.numVx;
            obj.isFail = false;
        end
        
        function imXYZ_out = computeReconstruction(obj)
            % Update the reconstruction algorithm and attempt the
            % reconstruction
            
            update(obj);
            
            if ~obj.isFail
                if ~obj.ant.isEPC && obj.isMult2Mono
                    freehandCompensation(obj);
                end
                
                reconstruct(obj);
                imXYZ_out = obj.imXYZ;
            else
                imXYZ_out = single(zeros(obj.im.numX,obj.im.numY,obj.im.numZ));
            end
        end
        
        function reconstruct(obj)
            % Reconstruct the image using the 3-D Range Migration Algorithm
            
            % Compute Wavenumbers
            k = single(reshape(obj.k_vec,1,1,[]));
            L_x = obj.nFFTx * obj.xStep_m;
            dkX = 2*pi/L_x;
            kX = make_kX(obj,dkX,obj.nFFTx);
            
            L_y = obj.nFFTy * obj.yStep_m;
            dkY = 2*pi/L_y;
            kY = make_kX(obj,dkY,obj.nFFTy)';
            
            % sarData is of size (scanner.numY, scanner.numX, wav.Nk)
            % Zero-Pad Data: s(y,x,k)
            sarDataPadded = obj.sarData;
            sarDataPadded = padarray(sarDataPadded,[floor((obj.nFFTy-size(obj.sarData,1))/2) 0],0,'pre');
            sarDataPadded = padarray(sarDataPadded,[0 floor((obj.nFFTx-size(obj.sarData,2))/2)],0,'pre');
            clear sarData
            
            kZU = single(reshape(linspace(0,2*max(k) - 2*max(k)/obj.nFFTz,obj.nFFTz),1,1,[]));
            dkZU = kZU(2) - kZU(1);
            
            if obj.isGPU
                reset(gpuDevice)
                k = gpuArray(k);
                kX = gpuArray(kX);
                kY = gpuArray(kY);
                kZU = gpuArray(kZU);
                sarDataPadded = gpuArray(sarDataPadded);
            end
            
            kYU = repmat(kY,[1,obj.nFFTx,obj.nFFTz]);
            kXU = repmat(kX,[obj.nFFTy,1,obj.nFFTz]);
            kU = single(1/2 * sqrt(kX.^2 + kY.^2 + kZU.^2));
            kZ = single(sqrt((4 * k.^2 - kX.^2 - kY.^2) .* (4 * k.^2 > kX.^2 + kY.^2)));
            
            % Compute Focusing Filter
            focusingFilter = exp(-1j * kZ * obj.z0_m);
            if obj.isAmplitudeFactor
                focusingFilter = kZ .* focusingFilter;
            end
            focusingFilter(4 * k.^2 < kX.^2 + kY.^2) = 0;
            
            % Compute FFT across Y & X Dimensions: S(kY,kX,k)
            sarDataFFT = fftshift(fftshift(fft(fft(conj(sarDataPadded),obj.nFFTy,1),obj.nFFTx,2),1),2)/obj.nFFTx/obj.nFFTy;
            clear sarDataPadded sarData
            
            if obj.isGPU
                focusingFilter = gpuArray(focusingFilter);
            end
            
            % Stolt Interpolation
            try
                sarImageFFT = interpn(kY(:),kX(:),k(:), sarDataFFT .* focusingFilter ,kYU,kXU,kU,obj.stolt_method,0);
            catch
                sarImageFFT = zeros(size(kU));
                for indkY = 1:size(kU,1)
                    for indkX = 1:size(kU,2)
                        tempS = squeeze(sarDataFFT(indkY,indkX,:) .* focusingFilter(indkY,indkX,:));
                        kZTemp = squeeze(kZ(indkY,indkX,:));
                        [kZTemp_unique,~,ind_c] = uniquetol(kZTemp);
                        tempS = accumarray(ind_c,tempS);
                        if length(kZTemp_unique) > 2
                            sarImageFFT(indkY,indkX,:) = gather(interp1(kZTemp_unique,tempS,kZU,obj.stolt_method,0));
                        end
                    end
                end
            end
            clear sarDataFFT focusingFilter kY kX k kYU kXU kU kZ kZU
            
            if obj.isGPU
                sarImageFFT = gather(sarImageFFT);
                reset(gpuDevice);
                sarImageFFT = gpuArray(sarImageFFT);
            end
            
            % Recover Image by IFT: p(y,x,z)
            sarImageFFT = sarImageFFT .* blackmanharris(obj.nFFTx).' .* blackmanharris(obj.nFFTy);
            
            sarImage = single(ifftn(sarImageFFT));
            clear sarImageFFT focusingFilter
            
            % Reorient Image: p(x,y,z)
            sarImage = permute(sarImage,[2,1,3]);
            
            % Declare Spatial Vectors
            x_m_temp = make_x(obj,obj.xStep_m,obj.nFFTx);
            y_m_temp = make_x(obj,obj.yStep_m,obj.nFFTy);
            z_m_temp = single(2*pi / (dkZU * obj.nFFTz) * (1:obj.nFFTz));
            
            % Interpolate Image
            obj.imXYZ = imageInterp3D(obj,x_m_temp,y_m_temp,z_m_temp,sarImage);
        end
        
        function displayImage(obj)
            % Display the reconstructed x-y-z image
            displayImage3D(obj.im);
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