% uniform_Y_SAR_YZ_RMA is a reconstructor class that performs a 2-D Range
% Migration Algorithm image reconstruction. The synthetic aperture must
% span the y-dimension and the target must be a 2-D target in the y-z plane
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

classdef uniform_Y_SAR_YZ_RMA < handle
    properties
        sarData             % Computed beat signal

        nFFTy = 512         % Number of FFT points along the y-dimension, when using FFT-based reconstruction algorithms
        nFFTz = 512         % Number of FFT points along the z-dimension, when using FFT-based reconstruction algorithms

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
        function obj = uniform_Y_SAR_YZ_RMA(im)
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

            obj.nFFTy = obj.im.nFFTy;
            obj.nFFTz = obj.im.nFFTz;

            obj.y_m = obj.im.y_m;
            obj.z_m = obj.im.z_m;

            obj.sarData = obj.target.sarData;

            obj.isGPU = obj.im.isGPU;
            obj.isAmplitudeFactor = obj.target.isAmplitudeFactor;
            obj.isMult2Mono = obj.im.isMult2Mono;

            obj.zRef_m = obj.im.zRef_m;
            obj.k_vec = obj.wav.k;
            obj.z0_m = obj.ant.z0_m;
            obj.yStep_m = obj.scanner.yStep_m;
        end

        function verifyParameters(obj)
            % Verify the parameters allow for imaging

            kZU = single(reshape(linspace(0,2*max(obj.k_vec) - 2*max(obj.k_vec)/obj.nFFTz,obj.nFFTz),1,1,[]));
            dkZU = kZU(2) - kZU(1);
            y_m_temp = double(make_x(obj,obj.yStep_m,obj.nFFTy));
            z_m_temp = double(2*pi / (dkZU * obj.nFFTz) * (1:obj.nFFTz));

            if obj.im.isApp
                app = obj.im.app;
                app.YMinmEditField_im_5.Value = min(y_m_temp);
                app.YMaxmEditField_im_5.Value = max(y_m_temp);

                app.ZMinmEditField_im_3.Value = min(z_m_temp);
                app.ZMaxmEditField_im_3.Value = max(z_m_temp);
            end

            if max(abs(obj.y_m)) > max(abs(y_m_temp))
                showErrorMessage(obj.im,"yMax_m is too large for nFFTy. Decrease yMax_m or increase nFFTy","2D RMA Error")
                obj.isFail = true;
                return;
            end
            if max(obj.z_m) > max(z_m_temp)
                showErrorMessage(obj.im,"zMax_m is too large for nFFTz. Decrease zMax_m or increase nFFTz","2D RMA Error")
                obj.isFail = true;
                return;
            end
        end

        function verifyReconstruction(obj)
            % Verify the reconstruction can continue

            if obj.scanner.method ~= "Linear"
                showErrorMessage(obj.im,"Must use linear SAR scan along the Y-axis to perform Uniform 1-D SAR 2-D RMA image reconstruction method!","2D RMA Error");
                obj.isFail = true;
                return
            end

            % Ensure array is colinear
            if max(diff([obj.ant.tx.xy_m(:,1);obj.ant.rx.xy_m(:,1)])) > 8*sqrt(eps)
                showErrorMessage(obj.im,"MIMO array must be colinear. Please disable necessary elements.","2D RMA Error");
                obj.isFail = true;
                return
            end

            % Ensure virtual array is uniform
            if mean(diff(obj.ant.vx.xyz_m(:,2),2)) > 8*sqrt(eps)
                showErrorMessage(obj.im,"Virtual antenna array is nonuniform! Change antenna positions.","2D RMA Error");
                obj.isFail = true;
                return
            end

            % And sar step size is correct
            if obj.scanner.yStep_m - mean(diff(obj.ant.vx.xyz_m(:,2)))*obj.ant.vx.numVx > 8*sqrt(eps)
                showErrorMessage(obj.im,"SAR step size is incorrect!");
                obj.isFail = true;
                return
            end

            if ~(isAppEPC(obj.im.app) || obj.ant.isEPC) && ~obj.isMult2Mono
                % If using MIMO Array Ensure multistatic-to-monostatic
                % approximation is employed
                showErrorMessage(obj.im,"Must use multistatic-to-monostatic approximation to use uniform method!","2D RMA Error");
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
                    mult2mono(obj);
                end

                reconstruct(obj);
                imXYZ_out = obj.imXYZ;
            else
                imXYZ_out = single(zeros(obj.im.numY,obj.im.numZ));
            end
        end

        function reconstruct(obj)
            % Reconstruct the image using the 2-D Range Migration Algorithm

            % sarData is of size (sar.numY, sar.numX, fmcw.ADCSamples)
            % Zero-Pad Data: s(y,k)
            sarDataPadded = obj.sarData;
            sarDataPadded = padarray(sarDataPadded,[floor((obj.nFFTy-size(obj.sarData,1))/2) 0],0,'pre');
            clear sarData

            % Compute Wavenumbers
            k = single(reshape(obj.k_vec,1,[]));
            L_y = obj.nFFTy * obj.yStep_m;
            dkY = 2*pi/L_y;
            kY = make_kX(obj,dkY,obj.nFFTy)';

            kZU = single(reshape(linspace(0,2*max(k) - 2*max(k)/obj.nFFTz,obj.nFFTz),1,[]));
            dkZU = kZU(2) - kZU(1);

            if obj.isGPU
                k = gpuArray(k);
                kY = gpuArray(kY);
                kZU = gpuArray(kZU);
                sarDataPadded = gpuArray(sarDataPadded);
            end

            kYU = repmat(kY,[1,obj.nFFTz]);
            kU = single(1/2 * sqrt(kY.^2 + kZU.^2));
            kZ = single(sqrt((4 * k.^2 - kY.^2) .* (4 * k.^2 > kY.^2)));

            % Compute Focusing Filter
            focusingFilter = exp(-1j * kZ * obj.z0_m);
            if obj.isAmplitudeFactor
                focusingFilter = kZ .* focusingFilter;
            end
            focusingFilter(4 * k.^2 < kY.^2) = 0;

            % Compute FFT across Y Dimension: S(kY,k)
            sarDataFFT = fftshift(fft(conj(sarDataPadded),obj.nFFTy,1),1)/obj.nFFTy;
            clear sarDataPadded sarData

            if obj.isGPU
                focusingFilter = gpuArray(focusingFilter);
            end

            % Stolt Interpolation
            try
                sarImageFFT = interpn(kY(:),k(:), sarDataFFT .* focusingFilter ,kYU,kU,obj.stolt_method,0);
            catch
                sarImageFFT = zeros(size(kU));
                for indkY = 1:size(kU,1)
                    tempS = squeeze(sarDataFFT(indkY,:) .* focusingFilter(indkY,:));
                    kZTemp = squeeze(kZ(indkY,:));
                    [kZTemp_unique,~,ind_c] = uniquetol(kZTemp);
                    tempS = accumarray(ind_c,tempS);
                    if length(kZTemp_unique) > 2
                        sarImageFFT(indkY,:) = gather(interp1(kZTemp_unique,tempS,kZU,obj.stolt_method,0));
                    end
                end
            end
            clear sarDataFFT focusingFilter kY k kYU kU kZ kZU

            % Recover Image by IFT: p(y,z)
            sarImage = single(ifftn(sarImageFFT));
            clear sarImageFFT focusingFilter

            % Declare Spatial Vectors
            y_m_temp = make_x(obj,obj.yStep_m,obj.nFFTy);
            z_m_temp = single(2*pi / (dkZU * obj.nFFTz) * (1:obj.nFFTz));

            % Interpolate Image
            obj.imXYZ = imageInterp2D(obj,y_m_temp,z_m_temp,sarImage,"yz");
        end

        function displayImage(obj)
            % Display the reconstructed y-z image

            displayImage2D(obj.im,obj.im.y_m,obj.im.z_m,"y (m)","z (m)");
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