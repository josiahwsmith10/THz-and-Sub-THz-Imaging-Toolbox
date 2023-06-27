% uniform_Y_SAR_Y_FFT is a reconstructor class that performs a 1-D FFT
% method image reconstruction. The synthetic aperture must span the
% y-dimension and the target can be a 1-D or 2-D target in the y-z plane.
% The reconstructed image is a slice along the y-dimension at the
% z-coordinate zSlice_m
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

classdef uniform_Y_SAR_Y_FFT < handle
    properties
        sarData             % Computed beat signal

        nFFTy = 512         % Number of FFT points along the y-dimension, when using FFT-based reconstruction algorithms

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
        yStep_m = 8e-3      % Step size along the y-dimension to move the antenna array in meters

        im_method           % Interpolation method for image interpolation

        wav                 % A THzWaveformParameters object handle
        ant                 % A THzAntennaArray object handle
        scanner             % A THzScanner object handle
        target              % A THzTarget object handle
        im                  % A THzImageReconstruction object handle
    end

    methods
        function obj = uniform_Y_SAR_Y_FFT(im)
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

            obj.nFFTy = obj.im.nFFTy;

            obj.y_m = obj.im.y_m;

            obj.sarData = obj.target.sarData;

            obj.isGPU = obj.im.isGPU;
            obj.isAmplitudeFactor = obj.target.isAmplitudeFactor;
            obj.isMult2Mono = obj.im.isMult2Mono;

            obj.zSlice_m = obj.im.zSlice_m;
            obj.zRef_m = obj.im.zRef_m;
            obj.k_vec = obj.wav.k;
            obj.z0_m = obj.ant.z0_m;
            obj.yStep_m = obj.scanner.yStep_m;

            if obj.im.isApp
                obj.im.zSlice_m = obj.im.app.ZSliceEditField_3.Value;
            end
        end

        function verifyParameters(obj)
            % Verify the parameters allow for imaging

            y_m_temp = double(make_x(obj,obj.yStep_m,obj.nFFTy));

            if obj.im.isApp
                app = obj.im.app;
                app.YMinmEditField_im_4.Value = min(y_m_temp);
                app.YMaxmEditField_im_4.Value = max(y_m_temp);
            end

            if max(abs(obj.y_m)) > max(abs(y_m_temp))
                showErrorMessage(obj.im,"yMax_m is too large for nFFTy. Decrease yMax_m or increase nFFTy","1D FFT Error")
                obj.isFail = true;
                return;
            end
        end

        function verifyReconstruction(obj)
            % Verify the reconstruction can continue

            if obj.scanner.method ~= "Linear"
                showErrorMessage(obj.im,"Must use linear SAR scan along the Y-axis to perform Uniform 1-D SAR 1-D FFT image reconstruction method!","1D FFT Error");
                obj.isFail = true;
                return
            end

            % Ensure array is colinear
            if max(diff([obj.ant.tx.xy_m(:,1);obj.ant.rx.xy_m(:,1)])) > 8*sqrt(eps)
                showErrorMessage(obj.im,"MIMO array must be colinear. Please disable necessary elements.","1D FFT Error");
                obj.isFail = true;
                return
            end

            % Ensure virtual array is uniform
            if mean(diff(obj.ant.vx.xyz_m(:,2),2)) > sqrt(eps)
                showErrorMessage(obj.im,"Virtual antenna array is nonuniform! Change antenna positions.","1D FFT Error");
                obj.isFail = true;
                return
            end

            % And sar step size is correct
            if obj.scanner.yStep_m - mean(diff(obj.ant.vx.xyz_m(:,2)))*obj.ant.vx.numVx > 8*sqrt(eps)
                showErrorMessage(obj.im,"SAR step size is incorrect!","1D FFT Error");
                obj.isFail = true;
                return
            end

            if ~obj.ant.isEPC && ~obj.isMult2Mono
                % If using MIMO Array Ensure multistatic-to-monostatic
                % approximation is employed
                showErrorMessage(obj.im,"Must use multistatic-to-monostatic approximation to use uniform method!","1D FFT Error");
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
            % Reconstruct the image using the 1-D FFT Method

            % Start the Progress Bar sarData is of size (sar.numY,
            % fmcw.ADCSamples) Zero-Pad Data: s(y,k)
            sarDataPadded = obj.sarData;
            sarDataPadded = padarray(sarDataPadded,[floor((obj.nFFTy-size(obj.sarData,1))/2) 0],0,'pre');
            clear sarData

            % Compute Wavenumbers
            k = single(reshape(obj.k_vec,1,[]));
            L_y = obj.nFFTy * obj.yStep_m;
            dkY = 2*pi/L_y;
            kY = make_kX(obj,dkY,obj.nFFTy)';

            if obj.isGPU
                reset(gpuDevice)
                k = gpuArray(k);
                kY = gpuArray(kY);
                sarDataPadded = gpuArray(sarDataPadded);
            end

            kZ = single(sqrt((4 * k.^2 - kY.^2) .* (4 * k.^2 > kY.^2)));

            % Compute Focusing Filter
            focusingFilter = exp(-1j * kZ * (obj.zSlice_m - obj.z0_m));
            focusingFilter(4 * k.^2 < kY.^2) = 0;

            % Compute FFT across Y & X Dimensions: S(kY,k)
            sarDataFFT = fftshift(fft((sarDataPadded),obj.nFFTy,1),1)/obj.nFFTy;
            clear sarDataPadded sarData

            if obj.isGPU
                focusingFilter = gpuArray(focusingFilter);
            end

            sarImageFFT = sum(sarDataFFT .* focusingFilter,2);
            clear sarDataFFT focusingFilter kY k kZ

            % Recover Image by IFT: p(y)
            sarImage = single(ifftn(sarImageFFT));
            clear sarImageFFT focusingFilter

            % Declare Spatial Vectors
            y_m_temp = make_x(obj,obj.yStep_m,obj.nFFTy);

            obj.imXYZ = single(gather(interpn(y_m_temp(:),sarImage,obj.y_m(:),obj.im_method,0)));
            if obj.isGPU
                reset(gpuDevice);
            end
        end

        function displayImage(obj)
            % Display the y image

            displayImage1D(obj.im,obj.im.y_m,'y (m)');
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