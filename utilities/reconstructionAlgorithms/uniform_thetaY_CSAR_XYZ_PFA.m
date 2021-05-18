% uniform thetaY_CSAR_XYZ_PFA is a reconstructor class that performs
% 3-D Polar Formatting Algorithm image reconstruction. The synthetic
% aperture must span across the y and theta dimensions, forming a
% cylindrical aperture, and the target can be a 1-D, 2-D, or 3-D target
% in x-y-z space
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

classdef uniform_thetaY_CSAR_XYZ_PFA < handle
    properties
        sarData             % Computed beat signal
        
        nFFTx = 512         % Number of FFT points along the x-dimension, when using FFT-based reconstruction algorithms
        nFFTy = 512         % Number of FFT points along the y-dimension, when using FFT-based reconstruction algorithms
        nFFTz = 512         % Number of FFT points along the z-dimension, when using FFT-based reconstruction algorithms
        thetaUpsampleFactor = 1 % Upsample factor on the theta dimension of the SAR data operated in circular or cylindrical mode
        
        x_m                 % Reconstructed image x axis
        y_m                 % Reconstructed image y axis
        z_m                 % Reconstructed image z axis
        
        imXYZ               % Reconstructed image
        
        isGPU               % Boolean whether or not to use the GPU for image reconstruction
        isAmplitudeFactor   % Boolean whether or not to include the amplitude factor in the image reconstruction process
        isFail = false      % Boolean whether or not the reconstruction has failed
        isMult2Mono = false   % Boolean whether or not to use the multistatic-to-monostatic approximation
        
        zRef_m = 0.25       % z location of reference plane for multistatic-to-monostatic approximation
        theta_rad_vec       % Angular scan locations vector
        k_vec               % Instantaneous wavenumber vector
        R0_m                % Radius of scan
        xStep_m = 1e-3      % Step size along the x-dimension to move the antenna array in meters
        yStep_m = 8e-3      % Step size along the y-dimension to move the antenna array in meters
        
        wav                 % A THzWaveformParameters object handle
        ant                 % A THzAntennaArray object handle
        scanner             % A THzScanner object handle
        target              % A THzTarget object handle
        im                  % A THzImageReconstruction object handle
    end
    
    methods
        function obj = uniform_thetaY_CSAR_XYZ_PFA(im)
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
            obj.nFFTy = obj.im.nFFTy;
            obj.nFFTz = obj.im.nFFTz;
            obj.thetaUpsampleFactor = obj.im.thetaUpsampleFactor;
            
            obj.x_m = obj.im.x_m;
            obj.y_m = obj.im.y_m;
            obj.z_m = obj.im.z_m;
            
            obj.sarData = obj.target.sarData;
            
            obj.isGPU = obj.im.isGPU;
            obj.isAmplitudeFactor = obj.target.isAmplitudeFactor;
            obj.isMult2Mono = obj.im.isMult2Mono;
            
            obj.zRef_m = obj.im.zRef_m;
            obj.theta_rad_vec = obj.scanner.theta_rad;
            obj.k_vec = obj.wav.k;
            obj.R0_m = obj.ant.z0_m;
            obj.xStep_m = obj.scanner.xStep_m;
            obj.yStep_m = obj.scanner.yStep_m;
            
            if obj.im.isApp
                obj.thetaUpsampleFactor = obj.im.app.ThetaUpsampleFactorEditField.Value;
            end
        end
        
        function verifyParameters(obj)
            % Verify the parameters allow for imaging
            
            [x_m_temp,y_m_temp,z_m_temp] = getTempXYZ(obj);
            
            if obj.im.isApp
                app = obj.im.app;
                app.XMinmEditField_im_4.Value = min(x_m_temp);
                app.XMaxmEditField_im_4.Value = max(x_m_temp);
                
                app.YMinmEditField_im_6.Value = min(y_m_temp);
                app.YMaxmEditField_im_6.Value = max(y_m_temp);
                
                app.ZMinmEditField_im_4.Value = min(z_m_temp);
                app.ZMaxmEditField_im_4.Value = max(z_m_temp);
            end
            
            if max(abs(obj.x_m)) > max(abs(x_m_temp))
                showErrorMessage(obj.im,"xMax_m is too large for nFFTx. Decrease xMax_m or increase nFFTx","3D PFA Error")
                obj.isFail = true;
                return;
            end
            if obj.thetaUpsampleFactor*obj.scanner.numTheta > obj.nFFTx
                showErrorMessage(obj.im,"thetaUpsampleFactor is too large for nFFTx. Decrease thetaUpsampleFactor or increase nFFTx","3D PFA Error")
                obj.isFail = true;
                return;
            end
            if max(obj.y_m) > max(y_m_temp)
                showErrorMessage(obj.im,"yMax_m is too large for nFFTy. Decrease yMax_m or increase nFFTy","3D PFA Error")
                obj.isFail = true;
                return;
            end
            if max(obj.z_m) > max(z_m_temp)
                showErrorMessage(obj.im,"zMax_m is too large for nFFTz. Decrease zMax_m or increase nFFTz","3D PFA Error")
                obj.isFail = true;
                return;
            end
        end
        
        function verifyReconstruction(obj)
            % Verify the reconstruction can continue
            
            if obj.scanner.method ~= "Cylindrical"
                showErrorMessage(obj.im,"Must use 2-D θY Cylindrical CSAR scan to use 2-D CSAR 3-D PFA image reconstruction method!","3D PFA Error");
                obj.isFail = true;
                return
            end
            
            % Ensure array is colinear
            if max(diff([obj.ant.tx.xy_m(:,1);obj.ant.rx.xy_m(:,1)])) > 8*sqrt(eps)
                showErrorMessage(obj.im,"MIMO array must be colinear. Please disable necessary elements.","3D PFA Error");
                obj.isFail = true;
                return
            end
            
            % Ensure virtual array is uniform
            if mean(diff(obj.ant.vx.xyz_m(:,2),2)) > sqrt(eps)
                showErrorMessage(obj.im,"Virtual antenna array is nonuniform! Change antenna positions.","3D PFA Error");
                obj.isFail = true;
                return
            end
            
            % And sar step size is correct
            if obj.scanner.yStep_m - mean(diff(obj.ant.vx.xyz_m(:,2)))*obj.ant.vx.numVx > 8*sqrt(eps)
                showErrorMessage(obj.im,"SAR step size is incorrect!","3D PFA Error");
                obj.isFail = true;
                return
            end
            
            if ~(isAppEPC(obj.im.app) || obj.ant.isEPC) && ~obj.isMult2Mono
                % If using MIMO Array
                % Ensure multistatic-to-monostatic approximation is employed
                showErrorMessage(obj.im,"Must use multistatic-to-monostatic approximation to use uniform method!","3D PFA Error");
                obj.isFail = true;
                return
            end
            
            % Everything is okay to continue
            obj.yStep_m = obj.scanner.yStep_m/obj.ant.vx.numVx;
            obj.isFail = false;
        end
        
        function [x_m_temp,y_m_temp,z_m_temp] = getTempXYZ(obj)
            % Compute the axes of the recovered image for parameter
            % verification
            
            theta_rad = single(reshape(obj.theta_rad_vec,1,[]));
            
            % Compute Wavenumbers
            k = single(reshape(obj.k_vec,1,1,[]));
            
            L_y = obj.nFFTy * obj.yStep_m;
            dkY = 2*pi/L_y;
            kY = make_kX(obj,dkY,obj.nFFTy)';
            
            kR = single(sqrt(4 * k.^2 - kY.^2) .* (4 * k.^2 > kY.^2));
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
            y_m_temp = double(make_x(obj,obj.yStep_m,obj.nFFTy));
            z_m_temp = double(make_x(obj,2*pi/(dkZU*obj.nFFTz),obj.nFFTz));
        end
        
        function imXYZ_out = computeReconstruction(obj)
            % Update the reconstruction algorithm and attempt the
            % reconstruction
            
            update(obj);
            
            if ~obj.isFail
                if ~obj.ant.isEPC && obj.isMult2Mono
                    mult2mono(obj);
                end
                
                %                 reconstruct4segs(obj);
                reconstruct(obj);
                imXYZ_out = obj.imXYZ;
            else
                imXYZ_out = single(zeros(obj.im.numX,obj.im.numY,obj.im.numZ));
            end
        end
        
        %         function reconstruct4segs(obj)
        %             theta_rad_orig = obj.theta_rad_vec;
        %             sarData_orig = obj.sarData;
        %
        %             % 1st segment
        %             obj.theta_rad_vec = theta_rad_orig(1:end/4);
        %             obj.sarData = sarData_orig(:,1:end/4,:);
        %             reconstruct(obj);
        %             im1 = obj.imXYZ;
        %
        %             % 2nd segment
        %             obj.theta_rad_vec = theta_rad_orig(end/4+1:end/2);
        %             obj.sarData = sarData_orig(:,end/4+1:end/2,:);
        %             reconstruct(obj);
        %             im2 = obj.imXYZ;
        %
        %             % 3rd segment
        %             obj.theta_rad_vec = theta_rad_orig(end/2+1:end*3/4);
        %             obj.sarData = sarData_orig(:,end/2+1:end*3/4,:);
        %             reconstruct(obj);
        %             im3 = obj.imXYZ;
        %
        %             % 4th segment
        %             obj.theta_rad_vec = theta_rad_orig(end*3/4+1:end);
        %             obj.sarData = sarData_orig(:,end*3/4+1:end,:);
        %             reconstruct(obj);
        %             im4 = obj.imXYZ;
        %
        %             obj.imXYZ = im1 + im2 + im3 + im4;
        %         end
        
        function reconstruct(obj)
            % Reconstruct the image using the 3-D Polar Formatting
            % Algorithm
            
            % sarData is of size (sar.numY, sar.numTheta, fmcw.ADCSamples)
            % Zero-Pad Data: s(y,theta,k)
            sarDataPadded = obj.sarData;
            sarDataPadded = padarray(sarDataPadded,[floor((obj.nFFTy-size(obj.sarData,1))/2) 0],0,'pre');
            
            theta_rad = single(reshape(obj.theta_rad_vec,1,[]));
            
            % Compute Wavenumbers
            k = single(reshape(obj.k_vec,1,1,[]));
            
            L_y = obj.nFFTy * obj.yStep_m;
            dkY = 2*pi/L_y;
            kY = make_kX(obj,dkY,obj.nFFTy)';
            
            kR = single(sqrt(4 * k.^2 - kY.^2) .* (4 * k.^2 > kY.^2));
            kX = kR.*cos(theta_rad);
            kZ = kR.*sin(theta_rad);
            
            kXmax = max(kX,[],'all');
            kZmax = max(kZ,[],'all');
            kXmin = min(kX,[],'all');
            kZmin = min(kZ,[],'all');
            clear kX kZ
            
            kYU = repmat(kY,[1,obj.nFFTx,obj.nFFTz]);
            kXU = single(reshape(linspace(kXmin,kXmax,obj.nFFTx),1,[],1));
            kZU = single(reshape(linspace(kZmin,kZmax,obj.nFFTz),1,1,[]));
            dkXU = kXU(2) - kXU(1);
            dkZU = kZU(2) - kZU(1);
            
            % Declare Uniform Grid for Interpolation
            theta_radU = repmat(atan2(kZU,kXU),[obj.nFFTy,1,1]);
            kU = 1/2*sqrt(kXU.^2 + kZU.^2 + kY.^2);
            clear kXU kZU
            
            % Upsample data
            if obj.thetaUpsampleFactor > 1
                theta_radUp = single(reshape(linspace(theta_rad(1),theta_rad(end),length(theta_rad)*obj.thetaUpsampleFactor),1,[]));
                [Y,T,K] = ndgrid(single(1:size(sarDataPadded,1))',theta_radUp,reshape(single(1:size(sarDataPadded,3)),1,1,[]));
                sarDataPadded = interpn(single(1:size(sarDataPadded,1))',theta_rad(:),single(1:size(sarDataPadded,3))',sarDataPadded,Y,T,K);
                clear Y T K
                theta_rad = theta_radUp;
            end
            
            % Compute Azimuth Filter H(kY,Theta,k)
            azimuthFilterFFT = single(fft(exp(1j*kR*obj.R0_m.*cos(theta_rad)),[],2));
            clear kR
            
            % Compute Azimuth Filtered Data: p(kY,theta,k) = IFT_Theta[ S(kY,Theta,k) * H*(kY,Theta,k)]
            azimuthFiltered = ifft(fftshift(fft(fft(sarDataPadded,[],2),obj.nFFTy,1),1) .* conj(azimuthFilterFFT),[],2);
            clear azimuthFilterFFT sarDataPadded
            
            if obj.isGPU
                reset(gpuDevice)
                kY = gpuArray(kY);
                theta_rad = gpuArray(theta_rad);
                k = gpuArray(k);
                azimuthFiltered = gpuArray(azimuthFiltered);
                kYU = gpuArray(kYU);
                theta_radU = gpuArray(theta_radU);
                kU = gpuArray(kU);
            end
            
            % Stolt Interpolation
            sarImageFFT = interpn(kY(:),theta_rad(:),k(:), azimuthFiltered ,kYU,theta_radU,kU,'linear',0);
            clear azimuthFiltered kY theta_rad k azimuthfiltered kYU theta_radU kU
            
            if obj.isGPU
                sarImageFFT = gather(sarImageFFT);
                reset(gpuDevice);
                sarImageFFT = gpuArray(sarImageFFT);
            end
            
            % Recover Image by IFT: p(y,x,z)
            sarImage = single(ifftshift(ifftshift(ifftn(sarImageFFT),2),3));
            clear sarImageFFT
            
            % Reorient Image: p(x,y,z)
            sarImage = permute(sarImage,[2,1,3]);
            
            % Declare Spatial Vectors
            x_m_temp = make_x(obj,2*pi/(dkXU*obj.nFFTx),obj.nFFTx);
            y_m_temp = make_x(obj,obj.yStep_m,obj.nFFTy);
            z_m_temp = make_x(obj,2*pi/(dkZU*obj.nFFTz),obj.nFFTz);
            
            [X,Y,Z] = ndgrid(obj.x_m(:),obj.y_m(:),obj.z_m(:));
            obj.imXYZ = single(gather(interpn(x_m_temp(:),y_m_temp(:),z_m_temp(:),sarImage,X,Y,Z,"nearest",0)));
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