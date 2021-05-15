% uniform_theta_CSAR_XZ_PFA_app see uniform_theta_CSAR_XZ_PFA
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

classdef uniform_theta_CSAR_XZ_PFA_app < handle
    properties
        sarData
        
        nFFTx
        nFFTz
        thetaUpsampleFactor
        
        x_m
        z_m
        
        imXYZ
        
        isGPU
        isAmplitudeFactor
        isFail
        
        theta_rad_vec
        k_vec
        R0_m
        xStep_m
    end
    
    methods
        function obj = uniform_theta_CSAR_XZ_PFA_app(app,im)
            obj = getParameters(obj,app,im);
            obj = updateApp(obj,app);
        end
        
        function obj = update(obj,app,im)
            obj = getParameters(obj,app,im);
            obj = verifyParameters(obj,app);
            obj = verifyReconstruction(obj,app);
            obj = updateApp(obj,app);
        end
        
        function obj = getParameters(obj,app,im)
            obj.nFFTx = im.nFFTx;
            obj.nFFTz = im.nFFTz;
            obj.thetaUpsampleFactor = app.ThetaUpsampleFactorEditField_2.Value;
            
            obj.x_m = im.x_m;
            obj.z_m = im.z_m;
            
            obj.sarData = app.target.sarData;
            
            obj.isGPU = im.isGPU;
            obj.isAmplitudeFactor = app.target.isAmplitudeFactor;
            
            obj.theta_rad_vec = app.sar.theta_rad;
            obj.k_vec = app.fmcw.k;
            obj.R0_m = app.ant.tx.z0_m;
            obj.xStep_m = app.sar.xStep_m;
        end
        
        function obj = verifyParameters(obj,app)
            obj.isFail = false;
            
            [x_m_temp,z_m_temp] = getTempXZ(obj);
            
            if max(abs(obj.x_m)) > max(abs(x_m_temp))
                uiconfirm(app.UIFigure,"X Max (m) is too large for Number of X FFT Points. Decrease X Max (m) or increase Number of X FFT Points",'Imaging Parameter Error!',...
                    "Options",{'OK'},'Icon','warning');
                obj.isFail = true;
                return;
            end
            
            if obj.thetaUpsampleFactor*app.sar.numTheta > obj.nFFTx
                uiconfirm(app.UIFigure,"θ Upsample Factor is too large for Number of X FFT Points. Decrease θ Upsample Factor or increase Number of X FFT Points",'Imaging Parameter Error!',...
                    "Options",{'OK'},'Icon','warning');
                obj.isFail = true;
                return;
            end
            
            if max(obj.z_m) > max(z_m_temp)
                uiconfirm(app.UIFigure,"Z Max (m) is too large for Number of Z FFT Points. Decrease Z Max (m) or increase Number of Z FFT Points",'Imaging Parameter Error!',...
                    "Options",{'OK'},'Icon','warning');
                obj.isFail = true;
                return;
            end
        end
        
        function obj = verifyReconstruction(obj,app)
            if app.sar.method ~= "Circular"
                uiconfirm(app.UIFigure,"Must use 1-D θ Circular CSAR scan to use 1-D CSAR 2-D PFA image reconstruction method!",'SAR Scenario Error!',...
                    "Options",{'OK'},'Icon','warning');
                obj.isFail = true;
                return
            end
            
            % Ensure single element array
            if app.ant.tx.numTx ~= 1 || app.ant.rx.numRx ~= 1
                uiconfirm(app.UIFigure,"Array must have only 1 Tx and 1 Rx. Please disable necessary elements.",'Array Topology Error!',...
                    "Options",{'OK'},'Icon','warning');
                app.SARMethodDropDown.Value = "-";
                obj.isFail = true;
                return
            end
        end
        
        function obj = updateApp(obj,app)
            [x_m_temp,z_m_temp] = getTempXZ(obj);
            
            app.XMinmEditField_im_5.Value = min(x_m_temp);
            app.XMaxmEditField_im_5.Value = max(x_m_temp);
            
            app.ZMinmEditField_im_5.Value = min(z_m_temp);
            app.ZMaxmEditField_im_5.Value = max(z_m_temp);
        end
        
        function [x_m_temp,z_m_temp] = getTempXZ(obj)
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
        
        function [obj,imXYZ_out] = computeReconstruction(obj,app,im)
            obj = update(obj,app,im);
            
            if ~obj.isFail
                obj = reconstruct(obj,app);
            end
            imXYZ_out = obj.imXYZ;
        end
        
        function obj = reconstruct(obj,app)
            % Start the Progress Bar
            d = uiprogressdlg(app.UIFigure,'Title','Please Wait',...
                'Message','Reconstructing Image using Uniform PFA ');
            
            % sarData is of size (sar.numTheta, fmcw.ADCSamples)
            d.Value = 1/10;
            
            % Zero-Pad Data: s(theta,k)
            sarDataPadded = obj.sarData;
            d.Value = 2/10;
            
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
            d.Value = 3/10;
            
            % Upsample data
            if obj.thetaUpsampleFactor > 1
                theta_radUp = single(reshape(linspace(theta_rad(1),theta_rad(end),length(theta_rad)*obj.thetaUpsampleFactor),[],1));
                [T,K] = ndgrid(theta_radUp(:),reshape(single(1:size(sarDataPadded,2)),1,[]));
                sarDataPadded = interpn(theta_rad(:),single(1:size(sarDataPadded,2))',sarDataPadded,T,K);
                clear T K
                theta_rad = theta_radUp;
            end
            d.Value = 3.5/10;
            
            % Compute Azimuth Filter H(Theta,k)
            azimuthFilterFFT = single(fft(exp(1j*kR*obj.R0_m.*cos(theta_rad)),[],1));
            clear kR
            d.Value = 4/10;
            
            % Compute Azimuth Filtered Data: p(theta,k) = IFT_Theta[ S(Theta,k) * H*(Theta,k)]
            azimuthFiltered = ifft(fft(sarDataPadded,[],1) .* conj(azimuthFilterFFT),[],1);
            clear azimuthFilterFFT sarDataPadded
            d.Value = 5/10;
            
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
            d.Value = 6/10;
            
            % Recover Image by IFT: p(x,z)
            sarImage = single(abs(ifftshift(ifftshift(ifftn(sarImageFFT),1),2)));
            clear sarImageFFT
            d.Value = 7/10;
            
            % Declare Spatial Vectors
            x_m_temp = make_x(obj,2*pi/(dkXU*obj.nFFTx),obj.nFFTx);
            z_m_temp = make_x(obj,2*pi/(dkZU*obj.nFFTz),obj.nFFTz);
            d.Value = 9/10;
            
            [X,Z] = ndgrid(obj.x_m(:),obj.z_m(:));
            obj.imXYZ = single(gather(interpn(x_m_temp(:),z_m_temp(:),sarImage,X,Z,"nearest",0)));
            d.Value = 10/10;
        end
        
        function displayImage(obj,im)
            displayImage2D_app(im,im.x_m,im.z_m,"x (m)","z (m)");
        end
        
        function x = make_x(obj,xStep_m,nFFTx)
            x = xStep_m * (-(nFFTx-1)/2 : (nFFTx-1)/2);
            x = single(x);
        end
    end
end