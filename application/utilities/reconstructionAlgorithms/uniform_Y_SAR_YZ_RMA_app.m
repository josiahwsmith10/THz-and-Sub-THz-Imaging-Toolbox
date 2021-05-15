% uniform_Y_SAR_YZ_RMA_app see uniform_Y_SAR_YZ_RMA documentation
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

classdef uniform_Y_SAR_YZ_RMA_app < handle
    properties
        sarData
        
        nFFTy
        nFFTz
        
        y_m
        z_m
        
        imXYZ
        
        isGPU
        isAmplitudeFactor
        isFail
        
        k_vec
        z0_m
        yStep_m
    end
    
    methods
        function obj = uniform_Y_SAR_YZ_RMA_app(app,im)
            obj = update(obj,app,im);
        end
        
        function obj = update(obj,app,im)
            obj = getParameters(obj,app,im);
            obj = verifyParameters(obj,app);
            obj = verifyReconstruction(obj,app);
            obj = updateApp(obj,app);
        end
        
        function obj = getParameters(obj,app,im)
            obj.nFFTy = im.nFFTy;
            obj.nFFTz = im.nFFTz;
            
            obj.y_m = im.y_m;
            obj.z_m = im.z_m;
            
            obj.sarData = app.target.sarData;
            
            obj.isGPU = im.isGPU;
            obj.isAmplitudeFactor = app.target.isAmplitudeFactor;
            
            obj.k_vec = app.fmcw.k;
            obj.z0_m = app.ant.tx.z0_m;
            obj.yStep_m = app.sar.yStep_m;
        end
        
        function obj = verifyParameters(obj,app)
            obj.isFail = false;
            
            kZU = single(reshape(linspace(0,2*max(obj.k_vec) - 2*max(obj.k_vec)/obj.nFFTz,obj.nFFTz),1,1,[]));
            dkZU = kZU(2) - kZU(1);
            y_m_temp = make_x(obj,app.sar.yStep_m,obj.nFFTy);
            z_m_temp = single(2*pi / (dkZU * obj.nFFTz) * (1:obj.nFFTz));
            
            if max(abs(obj.y_m)) > max(abs(y_m_temp))
                uiconfirm(app.UIFigure,"X Max (m) is too large for Number of Y FFT Points. Decrease Y Max (m) or increase Number of Y FFT Points",'Imaging Parameter Error!',...
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
            if app.sar.method ~= "Linear"
                uiconfirm(app.UIFigure,"Must use linear SAR scan along the Y-axis to perform Uniform 1-D SAR 2-D RMA image reconstruction method!",'SAR Scenario Error!',...
                    "Options",{'OK'},'Icon','warning');
                obj.isFail = true;
                return
            end
            
            % Ensure array is colinear
            if max(diff([app.ant.tx.xy_m(:,1);app.ant.rx.xy_m(:,1)])) > 8*eps
                uiconfirm(app.UIFigure,"MIMO array must be colinear. Please disable necessary elements.",'Array Topology Error!',...
                    "Options",{'OK'},'Icon','warning');
                obj.isFail = true;
                return
            end
            
            % Ensure virtual array is uniform
            if mean(diff(app.ant.vx.xyz_m(:,2),2)) > eps
                uiconfirm(app.UIFigure,"Virtual antenna array is nonuniform! Change antenna positions.",'Array Topology Error!',...
                    "Options",{'OK'},'Icon','warning');
                obj.isFail = true;
                return
            end
            
            % And sar step size is correct
            if app.sar.yStep_m - mean(diff(app.ant.vx.xyz_m(:,2)))*app.ant.vx.numVx > 8*eps
                uiconfirm(app.UIFigure,"SAR step size is incorrect!",'SAR Scenario Error!',...
                    "Options",{'OK'},'Icon','warning');
                obj.isFail = true;
                return
            end
            if app.sar.isMIMO && ~app.isMult2MonoCheckBox.Value
                % If using MIMO Array
                % Ensure multistatic-to-monostatic approximation is employed
                uiconfirm(app.UIFigure,"SAR step size is incorrect!",'SAR Scenario Error!',...
                    "Options",{'OK'},'Icon','warning');
                obj.isFail = true;
                return
            end
            
            % Everything is okay to continue
            obj.yStep_m = app.sar.yStep_m/app.ant.vx.numVx;
        end
        
        function obj = updateApp(obj,app)
            kZU = single(reshape(linspace(0,2*max(obj.k_vec) - 2*max(obj.k_vec)/obj.nFFTz,obj.nFFTz),1,1,[]));
            dkZU = kZU(2) - kZU(1);
            y_m_temp = double(make_x(obj,obj.yStep_m,obj.nFFTy));
            z_m_temp = double(2*pi / (dkZU * obj.nFFTz) * (0:obj.nFFTz-1));
            
            app.YMinmEditField_im_5.Value = min(y_m_temp);
            app.YMaxmEditField_im_5.Value = max(y_m_temp);
            
            app.ZMinmEditField_im_3.Value = min(z_m_temp);
            app.ZMaxmEditField_im_3.Value = max(z_m_temp);
        end
        
        function [obj,imXYZ_out] = computeReconstruction(obj,app,im)
            obj = update(obj,app,im);
            if ~obj.isFail
                if app.sar.isMIMO && app.isMult2MonoCheckBox.Value
                    obj = mult2mono_app(obj,app);
                end
                
                obj = reconstruct(obj,app);
            end
            imXYZ_out = obj.imXYZ;
        end
        
        function obj = reconstruct(obj,app)
            % Start the Progress Bar
            d = uiprogressdlg(app.UIFigure,'Title','Please Wait',...
                'Message','Reconstructing Image using Uniform RMA ');
            
            % sarData is of size (sar.numY, sar.numX, fmcw.ADCSamples)
            d.Value = 1/10;
            
            % Zero-Pad Data: s(y,k)
            sarDataPadded = obj.sarData;
            sarDataPadded = padarray(sarDataPadded,[floor((obj.nFFTy-size(obj.sarData,1))/2) 0],0,'pre');
            clear sarData
            d.Value = 2/10;
            
            % Compute Wavenumbers
            k = single(reshape(obj.k_vec,1,[]));
            L_y = obj.nFFTy * obj.yStep_m;
            dkY = 2*pi/L_y;
            kY = make_kX(obj,dkY,obj.nFFTy)';
            
            kZU = single(reshape(linspace(0,2*max(k) - 2*max(k)/obj.nFFTz,obj.nFFTz),1,[]));
            dkZU = kZU(2) - kZU(1);
            
            if obj.isGPU
                reset(gpuDevice)
                k = gpuArray(k);
                kY = gpuArray(kY);
                kZU = gpuArray(kZU);
                sarDataPadded = gpuArray(sarDataPadded);
            end
            
            kYU = repmat(kY,[1,obj.nFFTz]);
            kU = single(1/2 * sqrt(kY.^2 + kZU.^2));
            kZ = single(sqrt((4 * k.^2 - kY.^2) .* (4 * k.^2 > kY.^2)));
            d.Value = 3/10;
            
            % Compute Focusing Filter
            focusingFilter = exp(-1j * kZ * obj.z0_m);
            if obj.isAmplitudeFactor
                focusingFilter = kZ .* focusingFilter;
            end
            focusingFilter(4 * k.^2 < kY.^2) = 0;
            d.Value = 4/10;
            
            % Compute FFT across Y & X Dimensions: S(kY,k)
            sarDataFFT = fftshift(fft(conj(sarDataPadded),obj.nFFTy,1),1)/obj.nFFTy;
            clear sarDataPadded sarData
            
            if obj.isGPU
                focusingFilter = gpuArray(focusingFilter);
            end
            d.Value = 5/10;
            
            % Stolt Interpolation
            try
                sarImageFFT = interpn(kY(:),k(:), sarDataFFT .* focusingFilter ,kYU,kU,'linear',0);
            catch
                sarImageFFT = zeros(size(kU));
                for indkY = 1:size(kU,1)
                    tempS = squeeze(sarDataFFT(indkY,:) .* focusingFilter(indkY,:));
                    kZTemp = squeeze(kZ(indkY,:));
                    [kZTemp_unique,~,ind_c] = uniquetol(kZTemp);
                    tempS = accumarray(ind_c,tempS);
                    if length(kZTemp_unique) > 2
                        sarImageFFT(indkY,:) = gather(interp1(kZTemp_unique,tempS,kZU,'nearest',0));
                    end
                end
            end
            clear sarDataFFT focusingFilter kY k kYU kU kZ kZU
            d.Value = 6/10;
            
            % Recover Image by IFT: p(y,z)
            sarImage = single(abs(ifftn(sarImageFFT)));
            clear sarImageFFT focusingFilter
            d.Value = 7/10;
            
            % Declare Spatial Vectors
            y_m_temp = make_x(obj,obj.yStep_m,obj.nFFTy);
            z_m_temp = single(2*pi / (dkZU * obj.nFFTz) * (1:obj.nFFTz));
            d.Value = 9/10;
            
            [Y,Z] = ndgrid(obj.y_m(:),obj.z_m(:));
            obj.imXYZ = single(gather(interpn(y_m_temp(:),z_m_temp(:),sarImage,Y,Z,'nearest',0)));
            d.Value = 10/10;
        end
        
        function displayImage(obj,im)
            displayImage2D_app(im,im.y_m,im.z_m,"y (m)","z (m)");
        end
        
        function x = make_x(obj,xStep_m,nFFTx)
            x = xStep_m * (-(nFFTx-1)/2 : (nFFTx-1)/2);
            x = single(x);
        end
        
        function kX = make_kX(obj,dkX,nFFTx)
            if mod(nFFTx,2)==0
                kX = dkX * ( -nFFTx/2 : nFFTx/2-1 );
            else
                kX = dkX * ( -(nFFTx-1)/2 : (nFFTx-1)/2 );
            end
            kX = single(kX);
        end
    end
end