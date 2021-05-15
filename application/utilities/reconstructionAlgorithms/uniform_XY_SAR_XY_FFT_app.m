% uniform_XY_SAR_XY_FFT_app see uniform_XY_SAR_XY_FFT documentation
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

classdef uniform_XY_SAR_XY_FFT_app < handle
    properties
        sarData
        
        nFFTx
        nFFTy
        
        x_m
        y_m
        
        imXYZ
        
        isGPU
        isAmplitudeFactor
        isFail
        
        k_vec
        z0_m
        z_target_m
        xStep_m
        yStep_m
    end
    
    methods
        function obj = uniform_XY_SAR_XY_FFT_app(app,im)
            obj = update(obj,app,im);
        end
        
        function obj = update(obj,app,im)
            obj = getParameters(obj,app,im);
            obj = verifyParameters(obj,app);
            obj = verifyReconstruction(obj,app);
            obj = updateApp(obj,app);
        end
        
        function obj = getParameters(obj,app,im)
            obj.nFFTx = im.nFFTx;
            obj.nFFTy = im.nFFTy;
            
            obj.x_m = im.x_m;
            obj.y_m = im.y_m;
            
            obj.sarData = app.target.sarData;
            
            obj.isGPU = im.isGPU;
            obj.isAmplitudeFactor = app.target.isAmplitudeFactor;
            
            obj.k_vec = app.fmcw.k;
            obj.z0_m = app.ant.tx.z0_m;
            obj.z_target_m = app.ZSliceEditField_2.Value;
            obj.xStep_m = app.sar.xStep_m;
            obj.yStep_m = app.sar.yStep_m;
        end
        
        function obj = verifyParameters(obj,app)
            obj.isFail = false;
            
            x_m_temp = make_x(obj,app.sar.xStep_m,obj.nFFTx);
            y_m_temp = make_x(obj,app.sar.yStep_m,obj.nFFTy);
            
            % Resize Image
            if max(abs(obj.x_m)) > max(abs(x_m_temp))
                uiconfirm(app.UIFigure,"X Max (m) is too large for Number of X FFT Points. Decrease X Max (m) or increase Number of X FFT Points",'Imaging Parameter Error!',...
                    "Options",{'OK'},'Icon','warning');
                obj.isFail = true;
                return;
            end
            if max(abs(obj.y_m)) > max(abs(y_m_temp))
                uiconfirm(app.UIFigure,"X Max (m) is too large for Number of Y FFT Points. Decrease Y Max (m) or increase Number of Y FFT Points",'Imaging Parameter Error!',...
                    "Options",{'OK'},'Icon','warning');
                obj.isFail = true;
                return;
            end
        end
        
        function obj = verifyReconstruction(obj,app)
            if app.sar.method ~= "Rectilinear"
                uiconfirm(app.UIFigure,"Must use 2-D XY SAR scan to use 2-D SAR 2-D FFT image reconstruction method!",'SAR Scenario Error!',...
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
            x_m_temp = double(make_x(obj,obj.xStep_m,obj.nFFTx));
            y_m_temp = double(make_x(obj,obj.yStep_m,obj.nFFTy));
            
            app.XMinmEditField_im_3.Value = min(x_m_temp);
            app.XMaxmEditField_im_3.Value = max(x_m_temp);
            
            app.YMinmEditField_im_3.Value = min(y_m_temp);
            app.YMaxmEditField_im_3.Value = max(y_m_temp);
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
                'Message','Reconstructing Image using Uniform FFT Method');
            
            % sarData is of size (sar.numY, sar.numX, fmcw.ADCSamples)
            d.Value = 1/10;
            
            % Zero-Pad Data: s(y,x,k)
            sarDataPadded = obj.sarData;
            sarDataPadded = padarray(sarDataPadded,[floor((obj.nFFTy-size(obj.sarData,1))/2) 0],0,'pre');
            sarDataPadded = padarray(sarDataPadded,[0 floor((obj.nFFTx-size(obj.sarData,2))/2)],0,'pre');
            clear sarData
            d.Value = 2/10;
            
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
            
            kZ = single(sqrt((4 * k.^2 - kX.^2 - kY.^2) .* (4 * k.^2 > kX.^2 + kY.^2)));
            d.Value = 3/10;
            
            % Compute Focusing Filter
            focusingFilter = exp(-1j * kZ * (obj.z0_m + obj.z_target_m));
            focusingFilter(4 * k.^2 < kX.^2 + kY.^2) = 0;
            d.Value = 4/10;
            
            % Compute FFT across Y & X Dimensions: S(kY,kX,k)
            sarDataFFT = fftshift(fftshift(fft(fft(sarDataPadded,obj.nFFTy,1),obj.nFFTx,2),1),2)/obj.nFFTx/obj.nFFTy;
            clear sarDataPadded sarData
            
            if obj.isGPU
                focusingFilter = gpuArray(focusingFilter);
            end
            d.Value = 5/10;
            
            sarImageFFT = sum(sarDataFFT .* focusingFilter,3);
            clear sarDataFFT focusingFilter kY kX k kZ
            
            d.Value = 6/10;
            
            % Recover Image by IFT: p(y,x,z)
            sarImage = single(abs(ifftn(sarImageFFT)));
            clear sarImageFFT focusingFilter
            d.Value = 7/10;
            
            % Reorient Image: p(x,y,z)
            sarImage = permute(sarImage,[2,1,3]);
            d.Value = 8/10;
            
            % Declare Spatial Vectors
            x_m_temp = make_x(obj,obj.xStep_m,obj.nFFTx);
            y_m_temp = make_x(obj,obj.yStep_m,obj.nFFTy);
            d.Value = 9/10;
            
            [X,Y] = ndgrid(obj.x_m(:),obj.y_m(:));
            obj.imXYZ = single(gather(interpn(x_m_temp(:),y_m_temp(:),sarImage,X,Y,"nearest",0)));
            d.Value = 10/10;
        end
        
        function displayImage(obj,im)
            displayImage2D_app(im,im.x_m,im.y_m,"x (m)","y (m)");
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