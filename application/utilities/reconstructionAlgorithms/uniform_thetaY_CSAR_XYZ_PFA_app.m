% uniform_thetaY_CSAR_XYZ_PFA_app see uniform_thetaY_CSAR_XYZ_PFA
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

classdef uniform_thetaY_CSAR_XYZ_PFA_app < handle
    properties
        sarData
        
        nFFTx
        nFFTy
        nFFTz
        thetaUpsampleFactor
        
        x_m
        y_m
        z_m
        
        imXYZ
        
        isGPU
        isAmplitudeFactor
        isFail
        
        theta_rad_vec
        k_vec
        R0_m
        xStep_m
        yStep_m
    end
    
    methods
        function obj = uniform_thetaY_CSAR_XYZ_PFA_app(app,im)
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
            obj.nFFTy = im.nFFTy;
            obj.nFFTz = im.nFFTz;
            obj.thetaUpsampleFactor = app.ThetaUpsampleFactorEditField.Value;
            
            obj.x_m = im.x_m;
            obj.y_m = im.y_m;
            obj.z_m = im.z_m;
            
            obj.sarData = app.target.sarData;
            
            obj.isGPU = im.isGPU;
            obj.isAmplitudeFactor = app.target.isAmplitudeFactor;
            
            obj.theta_rad_vec = app.sar.theta_rad;
            obj.k_vec = app.fmcw.k;
            obj.R0_m = app.ant.tx.z0_m;
            obj.xStep_m = app.sar.xStep_m;
            obj.yStep_m = app.sar.yStep_m;
        end
        
        function obj = verifyParameters(obj,app)
            obj.isFail = false;
            
            [x_m_temp,y_m_temp,z_m_temp] = getTempXYZ(obj);
            
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
            if app.sar.method ~= "Cylindrical"
                uiconfirm(app.UIFigure,"Must use 2-D θY Cylindrical CSAR scan to use 2-D CSAR 3-D PFA image reconstruction method!",'SAR Scenario Error!',...
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
            [x_m_temp,y_m_temp,z_m_temp] = getTempXYZ(obj);
            
            app.XMinmEditField_im_4.Value = min(x_m_temp);
            app.XMaxmEditField_im_4.Value = max(x_m_temp);
            
            app.YMinmEditField_im_6.Value = min(y_m_temp);
            app.YMaxmEditField_im_6.Value = max(y_m_temp);
            
            app.ZMinmEditField_im_4.Value = min(z_m_temp);
            app.ZMaxmEditField_im_4.Value = max(z_m_temp);
        end
        
        function [x_m_temp,y_m_temp,z_m_temp] = getTempXYZ(obj)
            theta_rad = single(reshape(obj.theta_rad_vec,1,[]));
            
            % Compute Wavenumbers
            k = single(reshape(obj.k_vec,1,1,[]));
            
            L_y = obj.nFFTy * obj.yStep_m;
            dkY = 2*pi/L_y;
            kY = make_kX(obj,dkY,obj.nFFTy)';
            
            kSy = 2*pi/obj.yStep_m;
            kY = linspace(-kSy/2,kSy/2,obj.nFFTy)';
            
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
        
        function [obj,imXYZ_out] = computeReconstruction(obj,app,im)
            obj = update(obj,app,im);
            
            if ~obj.isFail
                if app.sar.isMIMO && app.isMult2MonoCheckBox.Value
                    obj = mult2mono_app(obj,app);
                end
                
                %                 obj = reconstruct4segs(obj,app);
                obj = reconstruct(obj,app);
            end
            imXYZ_out = abs(obj.imXYZ);
        end
        
        %         function obj = reconstruct4segs(obj,app)
        %             theta_rad_orig = obj.theta_rad_vec;
        %             sarData_orig = obj.sarData;
        %
        %             % 1st segment
        %             obj.theta_rad_vec = theta_rad_orig(1:end/4);
        %             obj.sarData = sarData_orig(:,1:end/4,:);
        %             obj = reconstruct(obj,app);
        %             im1 = obj.imXYZ;
        %
        %             % 2nd segment
        %             obj.theta_rad_vec = theta_rad_orig(end/4+1:end/2);
        %             obj.sarData = sarData_orig(:,end/4+1:end/2,:);
        %             obj = reconstruct(obj,app);
        %             im2 = obj.imXYZ;
        %
        %             % 3rd segment
        %             obj.theta_rad_vec = theta_rad_orig(end/2+1:end*3/4);
        %             obj.sarData = sarData_orig(:,end/2+1:end*3/4,:);
        %             obj = reconstruct(obj,app);
        %             im3 = obj.imXYZ;
        %
        %             % 4th segment
        %             obj.theta_rad_vec = theta_rad_orig(end*3/4+1:end);
        %             obj.sarData = sarData_orig(:,end*3/4+1:end,:);
        %             obj = reconstruct(obj,app);
        %             im4 = obj.imXYZ;
        %
        %             obj.imXYZ = im1 + im2 + im3 + im4;
        %         end
        
        function obj = reconstruct(obj,app)
            % Start the Progress Bar
            d = uiprogressdlg(app.UIFigure,'Title','Please Wait',...
                'Message','Reconstructing Image using Uniform PFA ');
            
            % sarData is of size (sar.numY, sar.numTheta, fmcw.ADCSamples)
            d.Value = 1/10;
            
            % Zero-Pad Data: s(y,theta,k)
            sarDataPadded = obj.sarData;
            sarDataPadded = padarray(sarDataPadded,[floor((obj.nFFTy-size(obj.sarData,1))/2) 0],0,'pre');
            d.Value = 2/10;
            
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
            d.Value = 3/10;
            
            % Upsample data
            if obj.thetaUpsampleFactor > 1 || obj.nFFTz > length(k)
                theta_radUp = single(reshape(linspace(theta_rad(1),theta_rad(end),length(theta_rad)*obj.thetaUpsampleFactor),1,[]));
                kUp = single(reshape(linspace(min(k),max(k),obj.nFFTz),1,1,[]));
                kR = single(sqrt(4 * kUp.^2 - kY.^2) .* (4 * kUp.^2 > kY.^2));
                y = single(reshape(-(size(sarDataPadded,1)-1)/2:(size(sarDataPadded,1)-1)/2,[],1));
                [Yold,Told,Kold] = ndgrid(y,theta_rad,k);
                [Y,T,K] = ndgrid(y,theta_radUp,kUp);
                sarDataPadded = interpn(Yold,Told,Kold,sarDataPadded,Y,T,K);
                clear Y T K
                theta_rad = theta_radUp;
                k = kUp;
                clear theta_radUp kUp y Yold Told Kold
            end
            d.Value = 3.5/10;
            
            % Compute Azimuth Filter H(kY,Theta,k)
            azimuthFilterFFT = single(fft(exp(1j*kR*obj.R0_m.*cos(theta_rad)),[],2));
            clear kR
            d.Value = 4/10;
            
            % Compute Azimuth Filtered Data: p(kY,theta,k) = IFT_Theta[ S(kY,Theta,k) * H*(kY,Theta,k)]
            azimuthFiltered = ifft(fftshift(fft(fft(sarDataPadded,[],2),obj.nFFTy,1),1) .* conj(azimuthFilterFFT),[],2);
            clear azimuthFilterFFT sarDataPadded
            d.Value = 5/10;
            
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
            d.Value = 6/10;
            
            % Recover Image by IFT: p(y,x,z)
            sarImage = single(ifftshift(ifftshift(ifftn(sarImageFFT),2),3));
            clear sarImageFFT
            d.Value = 7/10;
            
            % Reorient Image: p(x,y,z)
            sarImage = permute(sarImage,[2,1,3]);
            d.Value = 8/10;
            
            % Declare Spatial Vectors
            x_m_temp = make_x(obj,2*pi/(dkXU*obj.nFFTx),obj.nFFTx);
            y_m_temp = make_x(obj,obj.yStep_m,obj.nFFTy);
            z_m_temp = make_x(obj,2*pi/(dkZU*obj.nFFTz),obj.nFFTz);
            d.Value = 9/10;
            
            [X,Y,Z] = ndgrid(obj.x_m(:),obj.y_m(:),obj.z_m(:));
            obj.imXYZ = single(gather(interpn(x_m_temp(:),y_m_temp(:),z_m_temp(:),sarImage,X,Y,Z,'nearest',0)));
            d.Value = 10/10;
        end
        
        function obj = displayImage(obj,im)
            displayImage3D_app(im);
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