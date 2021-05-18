% Attempts to compute the Backprojection Algorithm by iterating over the
% wavenumber and target voxel domains. Usually takes an excessive amount of
% time to compute
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

function mediumBPA(obj,k)

if obj.im.isApp
    d = uiprogressdlg(obj.im.app.UIFigure,'Title','Performing BPA',...
        'Message',"Estimated Time Remaining: 0:0:0","Cancelable","on");
end

try
    mediumBPA1(obj,k,d);
catch
    try
        mediumBPA2(obj,k,d);
    catch
        slowBPA(obj,k,d);
    end
end
end

function mediumBPA1(obj,k,d)
if obj.isGPU
    reset(gpuDevice)
end

obj.imXYZ = single(zeros(1,size(obj.target_xyz_m,1)));
numSAR = size(obj.tx_xyz_m,1);
tocs = single(zeros(1,numSAR));
for indSAR = 1:numSAR
    tic
    
    if obj.im.isApp && d.CancelRequested
        showErrorMessage(obj.im,"Image not computed!","User Canceled BPA")
        obj.isFail = true;
        return;
    end
    
    if ~obj.ant.isEPC
        Rt = pdist2(obj.tx_xyz_m(indSAR,:),obj.target_xyz_m);
        Rr = pdist2(obj.rx_xyz_m(indSAR,:),obj.target_xyz_m);
        amplitudeFactor = Rt .* Rr;
        R_T_plus_R_R = Rt + Rr;
    else
        R = pdist2(obj.vx_xyz_m(indSAR,:),obj.target_xyz_m);
        R_T_plus_R_R = 2*R;
        amplitudeFactor = R.^2;
    end
    
    if obj.isGPU
        R_T_plus_R_R = gpuArray(R_T_plus_R_R);
    end
    
    bpaKernel = gather(exp(-1j*k.*R_T_plus_R_R));
    if obj.isAmplitudeFactor
        bpaKernel = bpaKernel .* amplitudeFactor(indSAR);
    end
    obj.imXYZ = obj.imXYZ + sum(obj.sarData(indSAR,:,:) .* bpaKernel,3);
    % Update the progress dialog
    if ~obj.im.isSilent
        tocs(indSAR) = toc;
        if obj.im.isApp
            d.Value = indSAR/numSAR;
            d.Message = "Iteration " + indSAR + "/" + numSAR + ". Estimated Time Remaining: " + getEstTime(obj,tocs,indSAR,numSAR);
        else
            disp("Iteration " + indSAR + "/" + numSAR + ". Estimated Time Remaining: " + getEstTime(obj,tocs,indSAR,numSAR));
        end
    end
end
end

function mediumBPA2(obj,k,d)
if obj.isGPU
    reset(gpuDevice)
end

obj.imXYZ = single(zeros(1,size(obj.target_xyz_m,1)));
numTargetVoxels = size(obj.target_xyz_m,1);
tocs = single(zeros(1,numTargetVoxels));
for indTarget = 1:numTargetVoxels
    tic
    
    if obj.im.isApp && d.CancelRequested
        showErrorMessage(obj.im,"Image not computed!","User Canceled BPA")
        obj.isFail = true;
        return;
    end
    
    if ~obj.ant.isEPC
        Rt = pdist2(obj.tx_xyz_m,obj.target_xyz_m(indTarget,:));
        Rr = pdist2(obj.rx_xyz_m,obj.target_xyz_m(indTarget,:));
        amplitudeFactor = Rt .* Rr;
        R_T_plus_R_R = Rt + Rr;
    else
        R = pdist2(obj.vx_xyz_m,obj.target_xyz_m(indTarget,:));
        R_T_plus_R_R = 2*R;
        amplitudeFactor = R.^2;
    end
    
    if obj.isGPU
        R_T_plus_R_R = gpuArray(R_T_plus_R_R);
    end
    
    bpaKernel = gather(exp(-1j*k.*R_T_plus_R_R));
    if obj.isAmplitudeFactor
        bpaKernel = bpaKernel .* amplitudeFactor;
    end
    obj.imXYZ(indTarget) = sum(obj.sarData .* bpaKernel,'all');
    % Update the progress dialog
    if ~obj.im.isSilent
        tocs(indSAR) = toc;
        if obj.im.isApp
            d.Value = indTarget/numTargetVoxels;
            d.Message = "Iteration " + indTarget + "/" + numTargetVoxels + ". Estimated Time Remaining: " + getEstTime(obj,tocs,indTarget,numTargetVoxels);
        else
            disp("Iteration " + indTarget + "/" + numTargetVoxels + ". Estimated Time Remaining: " + getEstTime(obj,tocs,indTarget,numTargetVoxels));
        end
    end
end
end

function slowBPA(obj,k,d)
obj.imXYZ = single(zeros(1,size(obj.target_xyz_m,1)));
tocs = single(zeros(1,2^14));
count = 0;
for indTarget = 1:size(obj.target_xyz_m,1)
    for indK = 1:length(k)
        tic
        count = count + 1;
        
        if obj.im.isApp && d.CancelRequested
            showErrorMessage(obj.im,"Image not computed!","User Canceled BPA")
            obj.isFail = true;
            return;
        end
        
        if ~obj.ant.isEPC
            Rt = pdist2(obj.tx_xyz_m,obj.target_xyz_m(indTarget,:));
            Rr = pdist2(obj.rx_xyz_m,obj.target_xyz_m(indTarget,:));
            amplitudeFactor = Rt .* Rr;
            R_T_plus_R_R = Rt + Rr;
        else
            R = pdist2(obj.vx_xyz_m,obj.target_xyz_m(indTarget,:));
            R_T_plus_R_R = 2*R;
            amplitudeFactor = R.^2;
        end
        
        if obj.isGPU
            R_T_plus_R_R = gpuArray(R_T_plus_R_R);
        end
        
        bpaKernel = gather(exp(-1j*k(indK)*R_T_plus_R_R));
        if obj.isAmplitudeFactor
            bpaKernel = bpaKernel .* amplitudeFactor;
        end
        obj.imXYZ(indTarget) = obj.imXYZ(indTarget) + sum(obj.sarData(:,:,indK) .* bpaKernel,1);
        % Update the progress dialog
        if ~obj.im.isSilent
            tocs(mod(count,length(tocs))+1) = toc;
            if obj.im.isApp
                d.Value = count/(length(k)*size(obj.target_xyz_m,1));
                d.Message = "Iteration " + count + "/" + length(k)*size(obj.target_xyz_m,1) + ". Estimated Time Remaining: " + getEstTime(obj,tocs,indK,length(k)*size(obj.target_xyz_m,1));
            else
                disp("Iteration " + count + "/" + length(k)*size(obj.target_xyz_m,1) + ". Estimated Time Remaining: " + getEstTime(obj,tocs,indK,length(k)*size(obj.target_xyz_m,1)));
            end
        end
    end
end
end

function outstr = getEstTime(obj,tocs,currentInd,totalInd)
% Estimates the time remaining given the recent time-per-iteration values.
% Returns a string containing the number of hours, number of minutes, and
% number of seconds remaining, each separated by a ":"

avgtoc = mean(tocs(tocs~=0))*(totalInd - currentInd);

if avgtoc > obj.estTime.old && obj.estTime.count < 100
    obj.estTime.count = obj.estTime.count + 1;
    avgtoc = obj.estTime.old;
else
    obj.estTime.count = 0;
    obj.estTime.old = avgtoc;
end

hrrem = floor(avgtoc/3600);
avgtoc = avgtoc - floor(avgtoc/3600)*3600;
minrem = floor(avgtoc/60);
avgtoc = avgtoc - floor(avgtoc/60)*60;
secrem = round(avgtoc);

if hrrem < 10
    hrrem = "0" + hrrem;
end
if minrem < 10
    minrem = "0" + minrem;
end
if secrem < 10
    secrem = "0" + secrem;
end
outstr = hrrem + ":" + minrem + ":" + secrem;
end