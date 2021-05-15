% sarTarget object holds the properties and methods used in the FMCW MIMO-SAR scenario as specified by the user
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

classdef sarTarget < handle
    properties
        isGPU = true                % Boolean whether or not to use the GPU for beat signal computation
        isLong = false              % Boolean whether or not to use the long beat signal computation method
        numTargets = 0              % Number of target voxels
        xyz_m                       % Target x-y-z locations as a (numTargets)x3 array
        amp                         % Column vector of target reflectivities
        R                           % Radial distances from antennas to target (MIMO or EPC)
    end
    
    properties(SetObservable)
        isAmplitudeFactor = false   % Boolean whether or not to include the amplitude factor (path loss) in the beat signal computation
        isTable = false             % Boolean whether or not to include the table of targets
        isPNG = false               % Boolean whether or not to include the PNG file as a target
        isSTL = false               % Boolean whether or not to include the STL file as a target
        isRandomPoints = false      % Boolean whether or not to include the random points as a target
        isGPUVerified = false       % Boolean whether or not the GPU is has been verified
        
        % tableTarget - Array containing the x-y-z location of the targets
        % and their corresponding reflectivities in the form:
        %   |   x   |   y   |   z   |   r   |
        %   |   0   |   0   |   0   |   1   |
        % Where the example has a target located at (x,y,z) = (0,0,0) with
        % a reflectivity of 1
        tableTarget = [0 0 0 0]
        
        png = struct('isLoaded',false,'fileNameLoaded',"") % Structure containing the parameters of the PNG target
        stl = struct('isLoaded',false,'fileNameLoaded',"") % Structure containing the parameters of the STL target
        rp                          % Structure containing the parameters of the random point targets
        
        sarData                     % Computed beat signal
        
        fig = struct('f',[],'h',[]) % Structure containing the figure and handle used for showing the target
        fmcw                        % An fmcwChirpParameters object
        ant                         % A sarAntennaArray object
        sar                         % A sarScenario object
        
        isSilent = false            % Boolean whether or not to display updates to the user
    end
    
    methods
        function obj = sarTarget(fmcw,ant,sar)
            % Attaches the listener to the sarTarget object so
            % observable properties can be watched for changes and sets the
            % necessary parameters from the other SAR and FMCW type
            % objects. Then verifies whether the GPU can be used, if
            % necessary
            
            attachListener(obj)
            obj.fmcw = fmcw;
            obj.ant = ant;
            obj.sar = sar;
            
            verifyGPU(obj);
        end
        
        function getTarget(obj)
            % Gets the target parameters using the methods
            % (table/png/stl/randompoints) specified by the user, then
            % verifies the GPU can be used, if necessary
            
            obj.xyz_m = [];
            obj.amp = [];
            
            if obj.isTable
                getTargetTable(obj);
            end
            if obj.isPNG
                getPNGTarget(obj);
            end
            if obj.isSTL
                getSTLTarget(obj);
            end
            if obj.isRandomPoints
                getRandomTarget(obj);
            end
            
            obj.numTargets = size(obj.xyz_m,1);
            obj.xyz_m = single(obj.xyz_m);
            obj.amp = single(obj.amp(:).');
            
            verifyGPU(obj);
        end
        
        function getTargetTable(obj)
            % Gets the positions and amplitudes from the table of targets
            
            temp_table = single(obj.tableTarget);
            
            temp_xyz_m = temp_table(:,1:3);
            temp_amp = temp_table(:,4);
            obj.xyz_m = cat(1,obj.xyz_m,temp_xyz_m);
            obj.amp = cat(1,obj.amp,temp_amp);
        end
        
        function getPNGTarget(obj)
            % Gets the positions and amplitudes from the PNG target
            
            if ~obj.png.isLoaded || string(obj.png.fileNameLoaded) ~= string(obj.png.fileName)
                loadPNG(obj);
            end
            
            [numY,numX] = size(obj.png.tMat);
            xAxisT = obj.png.xStep_m * (-(numX-1)/2 : (numX-1)/2) + obj.png.xOffset_m;
            yAxisT = obj.png.yStep_m * (-(numY-1)/2 : (numY-1)/2) + obj.png.yOffset_m;
            zAxisT = obj.png.zOffset_m;
            
            [zT,xT,yT] = meshgrid(zAxisT,xAxisT,yAxisT);
            temp_xyz_m = reshape(permute([xT,yT,zT],[1 3 2]),[],3);
            
            indT = rot90(obj.png.tMat,-1)==true;
            temp_xyz_m = single(temp_xyz_m(indT,:));
            
            temp_pc = pointCloud(temp_xyz_m);
            temp_pc = pcdownsample(temp_pc,'random',1/obj.png.downsampleFactor);
            temp_xyz_m = temp_pc.Location;
            temp_amp = ones(size(temp_xyz_m,1),1)*obj.png.reflectivity;
            
            obj.xyz_m = cat(1,obj.xyz_m,temp_xyz_m);
            obj.amp = cat(1,obj.amp,temp_amp);
        end
        
        function loadPNG(obj)
            % Loads the PNG file and updates png.isLoaded and
            % png.fileNameLoaded
            
            try
                obj.png.tMat = imread(obj.png.fileName);
            catch
                obj.png.isLoaded = false;
                obj.isPNG = false;
                error("Invalid png.fileName value!");
            end
            obj.png.tMat = obj.png.tMat(:,:,1);
            obj.png.tMat(obj.png.tMat<64) = 0;
            obj.png.tMat(obj.png.tMat>0) = 1;
            obj.png.tMat = ~obj.png.tMat;
            obj.png.tMat = fliplr(obj.png.tMat);
            
            obj.png.isLoaded = true;
            obj.png.fileNameLoaded = obj.png.fileName;
        end
        
        function getSTLTarget(obj)
            % Gets the positions and amplitudes from the STL target
            
            if ~obj.stl.isLoaded || string(obj.stl.fileNameLoaded) ~= string(obj.stl.fileName)
                loadSTL(obj);
            end
            
            temp_xyz_m = single(zeros(size(obj.stl.v)));
            temp_xyz_m(:,1) = obj.stl.v(:,1);
            temp_xyz_m(:,2) = obj.stl.v(:,3);
            temp_xyz_m(:,3) = obj.stl.v(:,2);
            
            temp_xyz_m = temp_xyz_m + [obj.stl.xOffset_m,0,0];
            temp_xyz_m = temp_xyz_m + [0,obj.stl.yOffset_m,0];
            temp_xyz_m = temp_xyz_m + [0,0,obj.stl.zOffset_m];
            temp_xyz_m = temp_xyz_m(temp_xyz_m(:,3)<obj.stl.zCrop_m,:);
            
            temp_pc = pointCloud(temp_xyz_m);
            temp_pc = pcdownsample(temp_pc,'random',1/obj.stl.downsampleFactor);
            temp_xyz_m = temp_pc.Location;
            temp_amp = ones(size(temp_xyz_m,1),1)*obj.stl.reflectivity;
            
            obj.xyz_m = cat(1,obj.xyz_m,temp_xyz_m);
            obj.amp = cat(1,obj.amp,temp_amp);
        end
        
        function loadSTL(obj)
            % Loads the STL file and updates stl.isLoaded and
            % stl.fileNameLoaded
            
            try
                [~,obj.stl.v] = stlread2011(obj.stl.fileName);
            catch
                obj.stl.isLoaded = false;
                obj.isSTL = false;
                error("Invalid stl.fileName value!")
            end
            obj.stl.v = obj.stl.v*1e-3;
            obj.stl.isLoaded = true;
            obj.stl.fileNameLoaded = obj.stl.fileName;
        end
        
        function getRandomTarget(obj)
            % Generates the random targets given the parameters
            
            temp_x = obj.rp.xMin_m + (obj.rp.xMax_m - obj.rp.xMin_m)*rand(obj.rp.numTargets,1);
            temp_y = obj.rp.yMin_m + (obj.rp.yMax_m - obj.rp.yMin_m)*rand(obj.rp.numTargets,1);
            temp_z = obj.rp.zMin_m + (obj.rp.zMax_m - obj.rp.zMin_m)*rand(obj.rp.numTargets,1);
            
            temp_xyz_m = [temp_x,temp_y,temp_z];
            temp_amp = obj.rp.ampMin + (obj.rp.ampMax - obj.rp.ampMin)*rand(obj.rp.numTargets,1);
            
            obj.xyz_m = cat(1,obj.xyz_m,temp_xyz_m);
            obj.amp = cat(1,obj.amp,temp_amp);
        end
        
        function computeTarget(obj)
            % Computes the beat signal. First attempts the fast method then
            % the slow method, if the fast method fails
            
            if ~obj.ant.isEPC
                % Get distances
                try
                    obj.R = struct;
                    obj.R.tx = pdist2(obj.sar.rx.xyz_m,obj.xyz_m);
                    obj.R.rx = pdist2(obj.sar.tx.xyz_m,obj.xyz_m);
                    R_T_plus_R_R = obj.R.tx + obj.R.rx;
                    
                    % Amplitude Factor
                    if obj.isAmplitudeFactor
                        amplitudeFactor = obj.amp./(obj.R.tx .* obj.R.rx);
                    else
                        amplitudeFactor = 1;
                    end
                    
                    obj.R.tx = 0;
                    obj.R.rx = 0;
                    
                    obj.isLong = false;
                catch
                    obj.isLong = true;
                end
            else
                % Get distances
                try
                    obj.R = pdist2(obj.sar.vx.xyz_m,obj.xyz_m);
                    
                    % Get echo signal
                    R_T_plus_R_R = 2*obj.R;
                    if obj.isAmplitudeFactor
                        amplitudeFactor = obj.amp./(obj.R).^2;
                    else
                        amplitudeFactor = 1;
                    end
                    obj.isLong = false;
                catch
                    obj.isLong = true;
                end
            end
            
            if obj.isGPU && ~obj.isLong
                amplitudeFactor = gpuArray(amplitudeFactor);
                R_T_plus_R_R = gpuArray(R_T_plus_R_R);
            end
            if ~obj.isSilent
                % Create the progress dialog
                d = waitbar(0,'1','Name',' Generating Beat Signal...',...
                    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
            end
            
            obj.sarData = single(zeros(size(obj.sar.tx.xyz_m,1),obj.fmcw.ADCSamples));
            
            try
                if obj.isLong
                    error("oops");
                end
                
                % Fast method
                if ~obj.isSilent
                    tocs = zeros(1,obj.fmcw.ADCSamples);
                end
                for indK = 1:obj.fmcw.ADCSamples
                    if ~obj.isSilent
                        if getappdata(d,'canceling')
                            warning("Beat Signal not Computed!")
                            obj.sarData = single(zeros([obj.sar.sarSize,obj.fmcw.ADCSamples]));
                            delete(d);
                            return;
                        end
                        tic
                    end
                    temp = exp(1j*obj.fmcw.k(indK)*R_T_plus_R_R);
                    if obj.isAmplitudeFactor
                        temp = amplitudeFactor .* temp;
                    end
                    
                    obj.sarData(:,indK) = single(gather(sum(temp,2)));
                    if ~obj.isSilent
                        % Update the progress dialog
                        tocs(indK) = toc;
                        waitbar(indK/obj.fmcw.ADCSamples,d,"Generating Beat Signal. Estimated Time Remaining: " + getEstTime(obj,tocs,indK,obj.fmcw.ADCSamples));
                    end
                end
            catch
                try
                    % Fast method 2
                    if ~obj.isSilent
                        tocs = zeros(1,obj.numTargets);
                    end
                    for indTarget = 1:obj.numTargets
                        if ~obj.isSilent
                            if getappdata(d,'canceling')
                                warning("Beat Signal not Computed!")
                                obj.sarData = single(zeros([obj.sar.sarSize,obj.fmcw.ADCSamples]));
                                delete(d);
                                return;
                            end
                            tic
                        end
                        temp = exp(1j*obj.fmcw.k.*R_T_plus_R_R(:,indTarget));
                        if obj.isAmplitudeFactor
                            temp = amplitudeFactor .* temp;
                        end
                        
                        obj.sarData = obj.sarData + single(gather(temp));
                        % Update the progress dialog
                        if ~obj.isSilent
                            tocs(indTarget) = toc;
                            waitbar(indTarget/obj.numTargets,d,"Generating Beat Signal. Estimated Time Remaining: " + getEstTime(obj,tocs,indTarget,obj.numTargets));
                        end
                    end
                catch
                    try
                        R_T_plus_R_R = gather(R_T_plus_R_R);
                        computeTargetLarge(obj,R_T_plus_R_R,d);
                    catch
                        computeTargetSlow(obj,d);
                    end
                end
            end
            
            if ~obj.isSilent
                delete(d);
            end
            
            % Reshape echo signal
            obj.sarData = reshape(obj.sarData,[obj.sar.sarSize,obj.fmcw.ADCSamples]);
        end
        
        function computeTargetLarge(obj,R_T_plus_R_R,d)
            [R_unique,~,IC] = unique(R_T_plus_R_R,'stable');
            if obj.isGPU
                R_unique = gpuArray(R_unique);
            end
            % Remember R_T_plus_R_R = reshape(R_unique(IC),size(R_T_plus_R_R));
            
            sizeR = size(R_T_plus_R_R);
            R_T_plus_R_R = 0;
            
            try
                [temp,~,IC2] = unique(mod(R_unique.*obj.fmcw.k,2*pi),'stable');
                % Remember R_unique = reshape(temp(IC2),size(R_unique));
            catch
                R_unique = gather(R_unique);
                [temp,~,IC2] = unique(mod(R_unique.*obj.fmcw.k,2*pi),'stable');
                % Remember R_unique = reshape(temp(IC2),size(R_unique));
            end
            if obj.isGPU
                temp = gpuArray(temp);
            end
            temp2 = gather(exp(1j*temp));
            
            temp3 = reshape(temp2(IC2),[length(R_unique),obj.fmcw.ADCSamples]);
            
            if ~obj.isSilent
                tocs = zeros(1,obj.fmcw.ADCSamples);
            end
            for indK = 1:obj.fmcw.ADCSamples
                if ~obj.isSilent
                    if getappdata(d,'canceling')
                        warning("Beat Signal not Computed!")
                        obj.sarData = 0;
                        return;
                    end
                    tic
                end
                temp = reshape(temp3(IC,indK),sizeR);
                if obj.isAmplitudeFactor
                    temp = amplitudeFactor .* temp;
                end
                
                obj.sarData(:,indK) = single(gather(sum(temp,2)));
                % Update the progress dialog
                if ~obj.isSilent
                    tocs(indK) = toc;
                    waitbar(indK/obj.fmcw.ADCSamples,d,"Generating Beat Signal. Estimated Time Remaining: " + getEstTime(obj,tocs,indK,obj.fmcw.ADCSamples));
                end
            end
        end
        
        function d = computeTargetSlow(obj,d)
            % Always works method
            if ~obj.isSilent
                tocs = single(zeros(1,obj.fmcw.ADCSamples*obj.numTargets));
                count = 0;
            end
            for indSAR = 1:size(obj.sar.rx.xyz_m,1)
                if ~obj.isSilent
                    tic
                end
                if ~obj.ant.isEPC
                    obj.R = struct;
                    obj.R.tx = pdist2(obj.sar.tx.xyz_m(indSAR,:),obj.xyz_m);
                    obj.R.rx = pdist2(obj.sar.rx.xyz_m(indSAR,:),obj.xyz_m);
                    R_T_plus_R_R = obj.R.tx + obj.R.rx;
                    
                    % Amplitude Factor
                    if obj.isAmplitudeFactor
                        amplitudeFactor = obj.amp./(obj.R.tx .* obj.R.rx);
                    else
                        amplitudeFactor = 1;
                    end
                else
                    obj.R = pdist2(obj.sar.vx.xyz_m(indSAR,:),obj.xyz_m);
                    
                    % Get echo signal
                    R_T_plus_R_R = 2*obj.R;
                    if obj.isAmplitudeFactor
                        amplitudeFactor = obj.amp./(obj.R).^2;
                    else
                        amplitudeFactor = single(1);
                    end
                end
                
                if obj.isGPU
                    amplitudeFactor = gpuArray(amplitudeFactor);
                    R_T_plus_R_R = gpuArray(R_T_plus_R_R);
                end
                
                for indK = 1:obj.fmcw.ADCSamples
                    if ~obj.isSilent
                        if getappdata(d,'canceling')
                            warning("Beat Signal not Computed!")
                            obj.sarData = single(zeros([obj.sar.sarSize,obj.fmcw.ADCSamples]));
                            delete(d);
                            return;
                        end
                        count = count + 1;
                    end
                    temp = exp(1j*obj.fmcw.k(indK)*R_T_plus_R_R);
                    if obj.isAmplitudeFactor
                        temp = amplitudeFactor .* temp;
                    end
                    
                    obj.sarData(indSAR,indK) = single(gather(sum(temp,2)));
                    if ~obj.isSilent
                        % Update the progress dialog
                        tocs(count) = toc;
                        waitbar(count/(obj.fmcw.ADCSamples*size(obj.sar.rx.xyz_m,1)),d,"Generating Beat Signal Using Slow Method. Estimated Time Remaining: " + getEstTime(obj,tocs,count,obj.fmcw.ADCSamples*size(obj.sar.rx.xyz_m,1)));
                    end
                end
            end
        end
        
        % Plot/figure functions
        function initializeFigures(obj)
            % Initializes the figures
            
            closeFigures(obj);
            set(0,'DefaultFigureWindowStyle','docked')
            
            obj.fig.f = figure;
            obj.fig.h = handle(axes);
        end
        
        function closeFigures(obj)
            % Attempt to close the figures
            
            try
                close(obj.fig.f)
            catch
            end
        end
        
        function displayTarget(obj)
            % Display the target scenario
            
            if isempty(obj.xyz_m)
                return;
            end
            
            if isempty(obj.fig.f) || ~isvalid(obj.fig.h)
                initializeFigures(obj);
            end
            
            if obj.ant.isEPC
                displayTargetEPC(obj);
            else
                displayTargetMIMO(obj);
            end
        end
        
        function displayTargetMIMO(obj)
            % Plot the target with the MIMO-SAR scenario
            
            h = obj.fig.h;
            hold(h,'off')
            temp = obj.sar.tx.xyz_m;
            scatter3(h,temp(:,1),temp(:,3),temp(:,2),'.r')
            hold(h,'on')
            temp = obj.sar.rx.xyz_m;
            scatter3(h,temp(:,1),temp(:,3),temp(:,2),'.b')
            temp = obj.xyz_m;
            scatter3(h,temp(:,1),temp(:,3),temp(:,2),'.k')
            
            xlabel(h,"x (m)")
            temp1 = obj.sar.tx.xyz_m(:,1);
            temp2 = obj.sar.rx.xyz_m(:,1);
            temp3 = obj.xyz_m(:,1);
            xlim(h,[min([min(temp1),min(temp2),min(temp3)])-0.01,max([max(temp1),max(temp2),max(temp3)])+0.01])
            ylabel(h,"z (m)")
            temp1 = obj.sar.tx.xyz_m(:,3);
            temp2 = obj.sar.rx.xyz_m(:,3);
            temp3 = obj.xyz_m(:,3);
            ylim(h,[min([min(temp1),min(temp2),min(temp3)])-0.01,max([max(temp1),max(temp2),max(temp3)])+0.01])
            zlabel(h,"y (m)")
            temp1 = obj.sar.tx.xyz_m(:,2);
            temp2 = obj.sar.rx.xyz_m(:,2);
            temp3 = obj.xyz_m(:,2);
            zlim(h,[min([min(temp1),min(temp2),min(temp3)])-0.01,max([max(temp1),max(temp2),max(temp3)])+0.01])
            title(h,"MIMO Aperture Image Scenario")
            legend(h,"Tx","Rx","Target")
            
            view(h,3)
            daspect(h,[1 1 1])
        end
        
        function displayTargetEPC(obj)
            % Plot the target with the virtual array SAR scenario
            
            h = obj.fig.h;
            hold(h,'off')
            temp = obj.sar.vx.xyz_m;
            scatter3(h,temp(:,1),temp(:,3),temp(:,2),'.b')
            hold(h,'on')
            temp = obj.xyz_m;
            scatter3(h,temp(:,1),temp(:,3),temp(:,2),'.k')
            
            xlabel(h,"x (m)")
            temp1 = obj.sar.tx.xyz_m(:,1);
            temp2 = obj.sar.rx.xyz_m(:,1);
            temp3 = obj.xyz_m(:,1);
            xlim(h,[min([min(temp1),min(temp2),min(temp3)])-0.01,max([max(temp1),max(temp2),max(temp3)])+0.01])
            ylabel(h,"z (m)")
            temp1 = obj.sar.tx.xyz_m(:,3);
            temp2 = obj.sar.rx.xyz_m(:,3);
            temp3 = obj.xyz_m(:,3);
            ylim(h,[min([min(temp1),min(temp2),min(temp3)])-0.01,max([max(temp1),max(temp2),max(temp3)])+0.01])
            zlabel(h,"y (m)")
            temp1 = obj.sar.tx.xyz_m(:,2);
            temp2 = obj.sar.rx.xyz_m(:,2);
            temp3 = obj.xyz_m(:,2);
            zlim(h,[min([min(temp1),min(temp2),min(temp3)])-0.01,max([max(temp1),max(temp2),max(temp3)])+0.01])
            title(h,"Virtual Aperture Image Scenario")
            legend(h,"Vx","Target")
            
            view(h,3)
            daspect(h,[1 1 1])
        end
        
        function saveTarget(obj,saveName)
            % Save the sarTarget object to a file
            
            if exist(saveName + ".mat",'file')
                str = input('Are you sure you want to overwrite? Y/N: ','s');
                if str ~= 'Y'
                    warning("Target not saved!");
                    return;
                end
            end
            
            savedtarget = obj;
            savedtarget.fig.f = [];
            savedtarget.ant = [];
            savedtarget.sar = [];
            save(saveName,"savedtarget");
            disp("Target saved to: " + saveName);
        end
        
        function loadTarget(obj,loadName)
            % Loads a sartarget object from a file
            
            if ~exist(loadName + ".mat",'file')
                warning("No file called " + loadName + ".mat to load. Antenna array not loaded!");
                return;
            end
            
            load(loadName,"savedtarget");
            
            fieldlist = string(fieldnames(savedtarget));
            for indField = 1:length(fieldlist)
                if fieldlist(indField) ~= "fmcw" && fieldlist(indField) ~= "ant" && fieldlist(indField) ~= "sar" && fieldlist(indField) ~= "fig"
                    obj.(fieldlist(indField)) = savedtarget.(fieldlist(indField));
                end
            end
            
            
        end
        
        function verifyGPU(obj)
            % Verifies if the GPU can be used
            
            if obj.isGPU && ~obj.isGPUVerified
                try
                    reset(gpuDevice);
                catch
                    obj.isGPU = false;
                    warning("Unable to locate Nvidia GPU");
                    return;
                end
                obj.isGPUVerified = true;
            end
        end
        
        function outstr = getEstTime(obj,tocs,currentInd,totalInd)
            % Estimates the time until completion
            
            avgtoc = mean(tocs(1:currentInd))*(totalInd - currentInd);
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
        
        function attachListener(obj)
            % Attaches a listener to the object handle
            
            addlistener(obj,{'isAmplitudeFactor','isTable','isPNG','isSTL','isRandomPoints','tableTarget','png','stl','rp','fmcw','ant','sar'},'PostSet',@sarTarget.propChange);
        end
    end
    
    methods(Static)
        function propChange(metaProp,eventData)
            % Recomputes the target if the watched parameters are changed
            
            getTarget(eventData.AffectedObject);
        end
    end
end