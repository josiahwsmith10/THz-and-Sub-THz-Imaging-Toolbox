% THzTarget object holds the properties and methods used in the target as
% specified by the user
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

classdef THzTarget < handle
    properties
        isGPU = true                % Boolean whether or not to use the GPU for beat signal computation
        isLong = false              % Boolean whether or not to use the long beat signal computation method
        numTargets = 0              % Number of target voxels
        xyz_m = []                  % Target x-y-z locations as a (numTargets)x3 array
        amp = []                    % Column vector of target reflectivities
        R = []                      % Radial distances from antennas to target (MIMO or EPC)
        
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
        rp = struct                 % Structure containing the parameters of the random point targets
        
        sarData = []                % Computed beat signal
        
        fig = struct('f',[],'h',[]) % Structure containing the figure and handle used for showing the target
        wav                         % A THzWaveformParameters handle
        ant                         % A THzAntennaArray handle
        scanner                     % A THzScanner handle
        
        app = struct("UIFigure",[]) % App handle
        isApp                       % Boolean whether or not to use the app functionality
        
        isSilent = false            % Boolean whether or not to display updates to the user
    end
    
    methods
        function obj = THzTarget(wav,ant,scanner,app)
            % Verifies the GPU can be used, sets the wav, ant, scanner and
            % app properties
            
            obj.wav = wav;
            obj.ant = ant;
            obj.scanner = scanner;
            if nargin == 3
                obj.isApp = false;
            elseif nargin == 4
                obj.isApp = true;
                obj.app = app;
            else
                showErrorMessage("Must input wav, ant, and scanner at least!");
            end
            
            obj = VerifyGPU(obj);
        end
        
        function obj = Update(obj)
            % Update the target with the values in the app
            
            if obj.isApp
                obj.app.BeatSignalComputedLamp.Color = "red";
                obj.app.ImageReconstructionCompleteLamp.Color = "red";
            end
            obj = Get(obj);
            
            if isempty(obj.xyz_m)
                return;
            end
            
            Display(obj);
        end
        
        function obj = Get(obj)
            % Gets the target using the values in the app
            
            if obj.isApp
                obj.isAmplitudeFactor = obj.app.PathLossCheckBox.Value;
            end
            
            obj.xyz_m = [];
            obj.amp = [];
            
            if obj.isTable || (obj.isApp && obj.app.UseTableofTargetsCheckBox.Value)
                obj = getTargetTable(obj);
            end
            if obj.isPNG || (obj.isApp && obj.app.UsePNGFileCheckBox.Value)
                obj = getPNGTarget(obj);
            end
            if obj.isSTL || (obj.isApp && obj.app.UseSTLFileCheckBox.Value)
                obj = getSTLTarget(obj);
            end
            if obj.isRandomPoints
                obj = getRandomTarget(obj);
            end
            
            obj.numTargets = size(obj.xyz_m,1);
            obj.xyz_m = single(obj.xyz_m);
            obj.amp = single(obj.amp(:).');
            
            if obj.isApp
                obj.app.TotalNumTargetsEditField.Value = obj.numTargets;
            end
        end
        
        function obj = getTargetTable(obj)
            % Gets the positions and amplitudes from the table of targets
            
            if obj.isApp && isempty(obj.app.TargetTable.Data)
                return;
            end
            if obj.isApp
                temp_table = single(table2array(obj.app.TargetTable.Data));
            else
                temp_table = single(obj.tableTarget);
            end
            
            temp_xyz_m = temp_table(:,1:3);
            temp_amp = temp_table(:,4);
            obj.xyz_m = cat(1,obj.xyz_m,temp_xyz_m);
            obj.amp = cat(1,obj.amp,temp_amp);
        end
        
        function obj = getPNGTarget(obj)
            % Gets the positions and amplitudes from the PNG target
            
            if obj.isApp
                obj = getPNGParameters(obj);
            end
            
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
            
            if obj.isApp
                obj.app.NumTargetsEditField.Value = size(temp_xyz_m,1);
            end
        end
        
        function obj = getPNGParameters(obj)
            % Gets the parameters for the PNG target from the app
            
            if ~obj.isApp
                warning("Must be using app to call getPNGParameters")
            end
            
            obj.png.fileName = obj.app.PNGFileNameEditField.Value;
            obj.png.xStep_m = obj.app.XPixelSizemEditField.Value;
            obj.png.yStep_m = obj.app.YPixelSizemEditField.Value;
            obj.png.xOffset_m = obj.app.XOffsetmEditField.Value;
            obj.png.yOffset_m = obj.app.YOffsetmEditField.Value;
            obj.png.zOffset_m = obj.app.ZOffsetmEditField.Value;
            obj.png.reflectivity = obj.app.ReflectivityEditField.Value;
            obj.png.downsampleFactor = obj.app.DownsampleFactorEditField.Value;
        end
        
        function loadPNG(obj)
            % Loads the PNG file and updates png.isLoaded and
            % png.fileNameLoaded
            
            try
                obj.png.tMat = imread(obj.png.fileName);
            catch
                obj.png.isLoaded = false;
                obj.isPNG = false;
                showErrorMessage(obj,"Invalid png.fileName value!","File Load Error");
            end
            obj.png.tMat = obj.png.tMat(:,:,1);
            obj.png.tMat(obj.png.tMat<64) = 0;
            obj.png.tMat(obj.png.tMat>0) = 1;
            obj.png.tMat = ~obj.png.tMat;
            obj.png.tMat = fliplr(obj.png.tMat);
            
            obj.png.isLoaded = true;
            obj.png.fileNameLoaded = obj.png.fileName;
        end
        
        function obj = getSTLTarget(obj)
            % Gets the positions and amplitudes from the STL target
            
            if obj.isApp
                obj = getSTLParameters(obj);
            end
            
            if ~obj.stl.isLoaded || string(obj.stl.fileNameLoaded) ~= string(obj.stl.fileName)
                loadSTL(obj);
            end
            
            temp_xyz_m = zeros(size(obj.stl.v),'single');
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
            
            if obj.isApp
                obj.app.NumTargetsEditField_2.Value = size(temp_xyz_m,1);
            end
        end
        
        function obj = getSTLParameters(obj)
            % Gets the parameters for the STL target from the app
            
            if ~obj.isApp
                return;
            end
            
            obj.stl.fileName = obj.app.STLFileNameEditField.Value;
            obj.stl.zCrop_m = obj.app.ZCropmEditField.Value;
            obj.stl.xOffset_m = obj.app.XOffsetmEditField_2.Value;
            obj.stl.yOffset_m = obj.app.YOffsetmEditField_2.Value;
            obj.stl.zOffset_m = obj.app.ZOffsetmEditField_2.Value;
            obj.stl.reflectivity = obj.app.ReflectivityEditField_2.Value;
            obj.stl.downsampleFactor = obj.app.DownsampleFactorEditField_2.Value;
        end
        
        function obj = loadSTL(obj)
            % Loads the STL file and updates stl.isLoaded
            
            if obj.isApp
                obj.app.LoadSTLLamp.Color = "yellow";
            end
            
            obj.stl.fileName = obj.app.STLFileNameEditField.Value;
            
            drawnow
            try
                [~,obj.stl.v] = stlread2011("./saved/pngstl/" + obj.stl.fileName);
                obj.stl.fileNameLoaded = obj.stl.fileName;
            catch
                if obj.isApp
                    selection = uiconfirm(obj.app.UIFigure,'Would you liked to locate the STL file?','STL File Not Found',...
                        'Icon','warning');
                    if string(selection) == "Cancel"
                        warning("STL file not loaded!");
                        obj.stl.isLoaded = false;
                        return;
                    end
                end
                
                [filename,pathname] = uigetfile("./saved/pngstl/*.stl","Select Desired STL File");
                
                if filename == 0
                    warning("STL file not loaded!");
                    obj.stl.isLoaded = false;
                    return;
                end
                if obj.isApp
                    obj.app.STLFileNameEditField.Value = filename;
                end
                [~,obj.stl.v] = stlread2011(string(pathname) + string(filename));
                obj.stl.fileNameLoaded = filename;
            end
            obj.stl.v = obj.stl.v*1e-3;
            obj.stl.isLoaded = true;
            if obj.isApp
                obj.app.LoadSTLLamp.Color = "green";
            end
        end
        
        function obj = getRandomTarget(obj)
            % Generates the random targets given the parameters
            
            obj = getRandomParameters(obj);
            
            temp_x = obj.rp.xMin_m + (obj.rp.xMax_m - obj.rp.xMin_m)*rand(obj.rp.numTargets,1);
            temp_y = obj.rp.yMin_m + (obj.rp.yMax_m - obj.rp.yMin_m)*rand(obj.rp.numTargets,1);
            temp_z = obj.rp.zMin_m + (obj.rp.zMax_m - obj.rp.zMin_m)*rand(obj.rp.numTargets,1);
            
            temp_xyz_m = [temp_x,temp_y,temp_z];
            temp_amp = obj.rp.ampMin + (obj.rp.ampMax - obj.rp.ampMin)*rand(obj.rp.numTargets,1);
            
            obj.xyz_m = cat(1,obj.xyz_m,temp_xyz_m);
            obj.amp = cat(1,obj.amp,temp_amp);
            
            if obj.isApp
                obj.app.TargetTable.Data = array2table([temp_xyz_m,temp_amp]);
            end
        end
        
        function obj = getRandomParameters(obj)
            % Gets the parameters for the ramdom targets from the app
            
            if ~obj.isApp
                return;
            end
            
            obj.rp.numTargets = obj.app.NumTargetsEditField_3.Value;
            obj.rp.xMin_m = obj.app.XMinmEditField.Value;
            obj.rp.xMax_m = obj.app.XMaxmEditField.Value;
            obj.rp.yMin_m = obj.app.YMinmEditField.Value;
            obj.rp.yMax_m = obj.app.YMaxmEditField.Value;
            obj.rp.zMin_m = obj.app.ZMinmEditField.Value;
            obj.rp.zMax_m = obj.app.ZMaxmEditField.Value;
            obj.rp.ampMin = obj.app.ReflectivityMinEditField.Value;
            obj.rp.ampMax = obj.app.ReflectivityMaxEditField.Value;
        end
        
        function obj = Compute(obj)
            % Computes the beat signal. First attempts the fast method then
            % the large target method, if the fast method fails. If the
            % large target method fails, it attempts the slow method
            
            obj = Get(obj);
            
            obj.ant.ConfigurePattern();
            
            try
                [obj,W,w] = computeWeightsFast(obj);
            catch
                try
                    [obj,W,w] = computeWeightsSlow(obj);
                catch
                    W = 1;
                    w.isFail = true;
                end
            end
            
            if obj.isGPU
                reset(gpuDevice)
            end
            
            [obj,R_T_plus_R_R,amplitudeFactor] = computeDistances(obj);
            
            if obj.isGPU && ~obj.isLong
                amplitudeFactor = gpuArray(amplitudeFactor);
                R_T_plus_R_R = gpuArray(R_T_plus_R_R);
                
                if ismatrix(W)
                    W_temp = gpuArray(W);
                else
                    W_temp = 1;
                end
                
                if ndims(W) ~= 3
                    W = 1;
                end
            end
            
            if ~obj.isSilent
                % Create the progress dialog
                if obj.isApp
                    d = uiprogressdlg(obj.app.UIFigure,'Title','Generating Echo Signal',...
                        'Message',"Estimated Time Remaining: 0:0:0","Cancelable","on");
                else
                    d = waitbar(0,'1','Name',' Generating Beat Signal...',...
                        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
                end
            end
            
            try
                obj = computeTargetFast(obj,d,R_T_plus_R_R,amplitudeFactor,W_temp,W);
            catch
                try
                    obj = computeTargetFast2(obj,d,R_T_plus_R_R,amplitudeFactor,W_temp,W,w);
                catch
                    if obj.isSilent
                        obj.isSilent = false;
                        % Create the progress dialog
                        if obj.isApp
                            d = uiprogressdlg(obj.app.UIFigure,'Title','Generating Echo Signal',...
                                'Message',"Estimated Time Remaining: 0:0:0","Cancelable","on");
                        end
                    end
                    try
                        R_T_plus_R_R = gather(R_T_plus_R_R);
                        obj = computeTargetLarge(obj,R_T_plus_R_R,d);
                    catch
                        R_T_plus_R_R = [];
                        obj = computeTargetSlow(obj,d,W_temp,W);
                    end
                end
            end
            
            if ~obj.isSilent
                delete(d);
            end
            
            % Reshape echo signal
            obj.sarData = reshape(obj.sarData,[obj.scanner.sarSize,obj.wav.Nk]);
        end
        
        function obj = computeTargetFast(obj,d,R_T_plus_R_R,amplitudeFactor,W_temp,W)
            % Fast method - loops over k
            
            obj.sarData = zeros(size(obj.scanner.tx.xyz_m,1),obj.wav.Nk,'single');
            
            if ~obj.isSilent
                tocs = zeros(1,obj.wav.Nk,'single');
            end
            
            for indK = 1:obj.wav.Nk
                if ~obj.isSilent
                    if obj.isApp
                        if d.CancelRequested
                            warning("Beat Signal not Computed!")
                            obj.sarData = zeros([obj.scanner.sarSize,obj.wav.Nk],'single');
                            return;
                        end
                    else
                        if getappdata(d,'canceling')
                            warning("Beat Signal not Computed!")
                            obj.sarData = zeros([obj.scanner.sarSize,obj.wav.Nk]);
                            delete(d);
                            return;
                        end
                    end
                    tic
                end
                
                if ndims(W) == 3
                    W_temp = gpuArray(W(:,:,indK));
                end
                
                temp = exp(1j*obj.wav.k(indK)*R_T_plus_R_R);
                temp = amplitudeFactor .* temp;
                
                obj.sarData(:,indK) = single(gather(sum(W_temp .* temp,2)));
                if ~obj.isSilent
                    % Update the progress dialog
                    tocs(indK) = toc;
                    if obj.isApp
                        d.Value = indK/obj.wav.Nk;
                        d.Message = "Estimated Time Remaining: " + getEstTime(obj,tocs,indK,obj.wav.Nk);
                    else
                        waitbar(indK/obj.wav.Nk,d,"Generating Beat Signal. Estimated Time Remaining: " + getEstTime(obj,tocs,indK,obj.wav.Nk));
                    end
                end
            end
        end
        
        function obj = computeTargetFast2(obj,d,R_T_plus_R_R,amplitudeFactor,W_temp,W,w)
            % Fast method 2 - loops over target voxels
            
            obj.sarData = zeros(size(obj.scanner.tx.xyz_m,1),obj.wav.Nk,'single');
            
            if ~obj.isSilent
                tocs = zeros(1,obj.numTargets,'single');
            end
            
            for indTarget = 1:obj.numTargets
                if ~obj.isSilent
                    if obj.isApp
                        if d.CancelRequested
                            warning("Beat Signal not Computed!")
                            obj.sarData = zeros([obj.scanner.sarSize,obj.wav.Nk],'single');
                            return;
                        end
                    else
                        if getappdata(d,'canceling')
                            warning("Beat Signal not Computed!")
                            obj.sarData = zeros([obj.scanner.sarSize,obj.wav.Nk],'single');
                            delete(d);
                            return;
                        end
                    end
                    tic
                end
                
                if ndims(W) == 3
                    W_temp = squeeze(W(:,indTarget,:));
                    if obj.isGPU
                        W_temp = gpuArray(W_temp);
                    end
                elseif ~isscalar(W)
                    W_temp = squeeze(W(:,indTarget));
                    if obj.isGPU
                        W_temp = gpuArray(W_temp);
                    end
                end
                
                if w.isFail
                    [obj,W_temp] = computeWeightsDuring(obj);
                end
                
                temp = exp(1j*obj.wav.k.*R_T_plus_R_R(:,indTarget));
                temp = amplitudeFactor(:,indTarget) .* temp;
                
                W_temp2 = W_temp(:,indTarget,:);
                
                obj.sarData = obj.sarData + single(gather(W_temp2 .* temp));
                % Update the progress dialog
                if ~obj.isSilent
                    tocs(indTarget) = toc;
                    if obj.isApp
                        d.Value = indTarget/obj.app.target.numTargets;
                        d.Message = "Estimated Time Remaining: " + getEstTime(obj,tocs,indTarget,obj.app.target.numTargets);
                    else
                        waitbar(indTarget/obj.numTargets,d,"Generating Beat Signal. Estimated Time Remaining: " + getEstTime(obj,tocs,indTarget,obj.numTargets));
                    end
                end
            end
        end
        
        function obj = computeTargetSlow(obj,d,W_temp,W)
            % Compute sarData by looping over antenna elements and
            % frequency
            
            obj.sarData = zeros(size(obj.scanner.tx.xyz_m,1),obj.wav.Nk,'single');
            
            if obj.isApp
                d.Title = "Generating Echo Signal Using Slow Method";
            else
                disp("Generating Echo Signal Using Slow Method");
            end
            % Always works method
            tocs = zeros(1,obj.wav.Nk*obj.numTargets,'single');
            count = 0;
            for indSAR = 1:size(obj.scanner.rx.xyz_m,1)
                tic
                if ~(isAppEPC(obj.app) || obj.ant.isEPC)
                    obj.R = struct;
                    obj.R.tx = pdist2(obj.scanner.tx.xyz_m(indSAR,:),obj.xyz_m);
                    obj.R.rx = pdist2(obj.scanner.rx.xyz_m(indSAR,:),obj.xyz_m);
                    R_T_plus_R_R = obj.R.tx + obj.R.rx;
                    
                    % Amplitude Factor
                    if obj.isAmplitudeFactor
                        amplitudeFactor = obj.amp./(obj.R.tx .* obj.R.rx);
                    else
                        amplitudeFactor = 1;
                    end
                else
                    obj.R = pdist2(obj.scanner.vx.xyz_m(indSAR,:),obj.xyz_m);
                    
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
                
                for indK = 1:obj.wav.Nk
                    if obj.isApp && d.CancelRequested
                        warning("Beat Signal not Computed!")
                        obj.sarData = zeros([obj.scanner.sarSize,obj.wav.Nk],'single');
                        return;
                    end
                    
                    if ndims(W) == 3
                        W_temp = gpuArray(W(:,:,indK));
                    end
                    W_temp2 = W_temp(indSAR,:,:);
                    
                    count = count + 1;
                    temp = exp(1j*obj.wav.k(indK)*R_T_plus_R_R);
                    
                    if obj.isAmplitudeFactor
                        temp = amplitudeFactor .* temp;
                    end
                    
                    obj.sarData(indSAR,indK) = single(gather(sum(W_temp2 .* temp,2)));
                    % Update the progress dialog
                    tocs(count) = toc;
                    if obj.isApp
                        d.Value = count/(obj.wav.Nk*size(obj.scanner.rx.xyz_m,1));
                        d.Message = "Estimated Time Remaining: " + getEstTime(obj,tocs,count,obj.wav.Nk*size(obj.scanner.rx.xyz_m,1));
                    else
                        disp(count/(obj.wav.Nk*size(obj.scanner.rx.xyz_m,1))*100 + "% Done. Time Remaining: " + getEstTime(obj,tocs,count,obj.wav.Nk*size(obj.scanner.rx.xyz_m,1)));
                    end
                end
            end
        end
        
        function obj = computeTargetLarge(obj,d,R_T_plus_R_R,amplitudeFactor,W_temp,W)
            % Compute sarData for a large target using unique find method
            
            obj.sarData = zeros(size(obj.scanner.tx.xyz_m,1),obj.wav.Nk,'single');
            
            R_T_plus_R_R = gather(R_T_plus_R_R);
            
            if obj.isApp
                selection = uiconfirm(obj.app.UIFigure,'Would you like to use the large target computation method (CPU intensive, likely time intensive)?','Use Large Target Method?',...
                    'Icon','warning');
                if string(selection) == "Cancel"
                    UIFigureCloseRequest(obj.app, 0);
                    return;
                end
                d.Title = "Generating Echo Signal Using Large Target Method";
            end
            [R_unique,~,IC] = unique(R_T_plus_R_R,'stable');
            if obj.isGPU
                R_unique = gpuArray(R_unique);
            end
            % Remember R_T_plus_R_R = reshape(R_unique(IC),size(R_T_plus_R_R));
            
            sizeR = size(R_T_plus_R_R);
            R_T_plus_R_R = 0;
            
            try
                [temp,~,IC2] = unique(mod(R_unique.*obj.wav.k,2*pi),'stable');
                % Remember R_unique = reshape(temp(IC2),size(R_unique));
            catch
                R_unique = gather(R_unique);
                [temp,~,IC2] = unique(mod(R_unique.*obj.wav.k,2*pi),'stable');
                % Remember R_unique = reshape(temp(IC2),size(R_unique));
            end
            
            if obj.isGPU
                temp = gpuArray(temp);
            end
            
            temp2 = gather(exp(1j*temp));
            
            temp3 = reshape(temp2(IC2),[length(R_unique),obj.wav.Nk]);
            
            tocs = zeros(1,obj.wav.Nk,'single');
            for indK = 1:obj.wav.Nk
                if obj.isApp && d.CancelRequested
                    warning("Beat Signal not Computed!")
                    obj.sarData = zeros([obj.scanner.sarSize,obj.wav.Nk],'single');
                    return;
                end
                tic
                
                if ndims(W) == 3
                    W_temp = gpuArray(W(:,:,indK));
                end
                
                temp = reshape(temp3(IC,indK),sizeR);
                
                if obj.isAmplitudeFactor
                    temp = amplitudeFactor .* temp;
                end
                
                obj.sarData(:,indK) = single(gather(sum(W_temp .* temp,2)));
                % Update the progress dialog
                tocs(indK) = toc;
                if obj.isApp
                    d.Value = indK/obj.wav.Nk;
                    d.Message = "Estimated Time Remaining: " + getEstTime(obj,tocs,indK,obj.wav.Nk);
                else
                    disp(count/(obj.wav.Nk*size(obj.scanner.rx.xyz_m,1))*100 + "% Done. Time Remaining: " + getEstTime(obj,tocs,count,obj.wav.Nk*size(obj.scanner.rx.xyz_m,1)));
                end
            end
        end
        
        function [obj,R_T_plus_R_R,amplitudeFactor] = computeDistances(obj)
            % Computes the distances for the Tx and Rx elements
            
            if ~(isAppEPC(obj.app) || obj.ant.isEPC)
                % Get distances
                try
                    obj.R = struct;
                    obj.R.tx = pdist2(obj.scanner.tx.xyz_m,obj.xyz_m);
                    obj.R.rx = pdist2(obj.scanner.rx.xyz_m,obj.xyz_m);
                    R_T_plus_R_R = obj.R.tx + obj.R.rx;
                    
                    % Amplitude Factor
                    if obj.isAmplitudeFactor
                        amplitudeFactor = obj.amp./(obj.R.tx .* obj.R.rx);
                    else
                        amplitudeFactor = obj.amp;
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
                    obj.R = pdist2(obj.scanner.vx.xyz_m,obj.xyz_m);
                    
                    % Get echo signal
                    R_T_plus_R_R = 2*obj.R;
                    if obj.isAmplitudeFactor
                        amplitudeFactor = obj.amp./(obj.R).^2;
                    else
                        amplitudeFactor = obj.amp;
                    end
                    obj.isLong = false;
                catch
                    obj.isLong = true;
                end
            end
        end
        
        function [obj,W,w] = computeWeightsFast(obj)
            % Computes the weight matrix between each antenna element and
            % target voxel based on the antenna pattern
            
            % Determine size
            [az,el] = getAngles(obj.xyz_m(1,:),obj.scanner.rx.unique);
            gainRx = obj.ant.computePattern(az,el);
            
            if isscalar(gainRx)
                W = single(1);
                w.isFail = false;
                return;
            end
            
            if ~obj.isSilent
                % Create the progress dialog
                if obj.isApp
                    d = uiprogressdlg(obj.app.UIFigure,'Title','Computing Antenna Array Weights',...
                        'Message',"Estimated Time Remaining: 0:0:0","Cancelable","on");
                else
                    d = waitbar(0,'1','Name',' Computing Antenna Array Weights...',...
                        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
                end
                tocs = zeros(obj.numTargets,1,'single');
            end
            
            % Try all angles at once
            if ~(isAppEPC(obj.app) || obj.ant.isEPC)
                % Using MIMO elements
                numAntTx = size(obj.scanner.tx.unique.xyz_m,1);
                numAntRx = size(obj.scanner.rx.unique.xyz_m,1);
                azAllTx = zeros(numAntTx,obj.numTargets,1,'single');
                elAllTx = zeros(numAntTx,obj.numTargets,1,'single');
                azAllRx = zeros(numAntRx,obj.numTargets,1,'single');
                elAllRx = zeros(numAntRx,obj.numTargets,1,'single');
            else
                % Using EPC elements
                numAntVx = size(obj.scanner.tx.unique.xyz_m,1);
                azAllVx = zeros(numAntVx,obj.numTargets,1,'single');
                elAllVx = zeros(numAntVx,obj.numTargets,1,'single');
            end
            
            for indTarget = 1:obj.numTargets
                if ~obj.isSilent
                    if obj.isApp
                        if d.CancelRequested
                            warning("Array Weights not Computed!")
                            W = 1;
                            return;
                        end
                    else
                        if getappdata(d,'canceling')
                            warning("Array Weights not Computed!")
                            W = 1;
                            delete(d);
                            return;
                        end
                    end
                    tic
                end
                
                if ~(isAppEPC(obj.app) || obj.ant.isEPC)
                    % Using MIMO elements
                    [az,el] = getAngles(obj.xyz_m(indTarget,:),obj.scanner.tx.unique);
                    
                    azAllTx(:,indTarget) = az;
                    elAllTx(:,indTarget) = el;
                    
                    [az,el] = getAngles(obj.xyz_m(indTarget,:),obj.scanner.rx.unique);
                    
                    azAllRx(:,indTarget) = az;
                    elAllRx(:,indTarget) = el;
                else
                    % Using EPC elements
                    [az,el] = getAngles(obj.xyz_m(indTarget,:),obj.scanner.vx.unique);
                    
                    azAllVx(:,indTarget) = az;
                    elAllVx(:,indTarget) = el;
                end
                
                if ~obj.isSilent
                    % Update the progress dialog
                    tocs(indTarget) = toc;
                    if obj.isApp
                        d.Value = indTarget/obj.numTargets;
                        d.Message = "Estimated Time Remaining: " + getEstTime(obj,tocs,indTarget,obj.numTargets);
                    else
                        waitbar(indTarget/obj.numTargets,d,"Generating Beat Signal. Estimated Time Remaining: " + getEstTime(obj,tocs,indTarget,obj.numTargets));
                    end
                end
            end
            
            if ~obj.isSilent
                delete(d);
            end
            
            if ~(isAppEPC(obj.app) || obj.ant.isEPC)
                % Using MIMO elements
                gainTx = obj.ant.computePattern(azAllTx,elAllTx);
                azAllTx = 0;
                elAllTx = 0;
                gainRx = obj.ant.computePattern(azAllRx,elAllRx);
                azAllRx = 0;
                elAllRx = 0;
                
                W = gainTx(obj.scanner.tx.unique.IC,:,:) .* gainRx(obj.scanner.rx.unique.IC,:,:);
            else
                % Using EPC elements
                gainVx = obj.ant.computePattern(azAllVx,elAllVx);
                W = gainVx(obj.scanner.vx.unique.IC,:,:).^2;
            end
            
            w.isFail = false;
            
            function [az,el] = getAngles(a,b)
                % Will get the azimuth and elevation angles between the
                % vector a (1x3) and set of vectors b.xyz_m (Mx3)
                
                % Elevation
                u = a - b.xyz_m;
                u_inner_nz = sum(u .* b.nz,2) ./ vecnorm(u,2,2);
                el = 90 - acosd(u_inner_nz);
                
                % Azimuth
                u_proj = u - sum(u .* b.nz,2);
                u_proj_inner_nx = sum(u_proj .* b.nx,2) ./ vecnorm(u_proj,2,2);
                az = acosd(u_proj_inner_nx);
            end
        end
        
        function [obj,W,w] = computeWeightsSlow(obj)
            % Computes the weight matrix between each antenna element and
            % target voxel based on the antenna pattern
            
            % Determine size
            [az,el] = getAngles(obj.xyz_m(1,:),obj.scanner.rx.unique);
            gainRx = obj.ant.computePattern(az,el);
            
            if isscalar(gainRx)
                W = single(1);
                w.isFail = false;
                return;
            elseif iscolumn(gainRx)
                W = zeros(size(gainRx,1),obj.numTargets,'single');
            elseif ismatrix(gainRx)
                W = zeros(size(gainRx,1),obj.numTargets,obj.wav.Nk,'single');
            end
            
            if ~obj.isSilent
                % Create the progress dialog
                if obj.isApp
                    d = uiprogressdlg(obj.app.UIFigure,'Title','Computing Antenna Array Weights',...
                        'Message',"Estimated Time Remaining: 0:0:0","Cancelable","on");
                else
                    d = waitbar(0,'1','Name',' Computing Antenna Array Weights...',...
                        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
                end
                tocs = zeros(obj.numTargets,1,'single');
            end
            
            % Do one target voxel at a time
            
            for indTarget = 1:obj.numTargets
                if ~obj.isSilent
                    if obj.isApp
                        if d.CancelRequested
                            warning("Array Weights not Computed!")
                            W = 1;
                            return;
                        end
                    else
                        if getappdata(d,'canceling')
                            warning("Array Weights not Computed!")
                            W = 1;
                            delete(d);
                            return;
                        end
                    end
                    tic
                end
                
                if ~(isAppEPC(obj.app) || obj.ant.isEPC)
                    % Using MIMO elements
                    [az,el] = getAngles(obj.xyz_m(indTarget,:),obj.scanner.tx.unique);
                    gainTx = obj.ant.computePattern(az,el);
                    
                    [az,el] = getAngles(obj.xyz_m(indTarget,:),obj.scanner.rx.unique);
                    gainRx = obj.ant.computePattern(az,el);
                    
                    W(:,indTarget,:) = gainTx .* gainRx;
                else
                    % Using EPC elements
                    [az,el] = getAngles(obj.xyz_m(indTarget,:),obj.scanner.vx.unique);
                    gain = obj.ant.computePattern(az,el);
                    
                    W(:,indTarget,:) = gain.^2;
                end
                
                if ~obj.isSilent
                    % Update the progress dialog
                    tocs(indTarget) = toc;
                    if obj.isApp
                        d.Value = indTarget/obj.numTargets;
                        d.Message = "Estimated Time Remaining: " + getEstTime(obj,tocs,indTarget,obj.numTargets);
                    else
                        waitbar(indTarget/obj.numTargets,d,"Generating Beat Signal. Estimated Time Remaining: " + getEstTime(obj,tocs,indTarget,obj.numTargets));
                    end
                end
            end
            
            if ~obj.isSilent
                delete(d);
            end
            
            w.isFail = false;
            
            function [az,el] = getAngles(a,b)
                % Will get the azimuth and elevation angles between the
                % vector a (1x3) and set of vectors b (Mx3)
                
                % Elevation
                u = a - b.xyz_m;
                u_inner_nz = sum(u .* b.nz,2) ./ vecnorm(u,2,2);
                el = 90 - acosd(u_inner_nz);
                
                % Azimuth
                u_proj = u - sum(u .* b.nz,2);
                u_proj_inner_nx = sum(u_proj .* b.nx,2) ./ vecnorm(u_proj,2,2);
                az = acosd(u_proj_inner_nx);
                
                el = el(b.IC);
                az = az(b.IC);
            end
        end
        
        function obj = InitializeFigures(obj)
            % Initialize the figures
            
            CloseFigures(obj);
            set(0,'DefaultFigureWindowStyle','docked')
            
            % AntAxes
            obj.fig.f = figure;
            obj.fig.h = handle(axes);
        end
        
        function CloseFigures(obj)
            % Attempt to close the figures
            
            try
                close(obj.fig.f)
            catch
            end
        end
        
        function Display(obj)
            % Plots either the MIMO or EPC array
            
            if isempty(obj.fig.f) || ~isvalid(obj.fig.h)
                InitializeFigures(obj)
            end
            
            if isempty(obj.xyz_m)
                return;
            end
            
            if isAppEPC(obj.app) || obj.ant.isEPC
                DisplayEPC(obj);
            else
                DisplayMIMO(obj);
            end
        end
        
        function DisplayMIMO(obj)
            % Display the target with the MIMO array
            
            if isempty(obj.xyz_m)
                return;
            end
            
            h = obj.fig.h;
            hold(h,'off')
            temp = obj.scanner.tx.xyz_m;
            scatter3(h,temp(:,1),temp(:,3),temp(:,2),'.r')
            hold(h,'on')
            temp = obj.scanner.rx.xyz_m;
            scatter3(h,temp(:,1),temp(:,3),temp(:,2),'.b')
            temp = obj.xyz_m;
            scatter3(h,temp(:,1),temp(:,3),temp(:,2),'.k')
            
            xlabel(h,"x (m)","interpreter","latex")
            temp1 = obj.scanner.tx.xyz_m(:,1);
            temp2 = obj.scanner.rx.xyz_m(:,1);
            temp3 = obj.xyz_m(:,1);
            xlim(h,[min([min(temp1),min(temp2),min(temp3)])-0.01,max([max(temp1),max(temp2),max(temp3)])+0.01])
            ylabel(h,"z (m)","interpreter","latex")
            temp1 = obj.scanner.tx.xyz_m(:,3);
            temp2 = obj.scanner.rx.xyz_m(:,3);
            temp3 = obj.xyz_m(:,3);
            ylim(h,[min([min(temp1),min(temp2),min(temp3)])-0.01,max([max(temp1),max(temp2),max(temp3)])+0.01])
            zlabel(h,"y (m)","interpreter","latex")
            temp1 = obj.scanner.tx.xyz_m(:,2);
            temp2 = obj.scanner.rx.xyz_m(:,2);
            temp3 = obj.xyz_m(:,2);
            zlim(h,[min([min(temp1),min(temp2),min(temp3)])-0.01,max([max(temp1),max(temp2),max(temp3)])+0.01])
            title(h,"MIMO Aperture Image Scenario","interpreter","latex")
            legend(h,"Tx","Rx","Target","interpreter","latex")
            
            view(h,3)
            daspect(h,[1 1 1])
        end
        
        function DisplayEPC(obj)
            % Display the target with the EPC array
            
            h = obj.fig.h;
            hold(h,'off')
            temp = obj.scanner.vx.xyz_m;
            scatter3(h,temp(:,1),temp(:,3),temp(:,2),'.b')
            hold(h,'on')
            temp = obj.xyz_m;
            scatter3(h,temp(:,1),temp(:,3),temp(:,2),'.k')
            
            xlabel(h,"x (m)","interpreter","latex")
            temp1 = obj.scanner.tx.xyz_m(:,1);
            temp2 = obj.scanner.rx.xyz_m(:,1);
            temp3 = obj.xyz_m(:,1);
            xlim(h,[min([min(temp1),min(temp2),min(temp3)])-0.01,max([max(temp1),max(temp2),max(temp3)])+0.01])
            ylabel(h,"z (m)","interpreter","latex")
            temp1 = obj.scanner.tx.xyz_m(:,3);
            temp2 = obj.scanner.rx.xyz_m(:,3);
            temp3 = obj.xyz_m(:,3);
            ylim(h,[min([min(temp1),min(temp2),min(temp3)])-0.01,max([max(temp1),max(temp2),max(temp3)])+0.01])
            zlabel(h,"y (m)","interpreter","latex")
            temp1 = obj.scanner.tx.xyz_m(:,2);
            temp2 = obj.scanner.rx.xyz_m(:,2);
            temp3 = obj.xyz_m(:,2);
            zlim(h,[min([min(temp1),min(temp2),min(temp3)])-0.01,max([max(temp1),max(temp2),max(temp3)])+0.01])
            title(h,"Virtual Aperture Image Scenario","interpreter","latex")
            legend(h,"Vx","Target","interpreter","latex")
            
            view(h,3)
            daspect(h,[1 1 1])
        end
        
        function Load(obj,wav,ant,scanner)
            % Loads the target from a file. If the additional arguments are
            % present, sets obj.wav to wav, obj.ant to ant, and obj.scanner
            % to scanner so the object handles are synced
            
            if obj.isApp
                loadPathFull = "./saved/targets/" + obj.app.TargetLoadNameEditField.Value + ".mat";
            else
                loadPathFull = [];
            end
            if ~exist(loadPathFull,'file')
                [filename,pathname] = uigetfile("./saved/targets/*.mat","Select Desired File to Load");
                
                if filename == 0
                    showErrorMessage(obj,"Target file not loaded!","Load Error");
                    return;
                else
                    loadPathFull = string(pathname) + string(filename);
                end
            end
            
            load(loadPathFull,"savedtarget","savedwav","savedant","savedscanner")
            
            obj = getFields(obj,savedtarget,["app","fig","wav","ant","scanner","isApp"]);
            
            if nargin == 1
                showErrorMessage(obj,"Using savedwav, savedant, and savedscanner from file, might not be properly linked","Warning");
                
                if obj.isApp
                    obj.wav = THzWaveformParameters(obj.app);
                    obj.ant = THzAntennaArray(obj.wav,obj.app);
                    obj.scanner = THzScanner(obj.ant,obj.app);
                else
                    obj.wav = THzWaveformParameters();
                    obj.ant = THzAntennaArray(obj.wav);
                    obj.scanner = THzScanner(obj.ant);
                end
                
                % Load the savedwav
                obj.wav = getFields(obj.wav,savedwav,["app","isApp"]);
                obj.wav = obj.wav.Compute();
                obj.wav = obj.wav.Display();
                
                % Load the savedant
                obj.ant = getFields(obj.ant,savedant,["app","fig","wav","isApp"]);
                obj.ant.tx.xy_m = single(obj.ant.tableTx);
                obj.ant.rx.xy_m = single(obj.ant.tableRx);
                obj.app.TxTable.Data = array2table(obj.tx.xy_m);
                obj.app.RxTable.Data = array2table(obj.rx.xy_m);
                obj.ant = obj.ant.Compute();
                obj.ant.wav = obj.wav;
                
                % Load the savedscanner
                obj.scanner = getFields(obj.scanner,savedscanner,["app","fig","ant","isApp"]);
                obj.scanner = obj.scanner.Compute();
                if obj.isApp
                    obj.app.SARMethodDropDown.Value = obj.scanner.method;
                    
                    if obj.app.XStepSwitch.Value == "λ"
                        obj.app.XStepSizeEditField.Value = obj.scanner.xStep_m/obj.app.wav.lambda_m;
                    elseif obj.app.XStepSwitch.Value == "mm"
                        obj.app.XStepSizeEditField.Value = obj.scanner.xStep_m*1e3;
                    end
                    
                    if obj.app.YStepSwitch.Value == "λ"
                        obj.app.YStepSizeEditField.Value = obj.scanner.yStep_m/obj.app.wav.lambda_m;
                    elseif obj.app.YStepSwitch.Value == "mm"
                        obj.app.YStepSizeEditField.Value = obj.scanner.yStep_m*1e3;
                    end
                    
                    obj.app.ThetaSizedegEditField.Value = obj.scanner.thetaMax_deg;
                    obj.app.NumXStepsEditField.Value = obj.scanner.numX;
                    obj.app.NumYStepsEditField.Value = obj.scanner.numY;
                    obj.app.NumThetaStepsEditField.Value = obj.scanner.numTheta;
                    
                    obj.app.XSizemEditField.Value = (obj.scanner.numX-1) * obj.scanner.xStep_m;
                    obj.app.YSizemEditField.Value = (obj.scanner.numY-1) * obj.scanner.yStep_m;
                end
                obj.scanner.ant = obj.ant;
            elseif nargin == 4
                obj.wav = wav;
                obj.ant = ant;
                obj.ant.wav = obj.wav;
                obj.scanner = scanner;
                obj.scanner.ant = obj.ant;
            end
        end
        
        function Save(obj)
            % Saves the THzTarget object handle as a strcut for reloading
            
            if obj.isApp
                savePathFull = "./saved/targets/" + obj.app.TargetSaveNameEditField.Value + ".mat";
            else
                savePathFull = [];
            end
            if ~exist(savePathFull,'file')
                [filename,pathname] = uiputfile("./saved/targets/*.mat","Select Desired File Location + Name for Save");
                
                if filename == 0
                    showErrorMessage(obj,"Target file not saved!","Save Error");
                    return;
                else
                    savePathFull = string(pathname) + string(filename);
                end
            end
            
            savedtarget = struct();
            savedtarget = getFields(savedtarget,obj,["app","fig","wav","ant","scanner","isApp"]);
            
            savedwav = struct();
            savedwav = getFields(savedwav,obj.wav,["app","isApp"]);
            
            savedant = struct();
            savedant = getFields(savedant,obj.ant,["app","fig","wav","isApp"]);
            
            savedscanner = struct();
            savedscanner = getFields(savedscanner,obj.scanner,["app","fig","ant","isApp"]);
            
            save(savePathFull,"savedtarget","savedwav","savedant","savedscanner");
        end
        
        function obj = VerifyGPU(obj)
            % Verifies if the GPU can be used
            
            % If the checkbox in the app is unchecked
            if obj.isApp && ~obj.app.UseGPUCheckBox.Value
                obj.isGPU = false;
                return;
            end
            
            if obj.isGPU && ~obj.isGPUVerified
                try
                    reset(gpuDevice);
                    obj.isGPUVerified = true;
                catch
                    obj.isGPU = false;
                    obj.isGPUVerified = false;
                    showErrorMessage(obj,"Unable to locate Nvidia GPU","No GPU available");
                    if obj.isApp
                        obj.app.UseGPUCheckBox.Value = false;
                    end
                    return;
                end
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
    end
end