% sarTarget_app object holds the properties and methods
% used in the FMCW MIMO-SAR scenario as specified by the user
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

classdef sarTarget_app < handle
    properties
        isGPU = true                % Boolean whether or not to use the GPU for beat signal computation
        isMIMO                      % Boolean whether or not to use the MIMO physical element locations instead of the equivalent phase center virtual element locations
        
        isAmplitudeFactor = false   % Boolean whether or not to include the amplitude factor (path loss) in the beat signal computation
        
        isLong = false              % Boolean whether or not to use the long beat signal computation method
        numTargets = 0              % Number of target voxels
        xyz_m                       % Target x-y-z locations as a (numTargets)x3 array
        amp                         % Column vector of target reflectivities
        R                           % Radial distances from antennas to target (MIMO or EPC)
        
        png = struct('isLoaded',false,'fileNameLoaded',"") % Structure containing the parameters of the PNG target
        stl = struct('isLoaded',false,'fileNameLoaded',"") % Structure containing the parameters of the STL target
        rp                          % Structure containing the parameters of the random point targets
        
        sarData                     % Computed beat signal
        
        fig = struct('f',[],'h',[]) % Structure containing the figure and handle used for showing the target
    end
    
    methods
        function obj = sarTarget_app(app)
            % Verifies the GPU can be used, initializes the figures, then
            % updates the target with the values in the app
            
            obj = verifyGPU(obj,app);
            obj = initializeFigures(obj);
            obj = update(obj,app);
            
            obj.stl.isLoaded = false;
        end
        
        function obj = update(obj,app)
            % Update the target with the values in the app
            
            app.BeatSignalComputedLamp.Color = "red";
            app.ImageReconstructionCompleteLamp.Color = "red";
            obj = getTarget(obj,app);
            
            if isempty(obj.xyz_m)
                return;
            end
            
            if verifyMIMO(obj,app)
                displayTarget(obj,app);
            else
                displayVirtualTarget(obj,app);
            end
        end
        
        function obj = getTarget(obj,app)
            % Gets the target using the values in the app
            
            obj.isMIMO = verifyMIMO(obj,app);
            
            obj.isAmplitudeFactor = app.PathLossCheckBox.Value;
            
            obj.xyz_m = [];
            obj.amp = [];
            
            if app.UseTableofTargetsCheckBox.Value
                obj = getTargetTable(obj,app);
            end
            if app.UsePNGFileCheckBox.Value
                obj = getPNGTarget(obj,app);
            end
            if app.UseSTLFileCheckBox.Value
                obj = getSTLTarget(obj,app);
            end
            
            obj.numTargets = size(obj.xyz_m,1);
            obj.xyz_m = single(obj.xyz_m);
            obj.amp = single(obj.amp(:).');
            
            app.TotalNumTargetsEditField.Value = obj.numTargets;
        end
        
        function obj = getTargetTable(obj,app)
            % Gets the positions and amplitudes from the table of targets
            
            if isempty(app.TargetTable.Data)
                return;
            end
            temp_table = single(table2array(app.TargetTable.Data));
            
            temp_xyz_m = temp_table(:,1:3);
            temp_amp = temp_table(:,4);
            obj.xyz_m = cat(1,obj.xyz_m,temp_xyz_m);
            obj.amp = cat(1,obj.amp,temp_amp);
        end
        
        function obj = getPNGTarget(obj,app)
            % Gets the positions and amplitudes from the PNG target
            
            obj = getPNGParameters(obj,app);
            
            try
                tMat = imread("./saved/pngstl/" + obj.png.fileName);
            catch
                selection = uiconfirm(app.UIFigure,'Would you liked to locate the PNG file?','PNG File Not Found',...
                    'Icon','warning');
                if string(selection) == "Cancel"
                    warning("PNG file not loaded!");
                    return;
                end
                
                [filename,pathname] = uigetfile("./saved/pngstl/*.png","Select Desired PNG File");
                
                if filename == 0
                    warning("PNG file not loaded!");
                    return;
                end
                app.PNGFileNameEditField.Value = filename;
                tMat = imread(string(pathname) + string(filename));
            end
            tMat = tMat(:,:,1);
            tMat(tMat<64) = 0;
            tMat(tMat>0) = 1;
            tMat = ~tMat;
            tMat = fliplr(tMat);
            
            [numY,numX] = size(tMat);
            xAxisT = obj.png.xStep_m * (-(numX-1)/2 : (numX-1)/2) + obj.png.xOffset_m;
            yAxisT = obj.png.yStep_m * (-(numY-1)/2 : (numY-1)/2) + obj.png.yOffset_m;
            zAxisT = obj.png.zOffset_m;
            
            [zT,xT,yT] = meshgrid(zAxisT,xAxisT,yAxisT);
            temp_xyz_m = reshape(permute([xT,yT,zT],[1 3 2]),[],3);
            
            indT = rot90(tMat,-1)==true;
            temp_xyz_m = single(temp_xyz_m(indT,:));
            
            temp_pc = pointCloud(temp_xyz_m);
            temp_pc = pcdownsample(temp_pc,'random',1/obj.png.downsampleFactor);
            temp_xyz_m = temp_pc.Location;
            temp_amp = ones(size(temp_xyz_m,1),1)*obj.png.reflectivity;
            
            obj.xyz_m = cat(1,obj.xyz_m,temp_xyz_m);
            obj.amp = cat(1,obj.amp,temp_amp);
            
            app.NumTargetsEditField.Value = size(temp_xyz_m,1);
        end
        
        function obj = getPNGParameters(obj,app)
            % Gets the parameters for the PNG target from the app
            
            obj.png.fileName = app.PNGFileNameEditField.Value;
            obj.png.xStep_m = app.XPixelSizemEditField.Value;
            obj.png.yStep_m = app.YPixelSizemEditField.Value;
            obj.png.xOffset_m = app.XOffsetmEditField.Value;
            obj.png.yOffset_m = app.YOffsetmEditField.Value;
            obj.png.zOffset_m = app.ZOffsetmEditField.Value;
            obj.png.reflectivity = app.ReflectivityEditField.Value;
            obj.png.downsampleFactor = app.DownsampleFactorEditField.Value;
        end
        
        function obj = getSTLTarget(obj,app)
            % Gets the positions and amplitudes from the STL target
            
            if ~obj.stl.isLoaded
                selection = uiconfirm(app.UIFigure,'Would you liked to load the STL file?','STL File Not Yet Loaded',...
                    'Icon','warning');
                if string(selection) == "Cancel"
                    warning("STL file not loaded!");
                    return;
                end
                obj = loadSTL(obj,app);
            end
            
            obj = getSTLParameters(obj,app);
            
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
            
            app.NumTargetsEditField_2.Value = size(temp_xyz_m,1);
        end
        
        function obj = getSTLParameters(obj,app)
            % Gets the parameters for the STL target from the app
            
            obj.stl.fileName = app.STLFileNameEditField.Value;
            obj.stl.zCrop_m = app.ZCropmEditField.Value;
            obj.stl.xOffset_m = app.XOffsetmEditField_2.Value;
            obj.stl.yOffset_m = app.YOffsetmEditField_2.Value;
            obj.stl.zOffset_m = app.ZOffsetmEditField_2.Value;
            obj.stl.reflectivity = app.ReflectivityEditField_2.Value;
            obj.stl.downsampleFactor = app.DownsampleFactorEditField_2.Value;
        end
        
        function obj = loadSTL(obj,app)
            % Loads the STL file and updates stl.isLoaded
            
            app.LoadSTLLamp.Color = "yellow";
            drawnow
            obj = getSTLParameters(obj,app);
            try
                [~,obj.stl.v] = stlread2011("./saved/pngstl/" + obj.stl.fileName);
            catch
                selection = uiconfirm(app.UIFigure,'Would you liked to locate the STL file?','STL File Not Found',...
                    'Icon','warning');
                if string(selection) == "Cancel"
                    warning("STL file not loaded!");
                    obj.stl.isLoaded = false;
                    return;
                end
                
                [filename,pathname] = uigetfile("./saved/pngstl/*.stl","Select Desired STL File");
                
                if filename == 0
                    warning("STL file not loaded!");
                    obj.stl.isLoaded = false;
                    return;
                end
                app.STLFileNameEditField.Value = filename;
                [~,obj.stl.v] = stlread2011(string(pathname) + string(filename));
            end
            obj.stl.v = obj.stl.v*1e-3;
            obj.stl.isLoaded = true;
            app.LoadSTLLamp.Color = "green";
        end
        
        function obj = getRandomTarget(obj,app)
            % Generates the random targets given the parameters
            
            obj = getRandomParameters(obj,app);
            
            temp_x = obj.rp.xMin_m + (obj.rp.xMax_m - obj.rp.xMin_m)*rand(obj.rp.numTargets,1);
            temp_y = obj.rp.yMin_m + (obj.rp.yMax_m - obj.rp.yMin_m)*rand(obj.rp.numTargets,1);
            temp_z = obj.rp.zMin_m + (obj.rp.zMax_m - obj.rp.zMin_m)*rand(obj.rp.numTargets,1);
            
            temp_xyz_m = [temp_x,temp_y,temp_z];
            temp_amp = obj.rp.ampMin + (obj.rp.ampMax - obj.rp.ampMin)*rand(obj.rp.numTargets,1);
            
            app.TargetTable.Data = array2table([temp_xyz_m,temp_amp]);
        end
        
        function obj = getRandomParameters(obj,app)
            % Gets the parameters for the ramdom targets from the app
            
            obj.rp.numTargets = app.NumTargetsEditField_3.Value;
            obj.rp.xMin_m = app.XMinmEditField.Value;
            obj.rp.xMax_m = app.XMaxmEditField.Value;
            obj.rp.yMin_m = app.YMinmEditField.Value;
            obj.rp.yMax_m = app.YMaxmEditField.Value;
            obj.rp.zMin_m = app.ZMinmEditField.Value;
            obj.rp.zMax_m = app.ZMaxmEditField.Value;
            obj.rp.ampMin = app.ReflectivityMinEditField.Value;
            obj.rp.ampMax = app.ReflectivityMaxEditField.Value;
        end
        
        function obj = computeTarget(obj,app)
            % Computes the beat signal. First attempts the fast method then
            % the large target method, if the fast method fails. If the
            % large target method fails, it attempts the slow method
            
            obj = getTarget(obj,app);
            
            if obj.isGPU
                reset(gpuDevice)
            end
            
            if obj.isMIMO
                % Get distances
                try
                    obj.R = struct;
                    obj.R.tx = pdist2(app.sar.rx.xyz_m,obj.xyz_m);
                    obj.R.rx = pdist2(app.sar.tx.xyz_m,obj.xyz_m);
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
                    obj.R = pdist2(app.sar.vx.xyz_m,obj.xyz_m);
                    
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
            % Create the progress dialog
            d = uiprogressdlg(app.UIFigure,'Title','Generating Echo Signal',...
                'Message',"Estimated Time Remaining: 0:0:0","Cancelable","on");
            
            obj.sarData = single(zeros(size(app.sar.tx.xyz_m,1),app.fmcw.ADCSamples));
            
            try
                if obj.isLong
                    error("oops");
                end
                
                % Fast method
                tocs = zeros(1,app.fmcw.ADCSamples);
                for indK = 1:app.fmcw.ADCSamples
                    if d.CancelRequested
                        warning("Beat Signal not Computed!")
                        obj.sarData = single(zeros([app.sar.sarSize,app.fmcw.ADCSamples]));
                        return;
                    end
                    tic
                    temp = exp(1j*app.fmcw.k(indK)*R_T_plus_R_R);
                    if obj.isAmplitudeFactor
                        temp = amplitudeFactor .* temp;
                    end
                    
                    obj.sarData(:,indK) = single(gather(sum(temp,2)));
                    % Update the progress dialog
                    tocs(indK) = toc;
                    d.Value = indK/app.fmcw.ADCSamples;
                    d.Message = "Estimated Time Remaining: " + getEstTime(obj,tocs,indK,app.fmcw.ADCSamples);
                end
            catch
                try
                    % Fast method 2
                    tocs = zeros(1,app.target.numTargets);
                    for indTarget = 1:app.target.numTargets
                        if d.CancelRequested
                            warning("Beat Signal not Computed!")
                            obj.sarData = single(zeros([app.sar.sarSize,app.fmcw.ADCSamples]));
                            return;
                        end
                        tic
                        temp = exp(1j*app.fmcw.k.*R_T_plus_R_R(:,indTarget));
                        if obj.isAmplitudeFactor
                            temp = amplitudeFactor .* temp;
                        end
                        
                        obj.sarData = obj.sarData + single(gather(temp));
                        % Update the progress dialog
                        tocs(indTarget) = toc;
                        d.Value = indTarget/app.target.numTargets;
                        d.Message = "Estimated Time Remaining: " + getEstTime(obj,tocs,indTarget,app.target.numTargets);
                    end
                catch
                    try
                        obj = computeTargetLarge(obj,app,R_T_plus_R_R,d);
                    catch
                        obj = computeTargetSlow(obj,app,d);
                    end
                end
            end
            
            delete(d);
            
            % Reshape echo signal
            obj.sarData = reshape(obj.sarData,[app.sar.sarSize,app.fmcw.ADCSamples]);
        end
        
        function obj = computeTargetSlow(obj,app,d)
            d.Title = "Generating Echo Signal Using Slow Method";
            % Always works method
            tocs = single(zeros(1,app.fmcw.ADCSamples*obj.numTargets));
            count = 0;
            for indSAR = 1:size(app.sar.rx.xyz_m,1)
                tic
                if obj.isMIMO
                    obj.R = struct;
                    obj.R.tx = pdist2(app.sar.tx.xyz_m(indSAR,:),obj.xyz_m);
                    obj.R.rx = pdist2(app.sar.rx.xyz_m(indSAR,:),obj.xyz_m);
                    R_T_plus_R_R = obj.R.tx + obj.R.rx;
                    
                    % Amplitude Factor
                    if obj.isAmplitudeFactor
                        amplitudeFactor = obj.amp./(obj.R.tx .* obj.R.rx);
                    else
                        amplitudeFactor = 1;
                    end
                else
                    obj.R = pdist2(app.sar.vx.xyz_m(indSAR,:),obj.xyz_m);
                    
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
                
                for indK = 1:app.fmcw.ADCSamples
                    if d.CancelRequested
                        warning("Beat Signal not Computed!")
                        obj.sarData = single(zeros([app.sar.sarSize,app.fmcw.ADCSamples]));
                        return;
                    end
                    count = count + 1;
                    temp = exp(1j*app.fmcw.k(indK)*R_T_plus_R_R);
                    if obj.isAmplitudeFactor
                        temp = amplitudeFactor .* temp;
                    end
                    
                    obj.sarData(indSAR,indK) = single(gather(sum(temp,2)));
                    % Update the progress dialog
                    tocs(count) = toc;
                    d.Value = count/(app.fmcw.ADCSamples*size(app.sar.rx.xyz_m,1));
                    d.Message = "Estimated Time Remaining: " + getEstTime(obj,tocs,count,app.fmcw.ADCSamples*size(app.sar.rx.xyz_m,1));
                end
            end
        end
        
        function obj = computeTargetLarge(obj,app,R_T_plus_R_R)
            R_T_plus_R_R = gather(R_T_plus_R_R);
            
            d.Title = "Generating Echo Signal Using Large Target Method";
            [R_unique,~,IC] = unique(R_T_plus_R_R,'stable');
            if obj.isGPU
                R_unique = gpuArray(R_unique);
            end
            % Remember R_T_plus_R_R = reshape(R_unique(IC),size(R_T_plus_R_R));
            
            sizeR = size(R_T_plus_R_R);
            R_T_plus_R_R = 0;
            
            try
                [temp,~,IC2] = unique(mod(R_unique.*app.fmcw.k,2*pi),'stable');
                % Remember R_unique = reshape(temp(IC2),size(R_unique));
            catch
                R_unique = gather(R_unique);
                [temp,~,IC2] = unique(mod(R_unique.*app.fmcw.k,2*pi),'stable');
                % Remember R_unique = reshape(temp(IC2),size(R_unique));
            end
            if obj.isGPU
                temp = gpuArray(temp);
            end
            temp2 = gather(exp(1j*temp));
            
            temp3 = reshape(temp2(IC2),[length(R_unique),app.fmcw.ADCSamples]);
            
            tocs = zeros(1,app.fmcw.ADCSamples);
            for indK = 1:app.fmcw.ADCSamples
                if d.CancelRequested
                    warning("Beat Signal not Computed!")
                    obj.sarData = single(zeros([app.sar.sarSize,app.fmcw.ADCSamples]));
                    return;
                end
                tic
                temp = reshape(temp3(IC,indK),sizeR);
                if obj.isAmplitudeFactor
                    temp = amplitudeFactor .* temp;
                end
                
                obj.sarData(:,indK) = single(gather(sum(temp,2)));
                % Update the progress dialog
                tocs(indK) = toc;
                d.Value = indK/app.fmcw.ADCSamples;
                d.Message = "Estimated Time Remaining: " + getEstTime(obj,tocs,indK,app.fmcw.ADCSamples);
            end
        end
        
        function obj = initializeFigures(obj)
            % Initialize the figures
            
            closeFigures(obj);
            set(0,'DefaultFigureWindowStyle','docked')
            
            % AntAxes
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
        
        function displayTarget(obj,app)
            % Display the target with the MIMO-SAR scenario
            
            if isempty(obj.xyz_m)
                return;
            end
            
            h = obj.fig.h;
            hold(h,'off')
            temp = app.sar.tx.xyz_m;
            scatter3(h,temp(:,1),temp(:,3),temp(:,2),'.r')
            hold(h,'on')
            temp = app.sar.rx.xyz_m;
            scatter3(h,temp(:,1),temp(:,3),temp(:,2),'.b')
            temp = obj.xyz_m;
            scatter3(h,temp(:,1),temp(:,3),temp(:,2),'.k')
            
            xlabel(h,"x (m)")
            temp1 = app.sar.tx.xyz_m(:,1);
            temp2 = app.sar.rx.xyz_m(:,1);
            temp3 = obj.xyz_m(:,1);
            xlim(h,[min([min(temp1),min(temp2),min(temp3)])-0.01,max([max(temp1),max(temp2),max(temp3)])+0.01])
            ylabel(h,"z (m)")
            temp1 = app.sar.tx.xyz_m(:,3);
            temp2 = app.sar.rx.xyz_m(:,3);
            temp3 = obj.xyz_m(:,3);
            ylim(h,[min([min(temp1),min(temp2),min(temp3)])-0.01,max([max(temp1),max(temp2),max(temp3)])+0.01])
            zlabel(h,"y (m)")
            temp1 = app.sar.tx.xyz_m(:,2);
            temp2 = app.sar.rx.xyz_m(:,2);
            temp3 = obj.xyz_m(:,2);
            zlim(h,[min([min(temp1),min(temp2),min(temp3)])-0.01,max([max(temp1),max(temp2),max(temp3)])+0.01])
            title(h,"MIMO Aperture Image Scenario")
            legend(h,"Tx","Rx","Target")
            
            view(h,3)
            daspect(h,[1 1 1])
        end
        
        function displayVirtualTarget(obj,app)
            % Display the target with the virtual array SAR scenario
            
            h = obj.fig.h;
            hold(h,'off')
            temp = app.sar.vx.xyz_m;
            scatter3(h,temp(:,1),temp(:,3),temp(:,2),'.b')
            hold(h,'on')
            temp = obj.xyz_m;
            scatter3(h,temp(:,1),temp(:,3),temp(:,2),'.k')
            
            xlabel(h,"x (m)")
            temp1 = app.sar.tx.xyz_m(:,1);
            temp2 = app.sar.rx.xyz_m(:,1);
            temp3 = obj.xyz_m(:,1);
            xlim(h,[min([min(temp1),min(temp2),min(temp3)])-0.01,max([max(temp1),max(temp2),max(temp3)])+0.01])
            ylabel(h,"z (m)")
            temp1 = app.sar.tx.xyz_m(:,3);
            temp2 = app.sar.rx.xyz_m(:,3);
            temp3 = obj.xyz_m(:,3);
            ylim(h,[min([min(temp1),min(temp2),min(temp3)])-0.01,max([max(temp1),max(temp2),max(temp3)])+0.01])
            zlabel(h,"y (m)")
            temp1 = app.sar.tx.xyz_m(:,2);
            temp2 = app.sar.rx.xyz_m(:,2);
            temp3 = obj.xyz_m(:,2);
            zlim(h,[min([min(temp1),min(temp2),min(temp3)])-0.01,max([max(temp1),max(temp2),max(temp3)])+0.01])
            title(h,"Virtual Aperture Image Scenario")
            legend(h,"Vx","Target")
            
            view(h,3)
            daspect(h,[1 1 1])
        end
        
        function saveTarget(obj,app)
            % Save the sarTarget_app object
            
            savePathFull = "./saved/targets/" + app.TargetSaveNameEditField.Value + ".mat";
            if exist(savePathFull,'file')
                selection = uiconfirm(app.UIFigure,'Are you sure you want to overwrite?','Confirm Overwrite',...
                    'Icon','warning');
                if string(selection) == "Cancel"
                    warning("Target not saved!");
                    return;
                end
            end
            
            savedtarget = obj;
            savedtarget.fig.f = [];
            save(savePathFull,"savedtarget");
        end
        
        function obj = loadTarget(obj,app)
            % Load a sarTarget_app object from a file
            
            loadPathFull = "./saved/targets/" + app.TargetLoadNameEditField.Value + ".mat";
            if ~exist(loadPathFull,'file')
                uiconfirm(app.UIFigure,"No file called " + app.TargetLoadNameEditField.Value + ".mat to load",'Cannot Load',...
                    "Options",{'OK'},'Icon','warning');
                warning("Target not loaded!");
                return;
            end
            
            load(loadPathFull,"savedtarget");
            
            savedtarget.fig = obj.fig;
            obj = savedtarget;
            
            app.BeatSignalComputedLamp.Color = "red";
            app.ImageReconstructionCompleteLamp.Color = "red";
            
            if isempty(obj.xyz_m)
                return;
            end
            
            if verifyMIMO(obj,app)
                displayTarget(obj,app);
            else
                displayVirtualTarget(obj,app);
            end
        end
        
        function obj = verifyGPU(obj,app)
            % Verifies if the GPU can be used
            
            if app.UseGPUCheckBox.Value
                try
                    reset(gpuDevice);
                catch
                    obj.isGPU = false;
                    app.UseGPUCheckBox.Value = false;
                    uiconfirm(app.UIFigure,"Unable to locate Nvidia GPU",'No GPU Available',...
                        "Options",{'OK'},'Icon','warning');
                    return;
                end
                obj.isGPU = true;
            else
                obj.isGPU = false;
            end
        end
        
        function tf = verifyMIMO(obj,app)
            % Verifies if the user has specified MIMO or EPC in the app
            
            if app.MIMOSwitch.Value == "Use MIMO Array"
                tf = true;
            elseif app.MIMOSwitch.Value == "Use EPC Virtual Elements"
                tf = false;
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