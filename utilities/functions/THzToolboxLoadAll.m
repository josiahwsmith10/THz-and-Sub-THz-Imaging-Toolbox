% Loads all simulation objects from a file saved by THzToolboxSaveAll
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

function [wav,ant,scanner,target,im] = THzToolboxLoadAll(loadPathFull,isDisplay)
if nargin < 1
    loadPathFull = [];
end

if nargin < 2
    isDisplay = true;
end

if ~exist(loadPathFull,'file')
    [filename,pathname] = uiputfile("./*.mat","Select Desired File Location + Name for Save");
    
    if filename == 0
        warning("All scenario file not saved!");
        return;
    else
        loadPathFull = string(pathname) + string(filename);
    end
end

load(loadPathFull,"savedwav","savedant","savedscanner","savedtarget","savedreconstructor","savedim")

% Get wav fields from struct
wav = THzWaveformParameters();
wav = getFields(wav,savedwav,["app","isApp"]);

% Get ant fields from struct
ant = THzAntennaArray(wav);
ant = getFields(ant,savedant,["app","fig","wav","isApp"]);

% Get scanenr fields from struct
scanner = THzScanner(ant);
scanner = getFields(scanner,savedscanner,["app","fig","ant","isApp"]);

% Get target fields from struct
target = THzTarget(wav,ant,scanner);
target = getFields(target,savedtarget,["app","fig","wav","ant","scanner","isApp"]);

% Get im fields from struct
im = THzImageReconstruction(wav,ant,scanner,target);
im = getFields(im,savedim,["app","fig","wav","ant","scanner","target","isApp"]);
im.Get();
im.reconstructor = getFields(im.reconstructor,savedreconstructor,["wav","ant","scanner","target","im"]);

if isDisplay
    ant.Display();
    scanner.Display();
    target.Display();
    im.Display();
end
end