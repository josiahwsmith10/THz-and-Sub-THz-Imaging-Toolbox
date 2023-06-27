% Saves all simulation objects to a file to be loaded by THzToolboxLoadAll
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

function THzToolboxSaveAll(wav,ant,scanner,target,im,savePathFull)
if nargin < 6
    savePathFull = [];
end

if exist(savePathFull,'file')
    [filename,pathname] = uiputfile("./*.mat","Select Desired File Location + Name for Save");
    
    if filename == 0
        warning("All scenario file not saved!");
        return;
    else
        savePathFull = string(pathname) + string(filename);
    end
end

% Send wav fields to struct
savedwav = struct();
savedwav = getFields(savedwav,wav,["app","isApp"]);

% Send ant fields to struct
savedant = struct();
savedant = getFields(savedant,ant,["app","fig","wav","isApp"]);

% Send scanner fields to struct
savedscanner = struct();
savedscanner = getFields(savedscanner,scanner,["app","fig","ant","isApp"]);

% Send target fields to struct
savedtarget = struct();
savedtarget = getFields(savedtarget,target,["app","fig","wav","ant","scanner","isApp"]);

% Send reconstructor fields to struct
savedreconstructor = struct();
savedreconstructor = getFields(savedreconstructor,im.reconstructor,["wav","ant","scanner","target","im"]);

% Send im fields to struct
savedim = struct();
savedim = getFields(savedim,im,["app","fig","wav","ant","scanner","target","isApp","reconstructor"]);

save(savePathFull,"savedwav","savedant","savedscanner","savedtarget","savedreconstructor","savedim","-v7.3");

disp("Done saving: " + savePathFull)
end