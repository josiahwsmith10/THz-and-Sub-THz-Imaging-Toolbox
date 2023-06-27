% IMAGEINTERP2D     Interpolate the image to the region of interest
% specified by the user
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

function imXYZ = imageInterp3D(obj,x_temp_m,y_temp_m,z_temp_m,sarImage)

x_m = obj.x_m;
y_m = obj.y_m;
z_m = obj.z_m;

[X,Y,Z] = ndgrid(x_m(:),y_m(:),z_m(:));

if obj.im_method ~= "None" && obj.im_method ~= "none"
    % Use interpolation
    imXYZ = single(gather(interpn(x_temp_m(:),y_temp_m(:),z_temp_m(:),sarImage,X,Y,Z,obj.im.im_method,0)));
else
    % Crop image to desired size
    indX = x_temp_m >= y_m(1) & x_temp_m <= x_m(end);
    indY = y_temp_m >= x_m(1) & y_temp_m <= y_m(end);
    indZ = z_temp_m >= z_m(1) & z_temp_m <= z_m(end);

    imXYZ = single(gather(sarImage(indX,indY,indZ)));
    obj.im.x_m = x_temp_m(indX);
    obj.im.y_m = y_temp_m(indY);
    obj.im.z_m = z_temp_m(indZ);
    obj.im.numX = length(obj.im.x_m);
    obj.im.numY = length(obj.im.y_m);
    obj.im.numZ = length(obj.im.z_m);
end