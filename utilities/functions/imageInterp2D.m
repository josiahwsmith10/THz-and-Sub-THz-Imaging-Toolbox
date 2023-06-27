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

function imUV = imageInterp2D(obj,u_temp_m,v_temp_m,sarImage,xyzStr)

switch xyzStr
    case "xy"
        u_m = obj.x_m;
        v_m = obj.y_m;
    case "yz"
        u_m = obj.y_m;
        v_m = obj.z_m;
    case "xz"
        u_m = obj.x_m;
        v_m = obj.z_m;
end

[U,V] = ndgrid(u_m(:),v_m(:));

if obj.im_method ~= "None" && obj.im_method ~= "none"
    % Use interpolation
    imUV = single(gather(interpn(u_temp_m(:),v_temp_m(:),sarImage,U,V,obj.im.im_method,0)));
else
    % Crop image to desired size
    indU = u_temp_m >= u_m(1) & u_temp_m <= u_m(end);
    indV = v_temp_m >= v_m(1) & v_temp_m <= v_m(end);

    imUV = single(gather(sarImage(indU,indV)));
    switch xyzStr
        case "xy"
            obj.im.x_m = u_temp_m(indU);
            obj.im.y_m = v_temp_m(indV);
            obj.im.numX = length(obj.im.x_m);
            obj.im.numY = length(obj.im.y_m);
        case "yz"
            obj.im.y_m = u_temp_m(indU);
            obj.im.z_m = v_temp_m(indV);
            obj.im.numY = length(obj.im.y_m);
            obj.im.numZ = length(obj.im.z_m);
        case "xz"
            obj.im.x_m = u_temp_m(indU);
            obj.im.z_m = v_temp_m(indV);
            obj.im.numX = length(obj.im.x_m);
            obj.im.numZ = length(obj.im.z_m);
    end
end