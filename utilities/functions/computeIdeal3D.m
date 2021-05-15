% Computes the 3-D ideal image along XYZ
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

function imIdealXYZ = computeIdeal3D(target,im,o_x,o_y,o_z)
imIdealXYZ = single(zeros(im.numX,im.numY,im.numZ));
x_m = single(reshape(im.x_m,[],1,1));
y_m = single(reshape(im.y_m,1,[],1));
z_m = single(reshape(im.z_m,1,1,[]));
xyz_m = target.xyz_m;
amp = target.amp;

if target.isGPU
    x_m = gpuArray(x_m);
    y_m = gpuArray(y_m);
    z_m = gpuArray(z_m);
    xyz_m = gpuArray(xyz_m);
end

for indTarget = 1:target.numTargets
    temp = single(exp(-(o_x)^(-2)*(x_m-xyz_m(indTarget,1)).^2 -...
    (o_y)^(-2)*(y_m-xyz_m(indTarget,2)).^2 -(o_z)^(-2)*(z_m-xyz_m(indTarget,3)).^2));
    tempMax = max(temp(:));
    if tempMax == 0
        tempMax = 1;
    end
    
    temp = temp*amp(indTarget)/tempMax;
    imIdealXYZ = imIdealXYZ + temp;
end
imIdealXYZ(imIdealXYZ>1) = 1;
imIdealXYZ(imIdealXYZ<0) = 0;

imIdealXYZ = gather(imIdealXYZ);
end