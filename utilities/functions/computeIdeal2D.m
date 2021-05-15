% Computes the 2-D ideal image along the dimensions specified
% by the user
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

function imIdeal2D = computeIdeal2D(target,im,dimStr,o_a,o_b)
switch dimStr
    case "XY"
        a_m = im.x_m;
        b_m = im.y_m;
        numA = im.numX;
        numB = im.numY;
        
        indA = 1;
        indB = 2;
    case "YZ"
        a_m = im.y_m;
        b_m = im.z_m;
        numA = im.numY;
        numB = im.numZ;
        
        indA = 2;
        indB = 3;
    case "XZ"
        a_m = im.x_m;
        b_m = im.z_m;
        numA = im.numX;
        numB = im.numZ;
        
        indA = 1;
        indB = 3;
    otherwise
        warning("dimStr must be ""XY"", ""YZ"", or ""YZ"". Cannont continue!");
        imIdeal2D = 0;
        return;
end

imIdeal2D = single(zeros(numA,numB));
a_m = single(reshape(a_m,[],1));
b_m = single(reshape(b_m,1,[]));
xyz_m = target.xyz_m;
amp = target.amp;

if target.isGPU
    a_m = gpuArray(a_m);
    b_m = gpuArray(b_m);
    xyz_m = gpuArray(xyz_m);
end

for indTarget = 1:target.numTargets
    temp = single(exp(-(o_a)^(-2)*(a_m-xyz_m(indTarget,indA)).^2 - ...
        (o_b)^(-2)*(b_m-xyz_m(indTarget,indB)).^2));
    tempMax = max(temp(:));
    if tempMax == 0
        tempMax = 1;
    end
    
    temp = temp*amp(indTarget)/tempMax;
    imIdeal2D = imIdeal2D + temp;
end
imIdeal2D(imIdeal2D>1) = 1;
imIdeal2D(imIdeal2D<0) = 0;

imIdeal2D = gather(imIdeal2D);
end