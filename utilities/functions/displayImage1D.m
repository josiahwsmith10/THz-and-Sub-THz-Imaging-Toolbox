% Displays the reconstructed 1-D image
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

function displayImage1D(im,x_m,xlab)
h = im.fig.h;

% Organize in meshgrid format
imgY = abs(im.imXYZ);

% Normalize Image
imgY = imgY/max(imgY(:));
imgY_dB = db(imgY);
clear imgXYZ imgZXY

plot(h,x_m,imgY_dB)
xlabel(xlab,'fontsize',im.fontSize)
ylabel(h,'dB','fontsize',im.fontSize)
xlim(h,[im.y_m(1),im.y_m(end)])
ylim(h,[im.dBMin 0])
title(h,"Reconstructed Image",'fontsize',im.fontSize)
end