% Displays the reconstructed 2-D image
%
% See also PLOTXYDB.
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

function displayImage2D(im,x_m,y_m,xlab,ylab)
h = im.fig.h;

% Organize in meshgrid format
imgYX = abs(im.imXYZ).';

% Normalize Image
imgYX = imgYX/max(imgYX(:));
imgYX_dB = db(imgYX);
clear imgXYZ imgZXY

mesh(h,x_m,y_m,imgYX_dB,'FaceColor','interp','EdgeColor','none')

colormap(h,'jet')
hc = colorbar(h);

view(h,2)
caxis(h,[im.dBMin 0])

ylabel(hc, 'dB','fontsize',im.fontSize)
xlabel(h,xlab,'fontsize',im.fontSize)
ylabel(h,ylab,'fontsize',im.fontSize)
xlim(h,[x_m(1),x_m(end)])
ylim(h,[y_m(1),y_m(end)])
title(h,"Reconstructed Image",'fontsize',im.fontSize)
end