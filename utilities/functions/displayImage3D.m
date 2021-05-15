% Displays the 3-D reconstructed image
%
% See also PLOTXYZDB.
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

function displayImage3D(im)
h = im.fig.h;
im.vSliceIndex = 1:length(im.z_m);

% Organize in meshgrid format
imgZXY = permute(abs(im.imXYZ),[3,1,2]);

U = reshape(im.x_m,1,[],1);
V = reshape(im.z_m,[],1,1);
W = reshape(im.y_m,1,1,[]);

[meshu,meshv,meshw] = meshgrid(U,V,W);

% Normalize Image
imgZXY = imgZXY/max(imgZXY(:));
imgZXY_dB = db(imgZXY);
clear imgXYZ imgZXY

imgZXY_dB(imgZXY_dB<im.dBMin) = -1e10;
imgZXY_dB(isnan(imgZXY_dB)) = -1e10;

hs = slice(h,meshu,meshv,meshw,imgZXY_dB,single([]),V(im.vSliceIndex),single([]));
set(hs,'FaceColor','interp','EdgeColor','none');
axis(h,'vis3d');

for kk=1:length(im.vSliceIndex)
    set(hs(kk),'AlphaData',squeeze(imgZXY_dB(kk+im.vSliceIndex(1)-1,:,:)),'FaceAlpha','interp');
end

colormap(h,'jet')
hc = colorbar(h);

view(h,3)
daspect(h,[1 1 1])
caxis(h,[im.dBMin 0])

ylabel(hc, 'dB','fontsize',im.fontSize)
xlabel(h,'x (m)','fontsize',im.fontSize)
ylabel(h,'z (m)','fontsize',im.fontSize)
zlabel(h,'y (m)','fontsize',im.fontSize)
xlim(h,[im.x_m(1),im.x_m(end)])
ylim(h,[im.z_m(1),im.z_m(end)])
zlim(h,[im.y_m(1),im.y_m(end)])
title(h,"Reconstructed Image",'fontsize',im.fontSize)
end