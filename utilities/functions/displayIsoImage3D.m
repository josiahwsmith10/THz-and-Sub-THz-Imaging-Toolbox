% Displays the 3-D reconstructed image as an iso surface with threshold of
% im.dBMin
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

function displayIsoImage3D(im)
clf(im.fig.f)
im.fig.h = handle(axes(im.fig.f));
h = im.fig.h;

% Organize in meshgrid format
imgZXY = permute(im.imXYZ,[3,1,2]);

U = reshape(im.x_m,1,[],1);
V = reshape(im.z_m,[],1,1);
W = reshape(im.y_m,1,1,[]);

[meshu,meshv,meshw] = meshgrid(U,V,W);

% Normalize Image
imgZXY = imgZXY/max(imgZXY(:));
imgZXY_dB = db(imgZXY);
clear imgXYZ imgZXY

imgZXY_dB(imgZXY_dB<im.dBMin) = im.dBMin;
imgZXY_dB(isnan(imgZXY_dB)) = im.dBMin;

% Create iso surface
[faces,vertices,colors] = isosurface(meshu,meshv,meshw,imgZXY_dB,im.dBMin,ones(size(meshu)));
patch('Vertices',vertices,'Faces',faces,'FaceVertexCData',colors,...
    'FaceColor','interp','EdgeColor','interp');

grid on
colormap copper
view(h,3)
daspect(h,[1 1 1])

xlabel(h,'x (m)','fontsize',im.fontSize,"interpreter","latex")
ylabel(h,'z (m)','fontsize',im.fontSize,"interpreter","latex")
zlabel(h,'y (m)','fontsize',im.fontSize,"interpreter","latex")
xlim(h,[im.xMin_m,im.xMax_m])
ylim(h,[im.zMin_m,im.zMax_m])
zlim(h,[im.yMin_m,im.yMax_m])
title(h,"Reconstructed Image",'fontsize',im.fontSize,"interpreter","latex")
h.FontSize = im.fontSize;

if im.isApp
    figure(im.fig.f)
    drawnow
    figure(im.app.UIFigure)
end
end

