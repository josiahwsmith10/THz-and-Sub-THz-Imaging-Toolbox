% ISOIMAGEXYZDB     Plots a 3-D point cloud image using an isoSurface where
% dBMin is the threshold dB
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

function isoImageXYZdB(imgXYZ,x,y,z,dBMin,titlestr,fontSize)
figure;
h = handle(axes);

% Organize in meshgrid format
imgZXY = permute(imgXYZ,[3,1,2]);

U = reshape(x,1,[],1);
V = reshape(z,[],1,1);
W = reshape(y,1,1,[]);

[meshu,meshv,meshw] = meshgrid(U,V,W);

% Normalize Image
imgZXY = imgZXY/max(imgZXY(:));
imgZXY_dB = db(imgZXY);
clear imgXYZ imgZXY

imgZXY_dB(imgZXY_dB<dBMin) = dBMin;
imgZXY_dB(isnan(imgZXY_dB)) = dBMin;

% Create iso surface
[faces,vertices,colors] = isosurface(meshu,meshv,meshw,imgZXY_dB,dBMin,ones(size(meshu)));
patch('Vertices',vertices,'Faces',faces,'FaceVertexCData',colors,...
    'FaceColor','interp','EdgeColor','interp');

grid on
colormap copper
view(h,3)
daspect(h,[1 1 1])

xlabel(h,'x (m)','fontsize',fontSize,"interpreter","latex")
ylabel(h,'z (m)','fontsize',fontSize,"interpreter","latex")
zlabel(h,'y (m)','fontsize',fontSize,"interpreter","latex")
xlim([x(1) x(end)]); 
zlim([y(1) y(end)]); 
ylim([z(1) z(end)]);
title(h,titlestr,'fontsize',fontSize,"interpreter","latex")

