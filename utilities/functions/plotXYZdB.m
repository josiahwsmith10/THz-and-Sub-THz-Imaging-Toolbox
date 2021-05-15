% PLOTXYZDB     Plots a 3-D point cloud image in dB
%   PLOTXYZDB(imgXYZ,X,Y,Z,zSliceIndex,dBMin,titlestr,fontSize) plots the
%   image in the 3-D array imgXYZ along the dimensions specified in X,Y,Z,
%   where zSliceIndex is the choice of the slices along the Z-dimension
%   specified by the user, dBMin is the minimum threshold in terms of dB,
%   titlestr is the title of the figure, and fontSize is the font size used
%   for the axes and various labels
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

function plotXYZdB(imgXYZ,X,Y,Z,zSliceIndex,dBMin,titlestr,fontSize)
if isempty(zSliceIndex)
    zSliceIndex = 1:length(Z);
end

imgXYZ = double(imgXYZ);
X = double(X);
Y = double(Y);
Z = double(Z);
zSliceIndex = double(zSliceIndex);

f = figure;
h = handle(axes);

% Organize in meshgrid format
imgZXY = permute(imgXYZ,[3,1,2]);

U = reshape(X,1,[],1);
V = reshape(Z,[],1,1);
W = reshape(Y,1,1,[]);

[meshu,meshv,meshw] = meshgrid(U,V,W);

% Normalize Image
imgZXY = imgZXY/max(imgZXY(:));
imgZXY_dB = db(imgZXY);
clear imgXYZ imgZXY

imgZXY_dB(imgZXY_dB<dBMin) = -1e10;
imgZXY_dB(isnan(imgZXY_dB)) = -1e10;

hs = slice(h,meshu,meshv,meshw,imgZXY_dB,[],V(zSliceIndex),[]);
set(hs,'FaceColor','interp','EdgeColor','none');
set(f,'PaperUnits','inches','PaperPosition',[0 0 4 3],'PaperSize',[4 3])
axis(h,'vis3d');

for kk=1:length(zSliceIndex)
    set(hs(kk),'AlphaData',squeeze(imgZXY_dB(kk+zSliceIndex(1)-1,:,:)),'FaceAlpha','interp');
end

colormap(h,'jet')
hc = colorbar(h);

view(h,3)
daspect(h,[1 1 1])
caxis(h,[dBMin 0])

ylabel(hc, 'dB','fontsize',fontSize)
xlabel(h,'x (m)','fontsize',fontSize)
ylabel(h,'z (m)','fontsize',fontSize)
zlabel(h,'y (m)','fontsize',fontSize)
xlim(h,[X(1),X(end)])
ylim(h,[Z(1),Z(end)])
zlim(h,[Y(1),Y(end)])
title(h,titlestr,'fontsize',fontSize)