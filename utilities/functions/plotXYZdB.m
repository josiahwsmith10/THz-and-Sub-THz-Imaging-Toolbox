% PLOTXYZDB     Plots a 3-D point cloud image in dB
%   PLOTXYZDB(imgXYZ,X,Y,Z,dBMin,titlestr,fontSize) plots the image in the
%   3-D array imgXYZ along the dimensions specified in X,Y,Z, where dBMin
%   is the minimum threshold in terms of dB, titlestr is the title of the
%   figure, and fontSize is the font size used for the axes and various
%   labels, and sliceAxis is the axis along which to perform the slicing.
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

function plotXYZdB(imgXYZ,X,Y,Z,dBMin,titlestr,fontSize,sliceAxis)
if nargin < 4
    X = 1:size(imXYZ,1);
    Y = 1:size(imXYZ,2);
    Z = 1:size(imXYZ,3);
end

if nargin < 5
    dBMin = -10;
end

if nargin < 6
    titlestr = "";
end

if nargin < 7
    fontSize = 15;
end

if nargin < 8
    sliceAxis = "z";
end

imgXYZ = double(imgXYZ);
X = double(X);
Y = double(Y);
Z = double(Z);

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

imgZXY_dB(imgZXY_dB<dBMin) = dBMin;
imgZXY_dB(isnan(imgZXY_dB)) = dBMin;

if sliceAxis == "x"
    hs = slice(h,meshu,meshv,meshw,imgZXY_dB,U,[],[]);
    for kk=1:length(X)
        set(hs(kk),'AlphaData',squeeze(imgZXY_dB(:,kk,:)),'FaceAlpha','interp');
    end
elseif sliceAxis == "y"
    hs = slice(h,meshu,meshv,meshw,imgZXY_dB,[],[],W);
    for kk=1:length(Y)
        set(hs(kk),'AlphaData',squeeze(imgZXY_dB(:,:,kk)),'FaceAlpha','interp');
    end
elseif sliceAxis == "z"
    hs = slice(h,meshu,meshv,meshw,imgZXY_dB,[],V,[]);
    for kk=1:length(Z)
        set(hs(kk),'AlphaData',squeeze(imgZXY_dB(kk,:,:)),'FaceAlpha','interp');
    end
end

set(hs,'FaceColor','interp','EdgeColor','none');
set(f,'PaperUnits','inches','PaperPosition',[0 0 4 3],'PaperSize',[4 3])
axis(h,'vis3d');

colormap(h,'jet')
hc = colorbar(h);

view(h,3)
daspect(h,[1 1 1])
caxis(h,[dBMin 0])

ylabel(hc, 'dB','fontsize',fontSize+5,"interpreter","latex")
xlabel(h,'x (m)','fontsize',fontSize,"interpreter","latex")
ylabel(h,'z (m)','fontsize',fontSize,"interpreter","latex")
zlabel(h,'y (m)','fontsize',fontSize,"interpreter","latex")
xlim(h,[X(1),X(end)])
ylim(h,[Z(1),Z(end)])
zlim(h,[Y(1),Y(end)])
title(h,titlestr,'fontsize',fontSize,"interpreter","latex")

h.FontSize = fontSize;