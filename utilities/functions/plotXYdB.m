% PLOTXYDB     Plots a 2-D heatmap image in dB
%   PLOTXYDB(imgXY,X,Y,dBMin,xlab,ylab,titlestr,fontSize) plots the image
%   in the 2-D array imgXY along the dimensions specified in X and Y where
%   dBMin is the minimum threshold in terms of dB, xlab is the label for
%   the x-axis, ylab is the label for the y-axis, titlestr is the title of
%   the figure, and fontSize is the font size used for the axes and various
%   labels
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

function plotXYdB(imgXY,X,Y,dBMin,xlab,ylab,titlestr,fontSize,h,isLinear)
% Organize in meshgrid format
imgYX = imgXY.';

if nargin < 9 || isempty(h)
    clf
    h = handle(axes);
end

if nargin < 10
    isLinear = false;
end

% Normalize Image
imgYX = imgYX/max(imgYX(:));
imgYX_dB = db(imgYX);
clear imgXYZ imgZXY

% Remove voxels that are too small
imgYX_dB(imgYX_dB<dBMin) = dBMin;
imgYX_dB(isnan(imgYX_dB)) = dBMin;

if isLinear
    imgYX_dB = 10.^(imgYX_dB/20);
end

% Use mesh to plot the 2-D image
mesh(X,Y,imgYX_dB,'FaceColor','interp','EdgeColor','none')
view(2)
colormap('jet')

% Set the axes and various labels
hc = colorbar;
caxis([dBMin 0])
ylabel(hc,'dB','fontsize',fontSize+5,"Interpreter","latex")

h.FontSize = fontSize;

if isLinear
    ticks = hc.Ticks;
    tickLabels = hc.TickLabels;
    caxis([10^(dBMin/20) 1])
    hc.Ticks = linspace(10^(dBMin/20),1,length(ticks));
    hc.TickLabels = tickLabels;
end

xlabel(xlab,'fontsize',fontSize,"Interpreter","latex")
ylabel(ylab,'fontsize',fontSize,"Interpreter","latex")
xlim([X(1),X(end)])
ylim([Y(1),Y(end)])
title(titlestr,'fontsize',fontSize,"Interpreter","latex")