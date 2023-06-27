% PLOTXYDB     Plots a 1-D image in dB
%   PLOTXYDB(imgX,X,dBMin,xlab,titlestr,fontSize) plots the image in the
%   1-D array imgX along the dimensions specified in X where dBMin is the
%   minimum threshold in terms of dB, xlab is the label for the x-axis,
%   titlestr is the title of the figure, and fontSize is the font size used
%   for the axes and various labels
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

function plotXdB(imgX,X,dBMin,xlab,titlestr,fontSize)
% Organize in meshgrid format
imgX = imgX.';

clf

h = handle(axes);

% Normalize Image
imgX = imgX/max(imgX(:));
imgX_dB = db(imgX);
clear imgXYZ imgZXY

% Remove voxels that are too small
imgX_dB(imgX_dB<dBMin) = dBMin;
imgX_dB(isnan(imgX_dB)) = dBMin;

% Use mesh to plot the 2-D image
plot(X,imgX_dB)
view(2)
colormap('jet')
caxis([dBMin 0])

% Set the axes and various labels
hc = colorbar;
ylabel(hc,'dB','fontsize',fontSize,"Interpreter","latex")

xlabel(xlab,'fontsize',fontSize,"Interpreter","latex")
ylabel("dB",'fontsize',fontSize,"Interpreter","latex")
xlim([X(1),X(end)])
ylim([dBMin,0])
title(titlestr,'fontsize',fontSize,"Interpreter","latex")

h.FontSize = fontSize;