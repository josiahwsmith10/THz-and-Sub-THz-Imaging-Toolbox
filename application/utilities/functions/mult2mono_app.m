% Computes the Y-dimension multistatic-to-monostatic approxmiation. Only
% works for MIMO arrays whose elements are colinear
%
% Also see mult2mono
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

function obj = mult2mono_app(obj,app)
% Y-dimension only multistatic-to-monostatic approximation
if numel(app.sar.sarSize) == 3
    k = reshape(obj.k_vec,1,1,[]);
elseif numel(app.sar.sarSize) == 4
    k = reshape(obj.k_vec,1,1,1,[]);
end

obj.sarData = reshape(obj.sarData,[],size(obj.sarData,3),size(obj.sarData,4),size(obj.sarData,5));

obj.sarData = obj.sarData .* exp(-1j*k .* app.ant.vx.dxy(:,2).^2 / (4*app.ZReferencemEditField.Value));

obj.sarData = reshape(obj.sarData,[],size(obj.sarData,3),size(obj.sarData,4));

end